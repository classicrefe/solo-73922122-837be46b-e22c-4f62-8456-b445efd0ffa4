"""Self-tests for the molecular consensus API.

Run with the project virtualenv:
    .venv/bin/python -I tests/test_consensus.py

The tests upload real paired FASTQ tasks through the FastAPI app, then call
POST /tasks/{id}/consensus and inspect the returned ZIP byte-for-byte.
"""

from __future__ import annotations

import asyncio
import hashlib
import io
import os
import shutil
import sys
import tempfile
import time
import uuid
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import httpx  # noqa: E402

from app import config, db, worker  # noqa: E402
from app.main import app  # noqa: E402

RESULTS: list[tuple[str, bool, str]] = []
LEN = 40


def record(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((name, ok, detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" - {detail}" if detail else ""))


def q(score: int, length: int = LEN) -> bytes:
    return bytes([score + 33]) * length


def rec(rid: str, seq: bytes, qual: bytes, side: int) -> bytes:
    return (
        b"@" + rid.encode() + b"/" + str(side).encode() + b"\n"
        + seq + b"\n+\n" + qual + b"\n"
    )


def fastq_bytes(records: list[bytes]) -> io.BytesIO:
    return io.BytesIO(b"".join(records))


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_fastq(data: bytes) -> list[tuple[str, bytes, bytes]]:
    lines = data.splitlines(keepends=True)
    assert len(lines) % 4 == 0, f"records do not span four lines: {len(lines)}"
    out = []
    for i in range(0, len(lines), 4):
        header = lines[i].rstrip(b"\r\n")
        out.append((header[1:].decode(), lines[i + 1].rstrip(b"\r\n"),
                    lines[i + 3].rstrip(b"\r\n")))
    return out


async def wait_terminal(client: httpx.AsyncClient, task_id: str,
                        timeout: float = 30.0) -> dict:
    deadline = time.time() + timeout
    while time.time() < deadline:
        body = (await client.get(f"/tasks/{task_id}")).json()
        if body["status"] in ("succeeded", "failed", "cancelled"):
            return body
        await asyncio.sleep(0.05)
    raise AssertionError(f"task {task_id} did not finish in {timeout}s")


async def upload(client: httpx.AsyncClient, r1: list[bytes],
                 r2: list[bytes]) -> str:
    r = await client.post(
        "/tasks",
        files={"r1": ("r1.fastq", fastq_bytes(r1), "application/octet-stream"),
               "r2": ("r2.fastq", fastq_bytes(r2), "application/octet-stream")},
    )
    assert r.status_code == 202, r.text
    body = await wait_terminal(client, r.json()["id"])
    assert body["status"] == "succeeded", body.get("message")
    return body["id"]


async def upload_pairs(
    client: httpx.AsyncClient,
    rows: list[tuple[str, bytes, bytes, bytes, bytes]],
) -> str:
    """rows: (id, R1 seq, R1 qual, R2 seq, R2 qual) in R1 order."""
    r1 = [rec(rid, s1, q1, 1) for rid, s1, q1, _s2, _q2 in rows]
    r2 = [rec(rid, s2, q2, 2) for rid, _s1, _q1, s2, q2 in rows]
    return await upload(client, r1, r2)


async def consensus(client: httpx.AsyncClient, task_id: str,
                    reads: list[dict]) -> httpx.Response:
    return await client.post(f"/tasks/{task_id}/consensus", json={"reads": reads})


def zip_reads(response: httpx.Response) -> tuple[list, list]:
    assert response.status_code == 200, response.text
    assert response.headers["content-type"].startswith("application/zip")
    zf = zipfile.ZipFile(io.BytesIO(response.content))
    assert sorted(zf.namelist()) == ["R1.fastq", "R2.fastq"], zf.namelist()
    return parse_fastq(zf.read("R1.fastq")), parse_fastq(zf.read("R2.fastq"))


async def main_async() -> int:
    tmp_root = Path(tempfile.mkdtemp(prefix="consensus-test-",
                                     dir=str(ROOT / "tests" / "output")))
    data_dir = tmp_root / "data"
    os.environ["FASTQ_DATA_DIR"] = str(data_dir)
    config.DATA_DIR = data_dir
    config.DB_PATH = data_dir / "tasks.db"
    config.TASKS_DIR = data_dir / "tasks"

    db._initialized = False
    db.init_db()
    worker.recover_on_startup()

    transport = httpx.ASGITransport(app=app)
    async with httpx.AsyncClient(transport=transport, base_url="http://test",
                                 timeout=60) as client:
        # ---------------------------------------------------------------
        # Main task: 8 reads across four families.
        #   start 10 / UMI AAAAAAAA : r1,r2,r3  (+ r4 AAAAAAAC merges in)
        #   start 20 / UMI CCCCCCCC : r5,r6
        #   start 30 / UMI GGGGGGGG : r7,r8
        # ---------------------------------------------------------------
        a40, c40 = b"A" * LEN, b"C" * LEN
        t40 = b"T" * LEN
        main_rows = [
            ("r1", a40, q(40), t40, q(40)),
            ("r2", a40, q(40), t40, q(40)),
            ("r3", c40[:1] + a40[1:], q(40), t40, q(40)),
            ("r4", a40, q(40), c40[:1] + t40[1:], q(40)),
            ("r5", a40, q(40), t40, q(40)),
            ("r6", c40, q(30), t40, q(30)),
            ("r7", a40, q(40), t40, q(40)),
            ("r8", c40[:1] + a40[1:], q(40), t40, q(30)),
        ]
        main_id = await upload_pairs(client, main_rows)
        out_dir = config.TASKS_DIR / main_id / "output"
        before = {p.name: sha(p) for p in out_dir.iterdir()}
        dir_listing = sorted(p.name for p in (config.TASKS_DIR / main_id).iterdir())

        main_meta = [
            {"read_id": "r1", "umi": "AAAAAAAA", "start": 10},
            {"read_id": "r2", "umi": "AAAAAAAA", "start": 10},
            {"read_id": "r3", "umi": "AAAAAAAA", "start": 10},
            {"read_id": "r4", "umi": "AAAAAAAC", "start": 10},
            {"read_id": "r5", "umi": "CCCCCCCC", "start": 20},
            {"read_id": "r6", "umi": "CCCCCCCC", "start": 20},
            {"read_id": "r7", "umi": "GGGGGGGG", "start": 30},
            {"read_id": "r8", "umi": "GGGGGGGG", "start": 30},
        ]
        r = await consensus(client, main_id, main_meta)
        r1_out, r2_out = zip_reads(r)
        record("consensus returns one ZIP with R1.fastq and R2.fastq",
               len(r.content) > 0)
        record("consensus output ids follow original R1 order",
               [rid for rid, _, _ in r1_out] == ["r1", "r5", "r7"],
               str([rid for rid, _, _ in r1_out]))
        record("merged family uses its earliest read id",
               [rid for rid, _, _ in r2_out] == ["r1", "r5", "r7"])
        # F1 R1: r1 A40,r2 A40,r3 C at pos0 q40,r4 A40 -> A wins by 80.
        record("R1 position vote: majority base wins with capped quality",
               r1_out[0][1] == a40 and r1_out[0][2] == q(40),
               f"{r1_out[0][1][:4]!r} {r1_out[0][2][:4]!r}")
        record("R2 vote merged family", r2_out[0][1] == t40 and r2_out[0][2] == q(40),
               f"{r2_out[0][1][:4]!r}")
        # F2 (r5,r6) R1: A(40) vs C(30) at every position -> A with
        # difference 10 -> '+', never capped.
        record("quality difference drives winner quality",
               r1_out[1][1] == a40 and r1_out[1][2] == b"+" * LEN,
               f"{r1_out[1][1][:2]!r} {r1_out[1][2][:2]!r}")
        # F3 (r7,r8) R1: pos0 A(40) vs C(40) tie -> N with '!'; rest A 80 -> I.
        record("exact tie becomes N and '!'",
               r1_out[2][1] == b"N" + a40[1:] and
               r1_out[2][2] == b"!" + q(40)[1:],
               f"{r1_out[2][1][:2]!r} {r1_out[2][2][:2]!r}")
        record("consensus sequences/qualities have equal length",
               all(len(s) == len(qq) for _, s, qq in r1_out + r2_out))
        record("task originals unchanged after consensus",
               {p.name: sha(p) for p in out_dir.iterdir()} == before)
        record("consensus writes no files into the task directory",
               sorted(p.name for p in (config.TASKS_DIR / main_id).iterdir())
               == dir_listing)

        # ---------------------------------------------------------------
        # Request validation: whole request rejected, originals untouched.
        # ---------------------------------------------------------------
        async def expect_fail(name: str, response: httpx.Response,
                              code: int, needle: str = "") -> None:
            ok = response.status_code == code
            detail = ""
            if ok and needle:
                try:
                    detail = response.json().get("detail", "")
                except ValueError:
                    detail = ""
                ok = needle in detail
            record(name, ok, f"{response.status_code} {detail or response.text[:120]}")
            record(f"{name}: originals unchanged",
                   {p.name: sha(p) for p in out_dir.iterdir()} == before)

        missing = [m for m in main_meta if m["read_id"] != "r8"]
        await expect_fail("metadata missing a read id -> 422",
                          await consensus(client, main_id, missing), 422, "missing")

        dup = main_meta + [main_meta[0]]
        await expect_fail("metadata duplicate read id -> 422",
                          await consensus(client, main_id, dup), 422, "duplicate")

        unknown = [dict(m) for m in main_meta]
        unknown[-1]["read_id"] = "nope"
        await expect_fail("metadata read id absent from task -> 422",
                          await consensus(client, main_id, unknown), 422,
                          "not present")

        for label, mutate, needle in [
            ("lowercase UMI", lambda m: m.update(umi="aaaaaaaA"), "invalid UMI"),
            ("7-base UMI", lambda m: m.update(umi="AAAAAAA"), "invalid UMI"),
            ("UMI with N", lambda m: m.update(umi="AAAAAAAN"), "invalid UMI"),
            ("negative start", lambda m: m.update(start=-1), "negative"),
        ]:
            bad = [dict(m) for m in main_meta]
            mutate(bad[0])
            await expect_fail(f"{label} -> 422",
                              await consensus(client, main_id, bad), 422, needle)

        bad_start = [dict(m) for m in main_meta]
        bad_start[0]["start"] = "not-a-number"
        await expect_fail("non-integer start -> 422",
                          await consensus(client, main_id, bad_start), 422)

        await expect_fail("empty metadata -> 422",
                          await consensus(client, main_id, []), 422, "empty")

        # not succeeded / unknown task
        orphan_id = (await client.post(
            "/tasks",
            files={"r1": ("r1.fastq", fastq_bytes([rec("x", a40, q(40), 1)]),
                          "application/octet-stream"),
                   "r2": ("r2.fastq", io.BytesIO(b""),
                          "application/octet-stream")},
        )).json()["id"]
        orphan = await wait_terminal(client, orphan_id)
        assert orphan["status"] == "failed"
        await expect_fail("failed task -> 409",
                          await consensus(client, orphan_id, []), 409, "succeeded")
        await expect_fail("unknown task -> 404",
                          await consensus(client, uuid.uuid4().hex, []), 404)

        # ---------------------------------------------------------------
        # Family-level length consistency (per end).
        # ---------------------------------------------------------------
        len_rows = [
            ("l1", a40, q(40), t40, q(40)),
            ("l2", a40 + b"A", q(40, 41), t40, q(40)),
        ]
        len_id = await upload_pairs(client, len_rows)
        len_meta = [
            {"read_id": "l1", "umi": "AAAAAAAA", "start": 5},
            {"read_id": "l2", "umi": "AAAAAAAA", "start": 5},
        ]
        r = await consensus(client, len_id, len_meta)
        record("inconsistent R1 length inside family -> 422",
               r.status_code == 422 and "inconsistent R1 lengths" in r.json()["detail"],
               f"{r.status_code} {r.text[:160]}")

        r2len_rows = [
            ("m1", a40, q(40), t40, q(40)),
            ("m2", a40, q(40), t40 + b"T", q(40, 41)),
        ]
        r2len_id = await upload_pairs(client, r2len_rows)
        r = await consensus(client, r2len_id, [
            {"read_id": "m1", "umi": "AAAAAAAA", "start": 5},
            {"read_id": "m2", "umi": "AAAAAAAA", "start": 5},
        ])
        record("inconsistent R2 length inside family -> 422",
               r.status_code == 422 and "inconsistent R2 lengths" in r.json()["detail"],
               f"{r.status_code} {r.text[:160]}")

        # ---------------------------------------------------------------
        # Per-read input constraints: 40-150 bases ACGT, Phred+33 0-40.
        # ---------------------------------------------------------------
        tiny = [
            ("s1", a40[:39], q(40, 39), t40[:39], q(40, 39)),
            ("s2", a40 + b"A", q(40, 41), t40, q(40)),
        ]
        tiny_id = await upload_pairs(client, tiny)
        r = await consensus(client, tiny_id, [
            {"read_id": "s1", "umi": "AAAAAAAA", "start": 1},
            {"read_id": "s2", "umi": "CCCCCCCC", "start": 2},
        ])
        record("read shorter than 40 -> 422",
               r.status_code == 422 and "outside 40-150" in r.json()["detail"],
               f"{r.status_code} {r.text[:160]}")

        long_rows = [
            ("t1", a40 * 3 + b"A" * 30, q(40, 150), t40 * 3 + b"T" * 30,
             q(40, 150)),
            ("t2", a40 * 4 + b"A" * 11, q(40, 171), t40 * 4 + b"T" * 11,
             q(40, 171)),
        ]
        long_id = await upload_pairs(client, long_rows)
        r = await consensus(client, long_id, [
            {"read_id": "t1", "umi": "AAAAAAAA", "start": 1},
            {"read_id": "t2", "umi": "TTTTTTTT", "start": 2},
        ])
        record("read longer than 150 -> 422",
               r.status_code == 422 and "outside 40-150" in r.json()["detail"],
               f"{r.status_code} {r.text[:160]}")

        longest_rows = [("t3", a40 * 3 + b"A" * 30, q(40, 150),
                         t40 * 3 + b"T" * 30, q(40, 150))]
        longest_id = await upload_pairs(client, longest_rows)
        r = await consensus(client, longest_id, [
            {"read_id": "t3", "umi": "AAAAAAAA", "start": 1},
        ])
        longest_r1, _ = zip_reads(r)
        record("150-base reads are accepted",
               len(longest_r1[0][1]) == 150 and longest_r1[0][2] == q(40, 150))

        n_rows = [("n1", b"N" * LEN, q(40), t40, q(40))]
        n_id = await upload_pairs(client, n_rows)
        r = await consensus(client, n_id, [
            {"read_id": "n1", "umi": "AAAAAAAA", "start": 1},
        ])
        record("non-ACGT base -> 422",
               r.status_code == 422 and "non-ACGT" in r.json()["detail"],
               f"{r.status_code} {r.text[:160]}")

        hiq_rows = [("h1", a40, q(41), t40, q(40))]
        hiq_id = await upload_pairs(client, hiq_rows)
        r = await consensus(client, hiq_id, [
            {"read_id": "h1", "umi": "AAAAAAAA", "start": 1},
        ])
        record("quality above Phred 40 -> 422",
               r.status_code == 422 and "Phred" in r.json()["detail"],
               f"{r.status_code} {r.text[:160]}")

        okq_rows = [("k1", a40, q(0), t40, q(40))]
        okq_id = await upload_pairs(client, okq_rows)
        r = await consensus(client, okq_id, [
            {"read_id": "k1", "umi": "AAAAAAAA", "start": 0},
        ])
        good_r1, _ = zip_reads(r)
        # A single observed base wins unconditionally; runner-up sum is 0, so
        # the quality is the read's own Phred 0 -> '!'.
        record("Phred 0 and start 0 are accepted",
               good_r1[0][1] == a40 and good_r1[0][2] == q(0))
        record("single read consensus keeps its bases at quality 0",
               good_r1[0][1] == a40 and good_r1[0][2] == b"!" * LEN)

        # A one-read family at high quality is preserved as-is.
        solo_rows = [("zo", a40, q(40), t40, q(40))]
        solo_id = await upload_pairs(client, solo_rows)
        solo_r1, solo_r2 = zip_reads(await consensus(client, solo_id, [
            {"read_id": "zo", "umi": "AAAAAAAA", "start": 9},
        ]))
        record("single-read family is copied through unchanged",
               solo_r1[0][1] == a40 and solo_r1[0][2] == q(40)
               and solo_r2[0][1] == t40 and solo_r2[0][2] == q(40))

        # ---------------------------------------------------------------
        # UMI merging rules.
        # ---------------------------------------------------------------
        merge_rows = []
        merge_meta = []

        def add(rid: str, umi: str, start: int, n: int) -> None:
            for i in range(n):
                read_id = f"{rid}{i}"
                merge_rows.append((read_id, a40, q(40), t40, q(40)))
                merge_meta.append({"read_id": read_id, "umi": umi, "start": start})

        add("a", "AAAAAAAA", 100, 2)   # big (2)
        add("b", "AAAAAAAT", 100, 1)   # 1 diff, size 1 -> merges into a
        add("s", "CCCCCCCG", 200, 1)   # ambiguous: two parents each >= 2
        add("p", "CCCCCCCC", 200, 2)
        add("q", "CCCCCCCA", 200, 2)
        add("u", "GGGGGGGT", 300, 1)   # its only parent is itself merging
        add("v", "GGGGGGGG", 300, 2)   # merges into w
        add("w", "GGGGGGGA", 300, 4)
        add("x", "TTTTTTTT", 400, 3)   # 3 < 2*2 -> no merge
        add("y", "TTTTTTTA", 400, 2)
        add("z", "ACGTACGT", 500, 1)   # unique 1-diff parent -> merges
        add("zb", "ACGTACGA", 500, 2)
        add("zc", "ACGTTCGA", 500, 2)  # 2 diffs from z, not a candidate

        merge_id = await upload_pairs(client, merge_rows)
        r = await consensus(client, merge_id, merge_meta)
        mr1, _mr2 = zip_reads(r)
        out_ids = [rid for rid, _, _ in mr1]
        record("unique family at least twice the size absorbs small family",
               "a0" in out_ids and "b0" not in out_ids, str(out_ids))
        record("ambiguous small family is not merged",
               sum(1 for rid in out_ids if rid.startswith("s")) == 1
               and "p0" in out_ids and "q0" in out_ids, str(out_ids))
        record("family awaiting its own merge is not a parent",
               "u0" in out_ids and "v0" in out_ids and "w0" not in out_ids,
               str(out_ids))
        record("parent smaller than twice the child -> no merge",
               "x0" in out_ids and "y0" in out_ids, str(out_ids))
        record("unique sufficiently large parent merges despite a 2-diff family",
               "z0" in out_ids and "zb0" not in out_ids and "zc0" in out_ids,
               str(out_ids))
        record("merge output stays in original R1 order",
               out_ids == ["a0", "s0", "p0", "q0", "u0", "v0",
                           "x0", "y0", "z0", "zc0"], str(out_ids))
        # Merged families: a+b (3 reads), v+w (6 reads), z+zb (3 reads); the
        # rest stay alone.  The representative is always the earliest read id
        # in the final family (v0 for v+w), even though w is the parent UMI.
        record("merged families vote across all members",
               len(out_ids) == 10 and all(seq == a40 and qual == q(40)
                                          for _, seq, qual in mr1),
               f"rows={len(out_ids)}")

        # ---------------------------------------------------------------
        # Task result integrity: missing mate file / damaged tail.
        # ---------------------------------------------------------------
        small_rows = [("e1", a40, q(40), t40, q(40))]
        small_meta = [{"read_id": "e1", "umi": "AAAAAAAA", "start": 1}]
        small_id = await upload_pairs(client, small_rows)
        out = config.TASKS_DIR / small_id / "output"
        r2_file = out / "R2.fastq"
        original = r2_file.read_bytes()
        r2_file.unlink()
        r = await consensus(client, small_id, small_meta)
        record("missing R2 result file -> 409",
               r.status_code == 409 and "missing" in r.json()["detail"],
               f"{r.status_code} {r.text[:160]}")
        r2_file.write_bytes(original)

        r1_file = out / "R1.fastq"
        r1_original = r1_file.read_bytes()
        r1_file.unlink()
        r = await consensus(client, small_id, small_meta)
        record("missing R1 result file -> 409", r.status_code == 409,
               f"{r.status_code} {r.text[:160]}")
        r1_file.write_bytes(r1_original)

        r2_file.write_bytes(original + b"@broken\nACGT\n")  # truncated tail
        r = await consensus(client, small_id, small_meta)
        record("damaged trailing record in R2 -> 409",
               r.status_code == 409 and "corrupt" in r.json()["detail"],
               f"{r.status_code} {r.text[:160]}")
        r2_file.write_bytes(original)

        # A huge unterminated line must be rejected without loading it whole.
        r2_file.write_bytes(original + b"@" + b"A" * (2 * 1024 * 1024))
        r = await consensus(client, small_id, small_meta)
        record("oversized damaged line in R2 -> 409",
               r.status_code == 409 and "exceeds" in r.json()["detail"],
               f"{r.status_code} {r.text[:160]}")
        r1_file.write_bytes(r1_original)  # R1 untouched by the reader
        r2_file.write_bytes(original)

        r = await consensus(client, small_id, small_meta)
        good_r1, good_r2 = zip_reads(r)
        record("result files usable again once restored",
               good_r1[0][0] == "e1" and good_r2[0][0] == "e1")

        # ---------------------------------------------------------------
        # Pair limit: more than 500 pairs.
        # ---------------------------------------------------------------
        big_rows = [(f"big{i}", a40, q(40), t40, q(40)) for i in range(501)]
        big_id = await upload_pairs(client, big_rows)
        big_meta = [{"read_id": f"big{i}", "umi": "AAAAAAAA", "start": 1}
                    for i in range(501)]
        r = await consensus(client, big_id, big_meta)
        record("more than 500 read pairs -> 422",
               r.status_code == 422 and "500" in r.json()["detail"],
               f"{r.status_code} {r.text[:160]}")

        # Exactly 500 pairs is accepted.
        limit_rows = [(f"lim{i}", a40, q(40), t40, q(40)) for i in range(500)]
        limit_id = await upload_pairs(client, limit_rows)
        limit_meta = [{"read_id": f"lim{i}", "umi": "AAAAAAAA", "start": 1}
                      for i in range(500)]
        r = await consensus(client, limit_id, limit_meta)
        limit_r1, _ = zip_reads(r)
        record("exactly 500 read pairs is accepted and voted as one family",
               len(limit_r1) == 1 and limit_r1[0][1] == a40
               and limit_r1[0][2] == q(40),
               f"rows={len(limit_r1)}")

    shutil.rmtree(tmp_root, ignore_errors=True)
    failed = [name for name, ok, _ in RESULTS if not ok]
    print("\n%d/%d checks passed" % (len(RESULTS) - len(failed), len(RESULTS)))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(asyncio.run(main_async()))
