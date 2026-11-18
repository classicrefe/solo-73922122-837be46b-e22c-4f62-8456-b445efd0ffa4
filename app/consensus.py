"""Molecular consensus construction for a succeeded paired FASTQ task.

The caller supplies, for every read pair of a succeeded task, its read id, an
8-base ACGT UMI and a non-negative reference start.  Read pairs with the same
``(start, umi)`` form an original family; families at the same start whose
UMIs differ by exactly one position may be merged using directional UMI rules
that are decided solely on pre-merge family sizes.  For every final family
both mates are voted position by position using Phred scores.

Nothing in this module writes to the task directory: the published R1/R2
files are opened read-only and the consensus result is assembled in memory as
one ZIP file.  Any validation problem rejects the whole request; no partial
result is ever produced.
"""

from __future__ import annotations

import io
import re
import sqlite3
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from . import config, db
from .fastq import FastqFormatError, _iter_raw_records, extract_read_id

MAX_PAIRS = 500
MIN_READ_LENGTH = 40
MAX_READ_LENGTH = 150
MAX_PHRED = 40
# Longest line the consensus reader will load.  Consensus reads are 40-150
# bases and read ids are short; headers may carry a description, so the bound
# is far above any realistic record but still keeps a damaged task file from
# forcing an unbounded allocation.  Longer lines reject the request as corrupt.
_MAX_LINE_BYTES = 8192

_UMI_RE = re.compile(r"[ACGT]{8}\Z")


class ConsensusError(Exception):
    """A request-level failure; the whole request is rejected."""

    def __init__(self, message: str, status_code: int = 422) -> None:
        super().__init__(message)
        self.message = message
        self.status_code = status_code


@dataclass(frozen=True)
class ReadMeta:
    read_id: str
    umi: str
    start: int


@dataclass
class ReadPair:
    read_id: str
    ordinal: int
    seq1: bytes
    qual1: bytes
    seq2: bytes
    qual2: bytes


def _result_path(task_id: str, relative: str) -> Path:
    """Resolve a published result file, refusing paths outside its task dir."""
    base = (config.TASKS_DIR / task_id).resolve()
    path = (base / relative).resolve()
    if base not in path.parents or not path.is_file():
        raise ConsensusError(
            f"task result file '{relative}' is missing; both paired-end "
            "files must be present",
            409,
        )
    return path


def _parse_task_side(path: Path, label: str) -> dict[str, tuple[bytes, bytes, int]]:
    """Read one published FASTQ and return read id -> (sequence, quality, ordinal).

    Any malformation anywhere in the file - including a truncated or damaged
    final record - rejects the whole consensus request.
    """
    records: dict[str, tuple[bytes, bytes, int]] = {}
    ordinal = 0
    try:
        for header, sequence, _plus, quality in _iter_raw_records(
            path, label, _MAX_LINE_BYTES
        ):
            if ordinal >= MAX_PAIRS:
                raise ConsensusError(
                    f"{label}: task has more than {MAX_PAIRS} read pairs, "
                    "which exceeds the consensus limit"
                )
            read_id = extract_read_id(header.rstrip(b"\r\n"))
            if read_id in records:
                raise ConsensusError(
                    f"{label}: duplicate read id '{read_id}' in the task result",
                    409,
                )
            records[read_id] = (
                sequence.rstrip(b"\r\n"),
                quality.rstrip(b"\r\n"),
                ordinal,
            )
            ordinal += 1
    except FastqFormatError as exc:
        raise ConsensusError(f"corrupt task result: {exc}", 409) from exc
    if not records:
        raise ConsensusError(f"{label}: task result contains no records", 409)
    return records


def _validate_read(meta: ReadMeta, pair: ReadPair) -> None:
    for label, seq, qual in (
        ("R1", pair.seq1, pair.qual1),
        ("R2", pair.seq2, pair.qual2),
    ):
        if not (MIN_READ_LENGTH <= len(seq) <= MAX_READ_LENGTH):
            raise ConsensusError(
                f"read '{meta.read_id}' {label} length {len(seq)} is outside "
                f"{MIN_READ_LENGTH}-{MAX_READ_LENGTH}"
            )
        if any(base not in b"ACGT" for base in seq):
            raise ConsensusError(
                f"read '{meta.read_id}' {label} contains a non-ACGT base"
            )
        if len(qual) != len(seq):
            raise ConsensusError(
                f"read '{meta.read_id}' {label} sequence and quality lengths differ"
            )
        for value in qual:
            score = value - 33
            if score < 0 or score > MAX_PHRED:
                raise ConsensusError(
                    f"read '{meta.read_id}' {label} has quality outside "
                    "Phred+33 (0-40)"
                )


def _prepare_metadata(
    entries: Sequence[ReadMeta],
    known_ids: set[str],
) -> list[ReadMeta]:
    """Validate metadata against the task's read ids (complete coverage)."""
    if len(entries) > MAX_PAIRS:
        raise ConsensusError(
            f"at most {MAX_PAIRS} read pairs may be submitted, got {len(entries)}"
        )
    if not entries:
        raise ConsensusError("reads metadata must not be empty")

    seen: set[str] = set()
    prepared: list[ReadMeta] = []
    for entry in entries:
        read_id = entry.read_id
        if not isinstance(read_id, str) or not read_id:
            raise ConsensusError("each metadata entry needs a non-empty read id")
        if read_id in seen:
            raise ConsensusError(f"duplicate read id in metadata: '{read_id}'")
        seen.add(read_id)
        if not isinstance(entry.umi, str) or not _UMI_RE.fullmatch(entry.umi):
            raise ConsensusError(
                f"read '{read_id}' has invalid UMI; expected 8 ACGT bases"
            )
        if isinstance(entry.start, bool) or not isinstance(entry.start, int):
            raise ConsensusError(
                f"read '{read_id}' has a non-integer reference start"
            )
        if entry.start < 0:
            raise ConsensusError(
                f"read '{read_id}' has a negative reference start"
            )
        if read_id not in known_ids:
            raise ConsensusError(
                f"metadata read id '{read_id}' is not present in the task"
            )
        prepared.append(entry)

    missing = known_ids - seen
    if missing:
        example = sorted(missing)[0]
        raise ConsensusError(
            f"metadata is missing {len(missing)} read id(s) from the task, "
            f"e.g. '{example}'; it must cover every read pair"
        )
    return prepared


def _hamming_one(a: str, b: str) -> bool:
    """True when two equal-length strings differ in exactly one position."""
    if len(a) != len(b):
        return False
    differences = 0
    for left, right in zip(a, b):
        if left != right:
            differences += 1
            if differences > 1:
                return False
    return differences == 1


def _resolve_families(
    families: dict[tuple[int, str], list[ReadPair]],
) -> dict[tuple[int, str], list[ReadPair]]:
    """Apply directional UMI merging and return final family read lists.

    Candidate parents are evaluated with pre-merge family sizes.  A family is
    merged only when exactly one same-start family has a UMI one substitution
    away and a size at least twice its own, and that parent is not itself
    awaiting a merge.  Every qualifying parent is strictly larger than its
    child, so no candidate cycle can exist; processing families from largest
    to smallest makes every parent's "awaiting a merge" answer final before it
    is consulted.
    """
    by_start: dict[int, list[tuple[int, str]]] = {}
    for start, umi in families:
        by_start.setdefault(start, []).append((start, umi))

    ordered = sorted(
        families, key=lambda key: (-len(families[key]), key[0], key[1])
    )
    parent: dict[tuple[int, str], tuple[int, str]] = {}
    pending: dict[tuple[int, str], bool] = {}
    for key in ordered:
        start, umi = key
        count = len(families[key])
        candidates = [
            other
            for other in by_start[start]
            if other != key
            and len(families[other]) >= 2 * count
            and _hamming_one(umi, other[1])
        ]
        if len(candidates) == 1 and not pending.get(candidates[0], False):
            parent[key] = candidates[0]
            pending[key] = True
        else:
            pending[key] = False

    def root(key: tuple[int, str]) -> tuple[int, str]:
        seen: set[tuple[int, str]] = set()
        while key in parent:
            if key in seen:  # defensive; impossible, sizes strictly increase
                break
            seen.add(key)
            key = parent[key]
        return key

    merged: dict[tuple[int, str], list[ReadPair]] = {}
    for key, reads in families.items():
        merged.setdefault(root(key), []).extend(reads)
    return merged

def _check_family_lengths(
    key: tuple[int, str], reads: Sequence[ReadPair]
) -> None:
    """Each end must have one common length inside a (possibly merged) family."""
    for label, attr in (("R1", "seq1"), ("R2", "seq2")):
        lengths = {len(getattr(pair, attr)) for pair in reads}
        if len(lengths) != 1:
            raise ConsensusError(
                f"family at start {key[0]} / UMI {key[1]} has inconsistent "
                f"{label} lengths within the family: {sorted(lengths)}"
            )


def _vote(
    reads: Sequence[ReadPair], seq_attr: str, qual_attr: str
) -> tuple[bytes, bytes]:
    """Vote one mate position by position using summed Phred scores.

    Only bases actually observed at a position receive votes.  The winning
    base's quality is the winner's total minus the runner-up's total (0 when
    no other base was observed), capped at 40.  An exact tie for the highest
    total yields ``N`` and ``!``.
    """
    length = len(getattr(reads[0], seq_attr))
    sequence = bytearray()
    quality = bytearray()
    for position in range(length):
        scores: dict[bytes, int] = {}
        for pair in reads:
            base = getattr(pair, seq_attr)[position : position + 1]
            scores[base] = scores.get(base, 0) + (
                getattr(pair, qual_attr)[position] - 33
            )
        best = max(scores.values())
        winners = [base for base, score in scores.items() if score == best]
        if len(winners) != 1:
            # Tie for the highest total: ambiguous, emit N with '!'.
            sequence += b"N"
            quality.append(ord("!"))
            continue
        winner = winners[0]
        second = max(
            (score for base, score in scores.items() if base != winner),
            default=0,
        )
        phred = min(MAX_PHRED, best - second)
        sequence += winner
        quality.append(phred + 33)
    return bytes(sequence), bytes(quality)


def _fastq_record(read_id: str, sequence: bytes, quality: bytes) -> bytes:
    return (
        b"@" + read_id.encode("latin-1") + b"\n"
        + sequence + b"\n+\n" + quality + b"\n"
    )


def build_consensus_zip(
    conn: sqlite3.Connection,
    task_id: str,
    entries: Sequence[ReadMeta],
) -> bytes:
    """Validate the request and return both consensus FASTQs as one ZIP."""
    task = db.get_task(conn, task_id)
    if task is None:
        raise ConsensusError("task not found", 404)
    if task["status"] != "succeeded":
        raise ConsensusError(
            f"task is {task['status']}; consensus requires a succeeded task", 409
        )

    r1_path = _result_path(task_id, task["result_r1"])
    r2_path = _result_path(task_id, task["result_r2"])

    r1 = _parse_task_side(r1_path, "R1")
    r2 = _parse_task_side(r2_path, "R2")
    if set(r1) != set(r2):
        missing = sorted(set(r1) ^ set(r2))[0]
        raise ConsensusError(
            f"task result files are incomplete: read '{missing}' is present on "
            "only one end",
            409,
        )
    if len(r1) > MAX_PAIRS:
        raise ConsensusError(
            f"task has more than {MAX_PAIRS} read pairs, which exceeds the "
            "consensus limit"
        )

    metadata = _prepare_metadata(entries, set(r1))

    pairs: dict[str, ReadPair] = {}
    for meta in metadata:
        seq1, qual1, ordinal = r1[meta.read_id]
        seq2, qual2, _ = r2[meta.read_id]
        pair = ReadPair(meta.read_id, ordinal, seq1, qual1, seq2, qual2)
        _validate_read(meta, pair)
        pairs[meta.read_id] = pair

    families: dict[tuple[int, str], list[ReadPair]] = {}
    for meta in metadata:
        families.setdefault((meta.start, meta.umi), []).append(pairs[meta.read_id])

    # Check the original families first, then the units produced by merging;
    # a family whose ends do not share one length cannot be voted.
    for key, reads in families.items():
        _check_family_lengths(key, reads)

    merged = _resolve_families(families)
    for key, reads in merged.items():
        _check_family_lengths(key, reads)

    output_rows: list[tuple[int, str, bytes, bytes, bytes, bytes]] = []
    for key, reads in merged.items():
        earliest = min(reads, key=lambda pair: pair.ordinal)
        seq1, qual1 = _vote(reads, "seq1", "qual1")
        seq2, qual2 = _vote(reads, "seq2", "qual2")
        output_rows.append(
            (earliest.ordinal, earliest.read_id, seq1, qual1, seq2, qual2)
        )
    output_rows.sort(key=lambda row: row[0])

    r1_buffer = bytearray()
    r2_buffer = bytearray()
    for _ordinal, read_id, seq1, qual1, seq2, qual2 in output_rows:
        r1_buffer += _fastq_record(read_id, seq1, qual1)
        r2_buffer += _fastq_record(read_id, seq2, qual2)

    archive = io.BytesIO()
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("R1.fastq", bytes(r1_buffer))
        zf.writestr("R2.fastq", bytes(r2_buffer))
    return archive.getvalue()
