#!/usr/bin/env bash
# End-to-end demo against a running service (default http://localhost:8000).
# Shows: upload -> progress polling -> download verification, plus cancellation
# and a rejected (malformed) upload.
set -euo pipefail

BASE_URL="${BASE_URL:-http://localhost:8000}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

need() { command -v "$1" >/dev/null || { echo "missing tool: $1" >&2; exit 1; }; }
need curl
need python3

echo "== health =="
curl -fsS "$BASE_URL/health"
echo

mk() { # file id seq side
  printf '@%s/%s desc\n%s\n+\n%s\n' "$2" "$4" "$3" "$(printf '~%.0s' $(seq ${#3}))" >> "$1"
}

echo "== generating shuffled R1/R2 =="
R1="$WORK/R1.fastq"; R2="$WORK/R2.fastq"
for id in alpha beta gamma delta epsilon; do
  mk "$R1" "$id" "ACGTACGTAC" 1
done
for id in epsilon gamma alpha delta beta; do   # deliberately shuffled
  mk "$R2" "$id" "TGCATGCATG" 2
done

echo "== upload =="
TASK=$(curl -fsS -X POST "$BASE_URL/tasks" \
  -F "r1=@$R1;filename=R1.fastq" \
  -F "r2=@$R2;filename=R2.fastq")
echo "$TASK"
ID=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["id"])' <<<"$TASK")

echo "== poll until terminal =="
while :; do
  BODY=$(curl -fsS "$BASE_URL/tasks/$ID")
  STATUS=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["status"])' <<<"$BODY")
  echo "$BODY"
  case "$STATUS" in
    succeeded|failed|cancelled) break ;;
  esac
  sleep 0.5
done

echo "== download and verify order =="
curl -fsS "$BASE_URL/tasks/$ID/download/r1" -o "$WORK/out_R1.fastq"
curl -fsS "$BASE_URL/tasks/$ID/download/r2" -o "$WORK/out_R2.fastq"
python3 - "$WORK/out_R1.fastq" "$WORK/out_R2.fastq" <<'PY'
import sys
ids = []
for path in sys.argv[1:]:
    got = []
    with open(path, "rb") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            got.append(header[1:].split()[0].decode().rsplit("/", 1)[0])
            fh.readline(); fh.readline(); fh.readline()
    ids.append(got)
print("R1 ids:", ids[0])
print("R2 ids:", ids[1])
assert ids[0] == ids[1] == ["alpha", "beta", "gamma", "delta", "epsilon"]
print("OK: paired in R1 order")
PY

echo "== malformed input is rejected as a whole task =="
printf '@x\nACGT\n+\n~~\n' > "$WORK/bad_R2.fastq"     # seq 4 != qual 2
BAD=$(curl -fsS -X POST "$BASE_URL/tasks" \
  -F "r1=@$R1;filename=R1.fastq" \
  -F "r2=@$WORK/bad_R2.fastq;filename=R2.fastq")
BAD_ID=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["id"])' <<<"$BAD")
sleep 1
curl -fsS "$BASE_URL/tasks/$BAD_ID"
echo
echo "(download attempt on failed task -> HTTP 409)"
curl -s -o /dev/null -w '%{http_code}\n' "$BASE_URL/tasks/$BAD_ID/download/r1"

echo "== cancellation returns immediately, task never publishes =="
BIG_R1="$WORK/big_R1.fastq"; BIG_R2="$WORK/big_R2.fastq"
python3 - "$BIG_R1" "$BIG_R2" <<'PY'
import sys
n = 200_000
with open(sys.argv[1], "wb") as f1, open(sys.argv[2], "wb") as f2:
    for i in range(n):
        f1.write(b"@r%09d/1\n" % i + b"A" * 40 + b"\n+\n" + b"~" * 40 + b"\n")
    for i in range(n - 1, -1, -1):
        f2.write(b"@r%09d/2\n" % i + b"T" * 40 + b"\n+\n" + b"!" * 40 + b"\n")
PY
CANCEL=$(curl -fsS -X POST "$BASE_URL/tasks" \
  -F "r1=@$BIG_R1;filename=big_R1.fastq" \
  -F "r2=@$BIG_R2;filename=big_R2.fastq")
CANCEL_ID=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["id"])' <<<"$CANCEL")
curl -fsS -X POST "$BASE_URL/tasks/$CANCEL_ID/cancel"
echo
sleep 2
curl -fsS "$BASE_URL/tasks/$CANCEL_ID"
echo

echo "== cancelling the already succeeded task -> HTTP 409 =="
curl -s -o /dev/null -w '%{http_code}\n' -X POST "$BASE_URL/tasks/$ID/cancel"

echo "== molecular consensus: UMI families -> ZIP with R1/R2 =="
# Build a consensus-valid task (40 bp reads, Phred 0-40) and merge the delta
# family (1 read) into gamma (2 reads; UMI one substitution apart, size at
# least twice).  epsilon stays a family of its own.
CONS_R1="$WORK/cons_R1.fastq"; CONS_R2="$WORK/cons_R2.fastq"
python3 - "$CONS_R1" "$CONS_R2" <<'PY'
import sys
seqs = {"gamma": b"A" * 40, "delta": b"G" * 40, "epsilon": b"T" * 40}
with open(sys.argv[1], "wb") as f1, open(sys.argv[2], "wb") as f2:
    for rid, seq in seqs.items():
        f1.write(b"@%s/1\n%s\n+\n%s\n" % (rid.encode(), seq, b"I" * 40))
        f2.write(b"@%s/2\n%s\n+\n%s\n" % (rid.encode(), b"C" * 40, b"I" * 40))
PY
CONS=$(curl -fsS -X POST "$BASE_URL/tasks" \
  -F "r1=@$CONS_R1;filename=R1.fastq" \
  -F "r2=@$CONS_R2;filename=R2.fastq")
CONS_ID=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["id"])' <<<"$CONS")
for _ in 1 2 3 4 5 6 7 8 9 10; do
  S=$(curl -fsS "$BASE_URL/tasks/$CONS_ID" |
      python3 -c 'import json,sys; print(json.load(sys.stdin)["status"])')
  [ "$S" = "succeeded" ] && break
  sleep 0.3
done
curl -fsS -X POST "$BASE_URL/tasks/$CONS_ID/consensus" \
  -H 'Content-Type: application/json' \
  -d '{"reads":[
        {"read_id":"gamma","umi":"CCCCCCCC","start":200},
        {"read_id":"delta","umi":"CCCCCCCA","start":200},
        {"read_id":"epsilon","umi":"GGGGGGGG","start":300}]}' \
  -o "$WORK/consensus.zip"
python3 - "$WORK/consensus.zip" <<'PY'
import sys, zipfile
zf = zipfile.ZipFile(sys.argv[1])
assert sorted(zf.namelist()) == ["R1.fastq", "R2.fastq"], zf.namelist()
for name in zf.namelist():
    print(name, "->", zf.read(name).decode().splitlines()[0])
print("OK: consensus ZIP has R1.fastq and R2.fastq")
PY

echo "== consensus request that omits a read id is rejected as a whole (422) =="
curl -s -X POST "$BASE_URL/tasks/$CONS_ID/consensus" \
  -H 'Content-Type: application/json' \
  -d '{"reads":[{"read_id":"gamma","umi":"CCCCCCCC","start":200}]}' \
  -o /dev/null -w '%{http_code}\n'

echo "ALL DEMO STEPS DONE"
