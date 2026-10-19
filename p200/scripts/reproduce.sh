#!/usr/bin/env bash
#
# Reproduce every computational claim in NOTES.md for Erdős problem #200.
#
# Stages:
#   1. fetch primary sources          (refs/fetch.sh; nothing is redistributed)
#   2. build the enumerator           (clang -O3)
#   3. verify the 26 known records    four checks, incl. an OEIS-independent source
#   4. assumption-free cross-check    ref_bruteforce.py, which shares no code with apsearch
#   5. exhaustive re-derivation       k = 8..20, the last to B = a(20) = 5.73e11
#   6. heuristic vs measurement       singular series, 30 count comparisons, KS calibration
#   7. sieve limits and tables        barrier.py, sieve_limits.py (~70 s, ~3.1 GB peak memory), analysis.py
#   8. final_check.py                 re-derives every claim from the artifacts; exit 0 = pass
#   9. lean/verify.sh                 axiom audit + statement fidelity for the formal part
#  10. lean-pnt/verify.sh             the same bound with the PNT hypothesis discharged
#
# Needs: clang, curl, python at ~/.venvs/main/bin/python (numpy, sympy).
# Stages 9 and 10 additionally need elan/lake; each is skipped with a notice if lake is
# absent. Their first runs download ~7.7 GB of Mathlib + formal-conjectures oleans and
# build ~7.4 GB of PrimeNumberTheoremAnd respectively.
# Runtime is dominated by stage 5.  On 10 otherwise-idle threads the k=20 run alone
# takes ~67 min and the whole stage several hours; k=16, despite a smaller bound, is
# also expensive because W_16 = 30030 leaves a large m-range.  Set FAST=1 to cap the
# exhaustive runs at k <= 15.  NT sets the thread count.  These figures degrade badly
# under competing load: a FAST=1 NT=3 run on a busy machine took 3.5 h.
#
# Usage:  bash scripts/reproduce.sh
#
# FAST=1 stops stage 5 at k=15.  Note that final_check.py scans every run in
# results/, so with FAST=1 it still reports the higher k from the committed runs
# rather than from ones this script just produced.

set -euo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
cd "$HERE"
PY="${PY:-$HOME/.venvs/main/bin/python}"
NT="${NT:-10}"
FAST="${FAST:-0}"

echo "== 1. primary sources"
sh refs/fetch.sh

echo "== 2. build"
cc -O3 -march=native -funroll-loops -pthread -o src/apsearch src/apsearch.c

echo "== 3. records"
"$PY" src/verify_records.py

echo "== 4. assumption-free cross-check"
"$PY" src/ref_bruteforce.py

echo "== 5. exhaustive re-derivation of a(k)"
mkdir -p results
run() { echo "   k=$1 B=$2"; ./src/apsearch "$1" "$2" "$NT" > "results/count_k$1.json"; }
run 8 100000000
run 9 100000000
run 10 200000000
run 11 1000000000
run 12 1000000000
run 13 3000000000
run 14 5000000000
run 15 8000000000
if [ "$FAST" = "0" ]; then
  run 16 8000000000
  run 17 20000000000
  run 18 30000000000
  run 19 150000000000
  # exactly one 20-AP at or below a(20) => no 21-AP below it either
  echo "   k=20 B=572945039351 (~67 min on 10 idle threads)"
  ./src/apsearch 20 572945039351 "$NT" > results/k20.json
fi

echo "== 6. heuristic vs measurement"
"$PY" src/heuristic.py
"$PY" src/compare_counts.py
"$PY" src/calibration.py

echo "== 7. sieve limits and tables"
"$PY" src/barrier.py
"$PY" src/sieve_limits.py
"$PY" src/analysis.py

echo "== 8. final check"
"$PY" src/final_check.py

echo "== 9. Lean formalisation"
if command -v lake >/dev/null 2>&1 || [ -x "$HOME/.elan/bin/lake" ]; then
  bash lean/verify.sh
else
  echo "   skipped: lake not found (see lean/README.md)"
fi

echo "== 10. Lean: unconditional bound via PrimeNumberTheoremAnd"
if ! command -v lake >/dev/null 2>&1 && [ ! -x "$HOME/.elan/bin/lake" ]; then
  echo "   skipped: lake not found (see lean-pnt/README.md)"
elif [ "$FAST" = "1" ]; then
  echo "   skipped in FAST mode: builds PrimeNumberTheoremAnd from source (~4 min)"
else
  bash lean-pnt/verify.sh
fi
echo
echo "reproduce.sh: all stages completed"
