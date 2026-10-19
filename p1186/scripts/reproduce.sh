#!/usr/bin/env bash
# Reproduce every artifact in p1186 from source.
#
#   bash scripts/reproduce.sh            full run   (many CPU-hours)
#   FAST=1 bash scripts/reproduce.sh     validation + literature only (~10 min)
#   NT=4 bash scripts/reproduce.sh       cap the number of concurrent searches
#
# Stages, in order:
#   0  environment
#   1  fetch third-party references          (never redistributed; see refs/)
#   2  build the independent C brute forces and the quadratic-residue generator
#   3  validate the exact evaluator against closed forms and a second evaluator
#   4  the wildcard chain (k=5, k=6) and the bases (k=4..8) -- runs even under FAST
#   5  reproduce published values, then run the searches
#   6  final_check.py -- re-derives every claim; exit 0 = pass
set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
PY="${PY:-$HOME/.venvs/main/bin/python}"
FAST="${FAST:-0}"
NT="${NT:-4}"
mkdir -p results

say() { printf '\n=== %s ===\n' "$*"; }
die() { printf 'FAILED: %s\n' "$*" >&2; exit 1; }

say "0. environment"
command -v cc >/dev/null || die "no C compiler"
[ -x "$PY" ] || die "no python at $PY (override with PY=...)"
"$PY" -c 'import numpy, sympy; print("python", __import__("sys").version.split()[0],
      "numpy", numpy.__version__, "sympy", sympy.__version__)' || die "missing numpy/sympy"
"$PY" -c 'import pysat' 2>/dev/null \
  && echo "python-sat present (src/apfree.py, src/safesat.py and src/bases.py enabled)" \
  || echo "python-sat ABSENT -- src/apfree.py, src/safesat.py and src/bases.py will be skipped (pip install python-sat)"

say "1. fetch references"
bash refs/fetch.sh || echo "  (some references are paywalled; see refs/README.md)"

say "2. build the independent C brute forces and the quadratic-residue generator"
cc -O3 -o src/verify_zp   src/verify_zp.c   || die "verify_zp"
cc -O3 -o src/verify_n    src/verify_n.c    || die "verify_n"
cc -O3 -o src/verify_two  src/verify_two.c  || die "verify_two"
cc -O3 -o src/verify_zm   src/verify_zm.c   || die "verify_zm"
cc -O3 -o src/verify_base src/verify_base.c || die "verify_base"
cc -O3 -o src/rabung      src/rabung.c      || die "rabung"
cc -O3 -o src/search      src/search.c -lm  || die "search"
echo "built: verify_zp verify_n verify_two verify_zm verify_base rabung search"

say "3. validate the exact evaluator"
"$PY" src/final_check.py --fast --no-brute >/dev/null \
  || die "evaluator validation failed -- do not trust anything below"
echo "evaluator, second evaluator and re-certification all pass"

say "4. the wildcard chain (the k=5 and k=6 improvements)"
# Cheap and central, so it runs even under FAST.  The coset probe carries its
# own positive control: it must rediscover the known k=5 structure
# (b=44, f=4, r=3) before its k=6 answer is believed.  Each chain level
# re-certifies against the recursion and aborts on mismatch.
chain() { echo "  + $*"; nice -n 5 "$PY" -u "$@" ; }
chain src/cosetprobe.py --k 6 --b 86 --t 2 --self-test \
    > results/cosetprobe_k6.log 2>&1 || die "coset probe (control may have failed)"
chain src/chain.py --k 5 --b 44 --f 4 --r 3 \
    > results/chain_k5.log 2>&1 || die "k=5 chain"
chain src/chain.py --k 6 --b 86 --f 2 --r 33 \
    --store results/wildcard_k6_deep.json > results/chain_k6.log 2>&1 || die "k=6 chain"
tail -4 results/chain_k5.log results/chain_k6.log

# The bases behind the recursion limits (NOTES.md sections 3-5): quadratic-residue words from src/rabung plus B_44 and B_86 from results/zm_k*.json, each re-checked before it is written. The self-test confirms that symmetry breaking changes no SAT answer and that its comparison can see a changed one.
if "$PY" -c 'import pysat' 2>/dev/null; then
  "$PY" src/bases.py || die "bases"
  "$PY" src/safesat.py --self-test || die "safesat.py self-test"
fi

if [ "$FAST" = "1" ]; then
  say "FAST=1 -- skipping the searches"
  "$PY" src/final_check.py --fast --no-brute
  exit $?
fi

say "5. searches (this is the long part)"
# Each search writes its own shard.  Shards are merged on read, never in place:
# concurrent writers that share one file silently clobber each other.
run() { echo "  + $*"; nice -n 5 "$PY" -u "$@" ; }

# 5a. the F_p block family -- bounds delta~_k
for k in 4 5 6 7; do
  run src/relax.py "$k" 2-150 --setting zp --restarts 8 --tabu 12000 \
      --store "results/zp_k${k}.json" > "results/relax_zp_k${k}.log" 2>&1 &
  while [ "$(jobs -r | wc -l)" -ge "$NT" ]; do wait -n; done
done
wait

# 5b. the {1..n} block family -- reproduces PRS08 at k=3
for k in 3 4 5; do
  run src/relax.py "$k" 2-120 --setting n --restarts 8 --tabu 12000 \
      --store "results/n_k${k}.json" > "results/relax_n_k${k}.log" 2>&1 &
  while [ "$(jobs -r | wc -l)" -ge "$NT" ]; do wait -n; done
done
wait

# 5c. the periodic family -- Lu-Peng's Lemma 1; this is where delta_5 improves
for k in 3 4 5 6; do
  run src/periodic.py "$k" 2-150 --restarts 40 --tabu 6000 \
      --store "results/zm_k${k}.json" > "results/zm_k${k}.log" 2>&1 &
  while [ "$(jobs -r | wc -l)" -ge "$NT" ]; do wait -n; done
done
wait

# 5d. large moduli, seeded by tiling the best small-modulus word.
#     Random restarts are useless past m ~ 100; tiling preserves Phi exactly,
#     so the tiled record is a starting point that already matches the record.
run src/bigm.py 5 --base 44  --mult 2,3,4,5,6 --tabu 30000 --restarts 6 \
    --store results/zm_k5_big.json > results/zm_k5_big.log 2>&1

# 5e. the two-scale family (periodic x block); contains both parents
run src/twoscale.py 5 --fine 44 --coarse 1-8 --restarts 14 --tabu 10000 \
    --store results/two_k5.json > results/two_k5.log 2>&1 &
run src/twoscale.py 4 --fine 2-12 --coarse 2-12 --restarts 14 --tabu 10000 \
    --store results/two_k4.json > results/two_k4.log 2>&1 &
wait

# 5f. SAT: turns "the search found none" into a proof that none exists
if "$PY" -c 'import pysat' 2>/dev/null; then
  run src/apfree.py 5 2-130 --budget 600 --store results/apfree_k5.json \
      > results/apfree_k5.log 2>&1
  run src/apfree.py 6 87-240 --budget 3000000 --store results/apfree_k6.json \
      > results/apfree_k6.log 2>&1
  # Every (b, f) whose b + f beats the best base, up to the cap and for k <= 5 beyond it (NOTES.md section 5). --store merges on each write, so the k=6 runs must stay sequential. b=221 needs symmetry breaking; nothing else does.
  run src/safesat.py 3 2-60 --above 4 --store results/safesat_k3.json \
      > results/safesat_k3.log 2>&1
  run src/safesat.py 4 2-150 --above 12 --store results/safesat_k4.json \
      > results/safesat_k4.log 2>&1
  run src/safesat.py 5 2-200 --above 48 --store results/safesat_k5.json \
      > results/safesat_k5.log 2>&1
  { run src/safesat.py 6 150-220 --above 228 --store results/safesat_k6.json
    run src/safesat.py 6 222-226 --above 228 --store results/safesat_k6.json
    run src/safesat.py 6 221 --above 228 --symbreak 24 --solver cadical195 \
        --store results/safesat_k6.json
  } > results/safesat_k6.log 2>&1
fi

say "6. final check"
"$PY" src/final_check.py
