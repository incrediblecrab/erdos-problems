#!/usr/bin/env bash
# Verify the Lean development for Erdős problem 200.
#
#   1. the library builds;
#   2. every listed theorem depends on exactly {propext, Classical.choice, Quot.sound}
#      -- no sorryAx anywhere in the transitive dependency graph;
#   3. the conclusion of `erdos_200_variants_upper` is the elaborated statement of
#      `Erdos200.erdos_200.variants.upper` from FormalConjectures, not merely similar.
#
# Exit 0 = all three hold.  Both checks in Erdos200/Audit.lean were tested against planted
# defects (a `sorry` in a dependency, and an equivalent-but-reassociated conclusion); each
# made this script exit 1.
set -uo pipefail
cd "$(dirname "$0")"

export PATH="$HOME/.elan/bin:$PATH"
command -v lake >/dev/null || { echo "FAIL: lake not on PATH (install elan)"; exit 1; }

EXPECTED_AXIOM_LINES=20
fail() { echo "FAIL: $*"; exit 1; }

echo "== lake build =="
lake build || fail "lake build exited non-zero"

echo "== re-elaborating Erdos200/Audit.lean =="
out=$(lake env lean Erdos200/Audit.lean 2>&1) || { echo "$out"; fail "Audit.lean did not elaborate"; }

grep -q 'sorryAx' <<<"$out" && { echo "$out"; fail "sorryAx reported"; }

n=$(grep -c 'depends on: \[propext, Classical.choice, Quot.sound\]' <<<"$out")
[ "$n" -eq "$EXPECTED_AXIOM_LINES" ] || {
  echo "$out"
  fail "expected $EXPECTED_AXIOM_LINES clean axiom reports, got $n"
}

grep -q "FormalConjectures' erdos_200.variants.upper" <<<"$out" \
  || { echo "$out"; fail "statement-fidelity check did not report success"; }

echo
echo "PASS: $n theorems machine-checked on {propext, Classical.choice, Quot.sound};"
echo "      erdos_200_variants_upper : PNT -> FormalConjectures' erdos_200.variants.upper."
