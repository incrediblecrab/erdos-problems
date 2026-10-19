#!/usr/bin/env bash
# Verify the PrimeNumberTheoremAnd build of the Erdős 200 upper bound.
#
#   1. the shared modules regenerate from ../lean (the canonical source);
#   2. the definitions copied from formal-conjectures into FCShim.lean are still
#      character-for-character identical to that project's current `main`;
#   3. the library builds;
#   4. every listed theorem -- including the imported Prime Number Theorem -- depends on
#      exactly {propext, Classical.choice, Quot.sound}, with no sorryAx.
#
# Point 4 is the one that matters: PrimeNumberTheoremAnd itself contains two `sorry`s, in
# Wiener.lean at lines 323 and 342. They are dead code, and this check is what proves the
# result does not reach them.
#
# Exit 0 = all four hold.
set -uo pipefail
cd "$(dirname "$0")"

export PATH="$HOME/.elan/bin:$PATH"
command -v lake >/dev/null || { echo "FAIL: lake not on PATH (install elan)"; exit 1; }
PYTHON=${PYTHON:-python3}

EXPECTED_AXIOM_LINES=4
fail() { echo "FAIL: $*"; exit 1; }

echo "== regenerating shared modules from ../lean =="
./scripts/sync.sh || fail "sync.sh exited non-zero"

echo "== checking FCShim.lean against formal-conjectures main =="
"$PYTHON" scripts/check_shim.py || fail "FCShim.lean is not a verbatim copy of upstream"

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

echo
echo "PASS: $n theorems machine-checked on {propext, Classical.choice, Quot.sound},"
echo "      including chebyshev_asymptotic, so"
echo "      erdos_200_variants_upper_unconditional carries no hypothesis."
