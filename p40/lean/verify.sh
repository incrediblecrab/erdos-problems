#!/usr/bin/env bash
# Checks the Lean development for Erdős problem 40. Exit 0 means all five steps passed; README.md lists what each step checks and the planted defects that made it fail. `./plants.py` reruns those plants.
set -uo pipefail
cd "$(dirname "$0")"

export PATH="$HOME/.elan/bin:$PATH"
command -v lake >/dev/null || { echo "FAIL: lake not on PATH (install elan)"; exit 1; }

EXPECTED_AXIOM_LINES=26
MODULES="Erdos40.Basic Erdos40.Equivalence Erdos40.Problem158"
fail() { echo "FAIL: $*"; exit 1; }

echo "== 0. formal-conjectures checkout =="
# The statement check in step 1 compares against formal-conjectures as built here, so that checkout must be the pinned commit with no local edits.
pin=$(sed -n 's/^rev = "\([0-9a-f]\{40\}\)"$/\1/p' lakefile.toml)
[ -n "$pin" ] || fail "no 40-hex rev in lakefile.toml"
fc=.lake/packages/formal_conjectures
if [ -d "$fc" ]; then
  head=$(git -C "$fc" rev-parse HEAD) || fail "cannot read the formal_conjectures commit"
  [ "$head" = "$pin" ] || fail "formal_conjectures is at $head, lakefile.toml pins $pin"
  dirty=$(git -C "$fc" status --porcelain --untracked-files=all) || fail "cannot read formal_conjectures status"
  [ -z "$dirty" ] || { echo "$dirty"; fail "formal_conjectures has local changes"; }
  echo "formal_conjectures at $pin, no local changes"
else
  echo "formal_conjectures not fetched yet; lake build fetches the pinned commit, and step 0 runs again after it"
fi

echo "== 1. lake build (includes Erdos40/Audit.lean) =="
lake build || fail "lake build exited non-zero"
if [ -z "${head:-}" ]; then
  head=$(git -C "$fc" rev-parse HEAD) || fail "cannot read the formal_conjectures commit"
  [ "$head" = "$pin" ] || fail "formal_conjectures is at $head, lakefile.toml pins $pin"
  [ -z "$(git -C "$fc" status --porcelain --untracked-files=all)" ] || fail "formal_conjectures has local changes"
fi

echo "== 2. re-elaborating Erdos40/Audit.lean =="
out=$(lake env lean Erdos40/Audit.lean 2>&1) || { echo "$out"; fail "Audit.lean did not elaborate"; }
n=$(grep -c 'depends on: \[propext, Classical.choice, Quot.sound\]' <<<"$out")
[ "$n" -eq "$EXPECTED_AXIOM_LINES" ] || { echo "$out"; fail "expected $EXPECTED_AXIOM_LINES clean axiom reports, got $n"; }
grep -q '^.*statement check passed: exists_iff, answerSet_nonempty_iff, strong_implies_erdos_28, strong_implies_erdos_158, strong_iff_forall_B2, weaker_iff_forall_B2 match FormalConjectures$' <<<"$out" || { echo "$out"; fail "statement check did not report success"; }
echo "$n axiom reports clean; statement check passed"

echo "== 3. axiom probe, compiled without the library and loading it at runtime =="
probe=$(lake env lean --run Probe.lean 2>&1); rc=$?
echo "$probe"
[ "$rc" -eq 0 ] || fail "probe exited $rc"
tail -n 1 <<<"$probe" | grep -q -E '^probe: [0-9]+ constants checked, 0 problems$' || fail "probe did not report success"

echo "== 4. kernel replay with leanchecker =="
checker="$(lean --print-prefix)/bin/leanchecker"
[ -x "$checker" ] || fail "no leanchecker at $checker"
for m in $MODULES; do
  # No module is nested under these names, so each call replays exactly one module.
  lake env "$checker" "$m" || fail "leanchecker rejected $m"
  echo "replayed $m"
done

echo
echo "PASS: $n theorems on {propext, Classical.choice, Quot.sound}; statements match FormalConjectures at $pin; probe clean; $MODULES replayed by the kernel."
