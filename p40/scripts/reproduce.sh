#!/usr/bin/env bash
#
# Reproduce every result in NOTES.md for Erdős problem #40, then run the gate.
#
# Usage: bash scripts/reproduce.sh [STAGE ...]. With no arguments every stage runs, in this order:
#   fetch   refs/fetch.sh downloads the sources; nothing is redistributed. Not fatal: offline, src/final_check.py reports the checks that need refs/ as skipped, and a skip is not a pass.
#   thmb    src/theorem_b_check.py, then its two planted defects, each of which must exit 1
#   greedy  src/greedy_check.py, then its three planted defects, each of which must exit 1
#   lean    lean/verify.sh under /usr/bin/time, logged to results/lean_verify.log with the date, the Lean release, the platform and the hashes of the files it checked
#   plants  lean/plants.py: verify.sh on an unmodified copy and on eleven planted defects, which takes about a quarter of an hour
#   final   src/final_check.py --write regenerates the tables in the documentation from results/; src/final_check.py then checks everything, logging to results/final_check.log; last, a planted wrong number must make it exit 1. It runs with --strict when the fetch stage succeeded in the same invocation, or when STRICT=--strict is set because refs/ was populated earlier.
#
# Needs bash, curl, a C compiler, poppler's pdftotext, elan (for lake and lean) and a python3 with numpy and mpmath. $PY overrides the interpreter; otherwise ~/.venvs/main/bin/python is used when present, and python3 from $PATH otherwise.
#
# Local paths are replaced by ~ in every log this writes.
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="$HOME/.elan/bin:$PATH"

PY="${PY:-}"
if [ -z "$PY" ]; then
  if [ -x "$HOME/.venvs/main/bin/python" ]; then
    PY="$HOME/.venvs/main/bin/python"
  else
    PY="$(command -v python3 || true)"
  fi
fi
if [ -z "$PY" ] || ! "$PY" -c 'import numpy, mpmath' >/dev/null 2>&1; then
  echo "need a python3 with numpy and mpmath; set PY=/path/to/python" >&2
  exit 1
fi

STAGES="${*:-fetch thmb greedy lean plants final}"
STRICT="${STRICT:-}"
mkdir -p results

# A planted defect must make its checker exit 1, not merely fail somehow.
expect_exit1() {
  local label=$1 rc
  shift
  set +e
  "$@" >/dev/null 2>&1
  rc=$?
  set -e
  if [ "$rc" -ne 1 ]; then
    echo "  plant $label exited $rc, expected 1" >&2
    exit 1
  fi
  echo "  plant $label: caught, exit 1"
}

for stage in $STAGES; do
  case "$stage" in
  fetch)
    echo "=== fetch: refs/fetch.sh ==="
    if sh refs/fetch.sh >results/fetch.log 2>&1; then
      STRICT="--strict"
      echo "  refs/ populated; results/fetch.log lists the hashes"
    else
      echo "  fetch failed (offline?); the checks that need refs/ will report skip"
    fi
    ;;
  thmb)
    echo "=== thmb: src/theorem_b_check.py ==="
    "$PY" src/theorem_b_check.py | tail -n 2
    expect_exit1 denominators "$PY" src/theorem_b_check.py --plant denominators
    expect_exit1 k "$PY" src/theorem_b_check.py --plant k
    ;;
  greedy)
    echo "=== greedy: src/greedy_check.py ==="
    "$PY" src/greedy_check.py | tail -n 1
    for p in nongreedy notb2 note; do
      expect_exit1 "$p" "$PY" src/greedy_check.py --plant "$p"
    done
    ;;
  lean)
    echo "=== lean: lean/verify.sh ==="
    if /usr/bin/time -l true >/dev/null 2>&1; then
      TIMER=(/usr/bin/time -l)
    elif /usr/bin/time -v true >/dev/null 2>&1; then
      TIMER=(/usr/bin/time -v)
    else
      TIMER=()
    fi
    tmp="$(mktemp)"
    {
      echo "# lean/verify.sh, run by scripts/reproduce.sh"
      echo "# date: $(date '+%Y-%m-%dT%H:%M:%S%z')"
      echo "# $(cd lean && lean --version)"
      echo "# platform: $(uname -sm)"
    } >"$tmp"
    set +e
    ${TIMER[@]+"${TIMER[@]}"} bash lean/verify.sh >>"$tmp" 2>&1
    rc=$?
    set -e
    echo "# verify.sh exit $rc" >>"$tmp"
    echo "# sha256 of the files it checked:" >>"$tmp"
    (cd lean && for f in *.lean Erdos40/*.lean verify.sh lakefile.toml lake-manifest.json lean-toolchain; do
      "$PY" -c 'import hashlib, sys; print("#", hashlib.sha256(open(sys.argv[1], "rb").read()).hexdigest(), sys.argv[1])' "$f"
    done) >>"$tmp"
    sed "s|$HOME|~|g" "$tmp" >results/lean_verify.log
    rm -f "$tmp"
    grep '^PASS:' results/lean_verify.log || { echo "  lean/verify.sh exited $rc; see results/lean_verify.log" >&2; exit 1; }
    ;;
  plants)
    echo "=== plants: lean/plants.py ==="
    "$PY" lean/plants.py
    ;;
  final)
    echo "=== final: src/final_check.py ==="
    "$PY" src/final_check.py --write >/dev/null || true
    "$PY" src/final_check.py $STRICT 2>&1 | sed "s|$HOME|~|g" | tee results/final_check.log
    expect_exit1 number "$PY" src/final_check.py --plant number
    ;;
  *)
    echo "unknown stage $stage; the stages are fetch thmb greedy lean plants final" >&2
    exit 2
    ;;
  esac
done
