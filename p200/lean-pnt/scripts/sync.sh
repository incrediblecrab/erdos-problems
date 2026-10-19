#!/usr/bin/env bash
# Regenerate the shared Erdos200 modules for the PrimeNumberTheoremAnd build.
#
# `../lean/Erdos200/` is the canonical source. This project differs only because it
# builds against PNT+'s Mathlib (Lean 4.32.2) instead of formal-conjectures' Mathlib
# (Lean 4.33.1). Exactly two mechanical substitutions are needed:
#
#   1. `FormalConjectures.ErdosProblems.«200»` -> `Erdos200.FCShim`
#      formal-conjectures cannot be a dependency here (conflicting Mathlib revision),
#      so its definitions are copied verbatim into FCShim.lean and checked by
#      scripts/check_shim.py.
#
#   2. `Set.mem_ofPred_eq` -> `Set.mem_setOf_eq`
#      Mathlib renamed this on 2026-07-09; the new name does not exist on 4.32.2.
#      Both are `rfl` lemmas with identical statements.
#
# Generated files are gitignored. Edit ../lean/Erdos200/ instead.
set -euo pipefail
cd "$(dirname "$0")/.."

SRC=../lean/Erdos200
MODULES=(Structure Bound Asymptotic)

for m in "${MODULES[@]}"; do
  {
    echo "/- GENERATED FILE -- DO NOT EDIT."
    echo "   Produced by scripts/sync.sh from $SRC/$m.lean; edit that file instead. -/"
    sed -e 's|^import FormalConjectures.ErdosProblems.«200»|import Erdos200.FCShim|' \
        -e 's|Set\.mem_ofPred_eq|Set.mem_setOf_eq|g' \
        "$SRC/$m.lean"
  } > "Erdos200/$m.lean"
  echo "  generated Erdos200/$m.lean from $SRC/$m.lean"
done

# Guard: no formal-conjectures *import* may survive into the generated tree.
if grep -rn "^import FormalConjectures" Erdos200/ ; then
  echo "FAIL: a FormalConjectures import survived sync" >&2
  exit 1
fi
echo "sync ok (${#MODULES[@]} modules)"
