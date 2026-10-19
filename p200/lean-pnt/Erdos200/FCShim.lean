/-
Copyright 2024 Google LLC. Licensed under the Apache License, Version 2.0.

VERBATIM EXCERPT -- NOT ORIGINAL WORK.

The declarations between the BEGIN/END markers below are copied character-for-character
from `google-deepmind/formal-conjectures` so that the Erdos 200 development can be built
against `PrimeNumberTheoremAnd`'s Mathlib (Lean 4.32.2), which cannot coexist in one Lake
project with `formal-conjectures`' Mathlib (Lean 4.33.1).

Sources:
  FormalConjecturesForMathlib/Combinatorics/AP/Basic.lean  (lines 44-45, 58-59, 122-123)
  FormalConjectures/ErdosProblems/200.lean                 (lines 35-36, 50-51)

`scripts/check_shim.py` re-downloads both files from `main` and fails if any copied line
below differs from the upstream text.

This file is Apache-2.0, not MIT: the excerpts are Google's work and are not relicensed
by this repository. See the "Licence" section of the top-level README.
-/
import Erdos200.Bound
import Mathlib.Data.ENat.Basic

open Filter Real

variable {α : Type*} [AddCommMonoid α]

-- BEGIN VERBATIM: FormalConjecturesForMathlib/Combinatorics/AP/Basic.lean
def Set.IsAPOfLengthWith (s : Set α) (l : ℕ∞) (a d : α) : Prop :=
  ENat.card s = l ∧ s = {a + n • d | (n : ℕ) (_ : n < l)}
-- END VERBATIM

-- BEGIN VERBATIM: FormalConjecturesForMathlib/Combinatorics/AP/Basic.lean
def Set.IsAPOfLength (s : Set α) (l : ℕ∞) : Prop :=
  ∃ a d : α, s.IsAPOfLengthWith l a d
-- END VERBATIM

namespace Set.IsAPOfLengthWith

variable {s : Set α} {l : ℕ∞} {a d : α}

-- BEGIN VERBATIM: FormalConjecturesForMathlib/Combinatorics/AP/Basic.lean
theorem card (h : s.IsAPOfLengthWith l a d) : ENat.card s = l := h.1
theorem eq (h : s.IsAPOfLengthWith l a d) : s = {a + n • d | (n : ℕ) (_ : n < l)} := h.2
-- END VERBATIM

end Set.IsAPOfLengthWith

namespace Erdos200

-- BEGIN VERBATIM: FormalConjectures/ErdosProblems/200.lean
noncomputable def longestPrimeArithmeticProgressions (n : ℕ) : ℕ :=
  sSup {(k : ℕ) | ∃ s ⊆ Set.Icc 1 n, s.IsAPOfLength k ∧ ∀ m ∈ s, m.Prime}
-- END VERBATIM

end Erdos200
