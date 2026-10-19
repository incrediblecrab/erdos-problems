/-
Copyright (c) 2026 Max Marquardt.
Released under the MIT licence; see LICENSE at the repository root.

Axiom audit for the `PrimeNumberTheoremAnd` build.

A Lean theorem can fail to mean what it appears to mean in two ways: its proof may rest
on `sorry`, or its statement may not be the intended one.  This file addresses the first.
`#assert_axioms` collects the axiom dependencies of a declaration with
`Lean.collectAxioms` and **fails the build** if `sorryAx` appears, so `lake build`
succeeding certifies that the listed theorems rest only on `propext`,
`Classical.choice` and `Quot.sound`.

This matters more here than in the sibling `../lean` project, because
`PrimeNumberTheoremAnd` does contain `sorry`: `Wiener.lean:323` (`prelim_decay_2`) and
`Wiener.lean:342` (`prelim_decay_3`).  Both are dead code -- nothing in the project
references them except `decay_alt`, which is itself referenced nowhere.  The assertion
on `chebyshev_asymptotic` below is what actually certifies that, rather than the
textual argument.

For the second failure mode, `scripts/check_shim.py` verifies that the definitions in
`FCShim.lean` are copied character-for-character from `formal-conjectures`.
-/
import Erdos200.Unconditional

open Lean Elab Command in
/-- `#assert_axioms foo` reports the axioms `foo` depends on, and errors if `sorryAx`
is among them. -/
elab "#assert_axioms " id:ident : command => do
  let n ← liftCoreM <| realizeGlobalConstNoOverload id
  let axs ← liftCoreM <| collectAxioms n
  if axs.contains ``sorryAx then
    throwError "{n} depends on sorryAx"
  logInfo m!"{n} depends on: {axs.toList}"

-- The imported Prime Number Theorem, which is what makes the result unconditional.
#assert_axioms chebyshev_asymptotic

namespace ErdosProblem200

-- The bridge from PNT+'s `θ ~[atTop] id` to the `PNT` hypothesis.
#assert_axioms pnt_of_primeNumberTheoremAnd

-- The two headline bounds.
#assert_axioms upper_unconditional
#assert_axioms erdos_200_variants_upper_unconditional

end ErdosProblem200
