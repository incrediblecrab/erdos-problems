/-
Copyright (c) 2026 Max Marquardt.
Released under the MIT licence; see LICENSE at the repository root.
-/
import Erdos200.Asymptotic

/-!
# Axiom audit

A declaration can compile and still be vacuous, either because it (or something it depends on)
used `sorry`, or because the statement is not the one intended.  This file checks the first
failure mode mechanically: `#assert_axioms` walks the full transitive dependency graph with
`Lean.collectAxioms` and **fails the build** if `sorryAx` appears.

`lake build` succeeding therefore certifies that every theorem listed below rests only on
`propext`, `Classical.choice` and `Quot.sound`, the three axioms of Lean's standard logic.

It does *not* certify the second failure mode.  For that, compare the statements against
`../../NOTES.md` by hand; `erdos_200_variants_upper` is stated to match
`Erdos200.erdos_200.variants.upper` in `FormalConjectures` verbatim except for the `PNT`
hypothesis.
-/

open Lean Elab Command in
/-- `#assert_axioms foo` reports the axioms `foo` depends on, and errors if `sorryAx`
is among them. -/
elab "#assert_axioms " id:ident : command => do
  let n ← liftCoreM <| realizeGlobalConstNoOverload id
  let axs ← liftCoreM <| collectAxioms n
  if axs.contains ``sorryAx then
    throwError "{n} depends on sorryAx"
  logInfo m!"{n} depends on: {axs.toList}"

namespace ErdosProblem200

-- Proposition 1 (`Structure.lean`)
#assert_axioms exists_lt_dvd
#assert_axioms le_of_not_dvd
#assert_axioms primorial_dvd
#assert_axioms structure_of_primeAP
#assert_axioms primorial_dvd_mul_iff

-- Proposition 2 (`Bound.lean`)
#assert_axioms mul_primorial_lt
#assert_axioms primorial_lt_two_mul
#assert_axioms theta_lt
#assert_axioms theta_lt_sharp

-- Asymptotics (`Asymptotic.lean`)
#assert_axioms exists_of_isAPOfLength
#assert_axioms longest_eq_zero
#assert_axioms longest_le_max
#assert_axioms log_le_two_mul_sqrt
#assert_axioms eventually_log_add_le
#assert_axioms eventually_longest_le
#assert_axioms exists_littleO_of_eventually_le
#assert_axioms theta_lower_log_two
#assert_axioms upper_unconditional
#assert_axioms theta_lower_of_pnt
#assert_axioms erdos_200_variants_upper

/-
### Statement fidelity

Compiling is not the same as proving the right thing.  This compares the elaborated
`Expr` of our conclusion, after stripping the `PNT` binder, against the elaborated `Expr` of
`Erdos200.erdos_200.variants.upper` in `FormalConjectures` (which is `sorry` there).  Any
difference in quantifier order, coercion, or the meaning of `=o` is a build failure.
-/
open Lean Elab Command in
run_cmd do
  let env ← getEnv
  let some fc := env.find? ``Erdos200.erdos_200.variants.upper
    | throwError "FormalConjectures' erdos_200.variants.upper not found"
  let some ours := env.find? ``ErdosProblem200.erdos_200_variants_upper
    | throwError "erdos_200_variants_upper not found"
  let .forallE _ hyp body _ := ours.type
    | throwError "expected erdos_200_variants_upper to start with a hypothesis binder"
  unless hyp.isConstOf ``ErdosProblem200.PNT do
    throwError "expected the hypothesis to be PNT, got {hyp}"
  unless body == fc.type do
    throwError "statement mismatch:\n  ours: {body}\n  theirs: {fc.type}"
  logInfo "erdos_200_variants_upper : PNT → (FormalConjectures' erdos_200.variants.upper)"

end ErdosProblem200
