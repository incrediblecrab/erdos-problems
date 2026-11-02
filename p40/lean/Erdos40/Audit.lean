/-
Copyright (c) 2026 Max Marquardt.
Released under the MIT licence; see LICENSE at the repository root.
-/
import Erdos40.Problem158
import Erdos40.TheoremB

/-!
# Axiom and statement audit

Two ways a compiled theorem can still fail to prove what it claims: it can rest on `sorry` or on an extra axiom, or its statement can differ from the one intended. This file checks both while the library builds, so `lake build` fails if either check does.

`#assert_axioms` walks the transitive dependency graph with `Lean.collectAxioms` and fails unless every axiom found is one of `propext`, `Classical.choice`, `Quot.sound`.

The `run_cmd` block below compares elaborated statements, as `Expr`s, against the declarations in `FormalConjectures` at the pinned commit. It checks the left side of `exists_iff` against the right side of `Erdos40.erdos_40.variants.weaker`, the set in `answerSet_nonempty_iff` against the set in `Erdos40.erdos_40`, and the conclusion of `strong_implies_erdos_28` against the type of `Erdos28.erdos_28`. For problem 158 it checks the conclusion of `strong_implies_erdos_158` against the right side of `Erdos158.erdos_158`, the right side of `strong_iff_forall_B2` with `g := 2` against that same right side less its hypothesis `A.Infinite`, and both sides of `weaker_iff_forall_B2` against those two anchors. For Theorem B it checks that `rpow_not_mem_answerSet` says, for every `ε > 0`, that `fun N ↦ (N : ℝ) ^ ε` is not in the set of `Erdos40.erdos_40`, and that `not_erdos40For_rpow` says that it lacks `Erdos40.Erdos40For`.

These checks run inside the environment they audit. `../verify.sh` repeats the axiom check from a separate process that loads the compiled modules at runtime, and replays them through the kernel with `leanchecker`.
-/

open Lean Elab Command in
/-- `#assert_axioms foo` logs the axioms `foo` depends on and fails unless they are among `propext`, `Classical.choice`, `Quot.sound`. -/
elab "#assert_axioms " id:ident : command => do
  let n ← liftCoreM <| realizeGlobalConstNoOverload id
  let axs ← liftCoreM <| collectAxioms n
  let bad := axs.filter fun a ↦ !([``propext, ``Classical.choice, ``Quot.sound].contains a)
  unless bad.isEmpty do
    throwError "{n} depends on non-standard axioms {bad.toList}"
  logInfo m!"{n} depends on: {axs.toList}"

namespace ErdosProblem40

#assert_axioms sumRep_le
#assert_axioms exists_bound_of_limsup_ne_top
#assert_axioms limsup_ne_top_of_bound
#assert_axioms sumRep_congr
#assert_axioms count_congr
#assert_axioms isClopen_of_local
#assert_axioms uniform_dip
#assert_axioms exists_of_strong
#assert_axioms exists_iff
#assert_axioms answerSet_nonempty_iff
#assert_axioms strong_implies_erdos_28
#assert_axioms Erdos40For.of_isBigO
#assert_axioms sumRep_univ
#assert_axioms univ_hypothesis_and_conclusion
#assert_axioms sumRep_singleton_one_le
#assert_axioms not_erdos40For_sqrt
#assert_axioms not_erdos40For_of_sqrt_isBigO
#assert_axioms sumRep_le_of_B2
#assert_axioms B2_of_sumRep_le
#assert_axioms count_sq_le
#assert_axioms rpow_neg_half
#assert_axioms tendsto_inv_sqrt
#assert_axioms liminf_eq_zero_of_bound
#assert_axioms strong_iff_forall_B2
#assert_axioms weaker_iff_forall_B2
#assert_axioms strong_implies_erdos_158
#assert_axioms TheoremB.measure_subset_le
#assert_axioms TheoremB.sum_powersetCard_prod_le
#assert_axioms TheoremB.measure_E_le
#assert_axioms TheoremB.sumRep_le_of_not_mem
#assert_axioms TheoremB.count_ge
#assert_axioms TheoremB.exists_good
#assert_axioms exists_dense_bounded
#assert_axioms not_erdos40For_rpow
#assert_axioms rpow_not_mem_answerSet
#assert_axioms not_erdos40For_of_rpow_isBigO

open Lean Elab Command Term Meta in
run_cmd do
  let env ← getEnv
  let get (n : Name) : CommandElabM ConstantInfo := do
    let some c := env.find? n | throwError "{n} not found"
    pure c
  -- exists_iff : (the right side of FC's weaker variant) ↔ Erdos40For (fun _ ↦ 1)
  let weaker ← get ``Erdos40.erdos_40.variants.weaker
  let ours ← get ``ErdosProblem40.exists_iff
  let some (_, fcRhs) := weaker.type.iff? | throwError "erdos_40.variants.weaker is not an ↔"
  let some (lhs, rhs) := ours.type.iff? | throwError "exists_iff is not an ↔"
  unless lhs == fcRhs do
    throwError "exists_iff: left side differs from FormalConjectures.\n  ours: {lhs}\n  theirs: {fcRhs}"
  let strong ← liftTermElabM do
    instantiateMVars (← elabTerm (← `(Erdos40.Erdos40For (fun _ : ℕ ↦ (1 : ℝ)))) none)
  unless rhs == strong do
    throwError "exists_iff: right side is not Erdos40For (fun _ ↦ 1).\n  ours: {rhs}"
  -- answerSet_nonempty_iff : (the set in FC's erdos_40).Nonempty ↔ Erdos40For (fun _ ↦ 1)
  let main ← get ``Erdos40.erdos_40
  let some (_, fcSet, _) := main.type.eq? | throwError "erdos_40 is not an equation"
  let ours2 ← get ``ErdosProblem40.answerSet_nonempty_iff
  let some (lhs2, rhs2) := ours2.type.iff? | throwError "answerSet_nonempty_iff is not an ↔"
  let expectedLhs2 ← liftTermElabM do
    instantiateMVars (← mkAppM ``Set.Nonempty #[fcSet])
  unless lhs2 == expectedLhs2 do
    throwError "answerSet_nonempty_iff: set differs from FormalConjectures.\n  ours: {lhs2}"
  unless rhs2 == strong do
    throwError "answerSet_nonempty_iff: right side is not Erdos40For (fun _ ↦ 1)."
  -- strong_implies_erdos_28 : Erdos40For (fun _ ↦ 1) → (the type of FC's erdos_28)
  let e28 ← get ``Erdos28.erdos_28
  let ours3 ← get ``ErdosProblem40.strong_implies_erdos_28
  let .forallE _ hyp body _ := ours3.type | throwError "strong_implies_erdos_28 has no hypothesis"
  unless hyp == strong do
    throwError "strong_implies_erdos_28: hypothesis is not Erdos40For (fun _ ↦ 1)."
  unless body == e28.type do
    throwError "strong_implies_erdos_28: conclusion differs from Erdos28.erdos_28."
  -- strong_implies_erdos_158 : Erdos40For (fun _ ↦ 1) → (the right side of FC's erdos_158)
  let e158 ← get ``Erdos158.erdos_158
  let some (_, fc158) := e158.type.iff? | throwError "erdos_158 is not an ↔"
  let ours4 ← get ``ErdosProblem40.strong_implies_erdos_158
  let .forallE _ hyp4 body4 _ := ours4.type | throwError "strong_implies_erdos_158 has no hypothesis"
  unless hyp4 == strong do
    throwError "strong_implies_erdos_158: hypothesis is not Erdos40For (fun _ ↦ 1)."
  unless !body4.hasLooseBVars && body4 == fc158 do
    throwError "strong_implies_erdos_158: conclusion differs from Erdos158.erdos_158.\n  ours: {body4}\n  theirs: {fc158}"
  -- strong_iff_forall_B2 : Erdos40For (fun _ ↦ 1) ↔ ∀ g A, B2 g A → L A. With g := 2 the right side must be FC's `∀ A, A.Infinite → B2 2 A → L A` without the hypothesis `A.Infinite`.
  let ours5 ← get ``ErdosProblem40.strong_iff_forall_B2
  let some (lhs5, rhs5) := ours5.type.iff? | throwError "strong_iff_forall_B2 is not an ↔"
  unless lhs5 == strong do
    throwError "strong_iff_forall_B2: left side is not Erdos40For (fun _ ↦ 1)."
  let .forallE _ (.const ``Nat []) gBody _ := rhs5 | throwError "strong_iff_forall_B2: right side does not start with ∀ g : ℕ"
  let .forallE _ tA oursA _ := gBody.instantiate1 (mkNatLit 2) | throwError "strong_iff_forall_B2: right side has no ∀ A"
  let .forallE _ tA' (.forallE _ _ fcA _) _ := fc158 | throwError "erdos_158: right side has another shape"
  unless tA == tA' && !fcA.hasLooseBVar 0 && oursA == fcA.lowerLooseBVars 1 1 do
    throwError "strong_iff_forall_B2: right side at g = 2 differs from Erdos158.erdos_158 without A.Infinite.\n  ours: {oursA}\n  theirs: {fcA}"
  -- weaker_iff_forall_B2 : (the right side of FC's weaker variant) ↔ (the right side of strong_iff_forall_B2)
  let ours6 ← get ``ErdosProblem40.weaker_iff_forall_B2
  let some (lhs6, rhs6) := ours6.type.iff? | throwError "weaker_iff_forall_B2 is not an ↔"
  unless lhs6 == fcRhs && rhs6 == rhs5 do
    throwError "weaker_iff_forall_B2: a side differs from erdos_40.variants.weaker or strong_iff_forall_B2."
  -- rpow_not_mem_answerSet : ∀ ε : ℝ, 0 < ε → (fun N ↦ (N : ℝ) ^ ε) ∉ (the set in FC's erdos_40)
  let ours7 ← get ``ErdosProblem40.rpow_not_mem_answerSet
  let expected7 ← liftTermElabM do
    let lam ← elabTerm (← `(fun S : Set (ℕ → ℝ) ↦ ∀ ε : ℝ, 0 < ε → (fun N : ℕ ↦ (N : ℝ) ^ ε) ∉ S)) none
    instantiateMVars (lam.beta #[fcSet])
  unless ours7.type == expected7 do
    throwError "rpow_not_mem_answerSet: statement differs from `∀ ε > 0, N ^ ε ∉` the set of Erdos40.erdos_40.\n  ours: {ours7.type}\n  expected: {expected7}"
  -- not_erdos40For_rpow : ∀ {ε : ℝ}, 0 < ε → ¬ Erdos40For (fun N ↦ (N : ℝ) ^ ε)
  let ours8 ← get ``ErdosProblem40.not_erdos40For_rpow
  let expected8 ← liftTermElabM do
    instantiateMVars (← elabTerm (← `(∀ {ε : ℝ}, 0 < ε → ¬ Erdos40.Erdos40For (fun N : ℕ ↦ (N : ℝ) ^ ε))) none)
  unless ours8.type == expected8 do
    throwError "not_erdos40For_rpow: statement differs from `∀ ε > 0, ¬ Erdos40For (N ^ ε)`.\n  ours: {ours8.type}"
  logInfo "statement check passed: exists_iff, answerSet_nonempty_iff, strong_implies_erdos_28, strong_implies_erdos_158, strong_iff_forall_B2, weaker_iff_forall_B2, rpow_not_mem_answerSet, not_erdos40For_rpow match FormalConjectures"

end ErdosProblem40
