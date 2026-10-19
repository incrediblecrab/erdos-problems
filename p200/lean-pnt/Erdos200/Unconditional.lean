/-
Copyright (c) 2026 Max Marquardt.
Released under the MIT licence; see LICENSE at the repository root.

Discharging the `PNT` hypothesis using `PrimeNumberTheoremAnd`.

`Erdos200.Asymptotic` proves the Erdos 200 upper bound with constant 1 under the
hypothesis `ErdosProblem200.PNT`, because Mathlib has no Prime Number Theorem.
`PrimeNumberTheoremAnd.Consequences.chebyshev_asymptotic` supplies exactly that fact,
and `#print axioms` certifies it is free of `sorryAx`.  Combining the two yields the
bound unconditionally.
-/
import Erdos200.Asymptotic
import PrimeNumberTheoremAnd.Consequences

namespace ErdosProblem200

open Filter Real Asymptotics Erdos200

/-- `PrimeNumberTheoremAnd`'s `chebyshev_asymptotic : θ ~[atTop] id` is definitionally the
hypothesis `PNT` assumed in `Erdos200.Asymptotic`. -/
theorem pnt_of_primeNumberTheoremAnd : PNT := by
  have hz : ∀ᶠ x : ℝ in atTop, id x ≠ 0 :=
    (eventually_gt_atTop (0 : ℝ)).mono fun x hx ↦ ne_of_gt hx
  have h := (isEquivalent_iff_tendsto_one hz).mp chebyshev_asymptotic
  show Tendsto (fun x : ℝ ↦ Chebyshev.theta x / x) atTop (nhds 1)
  exact h

/-- **The Erdos 200 upper bound, unconditionally.**

This is the statement left as `sorry` at `Erdos200.erdos_200.variants.upper` in
`formal-conjectures`, with no remaining hypothesis. -/
theorem erdos_200_variants_upper_unconditional :
    ∃ (o : ℕ → ℝ) (_ : o =o[atTop] (1 : ℕ → ℝ)),
      ∀ n, longestPrimeArithmeticProgressions n ≤ (1 + o n) * Real.log n :=
  erdos_200_variants_upper pnt_of_primeNumberTheoremAnd

end ErdosProblem200
