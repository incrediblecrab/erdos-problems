/-
Copyright (c) 2026 Max Marquardt.
Released under the MIT licence; see LICENSE at the repository root.
-/
import Erdos40.Equivalence
import FormalConjectures.ErdosProblems.«158»

/-!
# The strong form of problem 40 is problem 158 for every `g`

`Erdos158.B2 g A` (from `FormalConjectures`) says that for every `n` at most `g` pairs `(a, a')` with `a + a' = n`, `a ≤ a'` and `a, a' ∈ A` exist. Problem 158 asks whether every infinite `B2 2` set has `liminf |A ∩ {1, …, N}| / N^{1/2} = 0`; `Erdos158.erdos_158` states this with `|A ∩ [0, N)| N^{-1/2}`, which changes the count by at most one.

* `sumRep_le_of_B2`, `B2_of_sumRep_le`: `A` is `B2 g` for some `g` if and only if `sumRep A` is bounded, since `sumRep` counts ordered pairs and swapping coordinates maps the pairs with `a > a'` injectively into those with `a ≤ a'`.
* `count_sq_le`: if `sumRep A ≤ C` then `|A ∩ [0, N)|² ≤ 2 N C`, so `|A ∩ [0, N)| N^{-1/2}` stays bounded.
* `strong_iff_forall_B2`: the strong form `Erdos40For (fun _ ↦ 1)` holds if and only if every `B2 g` set, for every `g`, has `liminf |A ∩ [0, N)| N^{-1/2} = 0`.
* `weaker_iff_forall_B2`: with `exists_iff`, the same holds for the existence question of problem 40.
* `strong_implies_erdos_158`: the case `g = 2`, stated as the right side of `Erdos158.erdos_158`.

`Erdos158.erdos_158` counts `A ∩ Set.Iio N` and `Erdos40For` counts `A ∩ Set.Icc 1 N`. The two counts differ by at most one, and the proofs absorb the difference. In `ℝ`, `liminf` of a sequence tending to `+∞` is `0` by convention. `count_sq_le` is what excludes that case in the converse direction, and the forward direction shows directly that `0` is the greatest eventual lower bound.
-/

namespace ErdosProblem40

open Filter Set AdditiveCombinatorics Asymptotics Erdos40

/-- A `B2 g` set has `sumRep A n ≤ 2 g`: `sumRep` counts ordered pairs, and swapping coordinates maps the pairs with `a > a'` injectively into those with `a ≤ a'`. -/
theorem sumRep_le_of_B2 {g : ℕ} {A : Set ℕ} (hA : Erdos158.B2 g A) (n : ℕ) : sumRep A n ≤ 2 * g := by
  classical
  rw [sumRep_def]
  set s := (Finset.HasAntidiagonal.antidiagonal n).filter (fun p : ℕ × ℕ ↦ p.1 ∈ A ∧ p.2 ∈ A) with hs
  have hsplit := Finset.card_filter_add_card_filter_not (s := s) (fun p : ℕ × ℕ ↦ p.1 ≤ p.2)
  have h1 : (s.filter (fun p : ℕ × ℕ ↦ p.1 ≤ p.2)).card ≤ g := by
    have hB := hA n
    have hset : ((s.filter (fun p : ℕ × ℕ ↦ p.1 ≤ p.2) : Finset (ℕ × ℕ)) : Set (ℕ × ℕ)) = {x : ℕ × ℕ | x.1 + x.2 = n ∧ x.1 ≤ x.2 ∧ x.1 ∈ A ∧ x.2 ∈ A} := by
      ext x
      simp only [hs, Finset.coe_filter, Finset.mem_filter, Finset.HasAntidiagonal.mem_antidiagonal, Set.mem_ofPred_eq]
      tauto
    rw [← hset, Set.encard_coe_eq_coe_finsetCard] at hB
    exact_mod_cast hB
  have h2 : (s.filter (fun p : ℕ × ℕ ↦ ¬ p.1 ≤ p.2)).card ≤ (s.filter (fun p : ℕ × ℕ ↦ p.1 ≤ p.2)).card := by
    apply Finset.card_le_card_of_injOn Prod.swap
    · intro p hp
      simp only [hs, Finset.coe_filter, Finset.mem_filter, Finset.HasAntidiagonal.mem_antidiagonal, Set.mem_ofPred_eq] at hp ⊢
      refine ⟨⟨?_, hp.1.2.2, hp.1.2.1⟩, ?_⟩
      · simp only [Prod.fst_swap, Prod.snd_swap]; omega
      · simp only [Prod.fst_swap, Prod.snd_swap]; omega
    · intro p _ q _ h
      exact Prod.swap_injective h
  omega

/-- A set whose `sumRep` is bounded by `C` is `B2 C`: the pairs with `a ≤ a'` are among the ordered pairs that `sumRep` counts. -/
theorem B2_of_sumRep_le {A : Set ℕ} {C : ℕ} (h : ∀ n, sumRep A n ≤ C) : Erdos158.B2 C A := by
  classical
  intro n
  have hsub : {x : ℕ × ℕ | x.1 + x.2 = n ∧ x.1 ≤ x.2 ∧ x.1 ∈ A ∧ x.2 ∈ A} ⊆ (((Finset.HasAntidiagonal.antidiagonal n).filter (fun p : ℕ × ℕ ↦ p.1 ∈ A ∧ p.2 ∈ A) : Finset (ℕ × ℕ)) : Set (ℕ × ℕ)) := by
    intro x hx
    simp only [Finset.coe_filter, Finset.HasAntidiagonal.mem_antidiagonal, Set.mem_ofPred_eq] at hx ⊢
    exact ⟨hx.1, hx.2.2⟩
  refine (Set.encard_mono hsub).trans ?_
  rw [Set.encard_coe_eq_coe_finsetCard]
  have := h n
  rw [sumRep_def] at this
  exact_mod_cast this

/-- If `sumRep A ≤ C` then `|A ∩ [0, N)|² ≤ 2 N C`: every pair from `A ∩ [0, N)` is an ordered representation of some `n < 2 N`. -/
theorem count_sq_le {A : Set ℕ} {C : ℕ} (h : ∀ n, sumRep A n ≤ C) (N : ℕ) : (A ∩ Iio N).ncard ^ 2 ≤ 2 * N * C := by
  classical
  set s := (Finset.range N).filter (· ∈ A) with hs
  have hcount : (A ∩ Iio N).ncard = s.card := by
    rw [← Set.ncard_coe_finset]
    congr 1
    ext x
    simp [hs, and_comm]
  have hsub : s ×ˢ s ⊆ (Finset.range (2 * N)).biUnion (fun n ↦ (Finset.HasAntidiagonal.antidiagonal n).filter (fun p : ℕ × ℕ ↦ p.1 ∈ A ∧ p.2 ∈ A)) := by
    intro p hp
    simp only [hs, Finset.mem_product, Finset.mem_filter, Finset.mem_range] at hp
    simp only [Finset.mem_biUnion, Finset.mem_range, Finset.mem_filter, Finset.HasAntidiagonal.mem_antidiagonal]
    exact ⟨p.1 + p.2, by omega, rfl, hp.1.2, hp.2.2⟩
  rw [hcount, sq, ← Finset.card_product]
  calc (s ×ˢ s).card ≤ _ := Finset.card_le_card hsub
    _ ≤ ∑ n ∈ Finset.range (2 * N), ((Finset.HasAntidiagonal.antidiagonal n).filter (fun p : ℕ × ℕ ↦ p.1 ∈ A ∧ p.2 ∈ A)).card := Finset.card_biUnion_le
    _ = ∑ n ∈ Finset.range (2 * N), sumRep A n := by simp only [sumRep_def]
    _ ≤ ∑ _n ∈ Finset.range (2 * N), C := Finset.sum_le_sum fun n _ ↦ h n
    _ = 2 * N * C := by simp

theorem rpow_neg_half (N : ℕ) : (N : ℝ) ^ (- 1 / 2 : ℝ) = (√(N : ℝ))⁻¹ := by
  rw [show (- 1 / 2 : ℝ) = -(1 / 2) by ring, Real.rpow_neg (Nat.cast_nonneg _), ← Real.sqrt_eq_rpow]

theorem tendsto_inv_sqrt : Tendsto (fun N : ℕ ↦ (√(N : ℝ))⁻¹) atTop (nhds 0) :=
  (Real.tendsto_sqrt_atTop.comp tendsto_natCast_atTop_atTop).inv_tendsto_atTop

/-- Assuming the strong form, every `A` whose `sumRep` is bounded has `liminf |A ∩ [0, N)| N^{-1/2} = 0`. -/
theorem liminf_eq_zero_of_bound (hS : Erdos40For (fun _ ↦ 1)) {A : Set ℕ} {C : ℕ} (hC : ∀ n, sumRep A n ≤ C) :
    liminf (fun N : ℕ => (A ∩ .Iio N).ncard * (N : ℝ) ^ (- 1 / 2 : ℝ)) atTop = 0 := by
  set f : ℕ → ℝ := fun N ↦ (A ∩ .Iio N).ncard * (N : ℝ) ^ (- 1 / 2 : ℝ) with hf
  have hnot : ¬ (fun N : ℕ ↦ √N / (fun _ ↦ (1 : ℝ)) N) =O[atTop] (fun N ↦ ((A ∩ .Icc 1 N).ncard : ℝ)) :=
    fun h ↦ limsup_ne_top_of_bound hC (hS A h)
  have hfreq : ∀ a > 0, ∃ᶠ N in atTop, f N < a := by
    intro a ha
    rw [isBigO_iff] at hnot
    push Not at hnot
    have hc := hnot (2 / a)
    have hev : ∀ᶠ N : ℕ in atTop, (√(N : ℝ))⁻¹ < a / 2 := tendsto_inv_sqrt.eventually (gt_mem_nhds (by linarith))
    refine (hc.and_eventually (hev.and (eventually_ge_atTop 1))).mono ?_
    rintro N ⟨hN, hsmall, hN1⟩
    have hsq : 0 < √(N : ℝ) := Real.sqrt_pos.2 (by exact_mod_cast hN1)
    rw [div_one, Real.norm_of_nonneg (Real.sqrt_nonneg _), Real.norm_of_nonneg (Nat.cast_nonneg _)] at hN
    have hcount : ((A ∩ .Iio N).ncard : ℝ) ≤ (A ∩ .Icc 1 N).ncard + 1 := by
      have hsub : A ∩ .Iio N ⊆ insert 0 (A ∩ .Icc 1 N) := by
        intro x ⟨hxA, hxN⟩
        rcases Nat.eq_zero_or_pos x with rfl | hx
        · exact Set.mem_insert _ _
        · exact Set.mem_insert_of_mem _ ⟨hxA, hx, (Set.mem_Iio.1 hxN).le⟩
      have hfin : (insert 0 (A ∩ .Icc 1 N)).Finite := ((Set.finite_Icc 1 N).inter_of_right A).insert 0
      have := (Set.ncard_le_ncard hsub hfin).trans (Set.ncard_insert_le 0 (A ∩ .Icc 1 N))
      exact_mod_cast this
    simp only [hf]
    rw [rpow_neg_half]
    have e1 : ((A ∩ .Icc 1 N).ncard : ℝ) / √(N : ℝ) < a / 2 := by
      have : 2 / a * ((A ∩ .Icc 1 N).ncard : ℝ) < √(N : ℝ) := hN
      rw [div_mul_eq_mul_div, div_lt_iff₀ ha] at this
      rw [div_lt_iff₀ hsq]
      linarith
    calc ((A ∩ .Iio N).ncard : ℝ) * (√(N : ℝ))⁻¹
        ≤ (((A ∩ .Icc 1 N).ncard : ℝ) + 1) * (√(N : ℝ))⁻¹ :=
          mul_le_mul_of_nonneg_right hcount (inv_nonneg.2 hsq.le)
      _ = ((A ∩ .Icc 1 N).ncard : ℝ) / √(N : ℝ) + (√(N : ℝ))⁻¹ := by ring
      _ < a / 2 + a / 2 := add_lt_add e1 hsmall
      _ = a := by ring
  have h0 : (0 : ℝ) ∈ {a : ℝ | ∀ᶠ N in atTop, a ≤ f N} :=
    Eventually.of_forall fun N ↦ mul_nonneg (Nat.cast_nonneg _) (Real.rpow_nonneg (Nat.cast_nonneg _) _)
  have hle : ∀ a ∈ {a : ℝ | ∀ᶠ N in atTop, a ≤ f N}, a ≤ 0 := by
    intro a ha
    by_contra hpos
    push Not at hpos
    obtain ⟨N, hN1, hN2⟩ := ((hfreq a hpos).and_eventually ha).exists
    linarith
  rw [liminf_eq]
  exact le_antisymm (csSup_le ⟨0, h0⟩ hle) (le_csSup ⟨0, hle⟩ h0)

/-- The strong form of problem 40 holds if and only if, for every `g`, every `B2 g` set has `liminf |A ∩ [0, N)| N^{-1/2} = 0`. For `g = 1` the right side is Erdős's theorem on Sidon sets; for `g = 2` it is problem 158. -/
theorem strong_iff_forall_B2 : Erdos40For (fun _ ↦ 1) ↔ ∀ (g : ℕ) (A : Set ℕ), Erdos158.B2 g A →
    liminf (fun N : ℕ => (A ∩ .Iio N).ncard * (N : ℝ) ^ (- 1 / 2 : ℝ)) atTop = 0 := by
  refine ⟨fun hS g A hB ↦ liminf_eq_zero_of_bound hS (sumRep_le_of_B2 hB), fun H A hA ↦ ?_⟩
  by_contra hne
  obtain ⟨C, hC⟩ := exists_bound_of_limsup_ne_top hne
  have h0 := H C A (B2_of_sumRep_le hC)
  set f : ℕ → ℝ := fun N ↦ (A ∩ .Iio N).ncard * (N : ℝ) ^ (- 1 / 2 : ℝ) with hf
  obtain ⟨c, hc, hcA⟩ := hA.exists_pos
  rw [isBigOWith_iff] at hcA
  have hlow : ∀ᶠ N : ℕ in atTop, 1 / (2 * c) ≤ f N := by
    have hev : ∀ᶠ N : ℕ in atTop, (√(N : ℝ))⁻¹ ≤ 1 / (2 * c) := tendsto_inv_sqrt.eventually (ge_mem_nhds (by positivity))
    filter_upwards [hcA, hev, eventually_ge_atTop 1] with N hN hsmall hN1
    have hsq : 0 < √(N : ℝ) := Real.sqrt_pos.2 (by exact_mod_cast hN1)
    rw [div_one, Real.norm_of_nonneg (Real.sqrt_nonneg _), Real.norm_of_nonneg (Nat.cast_nonneg _)] at hN
    have hcount : ((A ∩ .Icc 1 N).ncard : ℝ) ≤ (A ∩ .Iio N).ncard + 1 := by
      have hsub : A ∩ .Icc 1 N ⊆ insert N (A ∩ .Iio N) := by
        intro x ⟨hxA, hx1, hxN⟩
        rcases eq_or_lt_of_le hxN with rfl | hlt
        · exact Set.mem_insert _ _
        · exact Set.mem_insert_of_mem _ ⟨hxA, hlt⟩
      have hfin : (insert N (A ∩ .Iio N)).Finite := ((Set.finite_Iio N).inter_of_right A).insert N
      have := (Set.ncard_le_ncard hsub hfin).trans (Set.ncard_insert_le N (A ∩ .Iio N))
      exact_mod_cast this
    have e1 : √(N : ℝ) / c ≤ (A ∩ .Icc 1 N).ncard := by
      rw [div_le_iff₀ hc]
      linarith
    simp only [hf]
    rw [rpow_neg_half]
    calc 1 / (2 * c) = 1 / c - 1 / (2 * c) := by field_simp; ring
      _ ≤ 1 / c - (√(N : ℝ))⁻¹ := by linarith
      _ = (√(N : ℝ) / c - 1) * (√(N : ℝ))⁻¹ := by field_simp
      _ ≤ ((A ∩ .Iio N).ncard : ℝ) * (√(N : ℝ))⁻¹ := mul_le_mul_of_nonneg_right (by linarith) (inv_nonneg.2 hsq.le)
  have hup : ∀ N : ℕ, f N ≤ √(2 * C) := by
    intro N
    simp only [hf]
    rw [rpow_neg_half]
    rcases Nat.eq_zero_or_pos N with rfl | hN
    · simp
    have hsq : 0 < √(N : ℝ) := Real.sqrt_pos.2 (by exact_mod_cast hN)
    have hc2 : ((A ∩ .Iio N).ncard : ℝ) ^ 2 ≤ 2 * C * N := by
      have := count_sq_le hC N
      have : (((A ∩ .Iio N).ncard ^ 2 : ℕ) : ℝ) ≤ ((2 * N * C : ℕ) : ℝ) := by exact_mod_cast this
      push_cast at this
      linarith
    have hle : ((A ∩ .Iio N).ncard : ℝ) ≤ √(2 * C) * √(N : ℝ) := by
      rw [← Real.sqrt_mul (by positivity)]
      exact Real.le_sqrt_of_sq_le hc2
    rw [← div_eq_mul_inv, div_le_iff₀ hsq]
    exact hle
  have hbdd : BddAbove {a : ℝ | ∀ᶠ N in atTop, a ≤ f N} := by
    refine ⟨√(2 * C), fun a ha ↦ ?_⟩
    obtain ⟨N, hN⟩ := ha.exists
    exact hN.trans (hup N)
  have hpos : 0 < 1 / (2 * c) := by positivity
  have := le_csSup hbdd hlow
  rw [liminf_eq] at h0
  linarith

/-- The existence question of problem 40, as the right side of `Erdos40.erdos_40.variants.weaker`, holds if and only if every `B2 g` set, for every `g`, has `liminf |A ∩ [0, N)| N^{-1/2} = 0`. -/
theorem weaker_iff_forall_B2 : (∃ g : ℕ → ℝ, Tendsto g atTop atTop ∧ Erdos40For g) ↔ ∀ (g : ℕ) (A : Set ℕ), Erdos158.B2 g A →
    liminf (fun N : ℕ => (A ∩ .Iio N).ncard * (N : ℝ) ^ (- 1 / 2 : ℝ)) atTop = 0 :=
  exists_iff.trans strong_iff_forall_B2

/-- The strong form of problem 40 implies the right side of `Erdos158.erdos_158`. -/
theorem strong_implies_erdos_158 (hS : Erdos40For (fun _ ↦ 1)) :
    ∀ A : Set ℕ, A.Infinite → Erdos158.B2 2 A →
      liminf (fun N : ℕ => (A ∩ .Iio N).ncard * (N : ℝ) ^ (- 1 / 2 : ℝ)) atTop = 0 :=
  fun A _ hB ↦ strong_iff_forall_B2.1 hS 2 A hB

end ErdosProblem40
