/-
Copyright (c) 2026 Max Marquardt.
Released under the MIT licence; see LICENSE at the repository root.
-/
import Erdos40.Basic

/-!
# The existence question of Erdős problem 40 is the "stronger conjecture" of Erdős problem 28

`Erdos40.Erdos40For g` (from `FormalConjectures`) says: every `A ⊆ ℕ` with `|A ∩ [1, N]| ≫ √N / g N` has `limsup sumRep A = ⊤`. With `g = 1` this is the conjecture that the site's page for problem 28 calls "another stronger conjecture": the hypothesis `|A ∩ [1, N]| ≫ N^{1/2}` alone forces unbounded representation counts.

`exists_iff` proves that some `g → ∞` has the property if and only if `g = 1` has it. The forward direction is immediate, since `g ≥ 1` eventually. The converse is a compactness argument. For a fixed bound `C`, the sets `A` with `sumRep A ≤ C` form a closed, hence compact, subset of `ℕ → Bool`. So if every such `A` has `liminf |A ∩ [1, N]| / √N = 0`, the dips below `δ √N` happen inside one window `[N₀, M]` that works for all of them at once (`uniform_dip`). Chaining those windows with `δ → 0` and `C → ∞`, and letting `g` grow by one per window, gives a single `g → ∞`.

This settles nothing about which `g` have the property. Both sides of `exists_iff` are open, and each implies the Erdős–Turán conjecture (`strong_implies_erdos_28`).
-/

namespace ErdosProblem40

open Filter Set AdditiveCombinatorics Asymptotics Erdos40

/-- Assuming the strong form, for a bound `C`, a threshold `δ > 0` and a start `N₀`, one window `[N₀, M]` catches a dip `|A ∩ [1, N]| < δ √N` of every `A` whose representation function is bounded by `C`. -/
theorem uniform_dip (hS : Erdos40For (fun _ ↦ 1)) (C : ℕ) {δ : ℝ} (hδ : 0 < δ) (N₀ : ℕ) :
    ∃ M, N₀ ≤ M ∧ ∀ A : Set ℕ, (∀ n, sumRep A n ≤ C) →
      ∃ N, N₀ ≤ N ∧ N ≤ M ∧ ((A ∩ Icc 1 N).ncard : ℝ) < δ * √N := by
  by_contra hcon
  push Not at hcon
  let S : (ℕ → Bool) → Set ℕ := fun f ↦ {i | f i = true}
  let T : ℕ → Set (ℕ → Bool) := fun m ↦
    {f | (∀ n, sumRep (S f) n ≤ C) ∧
      ∀ N, N₀ ≤ N → N ≤ N₀ + m → δ * √N ≤ ((S f ∩ Icc 1 N).ncard : ℝ)}
  have hmono : ∀ m, T (m + 1) ⊆ T m := by
    rintro m f ⟨h1, h2⟩
    exact ⟨h1, fun N hN hNm ↦ h2 N hN (by omega)⟩
  have hne : ∀ m, (T m).Nonempty := by
    intro m
    classical
    obtain ⟨A, hA, hdense⟩ := hcon (N₀ + m) (by omega)
    have hSA : S (fun i ↦ decide (i ∈ A)) = A := by
      ext i
      simp [S]
    refine ⟨fun i ↦ decide (i ∈ A), ?_, ?_⟩
    · intro n
      rw [hSA]
      exact hA n
    · intro N hN hNm
      rw [hSA]
      exact hdense N hN hNm
  have hclosed : ∀ m, IsClosed (T m) := by
    intro m
    have e : T m = (⋂ n, {f : ℕ → Bool | sumRep (S f) n ≤ C}) ∩
        ⋂ N, {f : ℕ → Bool | N₀ ≤ N → N ≤ N₀ + m →
          δ * √N ≤ ((S f ∩ Icc 1 N).ncard : ℝ)} := by
      ext f
      simp only [T, mem_ofPred_eq, mem_inter_iff, mem_iInter]
    rw [e]
    refine (isClosed_iInter fun n ↦ ?_).inter (isClosed_iInter fun N ↦ ?_)
    · refine (isClopen_of_local n _ fun f f' hff' hf ↦ ?_).isClosed
      rwa [← sumRep_congr (A := S f) (B := S f') (n := n) fun i hi ↦ by simp [S, hff' i hi]]
    · refine (isClopen_of_local N _ fun f f' hff' hf ↦ ?_).isClosed
      intro h1 h2
      rw [← count_congr (A := S f) (B := S f') (N := N) fun i hi ↦ by simp [S, hff' i hi]]
      exact hf h1 h2
  obtain ⟨f, hf⟩ := IsCompact.nonempty_iInter_of_sequence_nonempty_isCompact_isClosed T hmono hne
    (hclosed 0).isCompact hclosed
  rw [mem_iInter] at hf
  have hbound : ∀ n, sumRep (S f) n ≤ C := (hf 0).1
  have hdense : ∀ N, N₀ ≤ N → δ * √N ≤ ((S f ∩ Icc 1 N).ncard : ℝ) :=
    fun N hN ↦ (hf N).2 N hN (by omega)
  apply limsup_ne_top_of_bound hbound
  apply hS (S f)
  apply IsBigO.of_bound (1 / δ)
  filter_upwards [eventually_ge_atTop N₀] with N hN
  rw [div_one, Real.norm_of_nonneg (Real.sqrt_nonneg _), Real.norm_natCast, one_div,
    ← div_eq_inv_mul, le_div_iff₀ hδ, mul_comm]
  exact hdense N hN

/-- The converse direction of `exists_iff`: the strong form yields a single `g → ∞`. -/
theorem exists_of_strong (hS : Erdos40For (fun _ ↦ 1)) :
    ∃ g : ℕ → ℝ, Tendsto g atTop atTop ∧ Erdos40For g := by
  have key : ∀ k N₀ : ℕ, ∃ M, N₀ ≤ M ∧ ∀ A : Set ℕ, (∀ n, sumRep A n ≤ k) →
      ∃ N, N₀ ≤ N ∧ N ≤ M ∧ ((A ∩ Icc 1 N).ncard : ℝ) < 1 / ((k : ℝ) + 1) ^ 2 * √N :=
    fun k N₀ ↦ uniform_dip hS k (by positivity) N₀
  choose M hM₀ hM using key
  -- `s k` is where window `k` starts; window `k` is `[s k, M k (s k)]` and the next starts one later.
  let s : ℕ → ℕ := fun k ↦ Nat.rec (motive := fun _ ↦ ℕ) 1 (fun k sk ↦ M k sk + 1) k
  have s_succ : ∀ k, s (k + 1) = M k (s k) + 1 := fun k ↦ rfl
  have s_strict : StrictMono s := strictMono_nat_of_lt_succ fun k ↦ by
    rw [s_succ]
    have := hM₀ k (s k)
    omega
  have s_ge : ∀ k, k + 1 ≤ s k := by
    intro k
    induction k with
    | zero => exact le_refl 1
    | succ k ih =>
      rw [s_succ]
      have := hM₀ k (s k)
      omega
  -- `g N - 1` is the index of the window containing `N`.
  let j : ℕ → ℕ := fun N ↦ @Nat.findGreatest (fun i ↦ s i ≤ N) (fun i ↦ Nat.decLe (s i) N) N
  let g : ℕ → ℝ := fun N ↦ (j N : ℝ) + 1
  refine ⟨g, ?_, ?_⟩
  · rw [tendsto_atTop_atTop]
    intro b
    refine ⟨s ⌈b⌉₊, fun N hN ↦ ?_⟩
    have h1 : ⌈b⌉₊ ≤ j N :=
      @Nat.le_findGreatest ⌈b⌉₊ (fun i ↦ s i ≤ N) (fun i ↦ Nat.decLe (s i) N) N
        (by have := s_ge ⌈b⌉₊; omega) hN
    have h2 : b ≤ (⌈b⌉₊ : ℝ) := Nat.le_ceil b
    have h3 : ((⌈b⌉₊ : ℕ) : ℝ) ≤ (j N : ℝ) := by exact_mod_cast h1
    show b ≤ (j N : ℝ) + 1
    linarith
  · intro A hA
    by_contra hne
    obtain ⟨C, hC⟩ := exists_bound_of_limsup_ne_top hne
    obtain ⟨c, hc, hcw⟩ := hA.exists_pos
    rw [IsBigOWith] at hcw
    obtain ⟨N₁, hN₁⟩ := eventually_atTop.1 hcw
    let k : ℕ := max C (max ⌈c⌉₊ N₁)
    have hkC : C ≤ k := le_max_left _ _
    have hkc : c ≤ (k : ℝ) + 1 := by
      have : ⌈c⌉₊ ≤ k := le_trans (le_max_left _ _) (le_max_right _ _)
      have : (⌈c⌉₊ : ℝ) ≤ k := by exact_mod_cast this
      linarith [Nat.le_ceil c]
    have hkN : N₁ ≤ k := le_trans (le_max_right _ _) (le_max_right _ _)
    obtain ⟨N, hN1, hN2, hlt⟩ := hM k (s k) A (fun n ↦ (hC n).trans hkC)
    have hjN : j N ≤ k := by
      by_contra hj
      push Not at hj
      have hspec : s (j N) ≤ N :=
        @Nat.findGreatest_spec 0 (fun i ↦ s i ≤ N) (fun i ↦ Nat.decLe (s i) N) N (Nat.zero_le _)
          (by have := s_ge 0; have := s_strict.monotone (Nat.zero_le k); omega)
      have := s_strict.monotone (Nat.succ_le_of_lt hj)
      rw [s_succ] at this
      omega
    have hgN : g N ≤ (k : ℝ) + 1 := by
      show (j N : ℝ) + 1 ≤ (k : ℝ) + 1
      exact_mod_cast Nat.add_le_add_right hjN 1
    have hg1 : 1 ≤ g N := by
      show (1 : ℝ) ≤ (j N : ℝ) + 1
      linarith [(Nat.cast_nonneg (j N) : (0 : ℝ) ≤ j N)]
    have hNpos : 0 < N := by have := s_ge k; omega
    have hx : 0 < √(N : ℝ) := Real.sqrt_pos.2 (by exact_mod_cast hNpos)
    have hbd := hN₁ N (by have := s_ge k; omega)
    rw [Real.norm_of_nonneg (div_nonneg (Real.sqrt_nonneg _) (by linarith)), Real.norm_natCast,
      div_le_iff₀ (by linarith)] at hbd
    -- hbd : √N ≤ c * a * g N, hlt : a < √N / (k+1)^2, with a = |A ∩ [1, N]|
    set a : ℝ := ((A ∩ Icc 1 N).ncard : ℝ)
    set K : ℝ := (k : ℝ) + 1
    have hK : 1 ≤ K := by
      show (1 : ℝ) ≤ (k : ℝ) + 1
      linarith [(Nat.cast_nonneg k : (0 : ℝ) ≤ k)]
    have ha : 0 ≤ a := Nat.cast_nonneg _
    have hlt' : a * K ^ 2 < √(N : ℝ) := by
      rw [one_div, ← div_eq_inv_mul, lt_div_iff₀ (by positivity)] at hlt
      exact hlt
    have step1 : c * a * g N ≤ c * a * K := mul_le_mul_of_nonneg_left hgN (by positivity)
    have step2 : c * a * K ≤ K * a * K := by
      apply mul_le_mul_of_nonneg_right _ (by linarith)
      exact mul_le_mul_of_nonneg_right hkc ha
    have : K * a * K = a * K ^ 2 := by ring
    linarith

/-- **Theorem A.** Some `g → ∞` has the property of Erdős problem 40 if and only if `g = 1` has it, which is the strong form of the Erdős–Turán conjecture stated on the page for problem 28. The left side is the right side of `Erdos40.erdos_40.variants.weaker` in `FormalConjectures`. -/
theorem exists_iff :
    (∃ g : ℕ → ℝ, Tendsto g atTop atTop ∧ Erdos40For g) ↔ Erdos40For (fun _ ↦ 1) := by
  constructor
  · rintro ⟨g, hg, h⟩ A hA
    apply h A
    refine IsBigO.trans ?_ hA
    apply IsBigO.of_bound 1
    filter_upwards [hg.eventually_ge_atTop 1] with N hN
    rw [one_mul, Real.norm_of_nonneg (div_nonneg (Real.sqrt_nonneg _) (by linarith)),
      Real.norm_of_nonneg (div_nonneg (Real.sqrt_nonneg _) zero_le_one), div_one]
    exact div_le_self (Real.sqrt_nonneg _) hN
  · exact exists_of_strong

/-- Theorem A restated for the answer set in `Erdos40.erdos_40`: that set is nonempty exactly when the strong form holds. -/
theorem answerSet_nonempty_iff :
    {g : ℕ → ℝ | Tendsto g atTop atTop ∧ Erdos40For g}.Nonempty ↔ Erdos40For (fun _ ↦ 1) :=
  exists_iff

/-- The strong form implies the Erdős–Turán conjecture as stated in `Erdos28.erdos_28`. The work is `FormalConjectures`' own `erdos_40.variants.implies_erdos_28`; Theorem A supplies its `g`. -/
theorem strong_implies_erdos_28 (hS : Erdos40For (fun _ ↦ 1)) : type_of% Erdos28.erdos_28 := by
  obtain ⟨g, hg, h⟩ := exists_of_strong hS
  exact erdos_40.variants.implies_erdos_28 g hg h

/-- The answer set is closed downward under `=O`: if `g` has the property and `g' = O(g)` is eventually nonzero, then `g'` has it too. -/
theorem Erdos40For.of_isBigO {g g' : ℕ → ℝ} (h : Erdos40For g) (hg' : ∀ᶠ N in atTop, g' N ≠ 0)
    (hle : g' =O[atTop] g) : Erdos40For g' := by
  intro A hA
  apply h A
  refine IsBigO.trans ?_ hA
  obtain ⟨K, hK⟩ := hle.bound
  apply IsBigO.of_bound K
  filter_upwards [hK, hg'] with N h1 h2
  rw [norm_div, norm_div]
  have hg'pos : 0 < ‖g' N‖ := norm_pos_iff.2 h2
  have hgpos : 0 < ‖g N‖ := by
    by_contra hg0
    push Not at hg0
    have : ‖g N‖ = 0 := le_antisymm hg0 (norm_nonneg _)
    rw [this, mul_zero] at h1
    linarith
  have hKpos : 0 < K := by
    by_contra hK0
    push Not at hK0
    nlinarith
  rw [div_le_iff₀ hgpos, mul_div_assoc', div_mul_eq_mul_div, le_div_iff₀ hg'pos]
  have hs : 0 ≤ ‖(√(N : ℝ) : ℝ)‖ := norm_nonneg _
  nlinarith

theorem sumRep_univ (n : ℕ) : sumRep univ n = n + 1 := by
  classical
  rw [sumRep_def]
  simp

/-- The hypotheses can hold together: `ℕ` itself satisfies the density hypothesis for `g = 1`, and its representation function is unbounded. -/
theorem univ_hypothesis_and_conclusion :
    ((fun N : ℕ ↦ √N / (fun _ ↦ (1 : ℝ)) N) =O[atTop] fun N ↦ ((univ ∩ Icc 1 N).ncard : ℝ)) ∧
      limsup (fun N ↦ (sumRep univ N : ℕ∞)) atTop = ⊤ := by
  constructor
  · apply IsBigO.of_bound 1
    filter_upwards [eventually_ge_atTop 1] with N hN
    have hcount : ((univ : Set ℕ) ∩ Icc 1 N).ncard = N := by
      rw [univ_inter, Set.ncard_eq_toFinset_card', Set.toFinset_Icc, Nat.card_Icc]
      omega
    rw [hcount, div_one, one_mul, Real.norm_of_nonneg (Real.sqrt_nonneg _), Real.norm_natCast]
    have h1 : (1 : ℝ) ≤ N := by exact_mod_cast hN
    calc √(N : ℝ) ≤ √(N : ℝ) * √(N : ℝ) := le_mul_of_one_le_left (Real.sqrt_nonneg _)
            (Real.one_le_sqrt.2 h1)
      _ = N := Real.mul_self_sqrt (Nat.cast_nonneg _)
  · by_contra hne
    obtain ⟨C, hC⟩ := exists_bound_of_limsup_ne_top hne
    have := hC C
    rw [sumRep_univ] at this
    omega

theorem sumRep_singleton_one_le (n : ℕ) : sumRep ({1} : Set ℕ) n ≤ 1 := by
  classical
  rw [sumRep_def]
  apply Finset.card_le_one.2
  intro p hp q hq
  simp only [Finset.mem_filter, mem_singleton_iff] at hp hq
  exact Prod.ext (hp.2.1.trans hq.2.1.symm) (hp.2.2.trans hq.2.2.symm)

/-- The property is not vacuous: `g(N) = √N` does not have it. The witness is `A = {1}`. -/
theorem not_erdos40For_sqrt : ¬ Erdos40For (fun N ↦ √N) := by
  intro h
  apply limsup_ne_top_of_bound sumRep_singleton_one_le
  apply h
  apply IsBigO.of_bound 1
  filter_upwards [eventually_ge_atTop 1] with N hN
  have hcount : (({1} : Set ℕ) ∩ Icc 1 N).ncard = 1 := by
    have : ({1} : Set ℕ) ∩ Icc 1 N = {1} := by
      ext i
      simp only [mem_inter_iff, mem_singleton_iff, mem_Icc]
      constructor
      · exact fun h ↦ h.1
      · rintro rfl
        exact ⟨rfl, le_refl 1, hN⟩
    rw [this, ncard_singleton]
  have hx : 0 < √(N : ℝ) := Real.sqrt_pos.2 (by exact_mod_cast hN)
  rw [hcount, div_self hx.ne', Nat.cast_one, norm_one, mul_one]

/-- Any `g` with `√N = O(g N)` lacks the property, by `not_erdos40For_sqrt` and downward closure. -/
theorem not_erdos40For_of_sqrt_isBigO {g : ℕ → ℝ} (hg : (fun N : ℕ ↦ √N) =O[atTop] g) :
    ¬ Erdos40For g := by
  intro h
  refine not_erdos40For_sqrt (Erdos40For.of_isBigO h ?_ hg)
  filter_upwards [eventually_ge_atTop 1] with N hN
  exact (Real.sqrt_pos.2 (by exact_mod_cast hN)).ne'

end ErdosProblem40
