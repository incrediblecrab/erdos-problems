/-
Copyright (c) 2026 Max Marquardt.
Released under the MIT licence; see LICENSE at the repository root.
-/
import Erdos40.Equivalence

/-!
# Theorem B: no power of `N` is in the answer set

For every `ε > 0` some `A ⊆ ℕ` has bounded `sumRep A` and `|A ∩ [1, N]| ≫ √N / N^ε` (`exists_dense_bounded`), so `fun N ↦ (N : ℝ) ^ ε` lacks the property `Erdos40.Erdos40For` (`not_erdos40For_rpow`) and is not in the set of `Erdos40.erdos_40` (`rpow_not_mem_answerSet`). With `Erdos40For.of_isBigO` the same holds for every `g` with `N^ε = O(g)` (`not_erdos40For_of_rpow_isBigO`).

The proof is probabilistic, in the spirit of Erdős and Rényi, but with one random point per block instead of independent Bernoulli choices, so the counting function is deterministic. Fix `d ≥ 1` with `1 / (2 d + 2) ≤ ε`, and put `base = 2^(2d+2)` and `len j = 2^((d+2) j)`. Level `j` is `[base^j, base^(j+1))`; it is cut into consecutive blocks of length `len j`, and the random set `S ω` takes exactly one point, uniformly and independently, from each block (`ω (j, c) : Fin (len j)` is the offset in block `c`). Then:

* `count_ge`: the `2^(d J)` blocks numbered `2^(d J), …, 2^(d J + 1) - 1` of level `J` lie in `[1, base^(J+1))`, so for every outcome `|S ω ∩ [1, N]| ≥ 2^(d J)` once `N ≥ base^(J+1)`, which gives `N^(d/(2d+2)) ≤ 4^d |S ω ∩ [1, N]|` (`rpow_le_count`).
* `measure_subset_le`: a finite set `F` lies in `S ω` with probability at most `∏ x ∈ F, p x`, where `p x = 1 / len (level x)`: two points of one block are never both chosen, and distinct blocks are independent (`Measure.infinitePi_pi`).
* `measure_E_le`, `sum_powersetCard_prod_le`: the probability that `n` has at least `2 d + 2` representations `n = a + (n - a)` with `a < n - a` is at most `σ n ^ (2 d + 2)`, with `σ n = ∑_{a < n/2} p a p (n - a)`.
* `σ_le`, `level_sum_le`: for `n` in level `J + 1`, `σ n ≤ 2 base 2^d / 4^J`, and the bad events of level `J + 1` have total probability at most `base^2 (2 base 2^d)^(2d+2) / base^J`, a geometric series in `J`.
* `exists_good`: from level `J₀ + 1` on the bad events have total probability below `1`, so some outcome avoids all of them; it has `sumRep ≤ 2 (2 d + 1) + 1` beyond `base^(J₀+1)`, and `sumRep A n ≤ n + 1` below.

The construction is not optimised: it gives `|A ∩ [1, N]| ≫ N^(1/2 - 1/(2d+2))` with `sumRep` eventually at most `4 d + 3`.
-/

namespace ErdosProblem40

open Filter Set AdditiveCombinatorics Asymptotics Erdos40 MeasureTheory
open scoped NNReal ENNReal

namespace TheoremB

variable (d : ℕ)

/-- Level `j` is `[base^j, base^(j+1))`. -/
def base : ℕ := 2 ^ (2 * d + 2)

/-- Level `j` is cut into blocks of length `len j`. -/
def len (j : ℕ) : ℕ := 2 ^ ((d + 2) * j)

/-- The level of `x`: the `j` with `base^j ≤ x < base^(j+1)`, for `x ≥ 1`. -/
def level (x : ℕ) : ℕ := Nat.log (base d) x

/-- The block of `x`: its level and its block number within the level. -/
def key (x : ℕ) : ℕ × ℕ := (level d x, x / len d (level d x))

instance (j : ℕ) : NeZero (len d j) := ⟨by unfold len; positivity⟩

/-- An outcome chooses an offset in each block. -/
abbrev Ω := (i : ℕ × ℕ) → Fin (len d i.1)

/-- The offsets are independent and uniform. -/
noncomputable def P : Measure (Ω d) :=
  Measure.infinitePi fun i : ℕ × ℕ ↦ (PMF.uniformOfFintype (Fin (len d i.1))).toMeasure

instance : IsProbabilityMeasure (P d) := by unfold P; infer_instance

/-- The random set: the positive `x` whose offset in its block is the chosen one. -/
def S (ω : Ω d) : Set ℕ := {x | 1 ≤ x ∧ ((ω (key d x) : Fin _) : ℕ) = x % len d (level d x)}

/-- The probability that `x` is chosen. -/
noncomputable def p (x : ℕ) : ℝ≥0 := ((len d (level d x) : ℝ≥0))⁻¹

theorem one_lt_base : 1 < base d := Nat.one_lt_two_pow (by omega)

theorem coe_p (x : ℕ) : (p d x : ℝ≥0∞) = ((len d (level d x) : ℝ≥0∞))⁻¹ := by
  rw [p, ENNReal.coe_inv (by exact_mod_cast NeZero.ne _), ENNReal.coe_natCast]

theorem uniform_le (n : ℕ) [NeZero n] (t : Set (Fin n)) (ht : t.Subsingleton) :
    (PMF.uniformOfFintype (Fin n)).toMeasure t ≤ (n : ℝ≥0∞)⁻¹ := by
  rcases ht.eq_empty_or_singleton with rfl | ⟨a, rfl⟩
  · simp
  · rw [PMF.toMeasure_apply_singleton _ _ (MeasurableSet.singleton a), PMF.uniformOfFintype_apply,
      Fintype.card_fin]

theorem eq_of_key {x y : ℕ} (h : key d x = key d y)
    (hr : x % len d (level d x) = y % len d (level d y)) : x = y := by
  simp only [key, Prod.mk.injEq] at h
  obtain ⟨h1, h2⟩ := h
  rw [h1] at h2 hr
  calc x = len d (level d y) * (x / len d (level d y)) + x % len d (level d y) :=
        (Nat.div_add_mod _ _).symm
    _ = len d (level d y) * (y / len d (level d y)) + y % len d (level d y) := by rw [h2, hr]
    _ = y := Nat.div_add_mod _ _

/-- The key estimate: the probability that a finite set `F` lies in the random set is at most `∏ x ∈ F, p x`. Two points of one block cannot both be chosen, and points of distinct blocks are independent. -/
theorem measure_subset_le (F : Finset ℕ) :
    P d {ω | ∀ x ∈ F, x ∈ S d ω} ≤ ∏ x ∈ F, (p d x : ℝ≥0∞) := by
  classical
  let t : (i : ℕ × ℕ) → Set (Fin (len d i.1)) :=
    fun i ↦ {o | ∀ x ∈ F, key d x = i → (o : ℕ) = x % len d (level d x)}
  have hsub : {ω | ∀ x ∈ F, x ∈ S d ω} ⊆ Set.pi ((F.image (key d) : Finset _) : Set (ℕ × ℕ)) t := by
    intro ω hω i _ x hx hxi
    subst hxi
    exact (hω x hx).2
  calc P d {ω | ∀ x ∈ F, x ∈ S d ω}
      ≤ P d (Set.pi ((F.image (key d) : Finset _) : Set (ℕ × ℕ)) t) := measure_mono hsub
    _ = ∏ i ∈ F.image (key d), (PMF.uniformOfFintype (Fin (len d i.1))).toMeasure (t i) := by
        unfold P
        exact Measure.infinitePi_pi _ (fun _ _ ↦ MeasurableSet.of_discrete)
    _ ≤ ∏ i ∈ F.image (key d), ∏ x ∈ F with key d x = i, (p d x : ℝ≥0∞) := by
        apply Finset.prod_le_prod'
        intro i hi
        by_cases h2 : ∃ x ∈ F, ∃ y ∈ F, key d x = i ∧ key d y = i ∧ x ≠ y
        · obtain ⟨x, hx, y, hy, hxi, hyi, hne⟩ := h2
          have : t i = ∅ := by
            ext o
            simp only [Set.mem_empty_iff_false, iff_false]
            intro ho
            apply hne
            apply eq_of_key d (hxi.trans hyi.symm)
            rw [← ho x hx hxi, ← ho y hy hyi]
          rw [this, measure_empty]
          exact zero_le
        · push Not at h2
          obtain ⟨x, hx, hxi⟩ := Finset.mem_image.1 hi
          have hG : F.filter (fun y ↦ key d y = i) = {x} := by
            ext y
            simp only [Finset.mem_filter, Finset.mem_singleton]
            constructor
            · rintro ⟨hy, hyi⟩
              by_contra hne
              exact hne (h2 y hy x hx hyi hxi)
            · rintro rfl
              exact ⟨hx, hxi⟩
          rw [hG, Finset.prod_singleton]
          subst hxi
          rw [coe_p]
          apply uniform_le
          intro o ho o' ho'
          exact Fin.ext ((ho x hx rfl).trans (ho' x hx rfl).symm)
    _ = ∏ x ∈ F, (p d x : ℝ≥0∞) :=
        Finset.prod_fiberwise_of_maps_to (fun x hx ↦ Finset.mem_image_of_mem _ hx) _

/-- `∑_{|F| = k} ∏_{a ∈ F} q a ≤ (∑ q)^k`. -/
theorem sum_powersetCard_prod_le (s : Finset ℕ) (q : ℕ → ℝ≥0) (k : ℕ) :
    ∑ F ∈ s.powersetCard k, ∏ a ∈ F, q a ≤ (∑ a ∈ s, q a) ^ k := by
  classical
  induction s using Finset.induction_on generalizing k with
  | empty => cases k <;> simp
  | insert a s ha ih =>
    cases k with
    | zero => simp
    | succ k =>
      have hdisj : Disjoint (s.powersetCard (k + 1)) ((s.powersetCard k).image (insert a)) := by
        rw [Finset.disjoint_left]
        intro F hF hF'
        obtain ⟨G, -, rfl⟩ := Finset.mem_image.1 hF'
        exact ha ((Finset.mem_powersetCard.1 hF).1 (Finset.mem_insert_self a G))
      have hinj : Set.InjOn (insert a) (s.powersetCard k : Set (Finset ℕ)) := by
        intro F hF G hG h
        have hF' : a ∉ F := fun h' ↦ ha ((Finset.mem_powersetCard.1 hF).1 h')
        have hG' : a ∉ G := fun h' ↦ ha ((Finset.mem_powersetCard.1 hG).1 h')
        rw [← Finset.erase_insert hF', ← Finset.erase_insert hG', h]
      rw [Finset.powersetCard_succ_insert ha, Finset.sum_union hdisj, Finset.sum_image hinj,
        Finset.sum_insert ha]
      have h2 : ∑ F ∈ s.powersetCard k, ∏ x ∈ insert a F, q x =
          q a * ∑ F ∈ s.powersetCard k, ∏ x ∈ F, q x := by
        rw [Finset.mul_sum]
        apply Finset.sum_congr rfl
        intro F hF
        rw [Finset.prod_insert (fun h ↦ ha ((Finset.mem_powersetCard.1 hF).1 h))]
      rw [h2]
      set σ := ∑ x ∈ s, q x
      calc ∑ F ∈ s.powersetCard (k + 1), ∏ x ∈ F, q x + q a * ∑ F ∈ s.powersetCard k, ∏ x ∈ F, q x
          ≤ σ ^ (k + 1) + q a * σ ^ k := by gcongr <;> exact ih _
        _ = (q a + σ) * σ ^ k := by ring
        _ ≤ (q a + σ) * (q a + σ) ^ k := by gcongr; exact le_add_self
        _ = (q a + σ) ^ (k + 1) := by ring

/-- The candidates for the smaller element of a pair summing to `n`. -/
def Q (n : ℕ) : Finset ℕ := (Finset.range n).filter (fun a ↦ 1 ≤ a ∧ 2 * a < n)

/-- The bad event at `n`: at least `2 d + 2` pairs `a < n - a` of the random set sum to `n`. -/
def E (n : ℕ) : Set (Ω d) :=
  {ω | ∃ F ∈ (Q n).powersetCard (2 * d + 2), ∀ a ∈ F, a ∈ S d ω ∧ n - a ∈ S d ω}

/-- The expected number of representations `n = a + (n - a)` with `a < n - a`. -/
noncomputable def σ (n : ℕ) : ℝ≥0 := ∑ a ∈ Q n, p d a * p d (n - a)

theorem measure_E_le (n : ℕ) : P d (E d n) ≤ ((σ d n ^ (2 * d + 2) : ℝ≥0) : ℝ≥0∞) := by
  classical
  have hE : E d n ⊆ ⋃ F ∈ (Q n).powersetCard (2 * d + 2),
      {ω | ∀ x ∈ F ∪ F.image (n - ·), x ∈ S d ω} := by
    rintro ω ⟨F, hF, hFS⟩
    simp only [Set.mem_iUnion]
    refine ⟨F, hF, fun x hx ↦ ?_⟩
    rcases Finset.mem_union.1 hx with h | h
    · exact (hFS x h).1
    · obtain ⟨a, ha, rfl⟩ := Finset.mem_image.1 h
      exact (hFS a ha).2
  have hprod : ∀ F ∈ (Q n).powersetCard (2 * d + 2),
      ∏ x ∈ F ∪ F.image (n - ·), p d x = ∏ a ∈ F, p d a * p d (n - a) := by
    intro F hF
    have hFQ := (Finset.mem_powersetCard.1 hF).1
    have hQ : ∀ a ∈ F, 1 ≤ a ∧ 2 * a < n := fun a ha ↦ (Finset.mem_filter.1 (hFQ ha)).2
    have hdisj : Disjoint F (F.image (n - ·)) := by
      rw [Finset.disjoint_left]
      intro a ha ha'
      obtain ⟨b, hb, hba⟩ := Finset.mem_image.1 ha'
      have := hQ a ha
      have := hQ b hb
      omega
    have hinj : Set.InjOn (n - ·) (F : Set ℕ) := by
      intro a ha b hb h
      have := hQ a ha
      have := hQ b hb
      simp only at h
      omega
    rw [Finset.prod_union hdisj, Finset.prod_image hinj, Finset.prod_mul_distrib]
  calc P d (E d n)
      ≤ P d (⋃ F ∈ (Q n).powersetCard (2 * d + 2), {ω | ∀ x ∈ F ∪ F.image (n - ·), x ∈ S d ω}) :=
        measure_mono hE
    _ ≤ ∑ F ∈ (Q n).powersetCard (2 * d + 2), P d {ω | ∀ x ∈ F ∪ F.image (n - ·), x ∈ S d ω} :=
        measure_biUnion_finset_le _ _
    _ ≤ ∑ F ∈ (Q n).powersetCard (2 * d + 2), ((∏ a ∈ F, p d a * p d (n - a) : ℝ≥0) : ℝ≥0∞) := by
        apply Finset.sum_le_sum
        intro F hF
        refine (measure_subset_le d _).trans (le_of_eq ?_)
        rw [← hprod F hF, ENNReal.ofNNReal_finsetProd]
    _ = ((∑ F ∈ (Q n).powersetCard (2 * d + 2), ∏ a ∈ F, p d a * p d (n - a) : ℝ≥0) : ℝ≥0∞) :=
        (ENNReal.ofNNReal_finsetSum _ _).symm
    _ ≤ _ := ENNReal.coe_le_coe.2 (sum_powersetCard_prod_le _ _ _)

/-- Outside the bad event at `n`, `sumRep` at `n` is at most `2 (2 d + 1) + 1`. -/
theorem sumRep_le_of_not_mem (ω : Ω d) (n : ℕ) (h : ω ∉ E d n) :
    sumRep (S d ω) n ≤ 2 * (2 * d + 1) + 1 := by
  classical
  set T := (Q n).filter (fun a ↦ a ∈ S d ω ∧ n - a ∈ S d ω)
  have hT : T.card ≤ 2 * d + 1 := by
    by_contra hc
    push Not at hc
    obtain ⟨F, hFT, hF⟩ := Finset.exists_subset_card_eq (show 2 * d + 2 ≤ T.card by omega)
    exact h ⟨F, Finset.mem_powersetCard.2 ⟨hFT.trans (Finset.filter_subset _ _), hF⟩,
      fun a ha ↦ (Finset.mem_filter.1 (hFT ha)).2⟩
  rw [sumRep_def]
  calc _ ≤ (T.image (fun a ↦ (a, n - a)) ∪ T.image (fun a ↦ (n - a, a)) ∪ {(n / 2, n / 2)}).card := by
        apply Finset.card_le_card
        intro q hq
        rw [Finset.mem_filter, Finset.HasAntidiagonal.mem_antidiagonal] at hq
        obtain ⟨hsum, h1, h2⟩ := hq
        have h1' : 1 ≤ q.1 := h1.1
        have h2' : 1 ≤ q.2 := h2.1
        simp only [Finset.mem_union, Finset.mem_image, Finset.mem_singleton]
        by_cases ha : 2 * q.1 < n
        · left; left
          refine ⟨q.1, ?_, Prod.ext rfl (by simp only; omega)⟩
          refine Finset.mem_filter.2 ⟨Finset.mem_filter.2 ⟨Finset.mem_range.2 (by omega), h1', ha⟩, h1, ?_⟩
          rwa [show n - q.1 = q.2 by omega]
        · by_cases hb : 2 * q.2 < n
          · left; right
            refine ⟨q.2, ?_, Prod.ext (by simp only; omega) rfl⟩
            refine Finset.mem_filter.2 ⟨Finset.mem_filter.2 ⟨Finset.mem_range.2 (by omega), h2', hb⟩, h2, ?_⟩
            rwa [show n - q.2 = q.1 by omega]
          · right
            exact Prod.ext (by simp only; omega) (by simp only; omega)
    _ ≤ T.card + T.card + 1 := by
        refine (Finset.card_union_le _ _).trans ?_
        exact add_le_add ((Finset.card_union_le _ _).trans
          (add_le_add Finset.card_image_le Finset.card_image_le)) (Finset.card_singleton _).le
    _ ≤ _ := by omega


theorem level_eq {x J : ℕ} (h1 : base d ^ J ≤ x) (h2 : x < base d ^ (J + 1)) : level d x = J :=
  Nat.log_eq_of_pow_le_of_lt_pow h1 h2

theorem p_le_one (x : ℕ) : p d x ≤ 1 := by
  rw [p]
  apply inv_le_one_of_one_le₀
  exact_mod_cast Nat.one_le_iff_ne_zero.2 (NeZero.ne _)

theorem p_le_of_le_level {x J : ℕ} (h : J ≤ level d x) : p d x ≤ ((len d J : ℝ≥0))⁻¹ := by
  rw [p]
  apply inv_anti₀ (by exact_mod_cast Nat.pos_of_ne_zero (NeZero.ne _))
  unfold len
  exact_mod_cast Nat.pow_le_pow_right (by norm_num) (Nat.mul_le_mul_left _ h)

theorem base_pow (J : ℕ) : base d ^ J = 2 ^ (d * J) * len d J := by
  unfold base len
  rw [← pow_mul, ← pow_add]
  congr 1
  ring

theorem base_pow_succ (J : ℕ) : base d ^ (J + 1) = base d * 2 ^ (d * J) * len d J := by
  rw [pow_succ, base_pow]
  ring

/-- The expected number of points below `base^(J+1)` is at most `2 base 2^(d J)`. -/
theorem sum_p_le (hd : 1 ≤ d) (J : ℕ) :
    ∑ a ∈ Finset.Ico 1 (base d ^ (J + 1)), p d a ≤ ((2 * base d * 2 ^ (d * J) : ℕ) : ℝ≥0) := by
  have hb := one_lt_base d
  induction J with
  | zero =>
    calc _ ≤ ∑ a ∈ Finset.Ico 1 (base d ^ (0 + 1)), (1 : ℝ≥0) :=
          Finset.sum_le_sum fun a _ ↦ p_le_one d a
      _ = ((base d - 1 : ℕ) : ℝ≥0) := by simp
      _ ≤ _ := by
          gcongr
          simp only [mul_zero, pow_zero, mul_one]
          omega
  | succ J ih =>
    rw [← Finset.sum_Ico_consecutive _ (show 1 ≤ base d ^ (J + 1) from Nat.one_le_pow _ _ (by omega))
      (show base d ^ (J + 1) ≤ base d ^ (J + 1 + 1) from Nat.pow_le_pow_right (by omega) (by omega))]
    have hblock : ∑ a ∈ Finset.Ico (base d ^ (J + 1)) (base d ^ (J + 1 + 1)), p d a ≤
        ((base d * 2 ^ (d * (J + 1)) : ℕ) : ℝ≥0) := by
      calc _ = ∑ a ∈ Finset.Ico (base d ^ (J + 1)) (base d ^ (J + 1 + 1)), ((len d (J + 1) : ℝ≥0))⁻¹ := by
            apply Finset.sum_congr rfl
            intro a ha
            rw [Finset.mem_Ico] at ha
            rw [p, level_eq d ha.1 ha.2]
        _ ≤ ((base d ^ (J + 1 + 1) : ℕ) : ℝ≥0) * ((len d (J + 1) : ℝ≥0))⁻¹ := by
            rw [Finset.sum_const, nsmul_eq_mul]
            gcongr
            exact (Nat.card_Ico _ _).le.trans (Nat.sub_le _ _)
        _ = _ := by
            rw [base_pow_succ d (J + 1),
              show ((base d * 2 ^ (d * (J + 1)) * len d (J + 1) : ℕ) : ℝ≥0) =
                ((base d * 2 ^ (d * (J + 1)) : ℕ) : ℝ≥0) * (len d (J + 1) : ℝ≥0) from Nat.cast_mul _ _,
              mul_inv_cancel_right₀ (by exact_mod_cast NeZero.ne _)]
    calc _ ≤ ((2 * base d * 2 ^ (d * J) : ℕ) : ℝ≥0) + ((base d * 2 ^ (d * (J + 1)) : ℕ) : ℝ≥0) :=
          add_le_add (ih) hblock
      _ ≤ ((2 * base d * 2 ^ (d * (J + 1)) : ℕ) : ℝ≥0) := by
          rw [← Nat.cast_add]
          gcongr
          have h2d : 2 ≤ 2 ^ d := by
            calc 2 = 2 ^ 1 := rfl
              _ ≤ 2 ^ d := Nat.pow_le_pow_right (by norm_num) hd
          rw [mul_add, mul_one, pow_add]
          nlinarith [Nat.zero_le (base d * 2 ^ (d * J))]

/-- For `n` in level `J + 1`, `σ n ≤ 2 base 2^d / 4^J`. -/
theorem σ_le (hd : 1 ≤ d) {n J : ℕ} (h1 : base d ^ (J + 1) ≤ n) (h2 : n < base d ^ (J + 2)) :
    σ d n ≤ ((2 * base d * 2 ^ d : ℕ) : ℝ≥0) / ((4 ^ J : ℕ) : ℝ≥0) := by
  have hb := one_lt_base d
  have h2X : 2 * base d ^ J ≤ base d ^ (J + 1) := by
    rw [pow_succ, mul_comm]
    exact Nat.mul_le_mul_left _ hb
  have hp : ∀ a ∈ Q n, p d (n - a) ≤ ((len d J : ℝ≥0))⁻¹ := by
    intro a ha
    have := (Finset.mem_filter.1 ha).2
    apply p_le_of_le_level
    exact Nat.le_log_of_pow_le hb (by omega)
  have hQ : Q n ⊆ Finset.Ico 1 (base d ^ (J + 2)) := by
    intro a ha
    have := Finset.mem_filter.1 ha
    rw [Finset.mem_range] at this
    rw [Finset.mem_Ico]
    omega
  calc σ d n ≤ ∑ a ∈ Q n, p d a * ((len d J : ℝ≥0))⁻¹ :=
        Finset.sum_le_sum fun a ha ↦ by gcongr; exact hp a ha
    _ = (∑ a ∈ Q n, p d a) * ((len d J : ℝ≥0))⁻¹ := (Finset.sum_mul _ _ _).symm
    _ ≤ (∑ a ∈ Finset.Ico 1 (base d ^ (J + 2)), p d a) * ((len d J : ℝ≥0))⁻¹ := by
        exact mul_le_mul_of_nonneg_right (Finset.sum_le_sum_of_subset hQ) zero_le
    _ ≤ ((2 * base d * 2 ^ (d * (J + 1)) : ℕ) : ℝ≥0) * ((len d J : ℝ≥0))⁻¹ := by
        gcongr
        exact sum_p_le d hd (J + 1)
    _ = ((2 * base d * 2 ^ d : ℕ) : ℝ≥0) / ((4 ^ J : ℕ) : ℝ≥0) := by
        have hL : ((len d J : ℕ) : ℝ≥0) ≠ 0 := by exact_mod_cast NeZero.ne _
        have h4 : ((4 ^ J : ℕ) : ℝ≥0) ≠ 0 := by positivity
        rw [mul_inv_eq_iff_eq_mul₀ hL, div_mul_eq_mul_div, eq_div_iff h4]
        have : 2 * base d * 2 ^ (d * (J + 1)) * 4 ^ J = 2 * base d * 2 ^ d * len d J := by
          unfold len
          rw [show (4 : ℕ) = 2 ^ 2 from rfl, ← pow_mul]
          ring
        exact_mod_cast this

/-- The bad events of level `J + 1` have total probability at most `base^2 (2 base 2^d)^(2d+2) / base^J`. -/
theorem level_sum_le (hd : 1 ≤ d) (J : ℕ) :
    ∑ n ∈ Finset.Ico (base d ^ (J + 1)) (base d ^ (J + 2)), P d (E d n) ≤
      (((base d ^ 2 * (2 * base d * 2 ^ d) ^ (2 * d + 2) : ℕ) : ℝ≥0) / ((base d ^ J : ℕ) : ℝ≥0) : ℝ≥0) := by
  set K : ℕ := 2 * base d * 2 ^ d
  have hterm : ∀ n ∈ Finset.Ico (base d ^ (J + 1)) (base d ^ (J + 2)),
      σ d n ^ (2 * d + 2) ≤ (((K : ℝ≥0) / ((4 ^ J : ℕ) : ℝ≥0)) ^ (2 * d + 2)) := by
    intro n hn
    rw [Finset.mem_Ico] at hn
    exact pow_le_pow_left₀ zero_le (σ_le d hd hn.1 hn.2) _
  calc ∑ n ∈ Finset.Ico (base d ^ (J + 1)) (base d ^ (J + 2)), P d (E d n)
      ≤ ∑ n ∈ Finset.Ico (base d ^ (J + 1)) (base d ^ (J + 2)), ((σ d n ^ (2 * d + 2) : ℝ≥0) : ℝ≥0∞) :=
        Finset.sum_le_sum fun n _ ↦ measure_E_le d n
    _ = ((∑ n ∈ Finset.Ico (base d ^ (J + 1)) (base d ^ (J + 2)), σ d n ^ (2 * d + 2) : ℝ≥0) : ℝ≥0∞) :=
        (ENNReal.ofNNReal_finsetSum _ _).symm
    _ ≤ _ := by
        apply ENNReal.coe_le_coe.2
        calc _ ≤ ∑ n ∈ Finset.Ico (base d ^ (J + 1)) (base d ^ (J + 2)),
              (((K : ℝ≥0) / ((4 ^ J : ℕ) : ℝ≥0)) ^ (2 * d + 2)) := Finset.sum_le_sum hterm
          _ ≤ ((base d ^ (J + 2) : ℕ) : ℝ≥0) * (((K : ℝ≥0) / ((4 ^ J : ℕ) : ℝ≥0)) ^ (2 * d + 2)) := by
              rw [Finset.sum_const, nsmul_eq_mul]
              gcongr
              exact (Nat.card_Ico _ _).le.trans (Nat.sub_le _ _)
          _ = (((base d ^ 2 * K ^ (2 * d + 2) : ℕ) : ℝ≥0) / ((base d ^ J : ℕ) : ℝ≥0) : ℝ≥0) := by
              have h4 : (4 ^ J) ^ (2 * d + 2) = base d ^ (2 * J) := by
                unfold base
                rw [show (4 : ℕ) = 2 ^ 2 from rfl, ← pow_mul, ← pow_mul, ← pow_mul]
                ring_nf
              have hb0 : ((base d ^ J : ℕ) : ℝ≥0) ≠ 0 := by
                have := one_lt_base d
                positivity
              have h40 : (((4 ^ J : ℕ) : ℝ≥0)) ^ (2 * d + 2) ≠ 0 := by positivity
              rw [div_pow, mul_div_assoc', div_eq_div_iff h40 hb0]
              rw [← Nat.cast_pow (4 ^ J), h4]
              push_cast
              ring

/-- Every outcome has exactly one point in each block, so the count is deterministic: `2^(d J)` blocks of level `J` lie below `base^(J+1)`. -/
theorem count_ge (ω : Ω d) (J N : ℕ) (hN : base d ^ (J + 1) ≤ N) :
    2 ^ (d * J) ≤ (S d ω ∩ Icc 1 N).ncard := by
  classical
  have hb := one_lt_base d
  have hL : 0 < len d J := Nat.pos_of_ne_zero (NeZero.ne _)
  let f : ℕ → ℕ := fun c ↦ c * len d J + ((ω (J, c) : Fin _) : ℕ)
  have hlt : ∀ c, ((ω (J, c) : Fin _) : ℕ) < len d J := fun c ↦ (ω (J, c)).isLt
  have hdiv : ∀ c, f c / len d J = c := by
    intro c
    simp only [f]
    rw [add_comm, Nat.add_mul_div_right _ _ hL, Nat.div_eq_of_lt (hlt c), zero_add]
  have hmod : ∀ c, f c % len d J = ((ω (J, c) : Fin _) : ℕ) := by
    intro c
    simp only [f]
    rw [add_comm, Nat.add_mul_mod_self_right, Nat.mod_eq_of_lt (hlt c)]
  have hmem : ∀ c ∈ Finset.Ico (2 ^ (d * J)) (2 * 2 ^ (d * J)), f c ∈ S d ω ∩ Icc 1 N := by
    intro c hc
    rw [Finset.mem_Ico] at hc
    have h1 : base d ^ J ≤ f c := by
      rw [base_pow]
      calc 2 ^ (d * J) * len d J ≤ c * len d J := Nat.mul_le_mul_right _ hc.1
        _ ≤ f c := Nat.le_add_right _ _
    have h2 : f c < base d ^ (J + 1) := by
      rw [base_pow_succ]
      calc f c < (c + 1) * len d J := by
            simp only [f]
            rw [add_mul, one_mul]
            exact Nat.add_lt_add_left (hlt c) _
        _ ≤ (2 * 2 ^ (d * J)) * len d J := Nat.mul_le_mul_right _ hc.2
        _ ≤ base d * 2 ^ (d * J) * len d J := Nat.mul_le_mul_right _ (Nat.mul_le_mul_right _ hb)
    have hlev : level d (f c) = J := level_eq d h1 h2
    have hkey : key d (f c) = (J, c) := by
      simp only [key, hlev, hdiv]
    have h1' : 1 ≤ f c := le_trans (Nat.one_le_pow _ _ (by omega)) h1
    refine ⟨⟨h1', ?_⟩, h1', (h2.trans_le hN).le⟩
    have := congrArg (fun i ↦ ((ω i : Fin _) : ℕ)) hkey
    simp only at this
    rw [hlev, hmod]
    exact this
  have hinj : Set.InjOn f (Finset.Ico (2 ^ (d * J)) (2 * 2 ^ (d * J)) : Set ℕ) := by
    intro c _ c' _ h
    rw [← hdiv c, ← hdiv c', h]
  calc 2 ^ (d * J) = (Finset.Ico (2 ^ (d * J)) (2 * 2 ^ (d * J))).card := by
        rw [Nat.card_Ico]
        omega
    _ = ((Finset.Ico (2 ^ (d * J)) (2 * 2 ^ (d * J))).image f).card :=
        (Finset.card_image_of_injOn hinj).symm
    _ = ((((Finset.Ico (2 ^ (d * J)) (2 * 2 ^ (d * J))).image f : Finset ℕ) : Set ℕ)).ncard :=
        (Set.ncard_coe_finset _).symm
    _ ≤ (S d ω ∩ Icc 1 N).ncard := by
        apply Set.ncard_le_ncard _ ((Set.finite_Icc 1 N).inter_of_right _)
        intro x hx
        rw [Finset.coe_image] at hx
        obtain ⟨c, hc, rfl⟩ := hx
        exact hmem c hc

/-- Some outcome has `sumRep` bounded: the bad events from level `J₀ + 1` on have total probability below `1`. -/
theorem exists_good (hd : 1 ≤ d) : ∃ ω : Ω d, ∃ C, ∀ n, sumRep (S d ω) n ≤ C := by
  have hb := one_lt_base d
  set K₂ : ℕ := base d ^ 2 * (2 * base d * 2 ^ d) ^ (2 * d + 2) with hK₂
  set J₀ : ℕ := 4 * K₂ with hJ₀def
  set BAD : Set (Ω d) :=
    ⋃ J, ⋃ n ∈ Finset.Ico (base d ^ (J + J₀ + 1)) (base d ^ (J + J₀ + 2)), E d n with hBAD_def
  have hJ₀ : 4 * K₂ ≤ base d ^ J₀ :=
    calc 4 * K₂ = J₀ := rfl
      _ ≤ 2 ^ J₀ := Nat.lt_two_pow_self.le
      _ ≤ base d ^ J₀ := Nat.pow_le_pow_left hb J₀
  set x : ℝ≥0 := ((K₂ : ℕ) : ℝ≥0) / ((base d ^ J₀ : ℕ) : ℝ≥0) with hx
  have hx2 : x * 2 < 1 := by
    have hpos : (0 : ℝ≥0) < ((base d ^ J₀ : ℕ) : ℝ≥0) := by positivity
    have h4 : x * 4 ≤ 1 := by
      rw [hx, div_mul_eq_mul_div, div_le_one hpos]
      exact_mod_cast (by omega : K₂ * 4 ≤ base d ^ J₀)
    have h4' : (x : ℝ) * 4 ≤ 1 := by exact_mod_cast h4
    have : (x : ℝ) * 2 < 1 := by nlinarith [x.2]
    exact_mod_cast this
  have hterm : ∀ J, ((((K₂ : ℕ) : ℝ≥0) / ((base d ^ (J + J₀) : ℕ) : ℝ≥0) : ℝ≥0) : ℝ≥0∞) ≤
      (x : ℝ≥0∞) * 2⁻¹ ^ J := by
    intro J
    have h2 : (2⁻¹ : ℝ≥0∞) = ((2⁻¹ : ℝ≥0) : ℝ≥0∞) := by
      rw [ENNReal.coe_inv two_ne_zero, ENNReal.coe_ofNat]
    rw [h2, ← ENNReal.coe_pow, ← ENNReal.coe_mul]
    apply ENNReal.coe_le_coe.2
    rw [hx, inv_pow, ← div_eq_mul_inv, div_div]
    apply div_le_div_of_nonneg_left zero_le (by positivity)
    have : base d ^ J₀ * 2 ^ J ≤ base d ^ (J + J₀) := by
      rw [pow_add, mul_comm]
      exact Nat.mul_le_mul_right _ (Nat.pow_le_pow_left hb J)
    exact_mod_cast this
  have hBAD : P d BAD < 1 := by
    calc P d BAD
        ≤ ∑' J, ∑ n ∈ Finset.Ico (base d ^ (J + J₀ + 1)) (base d ^ (J + J₀ + 2)), P d (E d n) :=
          (measure_iUnion_le _).trans (ENNReal.tsum_le_tsum fun J ↦ measure_biUnion_finset_le _ _)
      _ ≤ ∑' J, (x : ℝ≥0∞) * 2⁻¹ ^ J :=
          ENNReal.tsum_le_tsum fun J ↦ (level_sum_le d hd (J + J₀)).trans (hterm J)
      _ = (x : ℝ≥0∞) * 2 := by
          rw [ENNReal.tsum_mul_left, ENNReal.tsum_geometric, ENNReal.one_sub_inv_two, inv_inv]
      _ < 1 := by
          rw [← ENNReal.coe_ofNat, ← ENNReal.coe_mul]
          exact_mod_cast hx2
  obtain ⟨ω, hω⟩ : ∃ ω, ω ∉ BAD := by
    by_contra h
    push Not at h
    have : BAD = Set.univ := Set.eq_univ_of_forall h
    rw [this, measure_univ] at hBAD
    exact lt_irrefl _ hBAD
  refine ⟨ω, max (base d ^ (J₀ + 1)) (2 * (2 * d + 1) + 1), fun n ↦ ?_⟩
  by_cases hn : n < base d ^ (J₀ + 1)
  · exact (sumRep_le _ n).trans (le_max_of_le_left (by omega))
  · push Not at hn
    have hn0 : n ≠ 0 := by
      have := Nat.one_le_pow (J₀ + 1) (base d) (by omega)
      omega
    have hL1 : base d ^ Nat.log (base d) n ≤ n := Nat.pow_log_le_self _ hn0
    have hL2 : n < base d ^ (Nat.log (base d) n + 1) := Nat.lt_pow_succ_log_self hb n
    have hL : J₀ + 1 ≤ Nat.log (base d) n := Nat.le_log_of_pow_le hb hn
    have hE : ω ∉ E d n := by
      intro hE
      apply hω
      rw [hBAD_def]
      refine Set.mem_iUnion.2 ⟨Nat.log (base d) n - J₀ - 1, Set.mem_iUnion₂.2 ⟨n, Finset.mem_Ico.2 ⟨?_, ?_⟩, hE⟩⟩
      · rwa [show Nat.log (base d) n - J₀ - 1 + J₀ + 1 = Nat.log (base d) n by omega]
      · rwa [show Nat.log (base d) n - J₀ - 1 + J₀ + 2 = Nat.log (base d) n + 1 by omega]
    exact (sumRep_le_of_not_mem d ω n hE).trans (le_max_right _ _)

/-- The deterministic lower bound on the count: `N^(d/(2d+2)) ≤ 4^d |S ω ∩ [1, N]|` for `N ≥ base`. -/
theorem rpow_le_count (hd : 1 ≤ d) (ω : Ω d) {N : ℕ} (hN : base d ≤ N) :
    (N : ℝ) ^ ((d : ℝ) / (2 * d + 2)) ≤ 4 ^ d * ((S d ω ∩ Icc 1 N).ncard : ℝ) := by
  have hb := one_lt_base d
  have hN0 : N ≠ 0 := by omega
  have hL2 : N < base d ^ (Nat.log (base d) N + 1) := Nat.lt_pow_succ_log_self hb N
  have hL : 1 ≤ Nat.log (base d) N := Nat.le_log_of_pow_le hb (by simpa using hN)
  obtain ⟨M, hM⟩ : ∃ M, Nat.log (base d) N = M + 1 := ⟨Nat.log (base d) N - 1, by omega⟩
  have hL1 : base d ^ (M + 1) ≤ N := hM ▸ Nat.pow_log_le_self _ hN0
  rw [hM] at hL2
  have hcount := count_ge d ω M N hL1
  have hr : (0 : ℝ) ≤ (d : ℝ) / (2 * d + 2) := by positivity
  calc (N : ℝ) ^ ((d : ℝ) / (2 * d + 2))
      ≤ ((base d ^ (M + 1 + 1) : ℕ) : ℝ) ^ ((d : ℝ) / (2 * d + 2)) :=
        Real.rpow_le_rpow (Nat.cast_nonneg _) (by exact_mod_cast hL2.le) hr
    _ = (2 : ℝ) ^ (d * (M + 2)) := by
        rw [show base d ^ (M + 1 + 1) = 2 ^ ((2 * d + 2) * (M + 2)) by unfold base; rw [← pow_mul]]
        rw [Nat.cast_pow, Nat.cast_ofNat, ← Real.rpow_natCast (2 : ℝ) ((2 * d + 2) * (M + 2)),
          ← Real.rpow_mul (by norm_num), ← Real.rpow_natCast (2 : ℝ) (d * (M + 2))]
        congr 1
        push_cast
        field_simp
    _ = 4 ^ d * (2 ^ (d * M) : ℝ) := by
        rw [show (4 : ℝ) = 2 ^ 2 by norm_num, ← pow_mul, ← pow_add]
        congr 1
        ring
    _ ≤ 4 ^ d * ((S d ω ∩ Icc 1 N).ncard : ℝ) := by
        gcongr
        exact_mod_cast hcount

/-- Theorem B inside the construction, with `d = ⌈1/ε⌉₊ + 1`. -/
theorem exists_dense {ε : ℝ} (hε : 0 < ε) : ∃ A : Set ℕ, (∃ C, ∀ n, sumRep A n ≤ C) ∧
    (fun N : ℕ ↦ √(N : ℝ) / (N : ℝ) ^ ε) =O[atTop] (fun N ↦ ((A ∩ Icc 1 N).ncard : ℝ)) := by
  set d := ⌈1 / ε⌉₊ + 1 with hd_def
  have hd : 1 ≤ d := by omega
  have hεd : 1 / (2 * (d : ℝ) + 2) ≤ ε := by
    have h1 : 1 / ε ≤ (d : ℝ) := by
      have := Nat.le_ceil (1 / ε)
      rw [hd_def]
      push_cast
      linarith
    have h2 : 1 ≤ ε * d := by
      rw [div_le_iff₀ hε] at h1
      linarith
    rw [div_le_iff₀ (by positivity)]
    nlinarith
  obtain ⟨ω, C, hC⟩ := exists_good d hd
  refine ⟨S d ω, ⟨C, hC⟩, IsBigO.of_bound (4 ^ d) ?_⟩
  filter_upwards [eventually_ge_atTop (base d)] with N hN
  have hb := one_lt_base d
  have hN1 : (1 : ℝ) ≤ N := by exact_mod_cast (by omega : 1 ≤ N)
  have hNpos : (0 : ℝ) < N := by linarith
  rw [Real.norm_of_nonneg (by positivity), Real.norm_of_nonneg (by positivity)]
  calc √(N : ℝ) / (N : ℝ) ^ ε = (N : ℝ) ^ ((1 : ℝ) / 2 - ε) := by
        rw [Real.sqrt_eq_rpow, Real.rpow_sub hNpos]
    _ ≤ (N : ℝ) ^ ((d : ℝ) / (2 * d + 2)) := by
        apply Real.rpow_le_rpow_of_exponent_le hN1
        have : (d : ℝ) / (2 * d + 2) = 1 / 2 - 1 / (2 * d + 2) := by
          field_simp
          ring
        linarith
    _ ≤ 4 ^ d * ((S d ω ∩ Icc 1 N).ncard : ℝ) := rpow_le_count d hd ω hN

end TheoremB

/-- **Theorem B.** For every `ε > 0` some `A ⊆ ℕ` has bounded `sumRep` and `|A ∩ [1, N]| ≫ √N / N^ε`. -/
theorem exists_dense_bounded {ε : ℝ} (hε : 0 < ε) : ∃ A : Set ℕ, (∃ C, ∀ n, sumRep A n ≤ C) ∧
    (fun N : ℕ ↦ √(N : ℝ) / (N : ℝ) ^ ε) =O[atTop] (fun N ↦ ((A ∩ Icc 1 N).ncard : ℝ)) :=
  TheoremB.exists_dense hε

/-- No power `N^ε`, `ε > 0`, has the property of problem 40. -/
theorem not_erdos40For_rpow {ε : ℝ} (hε : 0 < ε) : ¬ Erdos40For (fun N ↦ (N : ℝ) ^ ε) := by
  intro h
  obtain ⟨A, ⟨C, hC⟩, hO⟩ := exists_dense_bounded hε
  exact limsup_ne_top_of_bound hC (h A hO)

/-- No power `N^ε`, `ε > 0`, lies in the set of `Erdos40.erdos_40`. -/
theorem rpow_not_mem_answerSet : ∀ ε : ℝ, 0 < ε →
    (fun N : ℕ ↦ (N : ℝ) ^ ε) ∉ {g : ℕ → ℝ | Tendsto g atTop atTop ∧ Erdos40For g} :=
  fun _ hε hg ↦ not_erdos40For_rpow hε hg.2

/-- The corollary: no `g` with `N^ε = O(g)`, that is with `g(N) ≥ c N^ε` for large `N`, has the property. -/
theorem not_erdos40For_of_rpow_isBigO {g : ℕ → ℝ} {ε : ℝ} (hε : 0 < ε)
    (hg : (fun N : ℕ ↦ (N : ℝ) ^ ε) =O[atTop] g) : ¬ Erdos40For g := fun h ↦
  not_erdos40For_rpow hε (Erdos40For.of_isBigO h (by
    filter_upwards [eventually_ge_atTop 1] with N hN
    exact (Real.rpow_pos_of_pos (by exact_mod_cast hN) _).ne') hg)

end ErdosProblem40
