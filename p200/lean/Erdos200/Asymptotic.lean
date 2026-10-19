/-
Copyright (c) 2026 Max Marquardt.
Released under the MIT licence; see LICENSE at the repository root.
-/
import Erdos200.Bound
import FormalConjectures.ErdosProblems.«200»

/-!
# The asymptotic upper bound for Erdős problem 200

Write `L n` for `Erdos200.longestPrimeArithmeticProgressions n`, the length of the longest
arithmetic progression of primes inside `[1, n]`.  Proposition 2 gives `θ(L n) < log n + log 2`,
so any lower bound `θ(k) ≥ (c - o(1)) k` turns into `L n ≤ (1/c + o(1)) log n`.

* Unconditionally, Mathlib's `Chebyshev.theta_ge` supplies `c = log 2`, giving
  `L n ≤ (1/log 2 + o(1)) log n` with `1/log 2 = 1.4426…`.
* With the prime number theorem `θ(x)/x → 1` it supplies `c = 1`, giving the
  literature form `L n ≤ (1 + o(1)) log n`, which is
  `FormalConjectures`' `Erdos200.erdos_200.variants.upper`.

Mathlib at the revision pinned by this project does **not** contain the prime number theorem,
so `PNT` below is an explicit hypothesis, not a theorem.  Nothing here addresses the open
part of Erdős 200, which asks for `L n = o(log n)`.
-/

namespace ErdosProblem200

open Filter Real Erdos200

/- ### Extracting a concrete progression -/

/-- `Set.IsAPOfLength` is an existential statement about a set; this unpacks it into the
concrete data `a, d` used by `Erdos200.Bound`, including `d ≥ 1` (a constant "progression"
`d = 0` is a one-element set, so it cannot have length `≥ 3`). -/
theorem exists_of_isAPOfLength {n k : ℕ} (hk : 3 ≤ k) {s : Set ℕ}
    (hsub : s ⊆ Set.Icc 1 n) (hAP : s.IsAPOfLength (k : ℕ∞)) (hprime : ∀ m ∈ s, Nat.Prime m) :
    ∃ a d, 0 < d ∧ (∀ i < k, Nat.Prime (a + i * d)) ∧ a + (k - 1) * d ≤ n := by
  obtain ⟨a, d, hcard, hset⟩ := hAP
  have hmem : ∀ i < k, a + i * d ∈ s := by
    intro i hi
    rw [hset]
    exact ⟨i, by exact_mod_cast hi, by simp⟩
  have hd : 0 < d := by
    rcases Nat.eq_zero_or_pos d with rfl | hd
    · exfalso
      have hs : s = {a} := by
        rw [hset]
        ext x
        simp only [Set.mem_ofPred_eq, Set.mem_singleton_iff, smul_zero, add_zero]
        constructor
        · rintro ⟨i, -, rfl⟩; rfl
        · rintro rfl; exact ⟨0, by exact_mod_cast (show 0 < k by omega), rfl⟩
      rw [hs, show ENat.card ({a} : Set ℕ) = 1 by simp] at hcard
      have : k = 1 := by exact_mod_cast hcard.symm
      omega
    · exact hd
  exact ⟨a, d, hd, fun i hi ↦ hprime _ (hmem i hi), (hsub (hmem (k - 1) (by omega))).2⟩

/-- `L n = 0` for `n ≤ 1`: the only subset of `Set.Icc 1 n` consisting of primes is empty,
since `1` is not prime.  Needed because the `variants.upper` statement quantifies over *all*
`n`, including those where `log n = 0`. -/
theorem longest_eq_zero {n : ℕ} (hn : n ≤ 1) : longestPrimeArithmeticProgressions n = 0 := by
  rw [longestPrimeArithmeticProgressions]
  refine Nat.le_zero.mp (csSup_le' ?_)
  rintro k ⟨s, hsub, hAP, hprime⟩
  have hs : s = ∅ := by
    ext x
    simp only [Set.mem_empty_iff_false, iff_false]
    intro hx
    have hx' := hsub hx
    rw [Set.mem_Icc] at hx'
    exact Nat.not_prime_one ((show x = 1 by omega) ▸ hprime x hx)
  obtain ⟨a, d, hcard, -⟩ := hAP
  rw [hs, show ENat.card (∅ : Set ℕ) = 0 by simp] at hcard
  have : k = 0 := by exact_mod_cast hcard.symm
  omega

/- ### The finite-`n` bound -/

/-- If `c * k ≤ θ(k)` for all `k ≥ K`, then every admissible length is at most
`max K ⌊(log n + log 2)/c⌋₊`.  This is the only place the supremum is touched. -/
theorem longest_le_max {n K : ℕ} {c : ℝ} (hc : 0 < c) (hK3 : 3 ≤ K)
    (hK : ∀ k, K ≤ k → c * k ≤ Chebyshev.theta k) :
    longestPrimeArithmeticProgressions n ≤ max K ⌊(Real.log n + Real.log 2) / c⌋₊ := by
  rw [longestPrimeArithmeticProgressions]
  refine csSup_le' ?_
  rintro k ⟨s, hsub, hAP, hprime⟩
  by_cases hkK : k ≤ K
  · exact le_max_of_le_left hkK
  · obtain ⟨a, d, hd, hp, hle⟩ := exists_of_isAPOfLength (by omega : 3 ≤ k) hsub hAP hprime
    have h1 : Chebyshev.theta k < Real.log n + Real.log 2 := theta_lt (by omega) hd hp hle
    have h2 : c * k ≤ Chebyshev.theta k := hK k (by omega)
    refine le_max_of_le_right (Nat.le_floor ?_)
    rw [le_div_iff₀ hc]
    linarith

/- ### Elementary estimates for the error terms in `Chebyshev.theta_ge` -/

/-- `log x ≤ 2√x` for `x ≥ 1`. -/
theorem log_le_two_mul_sqrt {x : ℝ} (hx : 1 ≤ x) : Real.log x ≤ 2 * Real.sqrt x := by
  have hx0 : (0:ℝ) < x := by linarith
  have h := Real.log_le_sub_one_of_pos (Real.sqrt_pos.mpr hx0)
  rw [Real.log_sqrt hx0.le] at h
  linarith

/-- The error terms in `Chebyshev.theta_ge` are `o(k)`.  Stated with an explicit threshold
`k ≥ (11/δ)^4` rather than through `IsLittleO`; the proof runs on `u = k^{1/4}`, for which
`log (k+1) + 2√k log k ≤ 3u² + 8u³ ≤ 11u³` and `δ k = δ u⁴ ≥ 11 u³`. -/
theorem eventually_log_add_le {δ : ℝ} (hδ : 0 < δ) :
    ∀ᶠ k : ℕ in atTop,
      Real.log ((k : ℝ) + 1) + 2 * Real.sqrt k * Real.log k ≤ δ * k := by
  filter_upwards [eventually_ge_atTop (max 1 ⌈(11 / δ) ^ 4⌉₊)] with k hk
  have hk1 : 1 ≤ k := le_trans (le_max_left _ _) hk
  have hkR : (1:ℝ) ≤ (k:ℝ) := by exact_mod_cast hk1
  have hk0 : (0:ℝ) ≤ (k:ℝ) := by linarith
  set u : ℝ := Real.sqrt (Real.sqrt k) with hu
  have hu0 : 0 ≤ u := Real.sqrt_nonneg _
  have hu2 : u ^ 2 = Real.sqrt k := Real.sq_sqrt (Real.sqrt_nonneg _)
  have hu4 : u ^ 4 = (k : ℝ) := by
    rw [show u ^ 4 = (u ^ 2) ^ 2 by ring, hu2, Real.sq_sqrt hk0]
  have hu1 : 1 ≤ u := by
    rw [hu, show (1:ℝ) = Real.sqrt (Real.sqrt 1) by simp]
    exact Real.sqrt_le_sqrt (Real.sqrt_le_sqrt hkR)
  have hsk1 : (1:ℝ) ≤ Real.sqrt k := by
    rw [show (1:ℝ) = Real.sqrt 1 by simp]
    exact Real.sqrt_le_sqrt hkR
  -- `log (k+1) ≤ 2√(k+1) ≤ 3u²`
  have h1 : Real.log ((k:ℝ) + 1) ≤ 3 * u ^ 2 := by
    have ha : Real.log ((k:ℝ) + 1) ≤ 2 * Real.sqrt ((k:ℝ) + 1) :=
      log_le_two_mul_sqrt (by linarith)
    have hb : Real.sqrt ((k:ℝ) + 1) ≤ (3/2) * u ^ 2 := by
      rw [← Real.sqrt_sq (by positivity : (0:ℝ) ≤ (3/2) * u ^ 2)]
      exact Real.sqrt_le_sqrt (by nlinarith [hu4, hkR])
    linarith
  -- `log k = 2 log √k ≤ 4u`, hence `2√k log k ≤ 8u³`
  have hlogk : Real.log k ≤ 4 * u := by
    have h := log_le_two_mul_sqrt hsk1
    rw [Real.log_sqrt hk0] at h
    linarith
  have h2 : 2 * Real.sqrt k * Real.log k ≤ 8 * u ^ 3 := by
    have hs0 : (0:ℝ) ≤ 2 * Real.sqrt k := by positivity
    have h := mul_le_mul_of_nonneg_left hlogk hs0
    rw [← hu2] at h ⊢
    nlinarith [h]
  -- `k ≥ (11/δ)^4` means `u ≥ 11/δ`
  have hkge : ((11 / δ) ^ 4 : ℝ) ≤ (k:ℝ) := by
    calc ((11 / δ) ^ 4 : ℝ) ≤ (⌈(11 / δ) ^ 4⌉₊ : ℝ) := Nat.le_ceil _
      _ ≤ (k:ℝ) := by exact_mod_cast le_trans (le_max_right _ _) hk
  have h3 : 11 / δ ≤ u := by
    rw [hu, Real.le_sqrt (le_of_lt (div_pos (by norm_num) hδ)) (Real.sqrt_nonneg _),
      Real.le_sqrt (by positivity) hk0]
    nlinarith [hkge]
  -- assemble: `11u³ ≤ δu⁴ = δk`
  have hkey : 11 * u ^ 3 ≤ δ * u ^ 4 := by
    rw [div_le_iff₀ hδ] at h3
    nlinarith [pow_nonneg hu0 3]
  have husq : u ^ 2 ≤ u ^ 3 := by nlinarith [hu1, hu0]
  have hgoal : δ * (k:ℝ) = δ * u ^ 4 := by rw [hu4]
  linarith

/- ### From a Chebyshev lower bound to an asymptotic bound on `L` -/

/-- **The transfer lemma.**  If `θ(k) ≥ c' k` eventually for every `c' < c`, then
`L n ≤ (1/c + ε) log n` eventually, for every `ε > 0`. -/
theorem eventually_longest_le {c : ℝ} (hc : 0 < c)
    (hθ : ∀ c' : ℝ, c' < c → ∀ᶠ k : ℕ in atTop, c' * k ≤ Chebyshev.theta k)
    {ε : ℝ} (hε : 0 < ε) :
    ∀ᶠ n : ℕ in atTop,
      (longestPrimeArithmeticProgressions n : ℝ) ≤ (1 / c + ε) * Real.log n := by
  have hcne : c ≠ 0 := ne_of_gt hc
  have hinv : (0:ℝ) < 1 / c + ε / 2 := by positivity
  -- `c'` is chosen so that `1/c' = 1/c + ε/2`, leaving `ε/2` of slack for the `log 2`.
  set c' : ℝ := 1 / (1 / c + ε / 2) with hc'def
  have hc'pos : 0 < c' := by rw [hc'def]; positivity
  have hc'inv : 1 / c' = 1 / c + ε / 2 := by rw [hc'def, one_div_one_div]
  have hc'lt : c' < c := by
    rw [hc'def, div_lt_iff₀ hinv, show c * (1 / c + ε / 2) = 1 + c * ε / 2 by field_simp]
    nlinarith
  obtain ⟨K₀, hK₀⟩ := Filter.eventually_atTop.mp (hθ c' hc'lt)
  have hK : ∀ k, max K₀ 3 ≤ k → c' * k ≤ Chebyshev.theta k := fun k hk ↦
    hK₀ k (le_trans (le_max_left _ _) hk)
  have hlog2 : (0:ℝ) ≤ Real.log 2 := Real.log_nonneg (by norm_num)
  have hCε : (0:ℝ) < 1 / c + ε := by positivity
  have hlog : Tendsto (fun n : ℕ ↦ Real.log n) atTop atTop :=
    Real.tendsto_log_atTop.comp tendsto_natCast_atTop_atTop
  filter_upwards [eventually_ge_atTop 1,
    hlog.eventually_ge_atTop ((max K₀ 3 : ℕ) / (1 / c + ε)),
    hlog.eventually_ge_atTop ((1 / c + ε / 2) * Real.log 2 * (2 / ε))] with n hn1 hn2 hn3
  have hlogn : 0 ≤ Real.log n := Real.log_nonneg (by exact_mod_cast hn1)
  have harg : (0:ℝ) ≤ (Real.log n + Real.log 2) / c' :=
    div_nonneg (by linarith) hc'pos.le
  have hcast : (longestPrimeArithmeticProgressions n : ℝ) ≤
      max ((max K₀ 3 : ℕ) : ℝ) ((Real.log n + Real.log 2) / c') := by
    have h := (Nat.cast_le (α := ℝ)).mpr (longest_le_max (n := n) hc'pos (le_max_right _ _) hK)
    rw [Nat.cast_max] at h
    exact h.trans (max_le_max le_rfl (Nat.floor_le harg))
  refine hcast.trans (max_le ?_ ?_)
  · rw [div_le_iff₀ hCε] at hn2
    linarith
  · rw [show (Real.log n + Real.log 2) / c'
        = (1 / c + ε / 2) * (Real.log n + Real.log 2) by rw [div_eq_inv_mul, ← one_div, hc'inv]]
    have h4 : (1 / c + ε / 2) * Real.log 2 ≤ (ε / 2) * Real.log n := by
      calc (1 / c + ε / 2) * Real.log 2
          = (1 / c + ε / 2) * Real.log 2 * (2 / ε) * (ε / 2) := by field_simp
        _ ≤ Real.log n * (ε / 2) := mul_le_mul_of_nonneg_right hn3 (by positivity)
        _ = (ε / 2) * Real.log n := by ring
    nlinarith [h4]

/-- Repackaging an "`∀ ε > 0`, eventually" bound as the `(C + o(1)) log n` form used by
`FormalConjectures`.  The witness is `o n = max 0 (L n / log n - C)`. -/
theorem exists_littleO_of_eventually_le {C : ℝ}
    (h : ∀ ε : ℝ, 0 < ε → ∀ᶠ n : ℕ in atTop,
      (longestPrimeArithmeticProgressions n : ℝ) ≤ (C + ε) * Real.log n) :
    ∃ (o : ℕ → ℝ) (_ : o =o[atTop] (1 : ℕ → ℝ)),
      ∀ n, (longestPrimeArithmeticProgressions n : ℝ) ≤ (C + o n) * Real.log n := by
  classical
  refine ⟨fun n ↦ if 2 ≤ n then
    max 0 ((longestPrimeArithmeticProgressions n : ℝ) / Real.log n - C) else 0, ?_, ?_⟩
  · rw [Asymptotics.isLittleO_iff]
    intro ε hε
    simp only [Pi.one_apply, norm_one, mul_one]
    filter_upwards [h ε hε, eventually_ge_atTop 2] with n hn hn2
    have hlog : 0 < Real.log n := Real.log_pos (by exact_mod_cast hn2)
    rw [if_pos hn2, Real.norm_eq_abs, abs_of_nonneg (le_max_left _ _)]
    refine max_le hε.le ?_
    have : (longestPrimeArithmeticProgressions n : ℝ) / Real.log n ≤ C + ε := by
      rw [div_le_iff₀ hlog]; linarith
    linarith
  · intro n
    by_cases hn2 : 2 ≤ n
    · dsimp only
      rw [if_pos hn2]
      have hlog : 0 < Real.log n := Real.log_pos (by exact_mod_cast hn2)
      set x : ℝ := (longestPrimeArithmeticProgressions n : ℝ) with hx
      have hm : x / Real.log n ≤ C + max 0 (x / Real.log n - C) := by
        linarith [le_max_right (0:ℝ) (x / Real.log n - C)]
      calc x = (x / Real.log n) * Real.log n := by field_simp
        _ ≤ (C + max 0 (x / Real.log n - C)) * Real.log n :=
            mul_le_mul_of_nonneg_right hm hlog.le
    · have hlog0 : Real.log n = 0 := by
        rcases (by omega : n = 0 ∨ n = 1) with rfl | rfl <;> simp
      dsimp only
      rw [if_neg (by omega), longest_eq_zero (by omega), hlog0]
      simp

/- ### The unconditional bound -/

/-- Mathlib's `Chebyshev.theta_ge` gives `θ(k) ≥ c' k` eventually for every `c' < log 2`. -/
theorem theta_lower_log_two (c' : ℝ) (hc' : c' < Real.log 2) :
    ∀ᶠ k : ℕ in atTop, c' * k ≤ Chebyshev.theta k := by
  by_cases h : c' ≤ 0
  · filter_upwards with k
    have hk0 : (0:ℝ) ≤ (k:ℝ) := Nat.cast_nonneg k
    nlinarith [Chebyshev.theta_nonneg (k:ℝ)]
  · filter_upwards [eventually_log_add_le (show (0:ℝ) < Real.log 2 - c' by linarith)] with k hk
    nlinarith [Chebyshev.theta_ge k]

/--
**Unconditional, machine-checked:** `L n ≤ (1/log 2 + o(1)) log n`, i.e. `L n ≲ 1.4427 log n`.

The constant is `1/log 2` rather than the `1` of the literature purely because Mathlib at the
pinned revision proves only Chebyshev's `θ(x) ≥ x log 2 - O(√x log x)` and not the prime
number theorem.  See `erdos_200_variants_upper` for the `(1 + o(1))` form.
-/
theorem upper_unconditional :
    ∃ (o : ℕ → ℝ) (_ : o =o[atTop] (1 : ℕ → ℝ)),
      ∀ n, (longestPrimeArithmeticProgressions n : ℝ) ≤ (1 / Real.log 2 + o n) * Real.log n :=
  exists_littleO_of_eventually_le fun _ hε ↦
    eventually_longest_le (Real.log_pos one_lt_two) theta_lower_log_two hε

/- ### The literature bound, conditional on the prime number theorem -/

/-- The prime number theorem in Chebyshev form, `θ(x)/x → 1`.  Mathlib at the revision pinned
by this project does not prove this, so it appears as an explicit hypothesis. -/
def PNT : Prop := Tendsto (fun x : ℝ ↦ Chebyshev.theta x / x) atTop (nhds 1)

/-- `PNT` gives `θ(k) ≥ c' k` eventually for every `c' < 1`. -/
theorem theta_lower_of_pnt (h : PNT) (c' : ℝ) (hc' : c' < 1) :
    ∀ᶠ k : ℕ in atTop, c' * k ≤ Chebyshev.theta k := by
  filter_upwards [tendsto_natCast_atTop_atTop.eventually (h.eventually_const_lt hc'),
    eventually_ge_atTop 1] with k hk hk1
  rw [lt_div_iff₀ (show (0:ℝ) < k by exact_mod_cast hk1)] at hk
  linarith

/--
**Conditional on the prime number theorem:** `L n ≤ (1 + o(1)) log n`.

This is exactly the statement left as `sorry` at `Erdos200.erdos_200.variants.upper` in
`FormalConjectures`, discharged modulo `PNT`.  It is the *known* half of Erdős 200; the
problem itself asks whether `L n = o(log n)`, which this says nothing about.
-/
theorem erdos_200_variants_upper (h : PNT) :
    ∃ (o : ℕ → ℝ) (_ : o =o[atTop] (1 : ℕ → ℝ)),
      ∀ n, (longestPrimeArithmeticProgressions n : ℝ) ≤ (1 + o n) * Real.log n := by
  obtain ⟨o, ho, hbound⟩ := exists_littleO_of_eventually_le (C := 1 / (1:ℝ))
    fun _ hε ↦ eventually_longest_le one_pos (theta_lower_of_pnt h) hε
  exact ⟨o, ho, fun n ↦ by simpa using hbound n⟩

end ErdosProblem200
