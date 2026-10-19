/-
Copyright (c) 2026 Max Marquardt.
Released under the MIT licence; see LICENSE at the repository root.
-/
import Erdos200.Structure
import Mathlib.NumberTheory.Chebyshev

/-!
# The explicit primorial bound

From the structure theorem we get, for a `k`-term AP of primes inside `[1, N]` with `k ≥ 3`,

* `(k - 1) * primorial k < k * N`,
* hence `primorial k < 2 * N`,
* hence `θ(k) < log N + log 2`, and the sharper `θ(k) < log N + 1/(k-1)`.

This is Proposition 2 of `../../NOTES.md`.  The factor `k/(k-1)` is exactly the loss incurred
by Case B of the structure theorem; Case A alone would give `primorial k ≤ N`.
-/

namespace ErdosProblem200

open Real

/--
**Proposition 2** (integer form).  If `a, a+d, …, a+(k-1)d` are primes with `k ≥ 3`, `d ≥ 1`,
and `a + (k-1)d ≤ N`, then `(k-1) * primorial k < k * N`.
-/
theorem mul_primorial_lt {k a d N : ℕ} (hk : 3 ≤ k) (hd : 0 < d)
    (hprime : ∀ i < k, Nat.Prime (a + i * d)) (hN : a + (k - 1) * d ≤ N) :
    (k - 1) * primorial k < k * N := by
  have ha2 : 2 ≤ a := (by simpa using hprime 0 (by omega) : Nat.Prime a).two_le
  obtain ⟨m, rfl⟩ : ∃ m, k = m + 1 := ⟨k - 1, by omega⟩
  simp only [Nat.add_sub_cancel] at hN ⊢
  have hm : 2 ≤ m := by omega
  rcases structure_of_primeAP hk hd hprime with ⟨hak, hdvd⟩ | ⟨-, hak, hdvd⟩
  · -- Case A: `primorial k ≤ d` and `a ≥ k + 1 ≥ 1`, so `m * P ≤ m * d < N ≤ (m+1) * N`.
    have h1 : primorial (m + 1) ≤ d := Nat.le_of_dvd hd hdvd
    calc m * primorial (m + 1) ≤ m * d := Nat.mul_le_mul_left m h1
      _ < N := by nlinarith
      _ ≤ (m + 1) * N := Nat.le_mul_of_pos_left N (by omega)
  · -- Case B: `primorial k ≤ k * d` and `N ≥ (m+1) + m*d`.
    have h1 : primorial (m + 1) ≤ (m + 1) * d := Nat.le_of_dvd (by positivity) hdvd
    subst hak
    nlinarith

/-- `primorial k < 2 * N`: the crude form of Proposition 2, valid for every `k ≥ 3`. -/
theorem primorial_lt_two_mul {k a d N : ℕ} (hk : 3 ≤ k) (hd : 0 < d)
    (hprime : ∀ i < k, Nat.Prime (a + i * d)) (hN : a + (k - 1) * d ≤ N) :
    primorial k < 2 * N := by
  have key := mul_primorial_lt hk hd hprime hN
  have hP : 0 < primorial k := primorial_pos k
  -- `k ≤ 2 * (k - 1)` for `k ≥ 2`, so `k * P ≤ 2 * (k-1) * P < 2 * k * N`.
  have h2 : k * primorial k ≤ 2 * ((k - 1) * primorial k) := by
    have : k ≤ 2 * (k - 1) := by omega
    calc k * primorial k ≤ (2 * (k - 1)) * primorial k := Nat.mul_le_mul_right _ this
      _ = 2 * ((k - 1) * primorial k) := by ring
  have h3 : k * primorial k < k * (2 * N) := by
    calc k * primorial k ≤ 2 * ((k - 1) * primorial k) := h2
      _ < 2 * (k * N) := by omega
      _ = k * (2 * N) := by ring
  exact Nat.lt_of_mul_lt_mul_left h3

/-- Positivity of `N`, extracted for reuse. -/
private theorem pos_of_le {k a d N : ℕ} (hk : 3 ≤ k)
    (hprime : ∀ i < k, Nat.Prime (a + i * d)) (hN : a + (k - 1) * d ≤ N) : 0 < N := by
  have ha2 : 2 ≤ a := (by simpa using hprime 0 (by omega) : Nat.Prime a).two_le
  have : a ≤ N := le_trans (Nat.le_add_right a _) hN
  omega

/-- `θ(k) < log N + log 2`.  Chebyshev's `θ` of the *length* is bounded by the log of the range. -/
theorem theta_lt {k a d N : ℕ} (hk : 3 ≤ k) (hd : 0 < d)
    (hprime : ∀ i < k, Nat.Prime (a + i * d)) (hN : a + (k - 1) * d ≤ N) :
    Chebyshev.theta k < Real.log N + Real.log 2 := by
  have hNpos : 0 < N := pos_of_le hk hprime hN
  have key := primorial_lt_two_mul hk hd hprime hN
  have hP : (0:ℝ) < (primorial k : ℝ) := by exact_mod_cast primorial_pos k
  have hcast : (primorial k : ℝ) < 2 * N := by exact_mod_cast key
  have h := Real.log_lt_log hP hcast
  rw [Real.log_mul (by norm_num) (by positivity)] at h
  rw [Chebyshev.theta_eq_log_primorial, Nat.floor_natCast]
  linarith

/-- The sharp form `θ(k) < log N + 1/(k-1)`, which is what the `k/(k-1)` factor of
Proposition 2 actually costs. -/
theorem theta_lt_sharp {k a d N : ℕ} (hk : 3 ≤ k) (hd : 0 < d)
    (hprime : ∀ i < k, Nat.Prime (a + i * d)) (hN : a + (k - 1) * d ≤ N) :
    Chebyshev.theta k < Real.log N + 1 / ((k : ℝ) - 1) := by
  have hNpos : 0 < N := pos_of_le hk hprime hN
  have hNR : (0:ℝ) < N := by exact_mod_cast hNpos
  have hkR : (3:ℝ) ≤ k := by exact_mod_cast hk
  have hkm : (0:ℝ) < (k : ℝ) - 1 := by linarith
  have key := mul_primorial_lt hk hd hprime hN
  have key' : ((k : ℝ) - 1) * (primorial k : ℝ) < (k : ℝ) * N := by
    have h := (Nat.cast_lt (α := ℝ)).mpr key
    push_cast [Nat.cast_sub (show 1 ≤ k by omega)] at h
    linarith
  have hP : (0:ℝ) < (primorial k : ℝ) := by exact_mod_cast primorial_pos k
  have h1 : (primorial k : ℝ) < ((k : ℝ) / ((k : ℝ) - 1)) * N := by
    rw [div_mul_eq_mul_div, lt_div_iff₀ hkm]
    linarith
  have h2 := Real.log_lt_log hP h1
  rw [Real.log_mul (by positivity) (by positivity)] at h2
  have h3 : Real.log ((k : ℝ) / ((k : ℝ) - 1)) ≤ 1 / ((k : ℝ) - 1) := by
    have hpos : (0:ℝ) < (k : ℝ) / ((k : ℝ) - 1) := div_pos (by linarith) hkm
    have h := Real.log_le_sub_one_of_pos hpos
    have heq : (k : ℝ) / ((k : ℝ) - 1) - 1 = 1 / ((k : ℝ) - 1) := by field_simp; ring
    linarith
  rw [Chebyshev.theta_eq_log_primorial, Nat.floor_natCast]
  linarith

end ErdosProblem200
