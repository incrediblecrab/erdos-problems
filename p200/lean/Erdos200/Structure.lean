/-
Copyright (c) 2026 Max Marquardt.
Released under the MIT licence; see LICENSE at the repository root.
-/
import Mathlib.NumberTheory.Primorial
import Mathlib.Data.ZMod.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Algebra.Field.ZMod

/-!
# The structure of an arithmetic progression of primes

If `a, a + d, …, a + (k-1)d` are all prime with `k ≥ 3` and `d ≥ 1`, then exactly one of

* **(A)** `k < a` and `primorial k ∣ d`;
* **(B)** `k` is prime, `a = k`, and `(primorial k / k) ∣ d`.

Case (B) is not vacuous: `3, 5, 7` and `5, 11, 17, 23, 29` are of that shape, and they are
the extremal examples for `k = 3` and `k = 5`.  It is also the case that binds in
`Erdos200.Bound`, so dropping it does not merely lose sharpness — it makes the bound false.

This is Proposition 1 of `../../NOTES.md`.
-/

namespace ErdosProblem200

open Finset

/-- If `p` is prime and `p ∤ d`, then `p` divides one of the first `p` terms of the
progression `a, a + d, a + 2d, …`; the index is the representative of `-a/d` in `ZMod p`. -/
theorem exists_lt_dvd {p : ℕ} (hp : p.Prime) (a d : ℕ) (hpd : ¬p ∣ d) :
    ∃ i < p, p ∣ a + i * d := by
  have : Fact p.Prime := ⟨hp⟩
  have hd : (d : ZMod p) ≠ 0 := fun h ↦ hpd ((ZMod.natCast_eq_zero_iff d p).mp h)
  refine ⟨(-(a : ZMod p) / (d : ZMod p)).val, ZMod.val_lt _, ?_⟩
  rw [← ZMod.natCast_eq_zero_iff]
  push_cast [ZMod.natCast_val, ZMod.cast_id]
  field_simp
  ring

/-- The basic constraint on a progression of primes: a prime `p ≤ k` that does not divide the
common difference `d` divides some term, and a prime divisible by `p` *equals* `p`.  Hence the
progression starts at or below `p`. -/
theorem le_of_not_dvd {k a d p : ℕ} (hp : p.Prime) (hpk : p ≤ k) (hpd : ¬p ∣ d)
    (hprime : ∀ i < k, Nat.Prime (a + i * d)) : a ≤ p := by
  obtain ⟨i, hip, hdvd⟩ := exists_lt_dvd hp a d hpd
  have h : p = a + i * d :=
    (Nat.prime_dvd_prime_iff_eq hp (hprime i (hip.trans_le hpk))).mp hdvd
  omega

/-- If every prime `p ≤ k` divides `d`, then `primorial k ∣ d`. -/
theorem primorial_dvd {k d : ℕ} (h : ∀ p, p.Prime → p ≤ k → p ∣ d) : primorial k ∣ d := by
  rw [primorial]
  refine Finset.prod_primes_dvd _ (fun p hp ↦ ?_) (fun p hp ↦ ?_) <;>
    · simp only [Finset.mem_filter, Finset.mem_range] at hp
      first
        | exact hp.2.prime
        | exact h p hp.2 (by omega)

/--
**Proposition 1.**  Let `a, a + d, …, a + (k-1)d` be primes, `k ≥ 3`, `d ≥ 1`.  Then either

* `k < a` and `primorial k ∣ d`  (Case A), or
* `k` is prime, `a = k`, and `primorial k ∣ k * d`  (Case B).

The second divisibility is the division-free form of `(primorial k / k) ∣ d`; see
`primorial_dvd_mul_iff`.
-/
theorem structure_of_primeAP {k a d : ℕ} (hk : 3 ≤ k) (hd : 0 < d)
    (hprime : ∀ i < k, Nat.Prime (a + i * d)) :
    (k < a ∧ primorial k ∣ d) ∨ (Nat.Prime k ∧ a = k ∧ primorial k ∣ k * d) := by
  have ha : Nat.Prime a := by simpa using hprime 0 (by omega)
  have ha2 : 2 ≤ a := ha.two_le
  by_cases hak : k < a
  · refine Or.inl ⟨hak, primorial_dvd fun p hp hpk ↦ ?_⟩
    by_contra hpd
    have := le_of_not_dvd hp hpk hpd hprime
    omega
  · replace hak : a ≤ k := Nat.not_lt.mp hak
    -- `a ≤ k`, so `a` is a prime `≤ k`.  The term of index `a` is `a * (1 + d)`, which is not
    -- prime unless the index `a` is out of range, i.e. unless `a = k`.
    have hak' : a = k := by
      by_contra hne
      have hlt : a < k := lt_of_le_of_ne hak hne
      have h1 : Nat.Prime (a + a * d) := hprime a hlt
      have h2 : a = a + a * d := (Nat.prime_dvd_prime_iff_eq ha h1).mp ⟨1 + d, by ring⟩
      have h3 : 0 < a * d := Nat.mul_pos ha.pos hd
      omega
    refine Or.inr ⟨hak' ▸ ha, hak', primorial_dvd fun p hp hpk ↦ ?_⟩
    by_cases hpa : p = k
    · exact hpa ▸ ⟨d, rfl⟩
    · refine Dvd.dvd.mul_left ?_ k
      by_contra hpd
      have := le_of_not_dvd hp hpk hpd hprime
      omega

/-- For `k` prime, the Case B conclusion `primorial k ∣ k * d` is exactly
`(primorial k / k) ∣ d`, the form stated in `NOTES.md`. -/
theorem primorial_dvd_mul_iff {k d : ℕ} (hk : k.Prime) :
    primorial k ∣ k * d ↔ primorial k / k ∣ d := by
  obtain ⟨m, hm⟩ : k ∣ primorial k := hk.dvd_primorial
  rw [hm, Nat.mul_div_cancel_left _ hk.pos, Nat.mul_dvd_mul_iff_left hk.pos]

/-- Cases A and B are mutually exclusive. -/
theorem structure_exclusive {k a : ℕ} : ¬((k < a) ∧ (a = k)) := by omega

end ErdosProblem200
