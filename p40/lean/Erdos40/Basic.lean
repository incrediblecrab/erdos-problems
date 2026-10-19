/-
Copyright (c) 2026 Max Marquardt.
Released under the MIT licence; see LICENSE at the repository root.
-/
import FormalConjectures.ErdosProblems.«40»

/-!
# Basic facts about `sumRep` and the counting function

Everything here is about the definitions in `FormalConjectures`: `sumRep A n` is the number of ordered pairs `(a, b)` with `a + b = n` and `a, b ∈ A` (so `0 ∈ A` counts and `(a, a)` counts once), and the counting function is `(A ∩ Set.Icc 1 N).ncard`.

* `sumRep_le`: `sumRep A n ≤ n + 1`.
* `exists_bound_of_limsup_ne_top`, `limsup_ne_top_of_bound`: a finite `limsup` of `sumRep A` is the same as a uniform bound on it, because every single value is finite.
* `sumRep_congr`, `count_congr`: `sumRep A n` depends only on `A ∩ [0, n]`, and the count up to `N` only on `A ∩ [0, N]`.
* `isClopen_of_local`: in `ℕ → Bool` with the product topology, a property that reads only finitely many coordinates cuts out a clopen set.
-/

namespace ErdosProblem40

open Filter Set AdditiveCombinatorics

theorem sumRep_le (A : Set ℕ) (n : ℕ) : sumRep A n ≤ n + 1 := by
  classical
  rw [sumRep_def]
  exact (Finset.card_filter_le _ _).trans (Finset.Nat.card_antidiagonal n).le

theorem exists_bound_of_limsup_ne_top {A : Set ℕ}
    (h : limsup (fun N ↦ (sumRep A N : ℕ∞)) atTop ≠ ⊤) : ∃ C : ℕ, ∀ n, sumRep A n ≤ C := by
  set L := limsup (fun N ↦ (sumRep A N : ℕ∞)) atTop with hL
  have hlt : L < ((L.toNat + 1 : ℕ) : ℕ∞) := by
    conv_lhs => rw [← ENat.natCast_toNat h]
    exact_mod_cast Nat.lt_succ_self _
  obtain ⟨n₀, hn₀⟩ := eventually_atTop.1 (eventually_lt_of_limsup_lt hlt)
  refine ⟨max L.toNat n₀, fun n ↦ ?_⟩
  by_cases hn : n₀ ≤ n
  · have h' : sumRep A n < L.toNat + 1 := by exact_mod_cast hn₀ n hn
    omega
  · have := sumRep_le A n
    omega

theorem limsup_ne_top_of_bound {A : Set ℕ} {C : ℕ} (h : ∀ n, sumRep A n ≤ C) :
    limsup (fun N ↦ (sumRep A N : ℕ∞)) atTop ≠ ⊤ := by
  have hle : limsup (fun N ↦ (sumRep A N : ℕ∞)) atTop ≤ C :=
    limsup_le_of_le (h := Eventually.of_forall fun n ↦ by exact_mod_cast h n)
  exact ne_top_of_le_ne_top (ENat.natCast_ne_top C) hle

theorem sumRep_congr {A B : Set ℕ} {n : ℕ} (h : ∀ i ≤ n, (i ∈ A ↔ i ∈ B)) :
    sumRep A n = sumRep B n := by
  classical
  rw [sumRep_def, sumRep_def]
  congr 1
  apply Finset.filter_congr
  intro p hp
  rw [Finset.HasAntidiagonal.mem_antidiagonal] at hp
  rw [h p.1 (by omega), h p.2 (by omega)]

theorem count_congr {A B : Set ℕ} {N : ℕ} (h : ∀ i ≤ N, (i ∈ A ↔ i ∈ B)) :
    (A ∩ Icc 1 N).ncard = (B ∩ Icc 1 N).ncard := by
  congr 1
  ext i
  simp only [mem_inter_iff, mem_Icc]
  constructor
  · rintro ⟨hi, h1, h2⟩
    exact ⟨(h i h2).1 hi, h1, h2⟩
  · rintro ⟨hi, h1, h2⟩
    exact ⟨(h i h2).2 hi, h1, h2⟩

theorem isClopen_of_local (n : ℕ) (P : (ℕ → Bool) → Prop)
    (hP : ∀ f f' : ℕ → Bool, (∀ i ≤ n, f i = f' i) → P f → P f') : IsClopen {f | P f} := by
  let ρ : (ℕ → Bool) → (Fin (n + 1) → Bool) := fun f i ↦ f i
  have hρ : Continuous ρ := continuous_pi fun i ↦ continuous_apply (i : ℕ)
  have e : {f | P f} = ρ ⁻¹' (ρ '' {f | P f}) := by
    ext f
    constructor
    · intro hf
      exact ⟨f, hf, rfl⟩
    · rintro ⟨f', hf', hρf⟩
      refine hP f' f (fun i hi ↦ ?_) hf'
      simpa [ρ] using congrFun hρf ⟨i, by omega⟩
  rw [e]
  exact (isClopen_discrete _).preimage hρ

end ErdosProblem40
