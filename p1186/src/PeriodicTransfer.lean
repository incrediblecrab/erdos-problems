import Std

namespace PeriodicTransfer

def Mono (c : Nat → Bool) (m k a d : Nat) : Prop :=
  ∀ i : Fin k, c ((a + i.val * d) % m) = c (a % m)

def CyclicAPFree (c : Nat → Bool) (m k : Nat) : Prop :=
  ∀ a d : Fin m, d.val ≠ 0 → ¬ Mono c m k a.val d.val

theorem mono_iff_period_dvd
    (c : Nat → Bool) (m k a d : Nat) (hm : 0 < m)
    (hfree : CyclicAPFree c m k) :
    Mono c m k a d ↔ m ∣ d := by
  constructor
  · intro hmono
    by_cases hz : d % m = 0
    · exact Nat.dvd_of_mod_eq_zero hz
    · apply False.elim
      apply hfree ⟨a % m, Nat.mod_lt _ hm⟩ ⟨d % m, Nat.mod_lt _ hm⟩ hz
      intro i
      simpa [Nat.add_mod, Nat.mul_mod, Nat.mod_mod] using hmono i
  · intro hdiv i
    have hz : d % m = 0 := Nat.mod_eq_zero_of_dvd hdiv
    simp [Nat.add_mod, Nat.mul_mod, hz]

def fourColor (n : Nat) : Bool := decide (2 ≤ n)

theorem fourColor_ap_free : CyclicAPFree fourColor 4 3 := by
  unfold CyclicAPFree Mono fourColor
  decide

theorem fourColor_mono_iff (a d : Nat) :
    Mono fourColor 4 3 a d ↔ 4 ∣ d :=
  mono_iff_period_dvd fourColor 4 3 a d (by decide) fourColor_ap_free

end PeriodicTransfer
