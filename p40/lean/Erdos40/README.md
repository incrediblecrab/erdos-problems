# Erdos40

**Objective.** The Lean sources: theorems about the definitions of formal-conjectures at the commit pinned in `../lakefile.toml`, each using no axiom beyond `propext`, `Classical.choice` and `Quot.sound`.

**Inputs.** formal-conjectures' `Erdos40.Erdos40For`, `sumRep`, `Erdos28.erdos_28`, `Erdos158.B2` and `Erdos158.erdos_158`, imported rather than restated.

| file | contents |
|---|---|
| `Basic.lean` | `sumRep A n ≤ n + 1`; a finite `limsup` is a uniform bound; both depend only on an initial segment of `A`; a property that reads finitely many coordinates of `ℕ → Bool` is clopen |
| `Equivalence.lean` | Theorem A, `exists_iff`, with `answerSet_nonempty_iff`, `strong_implies_erdos_28`, `Erdos40For.of_isBigO` and `not_erdos40For_sqrt` |
| `Problem158.lean` | `strong_iff_forall_B2`, `weaker_iff_forall_B2` and `strong_implies_erdos_158`: the same question for `B2[g]` sets |
| `Audit.lean` | fails the build if a listed theorem uses another axiom or a statement differs, as an `Expr`, from the formal-conjectures declaration it should match |

`Audit.lean` runs inside the environment it checks, so `../Probe.lean` repeats its checks from a separate process; `../../NOTES.md` §8 describes both.
