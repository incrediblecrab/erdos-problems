# Erdős 200, upper bound, without the PNT hypothesis

The sibling project [`../lean`](../lean) proves

```
theorem erdos_200_variants_upper (h : PNT) :
    ∃ (o : ℕ → ℝ) (_ : o =o[atTop] (1 : ℕ → ℝ)),
      ∀ n, longestPrimeArithmeticProgressions n ≤ (1 + o n) * Real.log n
```

where `PNT : Tendsto (fun x ↦ Chebyshev.theta x / x) atTop (𝓝 1)` is an *assumed* hypothesis, because Mathlib has no Prime Number Theorem. Its unconditional bound therefore carries the weaker constant `1 / log 2 = 1.4427`.

This project removes the hypothesis. `PrimeNumberTheoremAnd` (PNT+) proves `chebyshev_asymptotic : θ ~[atTop] id`, which is exactly `PNT`, so

```
theorem erdos_200_variants_upper_unconditional :
    ∃ (o : ℕ → ℝ) (_ : o =o[atTop] (1 : ℕ → ℝ)),
      ∀ n, longestPrimeArithmeticProgressions n ≤ (1 + o n) * Real.log n
```

holds with no hypothesis. That is the statement `formal-conjectures` leaves as `sorry` at `Erdos200.erdos_200.variants.upper`.

## What this does and does not settle

It does **not** touch Erdős 200. The open question is whether the bound can be improved to `o(log N)`; see [`../NOTES.md`](../NOTES.md) §6 and §8 for why that is out of reach of sieve methods. This is the *known* half of the problem, formalised — a transcription of the classical argument, not new mathematics.

## Why PNT+ contains `sorry`, and why it does not matter here

`PrimeNumberTheoremAnd` has two `sorry`s, both in `Wiener.lean`:

| line | declaration | reached by |
|---|---|---|
| 323 | `prelim_decay_2` | nothing |
| 342 | `prelim_decay_3` | `decay_alt` only |

`decay_alt` is itself referenced nowhere in the project's 27 source files. So both are dead code. That is a textual argument, and textual arguments are not proofs — the `#assert_axioms chebyshev_asymptotic` line in `Erdos200/Audit.lean` is what actually settles it, by failing the build if `sorryAx` appears anywhere in the dependency graph. It reports:

```
'chebyshev_asymptotic' depends on axioms: [propext, Classical.choice, Quot.sound]
```

For contrast, `#print axioms prelim_decay_3` reports `[propext, sorryAx, Classical.choice, Quot.sound]`, confirming the check is not vacuous.

## Relationship to `../lean`

`../lean` is canonical. PNT+ pins Lean 4.32.2 and `formal-conjectures` pins 4.33.1, and two Mathlib revisions cannot coexist in one Lake project, so the shared modules are **generated** here by `scripts/sync.sh` with exactly two mechanical substitutions:

1. `import FormalConjectures.ErdosProblems.«200»` → `import Erdos200.FCShim`
2. `Set.mem_ofPred_eq` → `Set.mem_setOf_eq` (Mathlib renamed this on 2026-07-09; both are `rfl` lemmas with identical statements)

`Erdos200/{Structure,Bound,Asymptotic}.lean` are generated and gitignored. Edit `../lean/Erdos200/` instead.

Because `formal-conjectures` cannot be a dependency, the definitions it supplies are copied into `Erdos200/FCShim.lean`. `scripts/check_shim.py` re-downloads `formal-conjectures` `main` and fails if any copied line differs. This replaces the `Expr`-level fidelity check that `../lean/Erdos200/Audit.lean` performs; it is a weaker guarantee, and the difference is deliberate and stated rather than glossed over.

## Verifying

```bash
bash verify.sh
```

which regenerates the shared modules, checks `FCShim.lean` against upstream, builds, and asserts that all 4 audited declarations rest only on `{propext, Classical.choice, Quot.sound}`. Exit 0 means all four checks passed.

Both gates were tested against planted defects:

| planted defect | caught by | result |
|---|---|---|
| `n < l` → `n ≤ l` in the copied `IsAPOfLengthWith` | `check_shim.py` | exit 1, exact line flagged |
| two non-contiguous declarations merged into one verbatim block | `check_shim.py` | exit 1 |

## Environment

| component | pin |
|---|---|
| Lean | `leanprover/lean4:v4.32.2` |
| `PrimeNumberTheoremAnd` | `a5154676af9aa3095150ee410cdda80555aa0642` (2026-08-30) |
| Mathlib (via PNT+) | `905b95818eb32af7874a58b427f50c1711a5e96c` |

## Measured

- Port cost from the 4.33.1 sources: **one** identifier and one import line.
- Building PNT+ through `Consequences`: 3708 jobs, 3 min 17 s wall, warm Mathlib cache.
- Clean rebuild of the five Erdős 200 modules against PNT+: 3 min 57 s wall.
