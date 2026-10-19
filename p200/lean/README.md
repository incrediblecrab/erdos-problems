# `lean/` — machine-checked upper bound for Erdős #200

**What this is.** A Lean 4 + Mathlib formalisation of Propositions 1 and 2 of [`../NOTES.md`](../NOTES.md) and the asymptotic upper bound they imply. DeepMind's [FormalConjectures](https://github.com/google-deepmind/formal-conjectures) states Erdős 200 in `FormalConjectures/ErdosProblems/200.lean` and leaves **both** declarations as `sorry`. This discharges the second one, `Erdos200.erdos_200.variants.upper`, modulo an explicit prime-number-theorem hypothesis — which [`../lean-pnt/`](../lean-pnt/) then removes.

**What this is not.** Erdős 200 asks whether $L(N) = o(\log N)$. That remains open and is untouched here. Everything below is the *known* upper bound $L(N)\le(1+o(1))\log N$, which Erdős and Graham already attribute to the prime number theorem. The contribution is that it is now checked by a machine rather than asserted.

## Results

Write `L n` for `Erdos200.longestPrimeArithmeticProgressions n`.

| theorem | statement | file |
|---|---|---|
| `structure_of_primeAP` | $k\ge3$, $d\ge1$, all of $a,\dots,a+(k-1)d$ prime $\implies$ ($k<a$ and $W_k\mid d$) or ($k$ prime, $a=k$, $W_k\mid kd$) | `Structure.lean` |
| `primorial_dvd_mul_iff` | for $k$ prime, $W_k\mid kd \iff (W_k/k)\mid d$ — the form in `NOTES.md` | `Structure.lean` |
| `mul_primorial_lt` | $(k-1)\,W_k < kN$ | `Bound.lean` |
| `theta_lt_sharp` | $\theta(k) < \log N + 1/(k-1)$ | `Bound.lean` |
| `upper_unconditional` | $L(n) \le (1/\log 2 + o(1))\log n$ — **no hypotheses** | `Asymptotic.lean` |
| `erdos_200_variants_upper` | `PNT →` $L(n) \le (1 + o(1))\log n$ | `Asymptotic.lean` |

`1/\log 2 = 1.4426950…`, so the unconditional constant is **44 % worse** than the literature's $1$. That gap is not mathematical slack: Mathlib at the pinned revision proves only Chebyshev's `Chebyshev.theta_ge`,

```
∀ n : ℕ, n * log 2 - log (n + 1) - 2 * √n * log n ≤ Chebyshev.theta n
```

and does **not** contain the prime number theorem. So `PNT` here is a `def`, an assumed hypothesis — not a theorem, and not proved anywhere in this project:

```lean
def PNT : Prop := Tendsto (fun x : ℝ ↦ Chebyshev.theta x / x) atTop (nhds 1)
```

Both bounds come from one transfer lemma, `eventually_longest_le`: any $\theta(k)\ge(c-o(1))k$ yields $L(n)\le(1/c+o(1))\log n$. Chebyshev's $c=\log 2$ gives the first row, `PNT`'s $c=1$ the second.

**The hypothesis is discharged in [`../lean-pnt/`](../lean-pnt/).** Mathlib still has no PNT, but `PrimeNumberTheoremAnd` does, and that project proves `erdos_200_variants_upper_unconditional` — the same statement, constant $1$, no hypothesis. It is a separate Lake project only because PNT+ pins Lean `v4.32.2` while `formal-conjectures` pins `v4.33.1`, and one project cannot hold two Mathlib revisions. This project stays canonical: it is the one that can compare against `formal-conjectures`'s own `Expr`, and `lean-pnt/scripts/sync.sh` generates its copies of `Structure.lean`, `Bound.lean` and `Asymptotic.lean` from here.

## Verification

```bash
./verify.sh          # 98 s measured, warm cache and deps prebuilt; exit 0 = pass
```

It asserts three things:

1. `lake build` succeeds;
2. all 20 listed theorems depend on exactly `{propext, Classical.choice, Quot.sound}` — the `#assert_axioms` command in `Audit.lean` walks the transitive graph with `Lean.collectAxioms` and **errors** if `sorryAx` appears;
3. the elaborated conclusion of `erdos_200_variants_upper`, with the `PNT` binder stripped, is `Expr`-identical to `Erdos200.erdos_200.variants.upper` in FormalConjectures.

Check 3 exists because check 2 cannot see the interesting failure: a proof can be complete and still prove the wrong statement.

**Both checks were tested against planted defects**, since a passing checker that cannot fail proves nothing:

| planted defect | result |
|---|---|
| `sorry` in `primorial_lt_two_mul` | 6 downstream theorems reported `sorryAx`; `lake build` exit 1 |
| `sorry` in `primorial_dvd` | `verify.sh` exit 1 |
| conclusion reassociated to `log n * (1 + o n)` — equivalent, and still fully proved | `statement mismatch`; `lake build` exit 1 |

## Environment

| | |
|---|---|
| toolchain | `leanprover/lean4:v4.33.1` |
| formal-conjectures | `e5f428182a3ee32dde401eceb4a94ba0382e8434` |
| Mathlib | `0df444a360ea` (branch `v4.33.1`), via formal-conjectures |
| source | 673 lines across 5 files |
| clean rebuild of this library | 95.7 s wall / 24.6 s user (deps prebuilt) |
| `.lake/packages` on disk | 7.7 GB — gitignored |

`lake-manifest.json` is committed so the revisions above are reproducible. First build downloads the Mathlib olean cache automatically (`lake update` runs `cache get`).

## Files

| file | contents |
|---|---|
| `Erdos200/Structure.lean` | Proposition 1 and its three helper lemmas |
| `Erdos200/Bound.lean` | Proposition 2: the primorial and $\theta$ bounds |
| `Erdos200/Asymptotic.lean` | extraction from `Set.IsAPOfLength`, the $\theta\Rightarrow L$ transfer, both final theorems |
| `Erdos200/Audit.lean` | `#assert_axioms` and the statement-fidelity check |
| `verify.sh` | runs everything; exit 0 = pass |
