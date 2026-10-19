# erdos-1186

This repository studies Erdős problem 1186, monochromatic $k$-term arithmetic progressions in 2-colourings, by certifying explicit upper-bound constructions. The problem is open and cannot be resolved with a finite computation; no lower-bound work was attempted here.

**Objective:** push the *upper* bounds down by explicit construction, certify every value in exact arithmetic, and reproduce the published numbers before claiming to beat them.

**Inputs:** the erdosproblems.com entry; [LuPe12], [PRS08], [CCS07], [RL12], [HHLM07] and other sources listed in [`refs/`](refs/README.md). Run `refs/fetch.sh`; none are redistributed here.

**Files:**

- [`NOTES.md`](NOTES.md): normalisation, reproductions, the recursion theorem, family limits, negative results, validation and provenance caveats
- [`src/`](src/): search, certification, exact arithmetic and independent C cross-checks
- [`scripts/reproduce.sh`](scripts/reproduce.sh): full and fast reproduction entry point
- [`results/`](results/README.md): certified colourings, stored bases, search logs and `claims.json`
- [`src/openai_transfer.py`](src/openai_transfer.py) and [`results/openai-transfer.json`](results/openai-transfer.json): the October release's digit-product interface and a checked obstruction to automatic wildcard transfer
- [`src/PeriodicTransfer.lean`](src/PeriodicTransfer.lean) and [`src/check_periodic_lean.py`](src/check_periodic_lean.py): a core-only formal proof of the direct periodic bridge, checked with the existing shared verifier components and no Mathlib installation
- [`refs/`](refs/README.md): fetch script and source ledger

**Try it:** `python3 src/final_check.py` re-derives every claim; `FAST=1 bash scripts/reproduce.sh` skips the long searches.

## Problem statement

> Let $\delta_k$ be such that in any $2$-colouring of $\{1,\ldots,n\}$ there exist at least $(\delta_k+o(1))n^2$ many monochromatic $k$-term arithmetic progressions. Give reasonable bounds (or even an asymptotic formula) for $\delta_k$.
>
> […] It is easier to study this quantity if we replace $\{1,\ldots,n\}$ with a finite field $\mathbb{F}_p$; let this analogue be denoted by $\tilde{\delta}_k$.

**Status: open** — "this is open, and cannot be resolved with a finite computation." Not solved here. Erdős [Er80, p. 93]; quoted from [erdosproblems.com](https://www.erdosproblems.com/1186) (see [Attribution](#attribution)).

## October 2026 transfer check

OpenAI's [family 160](https://github.com/openai/math/blob/adc7f1241b42e322a6451854ab7e4b4c146bf78a/lean/docs/160.md) claims superexponential van der Waerden bounds. Its [digit-product proposition](https://github.com/openai/math/blob/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/Quantitative-Superexponential-Bounds-for-van-der-Waerden-Numbers-September-23-2026/build/sections/07-transfers.tex) takes a cyclic progression-free coloring and produces larger interval colorings. That input property is weaker than the freely recolorable proper-coset property required by this repository's wildcard recursion.

`src/openai_transfer.py` checks the distinction with the cyclic three-term-progression-free word `0011` on `Z/4Z`. Its digit products pass the finite interval checks, but all proper cosets fail the existing `cosetprobe.safe_coset` checker. The known quadratic-residue wildcard base on `Z/11Z` remains a positive control. The exact progression counts and rejected cosets are in `results/openai-transfer.json`.

This is a negative result about automatic wildcard transfer, not a refutation of OpenAI's construction. It does not block the simpler periodic construction already used in this repository. The follow-up below establishes that route; the release's asymptotic threshold is not treated as an explicit construction at the small lengths studied here.

```bash
python3 src/openai_transfer.py --output results/openai-transfer.json
python3 src/openai_transfer.py --plant coloring   # must print FAIL and exit 1
python3 src/openai_transfer.py --plant interface  # must print FAIL and exit 1
```

### Direct periodic transfer and a conditional asymptotic bound

For any cyclic $k$-AP-free two-coloring of modulus $m$, its periodic extension has an increasing monochromatic progression exactly when its difference is divisible by $m$. [`PeriodicTransfer.lean`](src/PeriodicTransfer.lean) proves this equivalence and a nonvacuous example using only Lean's `Std`. The existing shared Erdős axiom probe and an independent kernel replay accept its declarations; [`results/periodic-lean.json`](results/periodic-lean.json) records the source and verifier hashes and the existing toolchain used. No `.lake` directory, Mathlib build, or new toolchain is needed.

The exact finite count, including difference zero as the existing computational gate does, is derived in [the periodic-transfer argument](PERIODIC-TRANSFER.md). It recovers the known upper bound $\delta_k\le1/(2(k-1)m)$ without a safe coset. Applying this to the [published cyclic theorem](https://github.com/openai/math/blob/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/Quantitative-Superexponential-Bounds-for-van-der-Waerden-Numbers-September-23-2026/build/sections/06-perturbation.tex) gives, **conditional on that theorem**, $\delta_k\le1/(2(k-1)k^{k/100000})$ for all sufficiently large $k$, hence $\delta_k^{1/k}\to0$.

The local bridge is checked; OpenAI's large-coloring theorem is not independently verified here. The selected Comparator challenge states an interval van der Waerden bound, which alone is not the cyclic premise used in this deduction. No explicit threshold, improved small-$k$ record, asymptotic equality, or priority claim is established. The earlier safe-coset counterexample remains valid.

```bash
python3 -B -m unittest discover -s src -p 'test_openai_transfer.py' -v
python3 -B src/check_periodic_lean.py --check
```

The Lean checker reuses `../formal/scripts/verify_lean.py` and its toolchain pin, whose hashes are recorded; it refuses to install a missing toolchain. That verifier lives in a sibling `formal/` checkout that is not part of this repository, so on a clone without it, or without `elan`, `src/final_check.py` reports the Lean check as skipped and its summary line says the run is incomplete; the stored `results/periodic-lean.json` is then not re-verified. Its `--plant sorry`, `--plant axiom`, and `--plant claim` cases must fail. The numerical transfer's additional `--plant degenerate` and `--plant normalization` cases exercise the zero-difference convention and the factor of two.

## Findings

Comparing the optimal periodic colourings at $m=88$ and $m=132$ showed both reduce to one 44-bit base with a free coset of size 4. Generalising gives a **wildcard recursion**, proved here as a theorem ([`NOTES.md`](NOTES.md) §3.1). A *base* is a colouring of $\mathbb{Z}_b$ with a proper coset $F$ of size $f<b$ such that every $k$-AP of nonzero difference that avoids $F$ is non-monochromatic, and every one that meets both $F$ and its complement is non-constant outside $F$. Every base gives

$$m_k(\mathbb{Z}_{bt}) \le \frac{(b-f)+f^2 m_k(\mathbb{Z}_{ft})}{b^2}, \qquad \delta_k \le \frac{1}{2(k-1)(b+f)} .$$

This reproduces **both** published constants exactly — $b=11,f=1\Rightarrow c_4\le1/72$ and $b=37,f=1\Rightarrow c_5\le1/304$ — which is the evidence that it measures the published quantity. Lu and Peng's two bases turn out to be the quadratic-residue (QR) colourings of $\mathbb{Z}_{11}$ and $\mathbb{Z}_{37}$, the colourings behind Rabung's lower bounds for van der Waerden numbers [RL12]. Every QR colouring that passes Rabung's test is a base, and so is every good zipped one [HHLM07] (NOTES §4). Together with one base found by search, at $k=5$, they give:

- $\delta_5 \le \mathbf{1/384}$, from a 5-AP-free $\mathbb{Z}_{44}$ base with $f=4$: a **20.83 % improvement** on Lu–Peng's $1/304$.
- $\delta_6 \le \mathbf{1/2280}$, from the zipped QR colouring of $\mathbb{Z}_{226}$ ($f=2$), 2.6 times below this repo's previous $k=6$ value $1807/1590140$.
- $\delta_7 \le \mathbf{1/7416}$ and $\delta_8 \le \mathbf{1/23016}$, from the QR colouring of $\mathbb{Z}_{617}$ and the zipped one of $\mathbb{Z}_{1642}$. No published bound was found for $k \ge 6$.

Over $\mathbb{F}_p$, block colourings improve $\tilde\delta_4$ by 2.8 % and $\tilde\delta_5$ by 10.6 %.

**The family is exhausted for $k \le 6$, subject to the remaining SAT caveat.** A base can be recoloured on $F$ to be cyclic AP-free, so $(k-1)b < W(k;2)$, the van der Waerden number, and an exact SAT scan up to that cap finds no base with larger $b+f$. So $1/72$, $1/384$ and $1/2280$ are the best this family gives for $k = 4, 5, 6$ (NOTES §5). The UNSAT answers are now LRAT-checked by `lrat-check` for every $k \le 5$ scan instance and for 68 of the 72 $k=6$ UNSAT instances; the four remaining $k=6$ instances are $(203,29)$, $(217,31)$, and the two symmetry-broken $b=221$ cases $(221,13)$ and $(221,17)$. The cap still uses $W(5;2)$ and $W(6;2)$ as reported in [RL12].

**$\delta_3$ and $\delta_4$ were not improved.** $\delta_3$ reproduces [PRS08] exactly. $\delta_4 \le 1/72$ is Lu and Peng's bound, reproduced exactly and now a theorem, and NOTES §6 shows this family cannot beat it, consistent with [LuPe12, Conjecture 1]. Bases for $k = 9, \dots, 12$ are in NOTES §4; they are checked by `src/verify_base.c` but are not in `results/claims.json`.

## Certified bounds

All exact rationals, re-derived by `src/final_check.py` from the stored colourings and bases.

| | bound | decimal | where | vs. random | vs. published |
|---|---|---|---|---|---|
| $\delta_3$ | $117/2192$ | 0.05337591 | blocks, $m=548$ | 0.854× | $=117/2192$ [PRS08] — **not improved** |
| $\delta_4$ | $1/72$ | 0.01388889 | recursion limit, QR base $\mathbb{Z}_{11}$, $f=1$ | 0.667× | $=1/72$ [LuPe12] — **not improved** |
| $\delta_5$ | $\mathbf{1/384}$ | 0.00260417 | recursion limit, base $\mathbb{Z}_{44}$, $f=4$ | 0.333× | **beats $1/304$ by 20.83 %** |
| $\delta_6$ | $\mathbf{1/2280}$ | 0.00043860 | recursion limit, zipped QR base $\mathbb{Z}_{226}$, $f=2$ | 0.140× | no published bound found |
| $\delta_7$ | $\mathbf{1/7416}$ | 0.00013484 | recursion limit, QR base $\mathbb{Z}_{617}$, $f=1$ | 0.104× | no published bound found |
| $\delta_8$ | $\mathbf{1/23016}$ | 0.00004345 | recursion limit, zipped QR base $\mathbb{Z}_{1642}$, $f=2$ | 0.078× | no published bound found |
| $\tilde\delta_4$ | $\mathbf{159/2888}$ | 0.05505540 | $\mathbb{F}_p$ blocks, $m=152$ | 0.881× | **beats $17/300$ by 2.84 %** |
| $\tilde\delta_5$ | $\mathbf{1283/51984}$ | 0.02468067 | $\mathbb{F}_p$ blocks, $m=114$ | 0.790× | **beats $3629/131424$ by 10.62 %** |
| $\tilde\delta_6$ | $2821/216000$ | 0.01306019 | $\mathbb{F}_p$ blocks, $m=60$ | 0.836× | no published bound found |
| $\tilde\delta_7$ | $1517/230640$ | 0.00657735 | $\mathbb{F}_p$ blocks, $m=62$ | 0.842× | no published bound found |

The $\delta_4$ to $\delta_8$ rows are limits of the recursion theorem: each is an infimum over explicit colourings, and every finite level lies strictly above it. The gate re-verifies each base and re-counts the theorem's identity on explicit lifts. Finite levels counted directly sit just above the limits: $\delta_5 \le 13421/5153632$ at $m=5324$ is $1.6\times10^{-8}$ above $1/384$ and is gate-checked, and $\delta_6 \le 12657/28857940$ at $m=25538$ is $3.0\times10^{-10}$ above $1/2280$ (`results/verify_lift_k6_25538.log`, not gate-checked).

## How the numbers are trusted

Certification is separated from search, and no self-reported score is accepted.

- Every bound is recomputed in exact `Fraction` arithmetic (`src/exact.py`), and a stored value is never read as a result. A recursion limit is recomputed from its base's $b$ and $f$, after the base itself is re-verified.
- Headline values are confirmed by code sharing nothing with the search: `src/verify_zm.c` (direct double loop over $\mathbb{Z}_m$, no weight table), `src/verify_two.c` (brute force over $\{1,\dots,n\}$), and `src/verify_base.c` (the two base conditions; it agrees with the Python check on all 21 cases of a planted cross-check).
- The periodic-to-integer step is measured, not assumed: the boundary term is $1/(2n)$, so (measured $-$ claimed) $\times n$ must stay flat at $\approx0.50$. For the 6-AP-free $\mathbb{Z}_{226}$ word at $n=20000,60000,120000$ it is **0.506, 0.501, 0.501**.
- Published values are reproduced first — six reproductions, including Lu–Peng's $B_{74}$ with its exact $d$-breakdown (74 APs at $d=0$, 72 at $d=37$) and the six two-colour lower bounds of [RL12, Table 1].
- The gate is tested against planted defects, not merely observed to pass. Six defects each make it exit 1: two overclaims, a one-bit corruption of a stored word, and three corruptions of the $\mathbb{Z}_{226}$ base. `src/safesat.py --self-test` likewise fails when its symmetry group is replaced by random permutations (NOTES §10).
- The SAT solvers' UNSAT answers support only the statement that the family is exhausted; no bound depends on them. `src/prove_unsat.py` rebuilds the CNFs with `safesat.encode_cnf`, runs Homebrew CaDiCaL 3.0.1 with LRAT output, and checks with `lrat-check` from `drat-trim` commit `2e3b2dc0ecf938addbd779d42877b6ed69d9a985`. `results/lrat_check.json` records 1118 accepted proofs and 4 timeouts: all UNSAT instances for $k=3,4,5$ are checked, and $k=6$ has 68 checked and 4 unchecked. The two $b=221$ cases use symmetry-breaking clauses; LRAT checks the encoded CNF, while soundness of those extra clauses rests on the symmetry argument in NOTES §5.

## Try it

```bash
bash scripts/reproduce.sh          # full pipeline
FAST=1 bash scripts/reproduce.sh   # skip the long searches
python src/final_check.py          # re-derive every claim; exit 0 = pass
```

| path | contents |
|---|---|
| `src/exact.py` | exact rational evaluator — the certification path |
| `src/periodic.py`, `src/blocks.py`, `src/twoscale.py` | the three colouring families |
| `src/bigm.py` | tiling-seeded search at large $m$ |
| `src/cosetprobe.py` | finds the wildcard coset; exact safety test + positive control |
| `src/chain.py` | builds and certifies one level of the recursion |
| `src/rabung.c` | quadratic-residue and zipped colourings; reproduces [RL12, Table 1] |
| `src/bases.py` | assembles and re-checks the stored bases |
| `src/wildcard.py` | the $k=5$ restricted search |
| `src/apfree.py` | SAT search for cyclic AP-free moduli |
| `src/safesat.py` | exact SAT scan for bases, with symmetry breaking and a self-test |
| `src/verify_zm.c`, `src/verify_two.c`, `src/verify_zp.c`, `src/verify_base.c` | independent C cross-checks |
| `src/final_check.py` | re-derives every claim; exit 0 = pass |
| `results/claims.json` | the asserted bounds, machine-checked |
| `results/` | certified colourings and search logs |
| `NOTES.md` | normalisation, reproductions, the recursion theorem, the family's limits, negative results, validation |

Read [`NOTES.md`](NOTES.md) for the derivation, the failure of the $k=4$ family, and the provenance caveats.

## Attribution

The question quoted at the top is from [erdosproblems.com/1186](https://www.erdosproblems.com/1186), maintained by Thomas Bloom, whose companion repository [teorth/erdosproblems](https://github.com/teorth/erdosproblems) is licensed Apache 2.0. The attribution there is to Erdős, *A survey of problems in combinatorial number theory*, Ann. Discrete Math. **6** (1980), 89–115, at p. 93 [Er80].

Short quotations from [LuPe12] in `NOTES.md` are cited in place and used to fix the normalisation and to record which published values were reproduced. No third-party paper or web page is redistributed here: `refs/fetch.sh` retrieves them and `refs/README.md` records which claim each one supports.

## License

MIT, for this work only; the quoted problem is from [erdosproblems.com](https://www.erdosproblems.com/1186). See [`LICENSE`](LICENSE).
