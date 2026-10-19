# erdos-1162

This repository studies Erdős problem 1162, the number of subgroups of $S_n$ and their order, with exact computations, subgroup-family counts and conjectural asymptotics. The problem is open and is not solved here.

**Objective:** fix the ambiguity in the displayed statement, compute exactly to $n=512$ the rank-stratified family of subgroups that carries the whole first-order answer, identify limiting constants in closed form, and turn the vague second half of the question into a precise conjecture with the gap to the published theorem located exactly.

**Inputs:** the problem statement from [erdosproblems.com](https://www.erdosproblems.com/1162), Roney-Dougal–Tracey, OEIS A005432, A000638, A000085 and A006116, and the GroupNames table of transitive groups of degree up to 15. Run `refs/fetch.sh`; nothing is redistributed here.

**Files:**

- [`NOTES.md`](NOTES.md): proofs, heuristics, source checks, caveats and the numbered sections cited below
- [`src/`](src/README.md): exact counters, lattice enumerators, asymptotic checks and verification code
- [`scripts/`](scripts/README.md): the reproduction script
- [`data/`](data/README.md): exact lattice, rank-count and model data
- [`results/`](results/README.md): measured artifacts and logs
- [`refs/`](refs/README.md): fetch script and source ledger

**Try it:** `python3 src/final_check.py` rechecks the artifacts and compares the documented decimal numbers with fresh computation; `bash scripts/reproduce.sh` regenerates every artifact first.

## Problem statement

> Give an asymptotic formula for the number of subgroups of $S_n$. Is there a statistical theorem on their order? (Erdős–Turán, [Va99, 5.73].)

**Status: open.** Not solved here.

## Plain-language summary

Take $n$ numbered cards. A *shuffle* is any way of rearranging them, and $S_n$ is the collection of all of them — for 8 cards that is $8! = 40\,320$ shuffles.

Some sub-collections of shuffles are self-contained: do any two of them in a row and you land on something still in the collection, and anything you do can be undone without leaving it. Those are the **subgroups**. "Leave everything alone" on its own is one. All 40 320 together is another. The interesting ones sit in between — for instance, every shuffle that disturbs only the first three cards and leaves the other five alone.

Erdős and Turán asked two questions about them.

1. How many subgroups are there? For 8 cards, exactly 151 221. For 18 cards it is about $7.6\times10^{18}$, so actually listing them is hopeless well before you get there and the goal is a formula for roughly how big the number is. This half is now largely answered by other people: it is about $2^{n^2/16}$. For 100 cards that is a number with 189 digits.
2. How big is a typical subgroup? Pick one of those 151 221 at random. How many shuffles does it contain? Erdős wanted a "statistical theorem" — a rule for what to expect. This half is barely touched. One published theorem addresses it, from 2025.

This repository mostly works on question 2, in two rounds. The first counted one well-behaved family of subgroups exactly, out to 512 cards: its members are overwhelmingly concentrated at about $2^{n/4}$ shuffles, with a spread that stays bounded as $n$ grows. The second round proved that this family is a vanishing fraction of all subgroups, so it cannot be the answer. It then counted a larger family exactly, also out to 512 cards, whose typical member has about $2^{3n/8}$ shuffles; that limit is argued, not proved. Along the way it proved that there are at least about $2^{n^2/16}\,n^{7n/8}$ subgroups, up to a factor $2^{O(n)}$, which for large $n$ is more than the 2025 paper's lower bound. Whether that proof is new was checked against the 2025 paper only.

The catch, and it is a real one: the larger family is not known to be typical either. Nobody has proved that the subgroups outside it are rare, and at the sizes that can be checked directly, 6 to 18 cards, between 41.5 % and 59.4 % of all subgroups lie outside it and its companion family for odd $n$. For $4\le n\le8$ the average subgroup is nearer $2^{n/2}$ than $2^{3n/8}$. So $2^{3n/8}$ is a conjecture with some evidence behind it, not an answer, and the small-$n$ data are still drifting.

## Findings

- The logarithm in the displayed statement is **base 2** (§1) — the constant 1/16 is meaningless otherwise, and the problem page never says so.
- The full subgroup lattice of $S_n$ recomputed for $n \le 8$ (151 221 subgroups, 296 classes at $n=8$), matching A005432 and A000638 exactly.
- Elementary abelian 2-subgroups counted exactly to $n = 512$ — this family already carries the entire leading term $2^{n^2/16}$. Their rank concentrates at $n/4$ with a spread that does not grow: numerically the law tends to $\mathbb{P}(n/4+j) \propto 2^{-j^2}$, and the standard deviation of the exact rank distribution agrees with the closed form $0.849305961022908\ldots$ to **35 decimal places** at $n=512$. That is exact computation at each $n$; the limit law itself, and $e(n) \sim G_{\lfloor n/2\rfloor}$ times the number of 2-, 4-block partitions, are proved in NOTES §9.10 (October 19, 2026) with rates far weaker than the data show.
- A single maximal elementary abelian gives only $n^2/16 + 2.8820499654\ldots$ — a *bounded* excess, with the constant $\log_2(\vartheta/\eta)$ in closed form. So the $n\log n$ gain in [RoTr25] Theorem 1 comes entirely from multiplicity.
- The second-order term of $\log_2 e(n)$, the elementary abelian count, is $\tfrac34 n\log_2 n$. It was derived from orbit structure (orbits of size 4 beat size 2 in a tie at the leading constant), confirmed by **out-of-sample** extrapolation ($A_\infty = 0.75057$), and is now proved up to $O(n)$ (Proposition 9.2).
- [RoTr25] Theorem 6 — the only published result on the *order* — has one lossy step, located exactly in §7: Lemma 2.4 is applied with $\ell \approx n$, where the maximal elementary abelian subgroups have rank $n/2$. With $n/2$ the threshold equation has **discriminant exactly zero** and a double root at $1/4$, and Conjecture 7.1 is that Theorem 6 holds for both of its constants up to $1/4$. The first version of this bullet said Theorem 6 *"is not optimal"*; §7 does not show that. The double root turned out to be automatic, under the assumptions of §9.9, for any family that reaches $2^{n^2/16}$, and the claim that $1/4$ is optimal for ν is withdrawn.
- **Theorem 9.1:** $\log_2 f(n) \ge n^2/16 + \tfrac78 n\log_2 n - O(n)$ (§9.4), by counting the Frattini-closed subgroups that induce a conjugate of the extraspecial group 8T22 $= 2^{1+4}_+$ on each of $n/8$ blocks of 8 points. It improves the large-$n$ constant: $\liminf_n (\log_2 f(n) - n^2/16)/(n\log_2 n) \ge 7/8$, where the explicit constant of the [RoTr25] proof is at most $5/906$. Theorem 1 of [RoTr25] is stated for all $n > 1$, where $n = 2$ caps the constant at $3/8$, so Theorem 9.1 does not let one put $7/8$ into it. Novelty was checked against [RoTr25] only, which singles out the same group in its upper-bound proof.
- **Proposition 9.2:** $e(n)/f(n) \le 2^{-\frac18 n\log_2 n + O(n)} \to 0$ (§9.5). The elementary abelian family of the bullets above is a vanishing fraction of all subgroups, so Conjecture 7.2, which put a typical subgroup at $2^{n/4}$, is withdrawn, not disproved.
- **A larger family,** the Frattini-closed 2-subgroups with every orbit of size at most 8, counted exactly for every $n \le 512$ (§9.1, §9.6). Inside it the mean of $\log_2|H|/n$ is 0.3835 at $n = 512$ and falling toward $3/8$, and the sd grows like $n^{1/4}$. The limits come from a saddle-point heuristic, not a proof.
- **Conjecture 9.3** replaces 7.2: $\log_2 f(n) = n^2/16 + \tfrac78 n\log_2 n + O(n)$, and a random subgroup has $\log_2|H| = (\tfrac38 + o(1))n$ in probability (§9.8). It rests on three unproved premises, which together say that the subgroups outside this family and its odd-$n$ companion $S_3(n)$ are negligible.

Stated plainly: the first family is proved atypical, and the second is not known to be typical. For $6 \le n \le 18$ roughly half of all subgroups lie outside it and its companion, the ratios of $f(n)$ to the formula the premises would give are still drifting, and for $4 \le n \le 8$ the mean of $\log_2|H|/n$ over all subgroups is between 0.47 and 0.51, above both $1/4$ and $3/8$. The asymptotic statements proved here are Theorem 9.1 and Propositions 9.2 and 9.4. The limit laws of both families, and the approximation $e(n) \approx G_{n/2}M(n)$ behind the first, rest on exact computation at each $n$ and on heuristic arguments, not on proofs. See [`NOTES.md`](NOTES.md) §8 and §9.8.

## File reference

| path | contents |
|---|---|
| [`src/elemab.py`](src/elemab.py) | exact $a(n,k)$ by exponential formula + $q$-binomial inversion |
| [`src/verify_elemab.py`](src/verify_elemab.py) | brute-force permutation check, sharing no code with it |
| [`src/subgroups.c`](src/subgroups.c) | full subgroup lattice of $S_n$, $n\le8$ (C, conjugacy-class expansion) |
| [`src/asymptotics.py`](src/asymptotics.py) | rank distribution, Galois numbers and the excess fit |
| [`src/sd_precision.py`](src/sd_precision.py) | the rank sd against its closed form, exact at 60 dp |
| [`src/orbit_model.py`](src/orbit_model.py) | orbit-structure model for the $n\log n$ term; held-out extrapolation |
| [`src/census8.py`](src/census8.py) | every 2-subgroup of $S_n$, $n\le8$, by brute force on permutations; the transitive types of §9.1 |
| [`src/block_families.py`](src/block_families.py) | the Frattini-closed family of §9, exact to $n = 512$ |
| [`src/final_check.py`](src/final_check.py) | checks the artifacts and compares every decimal number in `NOTES.md` with a fresh computation, apart from the timings and the marked first-version values; exit 0 = pass |
| [`scripts/reproduce.sh`](scripts/reproduce.sh) | regenerates every artifact, then runs `final_check.py` |
| [`data/`](data/README.md) and [`results/`](results/README.md) | exact lattices, rank counts, measured artifacts and logs; each folder's README lists its files |
| [`refs/`](refs/README.md) | fetch script and source ledger |

## Conventions

All logarithms base 2. This repository needs Python 3 with numpy and mpmath (developed on 3.14.7), and poppler's `pdftotext` for the quotation checks; the header of `scripts/reproduce.sh` says which interpreter it picks and how `$PY` overrides it. The C enumerator needs a C99 compiler; it uses NEON intrinsics on arm64 and falls back to portable code elsewhere.

## Attribution

The question as displayed above is quoted from [erdosproblems.com](https://www.erdosproblems.com/1162), maintained by Thomas Bloom, whose companion repository [teorth/erdosproblems](https://github.com/teorth/erdosproblems) is licensed Apache 2.0. The original attribution there is to Erdős and Turán via [Va99] — *Some of Paul's favorite problems*, the booklet produced for the conference "Paul Erdős and his mathematics", Budapest, July 1999 — problem 5.73.

No third-party paper, OEIS page, or API response is redistributed here. `refs/fetch.sh` retrieves them; [`refs/README.md`](refs/README.md) records which claim each one supports.

## License

MIT, for this work only; the quoted problem is from [erdosproblems.com](https://www.erdosproblems.com/1162). See [`LICENSE`](LICENSE).
