# erdos-40

This repository studies Erdős problem 40, how dense a set must be to force unbounded representation counts, with Lean reductions, a written reconstruction of a negative result and a checked greedy-set computation. The problem is open and is not solved here.

**Objective:** pin down what the problem asks, in its formal-conjectures form; prove what can be proved about the answer set, with a machine check wherever the claim is formal; reconstruct the known negative result with explicit constants; and check the one computation posted on the closely related problem 158.

**Inputs:** the problem's page on erdosproblems.com; the pages or threads of problems 28, 39, 158 and 1191; formal-conjectures at the commit pinned in `lean/lakefile.toml`; Pliego; Cilleruelo; O'Bryant's Sidon bibliography; Zeraoulia's note; and OEIS A005282. `refs/fetch.sh` downloads them, and nothing is redistributed.

**Files:**

- [`NOTES.md`](NOTES.md): statement, formal version, proofs, literature, computation and validation record
- [`lean/`](lean/README.md): Lean 4 proofs against formal-conjectures definitions and planted-defect checks
- [`src/`](src/README.md): Theorem B constant check, greedy $B_2[g]$ computation and final gate
- [`scripts/`](scripts/README.md): reproduction script
- [`results/`](results/README.md): committed outputs and logs
- [`refs/`](refs/README.md): fetch script and source ledger

**Try it:** `python3 src/final_check.py` rechecks the results and Lean records; `bash scripts/reproduce.sh` reruns everything first.

## Problem statement

> For what functions $g(N)\to \infty$ is it true that $\lvert A\cap \{1,\ldots,N\}\rvert \gg \frac{N^{1/2}}{g(N)}$ implies $\limsup 1_A\ast 1_A(n)=\infty$?

**Status: open.** Not solved here. $A$ is a set of natural numbers, and $1_A\ast 1_A(n)$ is the number of ordered pairs $(a,b)$ of elements of $A$ with $a+b=n$.

## Plain-language summary

Pick a set $A$ of whole numbers and, for each $n$, count the ways to write $n=a+b$ with $a$ and $b$ both in $A$. If $A$ holds every number, $n$ can be written in about $n$ ways, so the counts grow without limit. If $A$ holds only the powers of two, no $n$ can be written in more than two ways, so the counts stay bounded.

A set whose counts stay bounded cannot be very dense. If no count is above $C$, then $A$ has at most $\sqrt{2CN}$ elements below $N$, because pairs of such elements have fewer than $2N$ possible sums. Erdős asked how close to $\sqrt N$ elements such a set can keep. If $A$ has at least a constant times $\sqrt N/g(N)$ elements up to every large $N$, where $g$ grows to infinity, must its counts be unbounded? And for which $g$? The case $g=1$, where $A$ has at least a constant times $\sqrt N$ elements, is an open strengthening of the Erdős–Turán conjecture, problem 28.

**What is in this repository.** Three results about the set $\mathcal{G}$ of functions $g$ that work, and a check of one published computation.

1. Some $g$ works if and only if $g=1$ works. So the first half of the question, whether $\mathcal{G}$ is empty, is exactly the $g=1$ question, and a positive answer would also answer problems 28 and 158. Proved in Lean.
2. That question in turn is whether every set with bounded counts has $A(N)/\sqrt N$ dipping toward $0$ along some sequence of $N$. For sets in which each number has at most two representations $n=a+b$ with $a\le b$, this is problem 158, which is open. Proved in Lean.
3. No $g$ that grows like a power of $N$ works. Random sets of the kind Erdős and Rényi used have bounded counts and more than a constant times $N^{1/2-\epsilon}$ elements up to $N$. This is a written proof, a reconstruction of their method with explicit constants, not machine-checked.
4. A computation posted in the problem 158 thread, on the greedy set of that kind, is reproduced and carried to five times as many elements. Its ratio $A(N)/\sqrt N$ keeps falling. That decides nothing either way.

**The catch.** None of this touches the open part. Whether $\mathcal{G}$ is empty is exactly as open as the $g=1$ question, and nothing here decides whether a slowly growing function such as $\log N$ is in $\mathcal{G}$. The first result is a reduction, not progress on either question.

## Findings

Labels as in `NOTES.md`: [Lean] is machine-checked, [written] is a written proof, [computation] is a finite computation.

* **Theorem A** [Lean]. $\mathcal{G}$ is nonempty if and only if the strong form holds, the statement that $\lvert A\cap\{1,\ldots,N\}\rvert\gg N^{1/2}$ already implies $\limsup 1_A\ast 1_A(n)=\infty$. So formal-conjectures' `erdos_40.variants.weaker` is equivalent to the strong form. The statements are checked, as Lean expressions, against formal-conjectures' `erdos_40` and `erdos_40.variants.weaker`, and the strong form implies its `erdos_28` (§2).
* **The $B_2[g]$ form** [Lean]. The strong form holds if and only if every $B_2[g]$ set, for every $g$, has $\liminf A(N)/\sqrt N=0$; at $g=2$ that is the right side of formal-conjectures' `erdos_158`, less the hypothesis that $A$ is infinite (§3).
* **Closure** [Lean]. $\mathcal{G}$ is closed downward under $O(\cdot)$, and no $g$ with $\sqrt N=O(g)$ is in it (§2).
* **Theorem B** [written]. No $g$ with $g(N)\ge cN^{\epsilon}$ for large $N$ is in $\mathcal{G}$. The proof follows Erdős and Rényi's probabilistic method with explicit constants, and a script checks the constants numerically (§4).
* **The greedy $B_2[2]$ set** [computation]. The numbers posted in the problem 158 thread and the table of the Zenodo note are right. Carried to the 10,000th term the ratio falls to 0.515, and the local growth exponent of the terms is above $2$ at every checkpoint (§6).
* **A conflict in the literature**, left unresolved: a sharper liminf for Sidon sets is reported as $\le c$ on the problem 1191 page, as finite in formal-conjectures, and as $=0$ by Pliego. Problem 40 does not depend on it (§3).

**Stated plainly.** Problem 40 is open. The machine-checked results are reductions: the existence half of problem 40 is the strong form of the Erdős–Turán conjecture, which in turn is a statement about $B_2[g]$ sets of which problem 158 is the case $g=2$. The exclusion of powers of $N$ is a written reconstruction, and the computation is illustrative. I did not look for Theorem A beyond the sources in `NOTES.md` §9, so no claim of novelty is made.

Read [`NOTES.md`](NOTES.md).

## File reference

| path | contents |
|---|---|
| `NOTES.md` | the statement and its formal version, the proofs, the literature, the computation, and how each claim was checked |
| `lean/` | Lean 4 proofs of Theorem A and the $B_2[g]$ reductions against formal-conjectures' definitions; `verify.sh`, and `plants.py`, which plants defects to test it |
| `src/theorem_b_check.py` | a numerical check of the constants in Theorem B |
| `src/greedy_b2g.c`, `src/greedy_check.py` | the greedy $B_2[g]$ computation, and the checks that do not share its code |
| `src/final_check.py` | the gate: rechecks the results and Lean records, the statement, the Lean excerpts, every quotation, and every numeral in the prose of the documentation; plants a defect for each check; exit 0 = pass |
| `scripts/reproduce.sh` | reruns everything, then `src/final_check.py` |
| `results/` | the committed outputs and logs; its README lists them |
| `refs/` | `fetch.sh`, which downloads the sources, and a README saying which claim each supports |

## Conventions

$A(N)=\lvert A\cap\{1,\ldots,N\}\rvert$, and $1_A\ast 1_A(n)$, written $r_A(n)$ in the notes, counts ordered pairs. Lean 4 at the release in `lean/lean-toolchain`, with Mathlib and formal-conjectures at the revisions in `lean/lake-manifest.json`. Python 3 with numpy and mpmath, poppler's `pdftotext`, and a C compiler; the header of `scripts/reproduce.sh` says which interpreter it picks and how `$PY` overrides it.

## Attribution

The question at the top is quoted from [erdosproblems.com/40](https://www.erdosproblems.com/40), maintained by Thomas Bloom, whose companion repository [teorth/erdosproblems](https://github.com/teorth/erdosproblems) is licensed Apache 2.0. The Lean statements quoted in `NOTES.md` come from [google-deepmind/formal-conjectures](https://github.com/google-deepmind/formal-conjectures), licensed Apache 2.0, which `lean/` imports rather than copies.

No third-party paper, page or file is redistributed here. `refs/fetch.sh` retrieves them, and `refs/README.md` records which claim each one supports.

## License

MIT, for this work only; the quoted problem is from [erdosproblems.com](https://www.erdosproblems.com/40). See [`LICENSE`](LICENSE).
