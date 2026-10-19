# erdos-1086

This repository studies Erdős problem 1086, "many triangles of one area," by reproducing and extending tractable computations around the open problem. It does not solve the problem; the status recorded here is open, with $n^2\log\log n\ll g(n)\ll n^{20/9}$.

**Objective:** sharpen the estimate where it is tractable, reproducing the published grid count first and recording caveats against overclaiming.

**Inputs:** the problem page and the sources listed in [`refs/`](refs/README.md), fetched locally by `refs/fetch.sh`; no third-party sources are redistributed.

**Files:**

- [`NOTES.md`](NOTES.md): proofs, caveats, source disagreements and the derivation behind the stated findings
- [`src/`](src/README.md): code and build steps for the grid counts, certificates and checks
- [`results/`](results/README.md): claims, data and logs, including exact grid counts and small-$n$ results
- [`refs/`](refs/README.md): source list and fetch script

**Try it:** `python3 src/final_check.py` rechecks the committed claims.

## Problem statement

> Let $g(n)$ be minimal such that any set of $n$ points in $\mathbb{R}^2$ contains the vertices of at most $g(n)$ many triangles with the same area. Estimate $g(n)$.

**Status: open**, with $n^2\log\log n\ll g(n)\ll n^{20/9}$. Not solved here.

## Findings

Proofs and caveats are in [`NOTES.md`](NOTES.md).

- The square grid gives $g(n)\ge(6e^\gamma/\pi^2-o(1))\,n^2\log\log n$, where $6e^\gamma/\pi^2=1.08276\ldots$, and no area does asymptotically better there. The proof is informal and unreviewed, and no novelty is claimed.
- $g(5)=7$ by hand and by certificate; $g(6)=12$, $g(7)=18$ and $g(8)=28$ computer-assisted. Exactly two six-point sets span 12 triangles of one area, the $2\times3$ grid and the affine regular hexagon, and no seven-point set spans 19.
- Exact grid counts up to $800\times800$ give $g(640{,}000)\ge911{,}639{,}642{,}036$.
- The problem page disagrees with its sources in three places.

## License

MIT, for this work only; the quoted problem is from [erdosproblems.com](https://www.erdosproblems.com/1086). See [`LICENSE`](LICENSE).
