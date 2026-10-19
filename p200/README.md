# erdos-200

This repository studies Erdős problem 200, the longest arithmetic progression of primes in $\{1,\ldots,N\}$, through computation, heuristic comparison, sieve-limit analysis and Lean 4 formalisation. The problem is open and is not solved here.

**Objective:** establish what is rigorously true, recompute $L(N)$ as far as compute allows, test the standard heuristic against measured counts, and quantify why sieve methods cannot close the gap.

**Inputs:** a local mirror of erdosproblems.com problem 200, OEIS A005115, A113827, A093364 and A133277, Luhn's records page, Granville's *Prime Number Patterns* and the Green–Tao note. Run `refs/fetch.sh`; none are redistributed here.

**Files:**

- [`NOTES.md`](NOTES.md): analysis, computations, corrections, sieve limits and citations
- [`src/`](src/README.md): AP search, verification, heuristic and analysis scripts
- [`lean/`](lean/README.md): Lean 4 + Mathlib formalisation with the PNT hypothesis and unconditional constant $1/\log2$
- [`lean-pnt/`](lean-pnt/README.md): Lean proof using `PrimeNumberTheoremAnd` to recover the constant $1$
- [`scripts/reproduce.sh`](scripts/reproduce.sh): computation reproduction entry point
- [`data/`](data/): verified record data
- [`results/`](results/README.md): measured artifacts and logs
- [`refs/`](refs/README.md): fetch script and source ledger

**Try it:** `python3 src/final_check.py` re-derives the documented claims; `bash lean/verify.sh` and `bash lean-pnt/verify.sh` check the two Lean developments.

## Problem statement

> Does it have length $o(\log N)$?

**Status: open.** Not solved here.

## Findings

$L(N)$ independently re-derived for every $N\le5.73\times10^{11}$. The Hardy–Littlewood model matches exhaustive counts to within 1 % over 30 tests; it predicts $L(N)\sim2\log N/\log\log N$, so the answer is conjecturally **yes** — but the decay is $2/\log\log N$, invisible to any computation. Selberg's sieve, the large sieve and Gallagher's larger sieve cannot prove $L(N)<c\log N$ for any fixed $c<1$: at $k=c\log N$ a proof must save a factor $e^{(2/c-1)k}$, and these sieves save $e^{o(k)}$ (NOTES §6, corrected September 22, 2026, after the earlier argument for it turned out to be invalid). The known upper bound is now **machine-checked in Lean 4**: $L(N)\le(1+o(1))\log N$ with **no hypotheses** — the statement DeepMind's FormalConjectures leaves as `sorry`. Mathlib has no prime number theorem, so `lean/` proves that form conditionally and gets $1/\log2=1.4427$ unconditionally; `lean-pnt/` discharges the hypothesis against `PrimeNumberTheoremAnd` and recovers the constant $1$.

Read [`NOTES.md`](NOTES.md), and [`lean/README.md`](lean/README.md) and [`lean-pnt/README.md`](lean-pnt/README.md) for the formal part.

## File reference

| path | contents |
|---|---|
| `src/apsearch.c` | exhaustive $k$-AP enumerator (C, pthreads) |
| `src/ref_bruteforce.py` | assumption-free reference; validates it |
| `src/verify_records.py` | 26 records, four independent checks |
| `src/heuristic.py` | Hardy–Littlewood singular series |
| `src/compare_counts.py`, `src/calibration.py` | model vs measurement |
| `src/barrier.py` | $B(k,N)=2^kk!\,C_k(N)$, the formula behind the withdrawn sieve argument (NOTES §6.2) |
| `src/sieve_limits.py` | what Selberg's sieve, the large sieve and the larger sieve can and cannot prove (NOTES §6), and the tracker survey (§8.1) |
| `src/analysis.py` | the $L(N)$ tables |
| `src/final_check.py` | re-derives every claim; exit 0 = pass |
| `lean/` | Lean 4 + Mathlib formalisation; `lean/verify.sh`, exit 0 = pass |
| `lean-pnt/` | same bound without the PNT hypothesis, via `PrimeNumberTheoremAnd`; `lean-pnt/verify.sh`, exit 0 = pass |
| `data/`, `results/` | verified records, measured artifacts; [`results/README.md`](results/README.md) says which script writes each file |

## License

MIT, for this work only; the quoted problem is from [erdosproblems.com](https://www.erdosproblems.com/200). See [`LICENSE`](LICENSE).

One file is **not** MIT and is not relicensed here: [`lean-pnt/Erdos200/FCShim.lean`](lean-pnt/Erdos200/FCShim.lean) reproduces four short declarations character-for-character from [`google-deepmind/formal-conjectures`](https://github.com/google-deepmind/formal-conjectures), which is Apache-2.0, © 2024 Google LLC. They are copied only because `PrimeNumberTheoremAnd` and `formal-conjectures` pin incompatible Lean versions; the excerpts are marked with `BEGIN VERBATIM` / `END VERBATIM` and `lean-pnt/scripts/check_shim.py` fails if any line drifts from upstream.

Nothing in `refs/` is redistributed — it holds copyrighted papers, fetched by [`refs/fetch.sh`](refs/fetch.sh). Mathlib, `formal-conjectures` and `PrimeNumberTheoremAnd` are dependencies fetched by Lake, all Apache-2.0, none vendored here beyond the excerpt above.
