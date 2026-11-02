# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) and GitHub Copilot CLI when working in this repository.

This repository studies Erdős problem 1186, monochromatic $k$-term arithmetic progressions in 2-colourings, by certifying explicit upper-bound constructions. The problem is open and cannot be settled by finite computation; no lower-bound work was attempted. The repository is public. `NOTES.md` is the derivation and validation record, `README.md` the summary, and `PERIODIC-TRANSFER.md` the October 2026 note on the direct periodic transfer, whose asymptotic statement is conditional on an upstream cyclic-colouring theorem.

## Commands

```bash
python src/final_check.py         # re-derives every claim; exit 0 = pass, about 4 to 5 minutes (260 s with --threads 4 on November 2, 2026)
FAST=1 bash scripts/reproduce.sh  # full pipeline without the long searches
bash scripts/reproduce.sh         # everything, including the long searches
python src/check_periodic_lean.py # checks src/PeriodicTransfer.lean (core Lean, no Mathlib) with the shared axiom probe in ../formal, a sibling checkout; without it the gate skips this check and says so
bash refs/fetch.sh                # download sources; nothing is redistributed; without refs/bgl14.pdf the gate skips the BGL14 Table 1 quotation and says so
```

Python is the shared `~/.venvs/main` interpreter. The C cross-checks (`src/verify_zm.c`, `src/verify_two.c`, `src/verify_zp.c`, `src/verify_zpx.c`, `src/verify_base.c`) and `src/rabung.c` need a C compiler. `src/circlef.c` and `src/twof.c` are floating-point evaluators for search only, built on first use into git-ignored `src/lib*.so`; nothing they return is recorded. `final_check.py --threads N` sets the C programs' thread count (default 8).

## How a bound becomes a claim

`src/exact.py` is the certification path: every asserted value is evaluated in exact rational arithmetic, never floating point. For $\mathbb{F}_p$ words whose weight table exceeds `exact.PHI_TABLE_CELLS` (it is $O(m^2)$: 4.7 GB at $m=1500$, $k=4$), the certificate is `exact.phi_bits`, the same formula without the table, which the gate checks against `exact.phi` wherever both run. The colourings come from three families (`src/periodic.py`, `src/blocks.py`, `src/twoscale.py`) plus the searches (`src/bigm.py`, `src/wildcard.py`, the $\mathbb{F}_p$ breakpoint search `src/circleopt.py` and its exact stationary-point companion `src/circlestat.py`, the $\{1,\dots,n\}$ free-breakpoint search `src/twofree.py`, SAT-based `src/apfree.py` and `src/safesat.py`), and `src/chain.py` certifies one level of the recursion. The asserted bounds are listed in `results/claims.json`, which `src/final_check.py` re-derives together with the certified colourings in `results/`; each is also cross-checked by an independent C program. A new bound therefore needs a certified colouring in `results/`, an entry in `results/claims.json`, and a passing gate, and it must reproduce the published numbers it claims to beat first ([RL12, Table 1] via `src/rabung.c`, [BGL14, Table 1] in check B). Where a claim has `beats`, its `percent` field is the figure quoted in the prose, and the gate recomputes it exactly.
