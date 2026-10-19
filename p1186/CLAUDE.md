# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) and GitHub Copilot CLI when working in this repository.

This repository studies Erdős problem 1186, monochromatic $k$-term arithmetic progressions in 2-colourings, by certifying explicit upper-bound constructions. The problem is open and cannot be settled by finite computation; no lower-bound work was attempted. The repository is public. `NOTES.md` is the derivation and validation record, `README.md` the summary, and `PERIODIC-TRANSFER.md` the October 2026 note on the direct periodic transfer, whose asymptotic statement is conditional on an upstream cyclic-colouring theorem.

## Commands

```bash
python src/final_check.py         # re-derives every claim; exit 0 = pass, about 2.5 minutes
FAST=1 bash scripts/reproduce.sh  # full pipeline without the long searches
bash scripts/reproduce.sh         # everything, including the long searches
python src/check_periodic_lean.py # checks src/PeriodicTransfer.lean (core Lean, no Mathlib) with the shared axiom probe in ../formal, a sibling checkout; without it the gate skips this check and says so
bash refs/fetch.sh                # download sources; nothing is redistributed
```

Python is the shared `~/.venvs/main` interpreter. The C cross-checks (`src/verify_zm.c`, `src/verify_two.c`, `src/verify_zp.c`, `src/verify_base.c`) and `src/rabung.c` need a C compiler.

## How a bound becomes a claim

`src/exact.py` is the certification path: every asserted value is evaluated in exact rational arithmetic, never floating point. The colourings come from three families (`src/periodic.py`, `src/blocks.py`, `src/twoscale.py`) plus the searches (`src/bigm.py`, `src/wildcard.py`, SAT-based `src/apfree.py` and `src/safesat.py`), and `src/chain.py` certifies one level of the recursion. The asserted bounds are listed in `results/claims.json`, which `src/final_check.py` re-derives together with the certified colourings in `results/`; each is also cross-checked by an independent C program. A new bound therefore needs a certified colouring in `results/`, an entry in `results/claims.json`, and a passing gate, and it must reproduce the published numbers it claims to beat ([RL12, Table 1] via `src/rabung.c`) first.
