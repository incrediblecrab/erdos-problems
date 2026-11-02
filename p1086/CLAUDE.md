# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) and GitHub Copilot CLI when working in this repository.

This repository studies Erdős problem 1086, "many triangles of one area": $g(n)$, the most triangles of a single area that $n$ points in the plane can determine, known to satisfy $n^2\log\log n\ll g(n)\ll n^{20/9}$. It reproduces the published grid count, then extends the tractable exact and small-$n$ computations. The problem is open and nothing here solves it. The repository is public. `NOTES.md` is the write-up; `src/README.md` describes each program.

## Commands

```bash
python3 src/final_check.py      # rechecks every committed claim; exit 0 = pass, 3.5 to 7 minutes on a shared machine
python3 src/brute_check.py      # brute-force cross-check of the fast counters
bash refs/fetch.sh              # download sources; nothing is redistributed
```

Python is the shared `~/.venvs/main` interpreter. The counting kernels are C (`src/grid_count.c`, `src/brute_count.c`, `src/small_search.c`, `src/box_scan.c`) and need a C compiler; `src/README.md` says how each is built and run. The exact small-$n$ results are `src/g5_certificate.py` (a Gröbner certificate that no 5 points span 8 triangles of one area), `src/g6_exact.py` with `src/g6_brute.py` ($g(6)=12$, by an exact classification and a scan sharing no code with it), `src/g7_no21.py` (no 7 points span 21 equal-area triangles), `src/g7_exact.py` with `src/g7_box_check.py` ($g(7)=18$ and hence $g(8)=28$, by a finite classification of the extremal 5- and 6-point sets and a box scan sharing no code with it), and `src/g9_exact.py` with `src/g9_box_check.py` ($g(9)=36$, by the same classification carried to lower thresholds, with Gröbner bases for one six-point case, and a checker sharing no code with it). `src/lower_search.py` holds the searches beyond the small boxes.

## Rules

Reproduce before extending: a new count or bound is only stated after the published values it depends on have been reproduced by code here. A number in `README.md` or `NOTES.md` must be produced by a committed script and rechecked by `src/final_check.py`; regenerate outputs in `results/` with their script rather than editing them. Keep finite computations described as evidence about small $n$, never as deciding the asymptotic order. `refs/` holds only `refs/README.md` and `refs/fetch.sh`; never commit fetched papers.
