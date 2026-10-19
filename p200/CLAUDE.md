# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) and GitHub Copilot CLI when working in this repository.

This repository studies Erdős problem 200, the longest arithmetic progression of primes in $\{1,\ldots,N\}$, written $L(N)$, through exhaustive search, a Hardy–Littlewood heuristic compared with measured counts, sieve-limit analysis and a Lean 4 bound. The problem is open and nothing here solves it. The repository is public. `NOTES.md` is the full write-up; `README.md` has the file reference.

## Commands

```bash
python3 src/final_check.py      # re-derives every claim; exit 0 = pass
bash scripts/reproduce.sh       # every stage, then the gate; the script header lists the stages
bash lean/verify.sh             # Lean + Mathlib bound with a PNT hypothesis; exit 0 = pass
bash lean-pnt/verify.sh         # the same bound without the hypothesis, via PrimeNumberTheoremAnd
bash refs/fetch.sh              # download sources; nothing is redistributed
```

Python is the shared `~/.venvs/main` interpreter. `src/apsearch.c` is the exhaustive $k$-AP enumerator (C with pthreads); `src/ref_bruteforce.py` is the assumption-free reference that validates it, and `src/verify_records.py` checks the 26 records four independent ways. `lean/` and `lean-pnt/` are separate Lake projects with their own pinned toolchains and manifests.

## Rules

`results/README.md` names the script that writes each file in `data/` and `results/`; regenerate a result with that script, never by editing the JSON. `src/final_check.py` re-derives the claims from those artifacts, so a number changed in `NOTES.md` must come with the recomputation that produces it. `src/sieve_limits.py` also reads the tracker database at `../tracker-and-solver/data/erdos.db` (outside this repository), read-only, when it is present. The sieve argument built on `src/barrier.py` was withdrawn (NOTES §6.2); keep it marked as withdrawn rather than reviving it. Keep every heuristic or computational statement labelled as such, and never present a finite search as settling the asymptotic question.
