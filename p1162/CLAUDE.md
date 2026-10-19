# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) and GitHub Copilot CLI when working in this repository.

This repository studies Erdős problem 1162, the number of subgroups of $S_n$ and their orders, with exact computations (the rank-stratified elementary-abelian family to $n=512$, full lattices for $n\le 8$), closed-form limiting constants, and a precise conjecture for the vague second half of the question. The problem is open and nothing here solves it. The repository is public. `NOTES.md` is the write-up; `README.md` has the file reference.

## Commands

```bash
python3 src/final_check.py      # checks artifacts and every decimal in NOTES.md against fresh computation; exit 0 = pass, about 1 minute
bash scripts/reproduce.sh       # regenerate every artifact, then the gate
bash refs/fetch.sh              # sources for the quotation checks; nothing is redistributed
```

`scripts/reproduce.sh` uses `$PY` if set, else `~/.venvs/main/bin/python`, else `python3`; it needs numpy and mpmath, a C99 compiler for `src/subgroups.c` (NEON intrinsics on arm64, portable fallback elsewhere), curl, and poppler's `pdftotext`. Without the fetched sources the quotation checks report "skip", not "pass", and the gate keeps the two apart.

## Rules

All logarithms are base 2. `src/final_check.py` recomputes every decimal number in `NOTES.md` except timings and values explicitly marked as first-version, so edit a number in the notes only together with the code that produces it. Each exact computation has an independent check that shares no code with it (`src/verify_elemab.py` for `src/elemab.py`, `src/census8.py` by brute force on permutations); keep that separation when adding results. `data/` and `results/` each have a README listing their files and the script that writes them; regenerate rather than hand-edit. Conjectural asymptotics stay labelled as conjectures.
