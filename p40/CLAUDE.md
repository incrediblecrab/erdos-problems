# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) and GitHub Copilot CLI when working in this repository.

This repository studies Erdős problem 40 (how dense a set must be to force unbounded representation counts $1_A\ast 1_A$). The problem is open and nothing here solves it; the results are Lean-checked reductions to the strong Erdős–Turán form and to $B_2[g]$ sets (problem 158 is $g=2$), a written reconstruction of the Erdős–Rényi exclusion with explicit constants, and a checked greedy $B_2[2]$ computation. The repository is public. `NOTES.md` is the full write-up; `README.md` is the summary and file reference.

## Commands

```bash
python3 src/final_check.py      # the gate: exit 0 = pass, about 10 s
bash scripts/reproduce.sh       # rerun every computation and the Lean check, then the gate
bash lean/verify.sh             # Lean build, statement check against formal-conjectures, axiom audit
python3 lean/plants.py          # plant defects and confirm verify.sh catches each
bash refs/fetch.sh              # download the sources the quotation checks need
```

Python is the shared `~/.venvs/main` interpreter (numpy, mpmath); `scripts/reproduce.sh` picks it and `$PY` overrides it. The C program `src/greedy_b2g.c` needs a C compiler. Lean, Mathlib and formal-conjectures are pinned in `lean/lean-toolchain`, `lean/lakefile.toml` and `lean/lake-manifest.json`; the statement check requires the formal-conjectures checkout to be that commit with no local edits.

## Rules the gate enforces

`src/final_check.py` rechecks the results and Lean records, the problem statement, the Lean excerpts quoted in the docs, every quotation against the fetched sources, and every numeral in the prose of the files in its `DOCS` list; it then plants a defect for each check and fails if any plant goes uncaught. So a number or quotation changed in `README.md`, `NOTES.md` or a folder README must match a fresh computation or source in the same change, and a new check needs a plant. It also reads the tracker database at `../tracker-and-solver/data/erdos.db` (outside this repository) when present.

Every claim carries one label, [Lean], [written], [source] or [computation], and the labels must stay honest: nothing is promoted to [Lean] without a machine check, and no finite computation is presented as deciding the problem. `refs/` is git-ignored apart from `refs/README.md` and `refs/fetch.sh`; never commit fetched papers or pages.
