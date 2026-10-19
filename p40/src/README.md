# `src/` — the computations and the gate

**Objective.** The scripts behind [`../NOTES.md`](../NOTES.md) §4, §6 and §8. Each exits 1 when one of its checks fails, and `--plant` gives each a defect it must catch. `../scripts/reproduce.sh` runs them all.

**Inputs.** `greedy_check.py` compiles `greedy_b2g.c` with the system C compiler and compares its output for $g=1$ with the OEIS b-file in `../refs/`. `final_check.py` reads `../results/`, the documentation, the sources in `../lean/` and, when they are present, the files in `../refs/`.

| file | contents |
|---|---|
| `greedy_b2g.c` | the greedy $B_2[g]$ set, which keeps a table of representation counts and tries each candidate against it (§6) |
| `greedy_check.py` | runs that program under a CPU-time limit and checks its output without sharing its code; writes `../results/greedy_b2g.json` and `../results/greedy_b2g_g2.txt` |
| `theorem_b_check.py` | the numerical check of the constants in Theorem B; writes `../results/theorem_b_check.json` (§4) |
| `final_check.py` | the gate of §8. `--write` regenerates the tables in the documentation, `--strict` counts any skip as a failure, and `--plant NAME` runs one planted defect, which must exit 1 |
