# src/

The programs behind `NOTES.md`. `scripts/reproduce.sh` runs them in order, and its header says which stage writes which file; `NOTES.md` §10 has the timings.

| file | what it does | `NOTES.md` |
| --- | --- | --- |
| `subgroups.c` | every subgroup of $S_n$, found class by class and then expanded to the whole class. `reproduce.sh` compiles it to `src/subgroups`, which is not committed, and runs it for $n \le 8$ into `data/lattice_n*.json` | §2.1 |
| `elemab.py` | the exact number $a(n,k)$ of elementary abelian 2-subgroups of $S_n$ of rank $k$, by the exponential formula and $q$-binomial inversion, and the Galois numbers $G_r$ | §2.2 |
| `verify_elemab.py` | the same counts by building the subgroups out of explicit permutations; it imports `elemab.py` only to compare the final numbers | §2.3 |
| `asymptotics.py` | the rank statistics, the fits of the excess and the Galois data, written to `data/elemab_stats.json`, `data/rank_counts_small.json` and `data/galois.json` | §2.2–§4 |
| `sd_precision.py` | the rank sd recomputed from the exact counts at 60 decimal places and compared with its closed form | §3 |
| `orbit_model.py` | the orbit-structure model of the second-order term, with its held-out extrapolation | §6 |
| `census8.py` | every 2-subgroup of $S_n$ for $n \le 8$, by brute force on permutations, with the transitive types and the Frattini-closed subgroups; exits 1 if a check fails | §9.1 |
| `block_families.py` | the Frattini-closed family: $B(n)$ and its order statistics to $n = 512$, the count $F_8(n)$ in the proof of Theorem 9.1, the family $S_3(n)$ for odd $n$, and the comparisons of §9.2, §9.3 and §9.8; exits 1 if a gate fails | §9 |
| `final_check.py` | the gate: checks the artifacts and compares the decimal numbers in `NOTES.md` with fresh values; exit 0 = pass, and `--coverage` lists the numbers it does not compare | §10 |

`final_check.py` imports `asymptotics.py`, `elemab.py`, `block_families.py` and `orbit_model.py`, so a change to any of them changes what it checks.
