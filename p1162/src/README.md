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
| `nilpotent8.py` | the 2-subgroup counts of `census8.py` against OEIS A218939: a nilpotent subgroup is $P \times Q$ with $P$ a 2-subgroup of the centraliser of its odd part $Q$, and every $Q \ne 1$ and the 2-subgroups of its centraliser are found by brute force on permutations, $n \le 8$; writes `data/nilpotent8.json`; exits 1 on a mismatch | §9.8 |
| `block_families.py` | the Frattini-closed family: $B(n)$ and its order statistics to $n = 512$, the count $F_8(n)$ in the proof of Theorem 9.1, the family $S_3(n)$ for odd $n$, and the comparisons of §9.2, §9.3 and §9.8; exits 1 if a gate fails | §9 |
| `two_subgroups.py` | every 2-subgroup of $S_n$ for $n \le 15$, Frattini-closed or not: the subgroups of a Sylow 2-subgroup of each orbit-type stabiliser that are transitive on every orbit, weighted by the number of Sylow subgroups containing them, read off as a count of maximal chains; writes `data/two_subgroups.json`; exits 1 if a check fails | §9.8 |
| `verify_two_subgroups.py` | the same counts on explicit permutations, with the Sylow subgroups tested one by one and, for $n \le 12$, by conjugating and counting distinct subgroups; it reads `data/two_subgroups.json` only to compare the final numbers and writes `results/verify_two_subgroups.json` | §9.8 |
| `census16.c` | the subgroups of the Sylow subgroup of one orbit type of up to 16 points, blocks of size 2 to 16, transitive on every block, with $d$ and the Sylow numbers $s_b$ of their projections from nests of blocks; `reproduce.sh` compiles it to `src/census16`, which is not committed; exits 2 if an invariant fails | §9.8 |
| `two_subgroups16.py` | every 2-subgroup of $S_{16}$ and $S_{17}$: runs `census16` on the ten orbit types of 16 points, weights them and adds the stored $t(k)$, $k \le 14$; writes `data/two_subgroups16.json`; exits 1 if a check fails | §9.8 |
| `verify_two_subgroups16.py` | the same counts without `census16.c`: the orbit types with orbits of at most 8 points by `count_type()` of `two_subgroups.py`, the type (16) by extending each subgroup of type $(8, 8)$ by an element that swaps the halves, with $s_{16}$ from the literal definition of a block, and the regular subgroups of $S_{16}$ without weights; writes `results/verify_two_subgroups16.json` | §9.8 |
| `final_check.py` | the gate: checks the artifacts and compares the decimal numbers in `NOTES.md` with fresh values; exit 0 = pass, and `--coverage` lists the numbers it does not compare | §10 |

`final_check.py` imports `asymptotics.py`, `elemab.py`, `block_families.py` and `orbit_model.py`, so a change to any of them changes what it checks. It reruns `census8.py`, `nilpotent8.py`, `block_families.py`, `two_subgroups.py --n-max 12` and `census16.c` on the orbit type (16), which it compiles with `cc` in a temporary directory, and reads the stored output of `verify_two_subgroups.py` and `verify_two_subgroups16.py` without rerunning them.
