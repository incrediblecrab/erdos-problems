# data/

Artifacts written by the programs in `src/`. `scripts/reproduce.sh` regenerates all of them, and `src/final_check.py` checks them (`NOTES.md` §10). `FAST=1` overwrites most of them with smaller runs, so run it in a copy of the repository if the committed files matter.

| file | written by | contents |
| --- | --- | --- |
| `lattice_n1.json` … `lattice_n8.json` | `src/subgroups` | the subgroups of $S_n$: the total and the number of conjugacy classes, histograms by order of all subgroups, of the abelian ones and of the classes, and the elementary abelian 2-subgroups by rank |
| `elemab_stats.json` | `asymptotics.py` | the rank statistics of the elementary abelian 2-subgroups for $n = 2, \dots, 64$ and every eighth $n$ up to 512, as doubles, and least-squares fits of the excess |
| `rank_counts_small.json` | `asymptotics.py` | the exact $a(n,k)$ for $n = 2, \dots, 40$, as decimal strings, on one line |
| `galois.json` | `asymptotics.py` | $\log_2 G_r$ and its excess over $r^2/4$ for $r = 2, \dots, 65$ and for $r = 2^p - 1, 2^p, 2^p + 1$ with $p = 7, 8, 9$ |
| `orbit_model.json` | `orbit_model.py` | orbit profiles, pure types, local coefficients and the held-out extrapolation |
| `sd_precision.json` | `sd_precision.py` | the exact rank sd at $n = 64, 128, 256, 512$ to 30 digits and its error against the closed form, computed at 60 decimal places |
| `census8.json` | `census8.py` | the 2-subgroups of $S_n$ for $n \le 8$: counts by order, the transitive types, and the Frattini-closed counts. The key `sylow_subgroups` is the number of subgroups of one Sylow 2-subgroup of $S_8$, not the number of Sylow subgroups |
| `block_families.json` | `block_families.py` | the §9 family: $\log_2 B(n)$ for $n \le 512$, its order statistics, $F_8$, $S_3(n)/B(n)$, the comparisons of §9.2 and §9.3, and the fractions of $f(n)$ in §9.8 for $n \le 18$ |

The statistics in `elemab_stats.json` are computed through floating-point logarithms and are much less precise than their doubles (`NOTES.md` §3); `sd_precision.json` has the sd computed exactly.
