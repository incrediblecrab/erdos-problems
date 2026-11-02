# src/

**Objective.** Compute and check every number in `results/claims.json`. **Inputs.** None beyond the problem statement; no script reads `refs/`. Python scripts run under `~/.venvs/main/bin/python` (numpy, sympy, z3) from any directory.

Build the C programs first (no warnings under Apple clang 21):

```bash
cd src && for p in grid_count brute_count small_search box_scan; do clang -O2 -Wall -Wextra -o $p $p.c -lpthread -lm; done
```

| file | what it does |
|---|---|
| `grid_count.c` | exact doubled-area histogram of an $a\times b$ grid, summing over triangle shapes |
| `brute_count.c`, `brute_check.py` | independent triple enumerations in C and Python |
| `lattice_table.py` | runs all 16 grids (about 13 min), writing `results/lattice.json` and `results/raw/` |
| `verify_lattice.py` | rechecks `lattice.json` against the raw histograms |
| `g5_certificate.py` | Gröbner certificate that no 5 points span 8 triangles of one area |
| `g6_exact.py`, `g6_brute.py` | $g(6)=12$: exact classification, and a scan sharing no code with it |
| `g7_no21.py` | rules out 21 equal-area triangles on 7 points by enumerating the surviving triple systems and checking their polynomial systems |
| `g7_exact.py` | $g(7)=18$: classifies the 5-point sets with $\ge6$ and the 6-point sets with $\ge11$ equal-area triangles by a finite enumeration, then shows no 7-point set has 19; `--write` stores the classes in `results/g7_classes.json` |
| `g7_box_check.py` | cross-checks that classification against an exhaustive scan of the $6\times6$ box, sharing no code with it |
| `g9_exact.py` | $g(9)=36$: classifies the 6-point sets with $\ge10$, 7-point sets with $\ge16$ and 8-point sets with $\ge25$ equal-area triangles (Gröbner bases for the Type II six-point case) and shows no 9-point set has 37; `--write` stores the classes in `results/g9_classes.json` |
| `box_scan.c`, `g9_box_check.py` | cross-check of `g9_exact.py` sharing no code with it: $6\times6$ box scans, the extension steps redone in integer arithmetic, and z3 on the Type II hypergraphs |
| `small_search.c`, `small_table.py` | exhaustive and annealing search for $n\le16$ (about 7 min on 3 threads), writing `results/small_n.json` |
| `lower_search.py` | strip, $7\times7$, large-box annealing and box-free exact beam searches for $n=10,\dots,16$ (about 40 min on 3 threads), writing `results/lower_search.json` |
| `recount.py` | recounts every stored small-$n$ configuration |
| `final_check.py` | re-derives every claim, planted defects included; exit 0 = pass (3.5 to 7 min on a shared machine) |
