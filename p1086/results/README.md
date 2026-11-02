# results/

**Objective.** Hold the values this work asserts and the data they are checked against. **Inputs.** Everything except `claims.json` is written by the scripts in `src/`; `claims.json` is authored, and `src/final_check.py` re-derives every value in it from the other files.

| file | contents |
|---|---|
| `claims.json` | the asserted values: exact $g(3),\dots,g(10)$, bounds for $n=11,\dots,16$, four grid lower bounds, the reproduced [DST09] count, and Theorems A and B, marked informal |
| `g7_classes.json` | written by `src/g7_exact.py --write`: every 5-point set with $\ge6$ and every 6-point set with $\ge11$ triangles of one area, up to area-preserving affine maps, and the (empty) list of 7-point sets with $\ge19$ |
| `g9_classes.json` | written by `src/g9_exact.py --write`: every 6-point set with $\ge10$, 7-point set with $\ge16$ and 8-point set with $\ge25$ triangles of one area, up to area-preserving affine maps, the Type II counts, and the (empty) list of 9-point sets with $\ge37$ |
| `g10_classes.json` | written by `src/g10_exact.py --write`: every 8-point set with $\ge24$ and 9-point set with $\ge35$ triangles of one area, up to area-preserving affine maps, the (empty) list of 10-point sets with $\ge49$, the largest line multiplicities, the number of 9-regular 8-point sets (0), the number of 10-point candidates, and the largest 10-point extension, 48, with the one configuration attaining it |
| `g10_lrat.json` | written by `src/g10_certify.py`: for each of the lemma's 22 CNFs, its SHA-256, the solver (CaDiCaL, or kissat with drat-trim for the gate's CNF) with its arguments and exit code, the proof's size, and `lrat-check`'s exit code and output; the drat-trim commit; and the corrupted-proof control |
| `constructions.json` | written by `src/constructions.py`: the best row and lattice-section construction for each $n=10,\dots,16$ |
| `lattice.json` | the 16 grid runs: size, SHA-256 of the raw histogram, total, best doubled area and its count, top five areas, fixed-area counts |
| `lower_search.json` | written by `src/lower_search.py`: every strip, annealing and beam run for $n=10,\dots,16$ and the best configuration each method found |
| `small_n.json` | the exhaustive and annealing runs for $n\le16$, and the best configuration for each $n$ with its points, doubled area and count |
| `raw/` | the full histograms; see [`raw/README.md`](raw/README.md) |
| `brute_count.log` | one-off runs in which `src/brute_count` reproduced `raw/grid_200x200.txt` and `raw/grid_3x3333.txt` byte for byte; not gate-checked |
