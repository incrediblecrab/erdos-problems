# results/

**Objective.** Hold the values this work asserts and the data they are checked against. **Inputs.** Everything except `claims.json` is written by the scripts in `src/`; `claims.json` is authored, and `src/final_check.py` re-derives every value in it from the other files.

| file | contents |
|---|---|
| `claims.json` | the asserted values: exact $g(3),\dots,g(8)$, bounds for $n=9,\dots,16$, four grid lower bounds, the reproduced [DST09] count, and Theorems A and B, marked informal |
| `g7_classes.json` | written by `src/g7_exact.py --write`: every 5-point set with $\ge6$ and every 6-point set with $\ge11$ triangles of one area, up to area-preserving affine maps, and the (empty) list of 7-point sets with $\ge19$ |
| `lattice.json` | the 16 grid runs: size, SHA-256 of the raw histogram, total, best doubled area and its count, top five areas, fixed-area counts |
| `small_n.json` | the exhaustive and annealing runs for $n\le16$, and the best configuration for each $n$ with its points, doubled area and count |
| `raw/` | the full histograms; see [`raw/README.md`](raw/README.md) |
| `brute_count.log` | one-off runs in which `src/brute_count` reproduced `raw/grid_200x200.txt` and `raw/grid_3x3333.txt` byte for byte; not gate-checked |
