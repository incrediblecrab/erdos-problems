# results/raw/

Full doubled-area histograms, one file per grid: 14 files, about 20 MB. They are git-ignored except for this README, and `python src/lattice_table.py` regenerates them.

Each `grid_AxB.txt` is `src/grid_count`'s output for the $A\times B$ integer grid: a `# a=A b=B n=AB dmax=DMAX` header, then one `D count` line for every doubled area $D$ with a nonzero count. `results/lattice.json` stores each file's SHA-256, and `src/verify_lattice.py` rechecks the files against it.
