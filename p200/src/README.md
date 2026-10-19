# src

The computations behind [`../NOTES.md`](../NOTES.md). The table in [`../README.md`](../README.md) says what each script is for, [`../results/README.md`](../results/README.md) says which file each one writes, and [`../scripts/reproduce.sh`](../scripts/reproduce.sh) runs them all in order.

The Python scripts use the shared environment at `~/.venvs/main/bin/python` (numpy, sympy) and find the repository root from their own path, so they can be run from anywhere, e.g. `~/.venvs/main/bin/python src/sieve_limits.py`. `final_check.py` needs only the standard library, so `python3 src/final_check.py` works too; it writes nothing and exits nonzero if any check fails.

`apsearch.c` is the exhaustive enumerator. Build it with `cc -O3 -march=native -funroll-loops -pthread -o src/apsearch src/apsearch.c`; the binary is gitignored. `run_counts.sh` reruns the fixed-bound counts of `reproduce.sh` stage 5 on their own.
