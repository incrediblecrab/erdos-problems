# scripts/

`reproduce.sh` regenerates every artifact in `data/` and `results/`, then runs `src/final_check.py` and exits with its status. Its header lists the stages, the tools they need, and where `FAST=1` stops them.

```sh
bash scripts/reproduce.sh          # the full run
FAST=1 bash scripts/reproduce.sh   # the same stages at smaller n
```

Both overwrite `data/` and `results/`; run `FAST=1` in a copy of the repository if the committed artifacts matter. The timings are in `NOTES.md` §10.
