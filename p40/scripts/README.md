# `scripts/` — reproduction

**Objective.** `reproduce.sh` reruns every computation behind [`../NOTES.md`](../NOTES.md), the Lean checks and their planted defects, and last the gate, `../src/final_check.py`. It stops at the first failure.

**Inputs.** Its header lists the stages and the tools they need: bash, curl, a C compiler, poppler's `pdftotext`, elan, and a Python 3 with numpy and mpmath, which `$PY` selects.

```bash
bash scripts/reproduce.sh              # every stage
bash scripts/reproduce.sh greedy final # only the named stages
```

Each stage overwrites what it writes: `fetch` the downloads in `../refs/`, `final` the generated tables in the documentation, and the others their files in `../results/`. Run it in a copy if the committed outputs matter. The `plants` stage is the slow one. When `../refs/` is empty the gate reports the checks that need it as skipped, and a skip is not a pass: such a run ends `PASS, incomplete`. `final` runs the gate with `--strict`, which fails on any skip, after a successful `fetch` stage in the same run, or when `STRICT=--strict` is set because `../refs/` was filled earlier.

| file | contents |
|---|---|
| `reproduce.sh` | the stages `fetch`, `thmb`, `greedy`, `lean`, `plants` and `final` |
