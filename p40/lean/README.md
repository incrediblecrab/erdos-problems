# `lean/` — Lean proofs for Erdős #40

**Objective.** Machine-check Theorem A and the $B_2[g]$ reductions of [`../NOTES.md`](../NOTES.md) §2 and §3 against formal-conjectures' own definitions, using no axiom beyond `propext`, `Classical.choice` and `Quot.sound`. Neither `erdos_40` nor its weaker variant is proved or disproved here, and Theorem B is not formalised.

**Inputs.** formal-conjectures at the commit pinned in `lakefile.toml`, whose `FormalConjectures/ErdosProblems/40.lean` states the problem; that commit fixes Mathlib. The Lean release is in `lean-toolchain`.

**Run.** `./verify.sh`; exit 0 = pass. `./plants.py` plants defects to test it. `../NOTES.md` §8 says what each step checks and records the last run of both.

| file | contents |
|---|---|
| `Erdos40.lean` | the library root |
| `Erdos40/` | the Lean sources; its README describes each file |
| `Probe.lean` | the axiom and statement probe, run in a separate process |
| `verify.sh` | the checks |
| `plants.py` | the planted-defect harness; writes `../results/lean_plants.json` and `.log` |
| `lakefile.toml`, `lake-manifest.json`, `lean-toolchain` | the pinned formal-conjectures commit, every package revision, and the Lean release |
