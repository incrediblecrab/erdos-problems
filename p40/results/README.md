# `results/` — the committed outputs

**Objective.** The outputs that [`../NOTES.md`](../NOTES.md) reports, committed so that `../src/final_check.py` can check the documentation against them without rerunning anything.

**Inputs.** The stages of `../scripts/reproduce.sh`, each of which overwrites its files.

| file | stage | contents |
|---|---|---|
| `greedy_b2g.json` | `greedy` | the checkpoints of §6, the posted numbers beside the computed ones, the checks run, and the hashes of the program, checker and list |
| `greedy_b2g_g2.txt` | `greedy` | the greedy $B_2[2]$ set as far as it was computed, one term per line |
| `theorem_b_check.json` | `thmb` | the tables of §4 |
| `lean_verify.log` | `lean` | the output of `lean/verify.sh`, headed by the date, the Lean release and the platform, and followed by its time and memory and the hashes of the files it checked |
| `lean_plants.json`, `lean_plants.log` | `plants` | each planted-defect case of §8, its verdict and its full output |
| `final_check.log` | `final` | the gate's report; its last line counts the checks passed, failed and skipped and the plants caught and skipped |

`fetch.log`, from the `fetch` stage, lists the hashes of the downloads and is not committed.
