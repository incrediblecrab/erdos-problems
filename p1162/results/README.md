# results/

Logs of the stages of `scripts/reproduce.sh`, committed as the record of the runs that produced the files in `data/`. Each stage overwrites its log.

| file | stage of `reproduce.sh` |
| --- | --- |
| `lattice_n8.log` | 3, `src/subgroups 8`; the logs for $n \le 7$ are not committed |
| `verify_elemab.log`, `verify_elemab.json` | 4, the brute-force cross-check, to $n = 10$ |
| `asymptotics_512.log` | 5, `asymptotics.py --n-max 512` |
| `orbit_model.log` | 6 |
| `sd_precision.log` | 6b |
| `census8.log`, `nilpotent8.log`, `block_families.log` | 6c |
| `two_subgroups.log` | 6d, `two_subgroups.py` to $n = 15$ |
| `verify_two_subgroups.log`, `verify_two_subgroups.json` | 6d, `verify_two_subgroups.py --n-max 15 --dedup-max 12`, the independent recount |
| `two_subgroups16.log` | 6e, `two_subgroups16.py`, which runs `src/census16` on the orbit types of 16 points |
| `verify_two_subgroups16.log`, `verify_two_subgroups16.json` | 6e, `verify_two_subgroups16.py`, the independent recount of $n = 16$ and 17; not rerun under `FAST=1` |
| `final_check.log` | 7, the gate; its last lines give the number of checks run, failed and skipped |

`fetch.log`, from stage 1, is not committed.
