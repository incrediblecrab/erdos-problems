# results/

Logs of the stages of `scripts/reproduce.sh`, committed as the record of the runs that produced the files in `data/`. Each stage overwrites its log.

| file | stage of `reproduce.sh` |
| --- | --- |
| `lattice_n8.log` | 3, `src/subgroups 8`; the logs for $n \le 7$ are not committed |
| `verify_elemab.log`, `verify_elemab.json` | 4, the brute-force cross-check, to $n = 10$ |
| `asymptotics_512.log` | 5, `asymptotics.py --n-max 512` |
| `orbit_model.log` | 6 |
| `sd_precision.log` | 6b |
| `census8.log`, `block_families.log` | 6c |
| `final_check.log` | 7, the gate; its last lines give the number of checks run, failed and skipped |

`fetch.log`, from stage 1, is not committed.
