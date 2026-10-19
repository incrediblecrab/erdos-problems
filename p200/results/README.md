# results

Measured artifacts behind [`../NOTES.md`](../NOTES.md), which [`../src/final_check.py`](../src/final_check.py) re-checks. The scripts are in [`../src/`](../src/), and [`../scripts/reproduce.sh`](../scripts/reproduce.sh) reruns them, except where the table says otherwise.

| file | written by | contents |
|---|---|---|
| `count_k8.json` … `count_k19.json` | `src/apsearch`, via `reproduce.sh` stage 5 or `src/run_counts.sh` | the number of $k$-APs of primes with last term up to a fixed bound, with a histogram by last term, for the Hardy–Littlewood count test |
| `counts.log` | `src/run_counts.sh`, which `reproduce.sh` does not run | apsearch's progress output for those runs |
| `k14.json` … `k19.json` | `src/apsearch`, run to the record $a(k)$; no committed script reruns them | exactly one $k$-AP with last term up to $a(k)$ in each, which re-derives $a(14)$ … $a(19)$ |
| `k20.json` | `src/apsearch`, via `reproduce.sh` stage 5 | the same for $a(20)=572945039351$ |
| `k18.log`, `k19.log`, `k20.log` | `src/apsearch`; not rewritten by any committed script | progress output of those runs |
| `reference.json` | `src/ref_bruteforce.py` | the assumption-free cross-check of apsearch |
| `count_comparison.json` | `src/compare_counts.py` | model against measured counts |
| `calibration.json` | `src/calibration.py` | a Kolmogorov–Smirnov test of the records against the model |
| `barrier.json` | `src/barrier.py` | $B(k,N)=2^kk!\,C_k(N)$, kept for NOTES §6.2 |
| `sieve_limits.json` | `src/sieve_limits.py` | the sieve limits of NOTES §6 and the tracker survey of §8.1 |
| `analysis.json` | `src/analysis.py` | the $L(N)$ tables |
| `thresholds.json` | no committed script | cached Hardy–Littlewood thresholds per $k$; `src/analysis.py` uses them if present and otherwise recomputes them with `heuristic.threshold_ak` |

The record table everything is checked against, `../data/records.json`, is written by `src/verify_records.py`.
