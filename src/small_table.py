"""Search for small-n lower bounds and write results/small_n.json.

Runs src/small_search (build: clang -O2 -o src/small_search src/small_search.c -lpthread -lm) exhaustively over the
5 x 5 box for n <= 10 and the 6 x 6 box for n <= 12, and by simulated annealing in the 10 x 10 box for n <= 16 with
several target areas and seeds, plus a longer pass in the 12 x 12 box for n >= 11. Keeps the best configuration per n, recounts it with recount.py, and adds the
averaging upper bound g(n) <= floor(g(n-1) n / (n-3)) started from g(7) = 18, which gives g(8) = 28.

    python src/small_table.py
"""
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from recount import classes

ROOT = Path(__file__).resolve().parent.parent
BIN = str(ROOT / "src" / "small_search")
EXHAUSTIVE = [(n, 5) for n in range(4, 11)] + [(n, 6) for n in range(4, 13)]
TARGETS = (1, 2, 3, 4, 6, 8, 12)
# (n, box, target D, iterations, seed); the second pass re-searches n >= 11 in a larger box with 4x the iterations
ANNEAL = ([(n, 10, D, 1_000_000, seed) for n in range(4, 17) for D in TARGETS for seed in (1, 2, 3)]
          + [(n, 12, D, 4_000_000, seed) for n in range(11, 17) for D in TARGETS for seed in (4, 5, 6)])
RESTARTS = 4
EXACT = {3: 1, 4: 4, 5: 7, 6: 12, 7: 18, 8: 28}


def run(args):
    out = subprocess.run([BIN, *map(str, args)], capture_output=True, text=True, check=True).stdout.split()
    assert out[0] == "best"
    pts = [tuple(map(int, p.split(","))) for p in out[3:]]
    return int(out[1]), int(out[2]), pts


def main():
    runs = []
    for n, k in EXHAUSTIVE:
        count, D, pts = run(["exhaustive", n, k, 10])
        runs.append({"n": n, "method": f"exhaustive {k}x{k}", "count": count, "D": D, "points": pts})
        print(f"exhaustive n={n} box {k}x{k}: {count}", flush=True)
    with ThreadPoolExecutor(10) as ex:
        res = list(ex.map(lambda a: run(["anneal", a[0], a[1], a[2], a[3], a[4], RESTARTS]), ANNEAL))
    for (n, k, D, iters, seed), (count, _, pts) in zip(ANNEAL, res):
        runs.append({"n": n, "method": f"anneal {k}x{k} {iters // 10**6}M D={D} seed={seed}", "count": count,
                     "D": D, "points": pts})
    for r in runs:  # every run, not only the winners, is recounted before it can be kept
        c = classes(r["points"])
        if len(set(r["points"])) != r["n"] or c[r["D"]] != r["count"]:
            sys.exit(f"recount disagrees with small_search: {r}")
    upper, best = dict(EXACT), []
    for n in range(4, 17):
        if n not in upper:
            upper[n] = upper[n - 1] * n // (n - 3)
        mine = [r for r in runs if r["n"] == n]
        top = max(r["count"] for r in mine)
        winner = min((r for r in mine if r["count"] == top), key=lambda r: (r["method"][0] != "e", r["D"]))
        best.append({**winner, "upper": upper[n], "exact": n in EXACT,
                     "methods_reaching_best": sorted({r["method"].split(" D=")[0] for r in mine if r["count"] == top}),
                     "best_by_method": {m: max(r["count"] for r in mine if r["method"].startswith(m))
                                        for m in sorted({r["method"].split(" D=")[0] for r in mine})}})
    out = {"exhaustive": [f"n={n} box {k}x{k}" for n, k in EXHAUSTIVE],
           "anneal": {"runs": len(ANNEAL), "restarts_per_run": RESTARTS, "targets": list(TARGETS),
                      "passes": ["n=4..16, box 10x10, 1M iterations, seeds 1-3",
                                 "n=11..16, box 12x12, 4M iterations, seeds 4-6"]},
           "best": best}
    (ROOT / "results" / "small_n.json").write_text(json.dumps(out, indent=1) + "\n")
    for b in best:
        print(f"n={b['n']:>2}  best {b['count']:>3} (D={b['D']})  upper {b['upper']:>3}  {b['best_by_method']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
