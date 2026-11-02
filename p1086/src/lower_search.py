"""Search beyond the boxes of src/small_table.py for lower bounds at n = 10, ..., 16; writes results/lower_search.json.

    python src/lower_search.py        about 40 minutes on 3 threads

Three methods, every kept configuration recounted with recount.py's exact integer count before it is written:

1. Exhaustive scans with src/box_scan of the 3 x 12, 4 x 9 and 7 x 7 boxes, at the threshold of the best count in
   results/small_n.json, so a run either lists configurations reaching that count or shows that none exceeds it.
   A triangular-lattice or affine-hexagon configuration is an affine image of a subset of Z^2, so these boxes and
   the square boxes already cover them; what the new boxes add is length (strips) and width (7 x 7).
2. Simulated annealing with src/small_search in the 16 x 16 and 20 x 20 boxes.
3. A box-free beam search in exact rationals: starting from the 15 configurations of C(6, 10) in
   results/g9_classes.json, each step adds every point lying on two non-parallel lines {X : det(q - p, X - p) = +-1}
   of the current set and keeps the BEAM best area-preserving affine classes by their number of triangles of doubled
   area 1. Coordinates are unbounded rationals, so this reaches configurations no fixed box contains.

None of this can show that a lower bound is optimal.
"""
import json
import math
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction as Fr
from pathlib import Path

from g7_exact import canonical, intersections, lines
from recount import classes

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"
THREADS = 3
STRIPS = [(10, 3, 12), (10, 4, 9), (10, 7, 7), (11, 3, 12), (11, 4, 9), (12, 3, 12), (12, 4, 9)]
ANNEAL = [(n, k, D, 2_000_000, seed) for n in range(10, 17) for k in (16, 20) for D in (1, 2, 3, 4, 6, 8, 12)
          for seed in (7, 8)]
BEAM, BEAM_TOP = 150, 12


def recounted(n, pts, method):
    """An entry for recount.py: the largest class of pts, at its smallest doubled area."""
    c = classes(pts)
    top = max(c.values())
    D = min(d for d, k in c.items() if k == top)
    if len(set(pts)) != n:
        raise RuntimeError(f"{method}: repeated points")
    return {"n": n, "method": method, "count": top, "D": D, "points": [list(p) for p in pts]}


def strip(n, a, b, threshold):
    r = subprocess.run([str(x) for x in (SRC / "box_scan", a, b, n, threshold, THREADS)], capture_output=True,
                       text=True, check=True)
    hits = [l.split() for l in r.stdout.splitlines() if l and not l.startswith("done")]
    subsets = int(r.stdout.split("done")[1])
    best = None
    for f in hits:
        pts = [tuple(map(int, s.split(","))) for s in f[2:]]
        e = recounted(n, pts, f"exhaustive {a}x{b}")
        if classes(pts)[int(f[0])] != int(f[1]):
            raise RuntimeError(f"box_scan count not reproduced: {f}")
        if best is None or e["count"] > best["count"]:
            best = e
    return {"method": f"exhaustive {a}x{b}", "n": n, "threshold": threshold, "subsets": subsets,
            "subsets_at_threshold": len(hits), "best": best["count"] if best else None}, best


def anneal(args):
    n, k, D, iters, seed = args
    out = subprocess.run([str(x) for x in (SRC / "small_search", "anneal", n, k, D, iters, seed, 4)],
                         capture_output=True, text=True, check=True).stdout.split()
    pts = [tuple(map(int, p.split(","))) for p in out[3:]]
    if classes(pts)[D] != int(out[1]):
        raise RuntimeError(f"small_search count not reproduced: {out}")
    return recounted(n, pts, f"anneal {k}x{k} {iters // 10**6}M D={D} seed={seed}")


def beam():
    g9 = json.loads((ROOT / "results" / "g9_classes.json").read_text())
    current = {tuple((Fr(x), Fr(y)) for x, y in e["points"]): e["edges"] for e in g9["six_points_at_least_10"]}
    out = []
    for n in range(7, BEAM_TOP + 1):
        scored = {}
        for base, e in current.items():
            ld = lines(base)
            for pt in intersections(ld):
                if pt in base:
                    continue
                deg = sum(m for (a, b, c), m in ld.items() if a * pt[0] + b * pt[1] == c)
                scored.setdefault(e + deg, []).append(base + (pt,))
        kept = {}
        for e in sorted(scored, reverse=True):
            for cfg in scored[e]:
                kept.setdefault(canonical(cfg), e)
            if len(kept) >= BEAM:
                break
        current = dict(sorted(kept.items(), key=lambda kv: (-kv[1], kv[0]))[:BEAM])
        top = max(current, key=lambda c: (current[c], c))
        den = 1
        for x, y in top:
            den = math.lcm(den, x.denominator, y.denominator)
        entry = recounted(n, [(int(x * den), int(y * den)) for x, y in top], f"beam {BEAM}")
        if entry["count"] < current[top]:
            raise RuntimeError("beam count not reproduced")
        out.append(entry)
        print(f"beam n={n}: best {current[top]}, {BEAM}th {min(current.values())}", flush=True)
    return out


def main():
    start = time.time()
    small = {e["n"]: e for e in json.loads((ROOT / "results" / "small_n.json").read_text())["best"]}
    runs, found = [], []
    for n, a, b in STRIPS:
        run, best = strip(n, a, b, small[n]["count"])
        runs.append(run)
        if best:
            found.append(best)
        print(f"{run}  ({time.time() - start:.0f} s)", flush=True)
    with ThreadPoolExecutor(THREADS) as ex:
        annealed = list(ex.map(anneal, ANNEAL))
    found += annealed
    runs += [{"method": e["method"], "n": e["n"], "best": e["count"]} for e in annealed]
    print(f"annealing: {len(annealed)} runs ({time.time() - start:.0f} s)", flush=True)
    beamed = beam()
    found += beamed
    runs += [{"method": e["method"], "n": e["n"], "best": e["count"]} for e in beamed]
    best = []
    for n in sorted({e["n"] for e in found}):
        top = max(found, key=lambda e: (e["n"] == n, e["count"]))
        best.append({**top, "previous": small[n]["count"], "improved": top["count"] > small[n]["count"]})
    out = {"strips": [f"n={n} box {a}x{b}" for n, a, b in STRIPS], "anneal_runs": len(ANNEAL), "beam_width": BEAM,
           "runs": runs, "best": best}
    (ROOT / "results" / "lower_search.json").write_text(json.dumps(out, indent=1) + "\n")
    for b in best:
        print(f"n={b['n']:>2}  best {b['count']:>3} (D={b['D']}, {b['method']})  previous {b['previous']}")
    print(f"elapsed {time.time() - start:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
