"""Structured constructions for n = 10, ..., 16, counted exactly; writes results/constructions.json.

    python src/constructions.py        about 30 seconds on 3 processes

Two families, each scanned exhaustively over the parameters below, every configuration counted at every doubled
area with recount.py's integer count:

1. Rows. r = 2, 3 or 4 parallel rows of consecutive integers at heights 0, ..., r - 1 (or 0, 1, 3 and 0, 1, 3, 4),
   with lengths summing to n and integer offsets 0, 1 or 2 between rows. These are the unions of parallel lines with
   arithmetic structure; the [ErPu71] construction is the case of floor(sqrt(log n)) rows, which is a single,
   collinear row for every n here (sqrt(log 16) < 2), so it gives nothing at this size.
2. Lattice sections. The n points of Z^2 smallest under a positive definite form a(x-c)^2 + b(x-c)(y-d) + c'(y-d)^2
   with integer 1 <= a, c' <= 5, |b| <= 5 and centre (c, d) in {0, 1/4, 1/2, 3/4}^2, skipping cuts that are not
   unique. Triangular and hexagonal patches are affine images of these, and affine maps preserve ratios of areas.

Any configuration with rational coordinates is an affine image of a subset of Z^2, so both families lie inside the
boxes and the beam search of lower_search.py; this scan adds their shapes, not new territory. It shows nothing about
optimality.
"""
import itertools
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from pathlib import Path

from recount import classes

ROOT = Path(__file__).resolve().parent.parent
PROCESSES = 3
HEIGHTS = {2: [(0, 1)], 3: [(0, 1, 2), (0, 1, 3)], 4: [(0, 1, 2, 3), (0, 1, 3, 4)]}


def top(pts):
    c = classes(pts)
    best = max(c.values())
    return best, min(d for d, k in c.items() if k == best)


def rows(n):
    best = None
    for r, hs in HEIGHTS.items():
        for lengths in itertools.product(range(1, n), repeat=r):
            if sum(lengths) != n:
                continue
            for h in hs:
                for offs in itertools.product(range(3), repeat=r - 1):
                    o = (0,) + offs
                    pts = [(o[i] + k, h[i]) for i in range(r) for k in range(lengths[i])]
                    count, D = top(pts)
                    if best is None or count > best[0]:
                        best = (count, D, pts, f"rows r={r} heights={h} lengths={lengths} offsets={o}")
    return best


def sections(n):
    best = None
    grid = [(x, y) for x in range(-8, 9) for y in range(-8, 9)]
    centres = [(Fraction(i, 4), Fraction(j, 4)) for i in range(4) for j in range(4)]
    for a, b, c in itertools.product(range(1, 6), range(-5, 6), range(1, 6)):
        if 4 * a * c - b * b <= 0:
            continue
        for cx, cy in centres:
            q = lambda p: a * (p[0] - cx) ** 2 + b * (p[0] - cx) * (p[1] - cy) + c * (p[1] - cy) ** 2
            pts = sorted(grid, key=q)
            if q(pts[n - 1]) == q(pts[n]):
                continue
            count, D = top(pts[:n])
            if best is None or count > best[0]:
                best = (count, D, pts[:n], f"section form=({a},{b},{c}) centre=({cx},{cy})")
    return best


def job(args):
    family, n = args
    return family, n, (rows if family == "rows" else sections)(n)


def main():
    small = {e["n"]: e["count"] for e in json.loads((ROOT / "results" / "small_n.json").read_text())["best"]}
    tasks = [(f, n) for n in range(10, 17) for f in ("rows", "sections")]
    with ProcessPoolExecutor(PROCESSES) as ex:
        res = list(ex.map(job, tasks))
    best = []
    for family, n, (count, D, pts, method) in res:
        best.append({"n": n, "method": method, "count": count, "D": D, "points": [list(p) for p in pts],
                     "table": small[n], "improved": count > small[n]})
        print(f"n={n:>2} {family:>8}: {count:>3} at D={D} (table {small[n]})  {method}", flush=True)
    (ROOT / "results" / "constructions.json").write_text(json.dumps({"best": best}, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
