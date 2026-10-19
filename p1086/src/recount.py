"""Independent recount of the small-n configurations: exact integer areas, no code shared with small_search.c.

    python src/recount.py results/small_n.json    exits 1 unless every stored count and area is reproduced
"""
import json
import sys
from collections import Counter
from itertools import combinations


def classes(pts):
    """Doubled area -> number of triangles, over the nondegenerate triangles of pts."""
    c = Counter()
    for (ax, ay), (bx, by), (cx, cy) in combinations(pts, 3):
        d = abs((bx - ax) * (cy - ay) - (by - ay) * (cx - ax))
        if d:
            c[d] += 1
    return c


def problem(entry):
    pts = [tuple(p) for p in entry["points"]]
    if len(pts) != entry["n"] or len(set(pts)) != len(pts):
        return "wrong number of points, or repeated points"
    c = classes(pts)
    if c[entry["D"]] != entry["count"] or max(c.values()) != entry["count"]:
        return f"stored {entry['count']} at D={entry['D']}, recount {c[entry['D']]}, largest class {max(c.values())}"
    return None


def main(path):
    entries = json.load(open(path))["best"]
    bad = [(e["n"], p) for e in entries if (p := problem(e))]
    for n, p in bad:
        print(f"n={n}: {p}")
    print(f"{len(entries)} configurations recounted, {len(bad)} mismatches")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
