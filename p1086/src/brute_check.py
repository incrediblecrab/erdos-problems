"""Independent check of grid_count: enumerate every triple of the a x b grid and histogram |det|.

Shares no code or method with grid_count.c (which sums over triangle shapes). Usage:
    python brute_check.py A B            print the histogram in grid_count's format
    python brute_check.py --compare A B FILE   exit 1 unless FILE (grid_count output) matches exactly
"""
import sys
from itertools import combinations

import numpy as np


def histogram(a, b):
    pts = np.array([(x, y) for x in range(a) for y in range(b)], dtype=np.int64)
    n = len(pts)
    hist = {}
    for i, j in combinations(range(n), 2):
        k = np.arange(j + 1, n)
        if k.size == 0:
            continue
        ux, uy = pts[j] - pts[i]
        w = pts[k] - pts[i]
        d = np.abs(ux * w[:, 1] - uy * w[:, 0])
        d = d[d > 0]
        vals, cnts = np.unique(d, return_counts=True)
        for v, c in zip(vals.tolist(), cnts.tolist()):
            hist[v] = hist.get(v, 0) + c
    return hist


def read_grid_count(path):
    out = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            d, c = line.split()
            out[int(d)] = int(c)
    return out


def main(argv):
    if argv[:1] == ["--compare"]:
        a, b, path = int(argv[1]), int(argv[2]), argv[3]
        want = histogram(a, b)
        got = read_grid_count(path)
        bad = sorted(set(want) | set(got))
        bad = [d for d in bad if want.get(d, 0) != got.get(d, 0)]
        total = sum(want.values())
        if bad:
            for d in bad[:10]:
                print(f"MISMATCH a={a} b={b} D={d}: brute={want.get(d, 0)} grid_count={got.get(d, 0)}")
            return 1
        print(f"OK a={a} b={b}: {len(want)} areas, {total} triangles agree")
        return 0
    a, b = int(argv[0]), int(argv[1])
    h = histogram(a, b)
    print(f"# a={a} b={b} n={a * b} brute")
    for d in sorted(h):
        print(d, h[d])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
