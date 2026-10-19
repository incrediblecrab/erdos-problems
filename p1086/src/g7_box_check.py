"""Cross-check src/g7_exact.py against an exhaustive scan of a lattice box, sharing no code with it.

    python src/g7_box_check.py [results/g7_classes.json]      exit 0 iff every check passes

Three checks, all exact integer arithmetic:

1. The four-point lemma, on every 4-subset of the 7 x 7 box and every doubled area: at most 4 triangles of one
   area; 4 only for a parallelogram; exactly 3 only for a triangle with its centroid.
2. Every 5-subset of the 6 x 6 box with at least 6 triangles of one area, and every 6-subset with at least 11, has
   the area profile of a configuration listed in results/g7_classes.json. The profile is the multiset of all
   doubled areas divided by the common area, an invariant of area-preserving affine maps, so a configuration the
   enumeration missed would show up as an unlisted profile unless it happened to share a profile with a listed one.
3. The listed configurations themselves have the edge counts the file states, and the two 12-triangle 6-point
   configurations are a 2 x 3 grid and a centrally symmetric hexagon.

The box scan is a finite search, so it can only refute the classification, never confirm it; the proof is the
written argument in NOTES.md section 4 together with the enumeration in g7_exact.py.
"""
import itertools
import json
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def area2(a, b, c):
    return abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))


def profile(points, unit):
    """Sorted doubled areas of all triangles, in units of the common doubled area."""
    return tuple(sorted(Fraction(area2(*t), unit) for t in itertools.combinations(points, 3)))


def is_parallelogram(pts):
    a, b, c, d = pts
    return any((p[0] + q[0], p[1] + q[1]) == (r[0] + s[0], r[1] + s[1])
               for (p, q), (r, s) in (((a, b), (c, d)), ((a, c), (b, d)), ((a, d), (b, c))))


def is_centroid_config(pts, unit):
    for g in pts:
        tri = [p for p in pts if p != g]
        if (3 * g[0], 3 * g[1]) == (sum(p[0] for p in tri), sum(p[1] for p in tri)) and area2(*tri) == 3 * unit:
            return True
    return False


def check_lemma(side):
    box = [(x, y) for x in range(side) for y in range(side)]
    bad = 0
    for pts in itertools.combinations(box, 4):
        counts = Counter(area2(*t) for t in itertools.combinations(pts, 3) if area2(*t))
        for unit, k in counts.items():
            if k > 4 or (k == 4 and not is_parallelogram(pts)) or (k == 3 and not is_centroid_config(pts, unit)):
                bad += 1
    return bad


def scan(side, size, threshold):
    box = [(x, y) for x in range(side) for y in range(side)]
    seen = {}
    for pts in itertools.combinations(box, size):
        counts = Counter(area2(*t) for t in itertools.combinations(pts, 3) if area2(*t))
        for unit, k in counts.items():
            if k >= threshold:
                seen.setdefault(profile(pts, unit), (k, pts))
    return seen


def main(argv):
    path = Path(argv[0]) if argv else ROOT / "results" / "g7_classes.json"
    listed = json.loads(path.read_text())
    ok = True

    bad = check_lemma(7)
    print(f"four-point lemma on the 7x7 box: {bad} violations")
    ok &= bad == 0

    def parse(entry):
        pts = [(Fraction(x), Fraction(y)) for x, y in entry["points"]]
        den = 1
        for x, y in pts:
            den = den * x.denominator // __import__("math").gcd(den, x.denominator)
            den = den * y.denominator // __import__("math").gcd(den, y.denominator)
        ipts = tuple((int(x * den), int(y * den)) for x, y in pts)
        return ipts, den * den

    for key, size, threshold in (("five_points_at_least_6", 5, 6), ("six_points_at_least_11", 6, 11)):
        profiles = {}
        for entry in listed[key]:
            ipts, unit = parse(entry)
            k = sum(1 for t in itertools.combinations(ipts, 3) if area2(*t) == unit)
            if k != entry["edges"] or max(Counter(area2(*t) for t in itertools.combinations(ipts, 3)
                                             if area2(*t)).values()) != k:
                print(f"{key}: stored edge count {entry['edges']} is wrong")
                ok = False
            profiles[profile(ipts, unit)] = k
        found = scan(6, size, threshold)
        unlisted = [v for p, v in found.items() if p not in profiles]
        print(f"{key}: {len(listed[key])} listed, {len(found)} distinct profiles with >= {threshold} in the 6x6 box, "
              f"{len(unlisted)} unlisted")
        for k, pts in unlisted:
            print("  unlisted:", k, pts)
        ok &= not unlisted
        if len(found) < len(profiles):
            print(f"  note: {len(profiles) - len(found)} listed configuration(s) do not fit the 6x6 box")

    twelve = [parse(e) for e in listed["six_points_at_least_11"] if e["edges"] == 12]
    kinds = []
    for ipts, unit in twelve:
        xs, ys = sorted({p[0] for p in ipts}), sorted({p[1] for p in ipts})
        cx, cy = sum(p[0] for p in ipts), sum(p[1] for p in ipts)
        symmetric = all((2 * cx - 6 * p[0], 2 * cy - 6 * p[1]) in {(6 * q[0], 6 * q[1]) for q in ipts} for p in ipts)
        rows = Counter(p[1] for p in ipts)
        grid = sorted(rows.values()) == [3, 3] and all(
            sorted(p[0] for p in ipts if p[1] == y) == [min(p[0] for p in ipts if p[1] == y) + i * d for i in range(3)]
            for y in rows for d in [sorted(p[0] for p in ipts if p[1] == y)[1] - min(p[0] for p in ipts if p[1] == y)])
        kinds.append("grid" if grid else ("symmetric hexagon" if symmetric else "other"))
    print(f"12-triangle six-point configurations: {kinds}")
    ok &= sorted(kinds) == ["grid", "symmetric hexagon"]
    ok &= not listed["seven_points_at_least_19"]
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
