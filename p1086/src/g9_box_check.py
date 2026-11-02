"""Cross-check src/g9_exact.py, sharing no code with it or with src/g7_exact.py.

    python src/g9_box_check.py [results/g9_classes.json]      exit 0 iff every check passes

Needs src/box_scan built (see README.md) and z3. Four checks:

1. Listed counts. Every configuration in the file, scaled to integer coordinates, has the stated number of triangles
   of the common doubled area, and no larger class of another area.
2. Box scan. src/box_scan lists every k-subset of the 6 x 6 box with at least t triangles of one doubled area, for
   (k, t) = (6, 10), (7, 16), (8, 25), (9, 37). Each must have the area profile (all doubled areas divided by the
   common one, an invariant of area-preserving affine maps) of a listed configuration, and none may have 9 points.
3. Re-derived extensions. From the listed (k-1)-point configurations (the 5-point ones from g7_classes.json), every
   point lying on two non-parallel lines {X : det(q - p, X - p) = +-unit} is added, in integer arithmetic written
   here; the profiles reaching the next threshold must be exactly the listed profiles of the next level, and no
   line may serve as many pairs as the new point needs. This repeats the enumeration with independent code, so it
   guards against implementation errors, not against an error in the argument.
4. Type II, by a different method. For each 3-uniform hypergraph on 6 vertices with 10 edges and all degrees 5
   (enumerated here), z3's nonlinear real arithmetic must find no six distinct real points whose triangles of
   doubled area 1 are exactly those edges. As a control, the hexagon's 12-edge hypergraph must be satisfiable.

Checks 2 and 4 are finite searches or solver calls that can refute the classification; the proof is the written
argument in NOTES.md section 4 together with the enumeration in g9_exact.py.
"""
import itertools
import json
import math
import subprocess
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

import z3

ROOT = Path(__file__).resolve().parent.parent
SCAN = ROOT / "src" / "box_scan"
LEVELS = (("six_points_at_least_10", 6, 10), ("seven_points_at_least_16", 7, 16),
          ("eight_points_at_least_25", 8, 25), ("nine_points_at_least_37", 9, 37))


def area2(a, b, c):
    return abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))


def profile(points, unit):
    return tuple(sorted(Fraction(area2(*t), unit) for t in itertools.combinations(points, 3)))


def to_integers(entry):
    """Scale a listed configuration (rational strings, common doubled area 1) to integer points; return them and
    the common doubled area after scaling."""
    pts = [(Fraction(x), Fraction(y)) for x, y in entry["points"]]
    den = 1
    for x, y in pts:
        den = math.lcm(den, x.denominator, y.denominator)
    return [(int(x * den), int(y * den)) for x, y in pts], den * den


def unit_count(points, unit):
    return sum(1 for t in itertools.combinations(points, 3) if area2(*t) == unit)


def largest_class(points):
    return max(Counter(area2(*t) for t in itertools.combinations(points, 3) if area2(*t)).values())


def candidate_points(points, unit):
    """Points on two non-parallel lines det(q - p, X - p) = +-unit, with the number of pairs each line serves."""
    rows = Counter()
    for p, q in itertools.combinations(points, 2):
        a, b = p[1] - q[1], q[0] - p[0]  # det(q - p, X - p) = a x + b y - (a p_x + b p_y)
        k = a * p[0] + b * p[1]
        for s in (unit, -unit):
            g = math.gcd(math.gcd(a, b), k + s)
            row = (a // g, b // g, (k + s) // g)
            if row[0] < 0 or (row[0] == 0 and row[1] < 0):
                row = (-row[0], -row[1], -row[2])
            rows[row] += 1
    pts = set()
    for (a1, b1, c1), (a2, b2, c2) in itertools.combinations(rows, 2):
        d = a1 * b2 - a2 * b1
        if d:
            pts.add((Fraction(c1 * b2 - c2 * b1, d), Fraction(a1 * c2 - a2 * c1, d)))
    return pts, max(rows.values())


def extend_profiles(bases, threshold):
    """Profiles of every base plus a candidate point with at least threshold triangles of the base's unit area."""
    out, worst_ok = set(), True
    for points, unit in bases:
        need = threshold - unit_count(points, unit)
        cands, worst = candidate_points(points, unit)
        worst_ok &= worst < need
        for x, y in cands:
            if (x, y) in points:
                continue
            den = math.lcm(x.denominator, y.denominator)
            scaled = [(px * den, py * den) for px, py in points] + [(int(x * den), int(y * den))]
            if unit_count(scaled, unit * den * den) >= threshold:
                out.add(profile(scaled, unit * den * den))
    return out, worst_ok


def regular_hypergraphs():
    triples = list(itertools.combinations(range(6), 3))
    reps = set()
    for chosen in itertools.combinations(triples, 10):
        if all(sum(v in t for t in chosen) == 5 for v in range(6)):
            reps.add(min(tuple(sorted(tuple(sorted(p[v] for v in t)) for t in chosen))
                         for p in itertools.permutations(range(6))))
    return sorted(reps)


def exact_edge_set_realisable(edge_set):
    """z3: are there six distinct real points whose doubled-area-1 triangles are exactly edge_set?"""
    xs = [z3.Real(f"x{i}") for i in range(6)]
    ys = [z3.Real(f"y{i}") for i in range(6)]
    s = z3.SolverFor("QF_NRA")
    a, b, c = edge_set[0]
    s.add(xs[a] == 0, ys[a] == 0, xs[b] == 1, ys[b] == 0, xs[c] == 0, ys[c] == 1)
    for i, j, k in itertools.combinations(range(6), 3):
        d = (xs[j] - xs[i]) * (ys[k] - ys[i]) - (ys[j] - ys[i]) * (xs[k] - xs[i])
        s.add(z3.Or(d == 1, d == -1) if (i, j, k) in edge_set else z3.And(d != 1, d != -1))
    for i, j in itertools.combinations(range(6), 2):
        s.add(z3.Or(xs[i] != xs[j], ys[i] != ys[j]))
    return s.check()


def main(argv):
    path = Path(argv[0]) if argv else ROOT / "results" / "g9_classes.json"
    listed = json.loads(path.read_text())
    g7 = json.loads((ROOT / "results" / "g7_classes.json").read_text())
    ok = True

    levels = {}
    for key, k, t in LEVELS:
        levels[key] = []
        for entry in listed[key]:
            pts, unit = to_integers(entry)
            n = unit_count(pts, unit)
            if len(set(pts)) != k or n != entry["edges"] or n < t or largest_class(pts) != n:
                print(f"{key}: stored count {entry['edges']} not reproduced (recount {n})")
                ok = False
            levels[key].append((pts, unit))
    print("listed counts: " + ", ".join(f"{key} {len(levels[key])}" for key, _, _ in LEVELS))

    for key, k, t in LEVELS:
        mine = {profile(p, u) for p, u in levels[key]}
        if len(mine) != len(levels[key]):
            print(f"  note: {key} has {len(levels[key]) - len(mine)} profile collision(s)")
        r = subprocess.run([str(SCAN), "6", "6", str(k), str(t), "3"], capture_output=True, text=True, check=True)
        lines = r.stdout.split("\n")
        found = {}
        for line in lines:
            if not line or line.startswith("done"):
                continue
            f = line.split()
            pts = [tuple(map(int, s.split(","))) for s in f[2:]]
            if int(f[1]) < t or unit_count(pts, int(f[0])) != int(f[1]):
                print(f"  box_scan line not reproduced: {line}")
                ok = False
            found.setdefault(profile(pts, int(f[0])), pts)
        unlisted = [pts for p, pts in found.items() if p not in mine]
        print(f"box 6x6, {k}-subsets with >= {t}: {len(found)} profiles, {len(unlisted)} unlisted, "
              f"{len(mine - set(found))} listed profile(s) not in the box")
        for pts in unlisted:
            print("  unlisted:", pts)
        ok &= not unlisted

    five = [to_integers(e) for e in g7["five_points_at_least_6"]]
    previous = five
    for key, k, t in LEVELS:
        got, worst_ok = extend_profiles(previous, t)
        mine = {profile(p, u) for p, u in levels[key]}
        same = got == mine if key != "six_points_at_least_10" else (
            got == mine and listed["type_ii"]["configurations_found_only_by_type_ii"] == 0)
        print(f"re-derived {k}-point extensions with >= {t}: {len(got)} profiles, equal to the listed ones: {same}; "
              f"every line serves fewer pairs than needed: {worst_ok}")
        ok &= same and worst_ok
        previous = levels[key]
    ok &= not listed["nine_points_at_least_37"]

    reps = regular_hypergraphs()
    results = [exact_edge_set_realisable(h) for h in reps]
    print(f"Type II by z3: {len(reps)} hypergraph classes, results {[str(r) for r in results]}")
    ok &= len(reps) == listed["type_ii"]["hypergraph_classes"] and all(r == z3.unsat for r in results)
    hexagon = [(-1, 0), (1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)]
    hex_edges = tuple(t for t in itertools.combinations(range(6), 3) if area2(*(hexagon[i] for i in t)) == 2)
    order = sorted(range(len(hex_edges)), key=lambda i: hex_edges[i])
    control = exact_edge_set_realisable(tuple(hex_edges[i] for i in order))
    print(f"control: the hexagon's {len(hex_edges)}-edge hypergraph is {control}")
    ok &= len(hex_edges) == 12 and control == z3.sat

    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
