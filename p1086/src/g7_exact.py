"""Determine g(7) exactly, by a finite enumeration that shares no code with g5_certificate.py, g6_exact.py or g7_no21.py.

    python src/g7_exact.py            prints the classification and the value of g(7); exit 0 iff the run is complete
    python src/g7_exact.py --write    also writes results/g7_classes.json (--write-to PATH writes elsewhere)

Everything is exact rational arithmetic. "Edge" means a triangle of doubled area exactly 1 (after an area-preserving
affine normalisation the common area can always be taken to be 1/2, doubled area 1). The argument has one written
lemma and three enumeration steps; NOTES.md section 4 gives the proofs.

Lemma (four points). Four distinct points span at most 4 edges. If they span 4 they form a parallelogram. If they
span exactly 3, three of them form a triangle of doubled area 3 and the fourth is its centroid.

Step 1. A 5-point set with at least 6 edges has a 4-subset with at least 3 edges (each triple lies in 2 of the five
4-subsets, so the 4-subsets carry 2 * 6 = 12 > 5 * 2 edges between them). By the lemma that 4-subset is a
parallelogram, normalised to the unit square, or a centroid configuration, normalised to (0,0), (1,0), (0,3) with
centroid (1/3, 1). The fifth point then lies on at least 2 (square) or 3 (centroid) of the 12 lines on which a
point forms an edge with a pair of the 4-subset. No line serves more than one pair, so the fifth point is an
intersection of two distinct lines: finitely many candidates, each counted exactly.

Step 2. A 6-point set with at least 11 edges has a 5-subset with at least 6 (3 * 11 = 33 > 6 * 5), so it is a Step 1
configuration plus a point on at least 11 - 7 = 4 of its 20 lines. The script checks that no line serves 4 pairs,
so the sixth point is again an intersection of two distinct lines.

Step 3. A 7-point set with at least 19 edges has a vertex of degree at most floor(57/7) = 8, so a 6-subset with at
least 11 edges: a Step 2 configuration plus a point on at least 19 - 12 = 7 of its 30 lines. Again no line serves
that many pairs, so the seventh point is an intersection of two distinct lines.

The enumeration therefore finds every 7-point set with at least 19 edges up to area-preserving affine maps. If
there are none, g(7) <= 18, and the stored 7-point configuration with 18 gives g(7) = 18. As controls, Step 1 must
rediscover g(5) = 7 with a single configuration and Step 2 must rediscover g(6) = 12 with the 2 x 3 grid among the
extremal configurations, matching src/g5_certificate.py and src/g6_exact.py, which use different methods. Step 2
also finds the configurations those scripts could not see: the 6-point sets with 12 edges none of whose 5-subsets
has 7.
"""
import itertools
import json
import sys
import time
from fractions import Fraction as Fr
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def det(p, q, r):
    return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])


def edges(points):
    return sum(1 for t in itertools.combinations(points, 3) if abs(det(*t)) == 1)


def lines(points):
    """The lines {X : det(q - p, X - p) = s}, s = +-1, over unordered pairs {p, q}, as (a, b, c) with ax + by = c,
    normalised so that the first nonzero of (a, b) is 1. Returns a dict line -> number of pairs it serves."""
    out = {}
    for p, q in itertools.combinations(points, 2):
        a, b = -(q[1] - p[1]), q[0] - p[0]
        c0 = -(q[0] - p[0]) * p[1] + (q[1] - p[1]) * p[0]
        for s in (1, -1):
            c = s - c0
            scale = a if a else b
            key = (Fr(a) / scale, Fr(b) / scale, Fr(c) / scale)
            out[key] = out.get(key, 0) + 1
    return out


def intersections(line_dict):
    seen = set()
    for (a1, b1, c1), (a2, b2, c2) in itertools.combinations(line_dict, 2):
        d = a1 * b2 - a2 * b1
        if d == 0:
            continue
        pt = ((c1 * b2 - c2 * b1) / d, (a1 * c2 - a2 * c1) / d)
        if pt not in seen:
            seen.add(pt)
            yield pt


def canonical(points):
    """Smallest image of the set under the affine maps sending an edge to (0,0), (1,0), (0,1). These maps have
    determinant +-1, so the canonical form is a complete invariant under area-preserving affine equivalence."""
    best = None
    for p, q, r in itertools.permutations(points, 3):
        d = det(p, q, r)
        if abs(d) != 1:
            continue
        ux, uy, vx, vy = q[0] - p[0], q[1] - p[1], r[0] - p[0], r[1] - p[1]
        image = tuple(sorted(((vy * (x - p[0]) - vx * (y - p[1])) / d, (-uy * (x - p[0]) + ux * (y - p[1])) / d)
                             for x, y in points))
        if best is None or image < best:
            best = image
    return best


def extend(bases, needed_total, label):
    """All sets base + {X} with at least needed_total edges, for each base, up to equivalence."""
    found = {}
    checked = 0
    for base in bases:
        need = needed_total - edges(base)
        ld = lines(base)
        worst = max(ld.values())
        if worst >= need:
            raise RuntimeError(f"{label}: a line serves {worst} pairs, so a point on one line could add {need}")
        for pt in intersections(ld):
            if pt in base:
                continue
            new = base + (pt,)
            checked += 1
            e = edges(new)
            if e >= needed_total:
                c = canonical(new)
                found[c] = e
    print(f"{label}: {checked} candidate points, {len(found)} configurations with >= {needed_total} edges")
    return found


def main(argv):
    start = time.time()
    square = ((Fr(0), Fr(0)), (Fr(1), Fr(0)), (Fr(1), Fr(1)), (Fr(0), Fr(1)))
    centroid = ((Fr(0), Fr(0)), (Fr(1), Fr(0)), (Fr(0), Fr(3)), (Fr(1, 3), Fr(1)))
    if edges(square) != 4 or edges(centroid) != 3:
        print("the two normalised 4-point bases do not have 4 and 3 edges")
        return 1
    if det(*centroid[:3]) != 3 or any(abs(det(*t)) != 1 for t in itertools.combinations(centroid, 3)
                                       if centroid[3] in t):
        print("the centroid base is not a triangle of doubled area 3 with its centroid")
        return 1

    five = {}
    for base, need in ((square, 2), (centroid, 3)):
        ld = lines(base)
        if max(ld.values()) >= need:
            print("a line of the 4-point base serves too many pairs")
            return 1
        for pt in intersections(ld):
            if pt in base:
                continue
            new = base + (pt,)
            e = edges(new)
            if e >= 6:
                five[canonical(new)] = e
    print(f"step 1: {len(five)} five-point configurations with >= 6 edges; edge counts "
          f"{sorted(five.values(), reverse=True)}")
    six = extend(list(five), 11, "step 2")
    seven = extend(list(six), 19, "step 3")

    ok = True
    g5 = max(five.values())
    ok &= g5 == 7 and sum(1 for e in five.values() if e == 7) == 1
    print(f"control: g(5) = {g5}, {sum(1 for e in five.values() if e == 7)} configuration(s) with 7")
    g6 = max(six.values())
    twelve = [c for c, e in six.items() if e == 12]
    grid = canonical(tuple((Fr(x), Fr(y)) for x in range(3) for y in range(2)))
    ok &= g6 == 12 and grid in twelve
    print(f"control: g(6) = {g6}, {len(twelve)} configuration(s) with 12, one of them the 2 x 3 grid: "
          f"{grid in twelve}; {sum(1 for e in six.values() if e == 11)} with 11")
    for c in twelve:
        print("  12-edge configuration:", [(str(x), str(y)) for x, y in c])
    ok &= not seven
    print(f"seven-point configurations with >= 19 edges: {len(seven)}")

    target = None
    if "--write" in argv:
        target = ROOT / "results" / "g7_classes.json"
    elif "--write-to" in argv:
        target = Path(argv[argv.index("--write-to") + 1])
    if target is not None:
        out = {
            "doubled_area": 1,
            "five_points_at_least_6": [{"edges": e, "points": [[str(x), str(y)] for x, y in c]}
                                       for c, e in sorted(five.items(), key=lambda kv: (-kv[1], kv[0]))],
            "six_points_at_least_11": [{"edges": e, "points": [[str(x), str(y)] for x, y in c]}
                                       for c, e in sorted(six.items(), key=lambda kv: (-kv[1], kv[0]))],
            "seven_points_at_least_19": [{"edges": e, "points": [[str(x), str(y)] for x, y in c]}
                                         for c, e in sorted(seven.items(), key=lambda kv: (-kv[1], kv[0]))],
        }
        target.write_text(json.dumps(out, indent=1) + "\n")
    print(f"elapsed {time.time() - start:.2f} s")
    if ok:
        print("g(7) <= 18; with the stored 18-triangle configuration, g(7) = 18")
        return 0
    print("the controls or the enumeration did not come out as expected")
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
