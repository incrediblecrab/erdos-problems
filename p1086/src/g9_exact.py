"""Determine g(9) exactly, by extending the finite classification of src/g7_exact.py down to lower thresholds.

    python src/g9_exact.py            prints the classification and the bound; exit 0 iff the run is complete
    python src/g9_exact.py --write    also writes results/g9_classes.json (--write-to PATH writes elsewhere)

Everything is exact. As in g7_exact.py, an "edge" is a triangle of doubled area exactly 1: one affine map of
determinant 1/D first scales the common doubled area D to 1, and after that equivalence and the canonical form use
only affine maps of determinant +-1. C(n, t) denotes the n-point sets with at least t edges, up to those maps.
The written steps are in NOTES.md section 4; this script carries out the finite parts.

Averaging. A set of n points with e edges has a point of degree at most floor(3e/n), so deleting it leaves at least
e - floor(3e/n) edges, and e - floor(3e/n) does not decrease as e grows. Hence for t = 16, 25, 37 every member of
C(7, 16), C(8, 25), C(9, 37) is a member of C(6, 10), C(7, 16), C(8, 25) respectively plus one point, and the new
point lies on at least t - c of the lines of a base with c edges. The script checks that no line of any base serves
that many pairs, so the new point is an intersection of two distinct lines: finitely many candidates, each counted
exactly (the `extend` step of g7_exact.py).

Six points with at least 10 edges. With e >= 11 edges, some 5-subset has at least 6 (3e > 30). With e = 10, either
some 5-subset has at least 6 (Type I) or every 5-subset has at most 5; the five-subset counts sum to 3 * 10 = 30,
so then every 5-subset has exactly 5 and every point has degree exactly 5 (Type II). Type I is the four
configurations of C(5, 6) from g7_exact.py plus a point on at least 10 - c of their lines, handled as above. For
Type II the script lists every 3-uniform hypergraph on 6 vertices with 10 edges and all degrees 5, up to
isomorphism, sends one edge to (0,0), (1,0), (0,1) and solves det = +-1 on the other nine edges, for every sign
vector, with Groebner bases over QQ. Every consistent system must be zero-dimensional; its number of distinct
solutions must equal the quotient dimension (so no complex solution is missed); every solution must satisfy the
system exactly on substitution; a solution is discarded only if a coordinate is provably non-real, and the run
aborts if realness cannot be decided; every real solution must be an exact rational, and its ten edges are rechecked
in rational arithmetic. Real solutions with six distinct points and at least 10 edges join C(6, 10).

If C(9, 37) is empty then g(9) <= 36, and the stored 9-point configuration with 36 gives g(9) = 36.
"""
import itertools
import json
import sys
import time
from fractions import Fraction as Fr
from pathlib import Path

import sympy as sp

from g6_exact import n_standard_monomials
from g7_exact import canonical, det, edges, extend, intersections, lines

ROOT = Path(__file__).resolve().parent.parent
# (n, threshold for n points, threshold the averaging step needs for n - 1 points)
STEPS = ((7, 16, 10), (8, 25, 16), (9, 37, 25))


def averaging_holds(n, t, base_t):
    """Every n-set with e >= t edges has an (n-1)-subset with >= base_t edges, for every e up to C(n, 3)."""
    return all(e - (3 * e) // n >= base_t for e in range(t, n * (n - 1) * (n - 2) // 6 + 1))


def type_ii_hypergraphs():
    """3-uniform hypergraphs on 6 vertices with 10 edges and every degree 5, one per isomorphism class."""
    triples = list(itertools.combinations(range(6), 3))
    found = set()
    for chosen in itertools.combinations(triples, 10):
        deg = [0] * 6
        for t in chosen:
            for v in t:
                deg[v] += 1
        if deg != [5] * 6:
            continue
        found.add(min(tuple(sorted(tuple(sorted(p[v] for v in t)) for t in chosen))
                      for p in itertools.permutations(range(6))))
    return sorted(found)


def type_ii(hypergraphs):
    gens = sp.symbols("x3 y3 x4 y4 x5 y5")
    pos = {0: (0, 0), 1: (1, 0), 2: (0, 1), 3: gens[0:2], 4: gens[2:4], 5: gens[4:6]}
    found, n_systems, n_consistent, n_solutions, dims = {}, 0, 0, 0, []
    for hi, h in enumerate(hypergraphs):
        base, rest = h[0], h[1:]
        if base != (0, 1, 2):
            raise RuntimeError(f"hypergraph {h} does not start with the edge (0, 1, 2)")
        dets = []
        for a, b, c in rest:
            p, q, r = pos[a], pos[b], pos[c]
            dets.append(sp.expand((q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])))
        for signs in itertools.product((1, -1), repeat=len(rest)):
            n_systems += 1
            eqs = [d - s for d, s in zip(dets, signs)]
            if sp.groebner(eqs, *gens, order="grevlex", domain=sp.QQ).exprs == [1]:
                continue
            n_consistent += 1
            dim = n_standard_monomials(eqs, gens)
            if dim is None:
                raise RuntimeError(f"positive-dimensional Type II system {h} {signs}")
            sols = sp.solve_poly_system(eqs, *gens)
            sols = {tuple(sp.sympify(v) for v in sol) for sol in (sols or [])}
            dims.append((hi, signs, dim, len(sols)))
            if len(sols) != dim:
                raise RuntimeError(f"{len(sols)} distinct solutions but quotient dimension {dim} for {h} {signs}")
            for sol in sols:
                n_solutions += 1
                # every solution, real or not, must satisfy all nine equations exactly
                subs = dict(zip(gens, sol))
                if any(sp.expand(e.subs(subs)) != 0 for e in eqs):
                    raise RuntimeError(f"solution {sol} does not satisfy the system {h} {signs}")
                if any(v.is_real is None for v in sol):
                    raise RuntimeError(f"cannot decide whether solution {sol} is real")
                if any(v.is_real is False for v in sol):
                    continue  # provably non-real
                if not all(isinstance(v, sp.Rational) for v in sol):
                    raise RuntimeError(f"irrational real Type II solution {sol}: not handled by this script")
                pts = ((Fr(0), Fr(0)), (Fr(1), Fr(0)), (Fr(0), Fr(1))) + tuple(
                    (Fr(int(sol[2 * i].p), int(sol[2 * i].q)), Fr(int(sol[2 * i + 1].p), int(sol[2 * i + 1].q)))
                    for i in range(3))
                # recheck all ten edges of the hypergraph in exact rationals: det = +-1, with the stated signs
                d = lambda t: det(pts[t[0]], pts[t[1]], pts[t[2]])
                if d(base) != 1 or any(d(t) != s for t, s in zip(rest, signs)):
                    raise RuntimeError(f"rational solution {pts} fails an edge of {h} {signs}")
                if len(set(pts)) == 6 and edges(pts) >= 10:
                    found[canonical(pts)] = edges(pts)
    return found, n_systems, n_consistent, n_solutions, dims


def as_json(classes):
    return [{"edges": e, "points": [[str(x), str(y)] for x, y in c]}
            for c, e in sorted(classes.items(), key=lambda kv: (-kv[1], kv[0]))]


def main(argv):
    start = time.time()
    ok = True
    g7 = json.loads((ROOT / "results" / "g7_classes.json").read_text())
    parse = lambda e: tuple((Fr(x), Fr(y)) for x, y in e["points"])
    five = {parse(e): e["edges"] for e in g7["five_points_at_least_6"]}
    six11 = {parse(e): e["edges"] for e in g7["six_points_at_least_11"]}
    if sorted(five.values()) != [6, 6, 6, 7] or any(edges(c) != e for c, e in five.items()):
        print("results/g7_classes.json does not hold the four 5-point configurations of g7_exact.py")
        return 1

    for n, t, base_t in STEPS:
        holds = averaging_holds(n, t, base_t)
        print(f"averaging: every {n}-set with >= {t} edges has an {n - 1}-subset with >= {base_t}: {holds}")
        ok &= holds
    # the six-point step: e >= 11 forces a 5-subset with >= 6; e = 10 is split into Type I and Type II
    ok &= all(3 * e > 6 * 5 for e in (11, 12))

    hypergraphs = type_ii_hypergraphs()
    t2, n_systems, n_consistent, n_solutions, dims = type_ii(hypergraphs)
    for hi, signs, dim, nsol in dims:
        print(f"Type II system: hypergraph {hi}, signs {signs}: quotient dimension {dim}, {nsol} distinct "
              f"solution(s), all satisfying the system exactly")
    print(f"Type II: {len(hypergraphs)} hypergraph classes, {n_systems} sign systems, {n_consistent} consistent, "
          f"{n_solutions} solutions, {len(t2)} configuration(s) with >= 10 edges, edge counts {sorted(t2.values())}")
    # the affine regular hexagon {+-u, +-v, +-(u - v)} has 12 triangles of doubled area 2; (x, y) -> (x/2, y)
    hexagon = canonical(tuple((Fr(x, 2), Fr(y)) for x, y in ((1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1))))
    ok &= t2 == {hexagon: 12}
    print(f"control: every Type II solution is the affine regular hexagon, with 12 edges: {t2 == {hexagon: 12}}")
    six = extend(list(five), 10, "six points, Type I")
    t2_new = {c: e for c, e in t2.items() if c not in six}
    six.update(t2)
    print(f"C(6, 10): {len(six)} configurations, edge counts {sorted(six.values(), reverse=True)}; "
          f"{len(t2_new)} found only by Type II")
    ok &= all(six.get(c) == e for c, e in six11.items()) and sum(1 for e in six.values() if e >= 11) == len(six11)
    print(f"control: C(6, 10) contains exactly the {len(six11)} configurations of C(6, 11) in g7_classes.json")

    seven = extend(list(six), 16, "C(7, 16)")
    grid7 = canonical(tuple((Fr(x), Fr(y)) for x, y in ((0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1), (1, 2))))
    ok &= max(seven.values()) == 18 and seven.get(grid7) == 18
    print(f"control: largest count in C(7, 16) is {max(seven.values())} (g(7) = 18), with "
          f"{sum(1 for e in seven.values() if e == 18)} configurations; 2 x 3 grid plus (1, 2) among them: "
          f"{seven.get(grid7) == 18}")

    eight = extend(list(seven), 25, "C(8, 25)")
    # the 3 x 3 grid minus its centre has doubled area 2; (x, y) -> (x/2, y) halves every area
    ring = canonical(tuple((Fr(x, 2), Fr(y)) for x in range(3) for y in range(3) if (x, y) != (1, 1)))
    top8 = [c for c, e in eight.items() if e == 28]
    ok &= max(eight.values()) == 28 and top8 == [ring]
    print(f"control: largest count in C(8, 25) is {max(eight.values())} (g(8) = 28), attained only by the 3 x 3 "
          f"grid minus its centre: {top8 == [ring]}")

    nine = extend(list(eight), 37, "C(9, 37)")
    best9 = 0
    for base in eight:
        for pt in intersections(lines(base)):
            if pt not in base:
                best9 = max(best9, edges(base + (pt,)))
    print(f"largest count among all intersection-point extensions of C(8, 25): {best9}")
    ok &= not nine

    # why the method stops here: bounding g(10) this way needs C(6, 9), which is infinite, because R (7 edges) has
    # lines serving two pairs; R plus any point of the line x = 2 has at least 9 edges, and exactly 9 unless the
    # point also lies on one of the finitely many other lines of R
    R = ((Fr(0), Fr(0)), (Fr(1), Fr(0)), (Fr(0), Fr(1)), (Fr(1), Fr(1)), (Fr(1), Fr(-1)))
    ld = lines(R)
    double = ld.get((Fr(1), Fr(0), Fr(2)), 0)
    special = {(c - 2 * a) / b for (a, b, c) in ld if b != 0}
    ts = [Fr(t, 7) for t in range(-50, 51)]
    counts = {t: edges(R + ((Fr(2), t),)) for t in ts}
    generic = all(counts[t] == 9 for t in ts if t not in special) and all(e >= 9 for e in counts.values())
    print(f"limit: the line x = 2 serves {double} pairs of R; R plus (2, t) has >= 9 edges at all {len(ts)} values "
          f"of t tried and exactly 9 off the {len(special)} other lines: {generic}, so C(6, 9) is infinite")
    ok &= double == 2 and generic and sum(1 for t in ts if t not in special) > 90

    target = None
    if "--write" in argv:
        target = ROOT / "results" / "g9_classes.json"
    elif "--write-to" in argv:
        target = Path(argv[argv.index("--write-to") + 1])
    if target is not None:
        out = {
            "doubled_area": 1,
            "type_ii": {"hypergraph_classes": len(hypergraphs), "sign_systems": n_systems,
                        "consistent": n_consistent, "configurations_found_only_by_type_ii": len(t2_new)},
            "six_points_at_least_10": as_json(six),
            "seven_points_at_least_16": as_json(seven),
            "eight_points_at_least_25": as_json(eight),
            "nine_points_at_least_37": as_json(nine),
            "largest_nine_point_extension": best9,
        }
        target.write_text(json.dumps(out, indent=1) + "\n")
    print(f"elapsed {time.time() - start:.2f} s")
    if ok:
        print("no 9-point set has 37 triangles of one area, so g(9) <= 36; with the stored configuration, g(9) = 36")
        return 0
    print("the controls or the enumeration did not come out as expected")
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
