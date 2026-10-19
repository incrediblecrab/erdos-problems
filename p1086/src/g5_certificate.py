"""Exact certificate that 5 points in the plane span at most 7 triangles of one nonzero area, i.e. g(5) <= 7.

Independent of the hand proof in NOTES.md.  Suppose 5 points span >= K triangles of doubled area s != 0.  An affine
map sends one of those triangles to (0,0), (1,0), (0,1); affine maps scale every signed doubled area by the same
nonzero factor, so afterwards s = 1 and the other two points are P3 = (x3, y3), P4 = (x4, y4).  At least K - 1 of
the other 9 triangles t have det_t = +1 or -1.  For every choice S of K - 1 of them and every sign vector, the
polynomial system {det_t = sign_t : t in S} is sent to a Groebner basis over QQ.  A basis equal to [1] means the
system has no solution even over C.

    python src/g5_certificate.py 8    must report 0 consistent systems: nothing spans 8, so g(5) <= 7
    python src/g5_certificate.py 7    the control: must find consistent systems, and a real solution that spans 7
"""
import itertools
import sys
from multiprocessing import Pool

import sympy as sp

x3, y3, x4, y4 = sp.symbols("x3 y3 x4 y4")
P = [(0, 0), (1, 0), (0, 1), (x3, y3), (x4, y4)]
TRIPLES = [t for t in itertools.combinations(range(5), 3) if t != (0, 1, 2)]


def det(t):
    (ax, ay), (bx, by), (cx, cy) = (P[i] for i in t)
    return sp.expand((bx - ax) * (cy - ay) - (by - ay) * (cx - ax))


DETS = {t: det(t) for t in TRIPLES}


def self_test(trials=25):
    """Check every det polynomial against a shoelace area computed from random rational points."""
    import random
    rng = random.Random(1086)
    for _ in range(trials):
        vals = {v: sp.Rational(rng.randint(-40, 40), rng.randint(1, 9)) for v in (x3, y3, x4, y4)}
        pts = [(sp.Integer(0), sp.Integer(0)), (sp.Integer(1), sp.Integer(0)), (sp.Integer(0), sp.Integer(1)),
               (vals[x3], vals[y3]), (vals[x4], vals[y4])]
        for t, poly in DETS.items():
            (ax, ay), (bx, by), (cx, cy) = (pts[i] for i in t)
            shoelace = ax * by - bx * ay + bx * cy - cx * by + cx * ay - ax * cy
            if poly.subs(vals) != shoelace:
                return False
    return True


def solve_case(case):
    S, signs = case
    eqs = [DETS[t] - e for t, e in zip(S, signs)]
    G = sp.groebner(eqs, x3, y3, x4, y4, order="lex", domain=sp.QQ)
    if list(G.exprs) == [1]:
        return None
    return (S, signs, [str(g) for g in G.exprs])


def count_equal(pts):
    """Doubled-area histogram of a concrete 5-point set, computed from scratch."""
    hist = {}
    for a, b, c in itertools.combinations(pts, 3):
        d = sp.simplify(sp.Abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])))
        if d != 0:
            hist[d] = hist.get(d, 0) + 1
    return hist


def real_points(S, signs):
    """Real solutions of a consistent system, when it has finitely many."""
    eqs = [DETS[t] - e for t, e in zip(S, signs)]
    sols = sp.solve(eqs, [x3, y3, x4, y4], dict=True)
    out = []
    for s in sols:
        if len(s) < 4:
            out.append(("positive-dimensional", s))
            continue
        vals = [s[v] for v in (x3, y3, x4, y4)]
        if all(sp.im(sp.nsimplify(v)) == 0 for v in vals):
            out.append(("real", vals))
    return out


def main(argv):
    K = int(argv[0]) if argv else 8
    if not self_test():
        print("self-test FAILED: a det polynomial disagrees with the shoelace formula")
        return 1
    cases = [(S, signs) for S in itertools.combinations(TRIPLES, K - 1)
             for signs in itertools.product((1, -1), repeat=K - 1)]
    with Pool(10) as pool:
        results = pool.map(solve_case, cases, chunksize=16)
    consistent = [r for r in results if r is not None]
    print(f"K={K}: {len(cases)} systems, {len(consistent)} consistent over C")
    if K == 8:
        return 0 if not consistent else 1
    best = 0
    for S, signs, _ in consistent[:200]:
        for kind, vals in real_points(S, signs):
            if kind != "real":
                continue
            pts = [(sp.Integer(0), sp.Integer(0)), (sp.Integer(1), sp.Integer(0)), (sp.Integer(0), sp.Integer(1)),
                   (vals[0], vals[1]), (vals[2], vals[3])]
            if len(set(pts)) < 5:
                continue
            best = max(best, count_equal(pts).get(1, 0))
    print(f"K={K}: best real distinct configuration found among the first 200 consistent systems spans {best} unit triangles")
    return 0 if best >= K else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
