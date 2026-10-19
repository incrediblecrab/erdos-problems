"""g(6) exactly, from the g(5) certificate.

1. Classify: every real 5-point set with 7 triangles of one nonzero area is affinely equivalent to
   R = {(0,0), (1,0), (0,1), (1,1), (1,-1)}.  The consistent K = 7 systems of g5_certificate are solved
   (solve and solve_poly_system must agree), and an affine bijection onto R is exhibited for every real solution.
2. Six points with t >= 13 equal-area triangles contain a 5-subset with 7 of them: each triangle lies in 3 of the
   6 five-subsets, and 3t >= 39 > 36 = 6 * 6.  After an affine map that 5-subset is R with doubled area 1, and at
   least 6 of the 10 triangles through the sixth point X have det = +-1.  Each such condition puts X on a line
   parallel to a side of R; no direction occurs among more than 4 pairs of R, so X lies on two non-parallel such
   lines.  Every intersection point is tried with exact rationals, so g(6) = max(12, best count found).

    python src/g6_exact.py
"""
import itertools
import sys
from fractions import Fraction

import sympy as sp

import g5_certificate as g5

R = [(0, 0), (1, 0), (0, 1), (1, 1), (1, -1)]


def det3(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def profile(pts):
    h = {}
    for a, b, c in itertools.combinations(pts, 3):
        d = abs(det3(a, b, c))
        h[d] = h.get(d, 0) + 1
    return h


def affine_to_R(pts):
    """An affine bijection pts -> R, found by trying every labelling, or None."""
    for perm in itertools.permutations(range(5)):
        tgt = [R[i] for i in perm]
        (p0, p1, p2), (q0, q1, q2) = pts[:3], tgt[:3]
        M = sp.Matrix([[p1[0] - p0[0], p2[0] - p0[0]], [p1[1] - p0[1], p2[1] - p0[1]]])
        if M.det() == 0:
            continue
        N = sp.Matrix([[q1[0] - q0[0], q2[0] - q0[0]], [q1[1] - q0[1], q2[1] - q0[1]]])
        A = N * M.inv()
        t = sp.Matrix(q0) - A * sp.Matrix(p0)
        if all(sp.simplify(A * sp.Matrix(p) + t - sp.Matrix(q)) == sp.zeros(2, 1) for p, q in zip(pts, tgt)):
            return perm
    return None


def n_standard_monomials(eqs, gens):
    """dim_QQ QQ[gens]/(eqs): the number of complex solutions counted with multiplicity, or None if infinite."""
    G = sp.groebner(eqs, *gens, order="grevlex", domain=sp.QQ)
    lead = [tuple(sp.Poly(g, *gens).LM(order="grevlex")) for g in G.exprs]
    k = len(gens)
    bounds = []
    for i in range(k):
        pure = [m[i] for m in lead if all(m[j] == 0 for j in range(k) if j != i)]
        if not pure:
            return None
        bounds.append(min(pure))
    return sum(1 for m in itertools.product(*(range(b) for b in bounds))
               if not any(all(m[j] >= l[j] for j in range(k)) for l in lead))


def classify():
    cases = [(S, s) for S in itertools.combinations(g5.TRIPLES, 6) for s in itertools.product((1, -1), repeat=6)]
    consistent = [r for r in map(g5.solve_case, cases) if r is not None]
    n_real = 0
    for S, signs, _ in consistent:
        eqs = [g5.DETS[t] - e for t, e in zip(S, signs)]
        gens = (g5.x3, g5.y3, g5.x4, g5.y4)
        sols = sp.solve(eqs, list(gens), dict=True)
        if any(len(d) < 4 for d in sols):
            raise SystemExit(f"positive-dimensional solution set for {S} {signs}")
        a = {tuple(d[v] for v in gens) for d in sols}
        if n_standard_monomials(eqs, gens) != len(a):
            raise SystemExit(f"solution count differs from the quotient dimension for {S} {signs}")
        b = set(sp.solve_poly_system(eqs, *gens) or [])
        if {tuple(map(sp.nsimplify, s)) for s in a} != {tuple(map(sp.nsimplify, s)) for s in b}:
            raise SystemExit(f"solve and solve_poly_system disagree on {S} {signs}")
        for sol in a:
            if len(sol) != 4 or not all(v.is_real for v in sol):
                continue
            pts = [(0, 0), (1, 0), (0, 1), (sol[0], sol[1]), (sol[2], sol[3])]
            if len(set(pts)) < 5:
                continue
            n_real += 1
            if affine_to_R(pts) is None:
                raise SystemExit(f"a 7-configuration not equivalent to R: {pts}")
    return len(consistent), n_real


def sixth_point():
    pairs = list(itertools.combinations(range(5), 2))
    # det(R_b - R_a, X - R_a) = (bx-ax)(y-ay) - (by-ay)(x-ax) = -dy*x + dx*y + (dy*ax - dx*ay)
    lines = []
    for a, b in pairs:
        (ax, ay), (bx, by) = R[a], R[b]
        dx, dy = bx - ax, by - ay
        for s in (1, -1):
            lines.append((-dy, dx, s - (dy * ax - dx * ay)))  # alpha*x + beta*y = gamma
    counts = {}
    for (a1, b1, c1), (a2, b2, c2) in itertools.combinations(lines, 2):
        D = a1 * b2 - a2 * b1
        if D == 0:
            continue
        X = (Fraction(c1 * b2 - c2 * b1, D), Fraction(a1 * c2 - a2 * c1, D))
        if X in [(Fraction(p[0]), Fraction(p[1])) for p in R]:
            continue
        counts[X] = profile([(Fraction(p[0]), Fraction(p[1])) for p in R] + [X]).get(1, 0)
    best = max(counts.values())
    return best, sorted(X for X, c in counts.items() if c == best), len(counts)


def main():
    if not g5.self_test():
        print("g5 self-test failed")
        return 1
    n_cons, n_real = classify()
    print(f"classification: {n_cons} consistent systems, {n_real} real 5-point solutions, all affinely equivalent to R")
    print(f"R profile (doubled area: count): {profile(R)}")
    dirs = {}
    for a, b in itertools.combinations(R, 2):
        v = (b[0] - a[0], b[1] - a[1])
        gg = abs(sp.gcd(v[0], v[1]))
        v = (v[0] // gg, v[1] // gg)
        if v < (0, 0):
            v = (-v[0], -v[1])
        dirs[v] = dirs.get(v, 0) + 1
    print(f"largest direction class among pairs of R: {max(dirs.values())} (must be < 6)")
    best, argmax, n_cand = sixth_point()
    print(f"{n_cand} candidate sixth points; the largest count of doubled-area-1 triangles in R + X is {best}, at")
    print("  " + " ".join(f"({x},{y})" for x, y in argmax))
    lower = max(v for d, v in profile([(x, y) for x in range(3) for y in range(2)]).items() if d != 0)
    print(f"2 x 3 grid: {lower} triangles of one nonzero area")
    print(f"g(6) = max({lower}, {best}) = {max(lower, best)}")
    return 0 if max(dirs.values()) < 6 else 1


if __name__ == "__main__":
    sys.exit(main())
