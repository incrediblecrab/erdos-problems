"""Independent recount for the sixth-point step of g6_exact.py; shares no code with it.

If 6 points span >= 13 triangles of one nonzero area, an affine map sends 5 of them to
R = {(0,0), (1,0), (0,1), (1,1), (1,-1)} with doubled area 1, and the sixth point X satisfies |det| = 1 for at least
6 pairs of R.  Each condition is a line c1*x + c2*y = c3 with |c1|, |c2| <= 2 and |c3| <= 4; two non-parallel ones
meet where Cramer's rule gives denominator det(v, w) for difference vectors v, w of R, and numerators of size at most
4*2 + 4*2 = 16.  This script lists those denominators, then scans every X in (1/L)Z^2 with |x|, |y| <= 16 for their
lcm L, counting doubled-area-1 triangles of R + X in integer arithmetic.

    python src/g6_brute.py    prints the largest count; g(6) = 12 needs it to be 12
"""
import itertools
import math
import sys
from fractions import Fraction

BASE = [(0, 0), (1, 0), (0, 1), (1, 1), (1, -1)]


def main():
    diffs = [(q[0] - p[0], q[1] - p[1]) for p, q in itertools.combinations(BASE, 2)]
    dens = {abs(v[0] * w[1] - v[1] * w[0]) for v in diffs for w in diffs} - {0}
    L = math.lcm(*dens)
    pts = [(L * x, L * y) for x, y in BASE]
    target = L * L
    best, where = 0, []
    for i in range(-16 * L, 16 * L + 1):
        for j in range(-16 * L, 16 * L + 1):
            if (i, j) in pts:
                continue
            s = pts + [(i, j)]
            c = sum(1 for a, b, e in itertools.combinations(s, 3)
                    if abs((b[0] - a[0]) * (e[1] - a[1]) - (b[1] - a[1]) * (e[0] - a[0])) == target)
            if c > best:
                best, where = c, [(i, j)]
            elif c == best:
                where.append((i, j))
    print(f"denominators {sorted(dens)}, L = {L}, {(32 * L + 1) ** 2 - 5} points scanned")
    print(f"largest count {best}, at")
    print("  " + " ".join(f"({Fraction(i, L)},{Fraction(j, L)})" for i, j in sorted(where)))
    return 0 if best == 12 else 1


if __name__ == "__main__":
    sys.exit(main())
