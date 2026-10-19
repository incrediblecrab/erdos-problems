#!/usr/bin/env python
"""The quadratic colouring family for the finite-field problem (Erdos #1186, k=4).

Motivation.  Block colourings x -> c[floor(m*x/p)] are pull-backs along the
*linear* map x -> x/p.  Four-term progressions are not controlled by linear
(Fourier) structure -- they need the quadratic, U^3 layer -- so the natural next
family is the pull-back along the *quadratic* map x -> x^2/p:

    colour(x) = G[ floor( m * (x^2 mod p)/p ) ].

Why this is exactly computable.  For a progression x+jd,

    (x+jd)^2 / p  =  gamma + 2 j delta + j^2 epsilon   (mod 1),

with gamma = x^2/p, delta = xd/p, epsilon = d^2/p.  For every nonzero integer
triple (h3,h4,h5) the Weyl sum  sum_{x,d} e((h3 x^2 + h4 xd + h5 d^2)/p)  is a
two-variable Gauss sum of size O(p), hence o(p^2), so (gamma,delta,epsilon)
equidistributes in (R/Z)^3 as p -> infinity.

Write y_j = gamma + 2 j delta + j^2 epsilon.  This is a quadratic in j, so it is
determined by y_0, y_1, y_2 through Lagrange interpolation at the nodes 0,1,2:

    y_j = L0(j) y_0 + L1(j) y_1 + L2(j) y_2,
    L0(j) = (j-1)(j-2)/2,  L1(j) = j(2-j),  L2(j) = j(j-1)/2.

For k = 4 the only dependent point is y_3 = y_0 - 3 y_1 + 3 y_2.  The map
(gamma,delta,epsilon) -> (y_0,y_1,y_2) is a surjective homomorphism of tori, so
(y_0,y_1,y_2) is uniform on (R/Z)^3 and

    Theta_4(G) = integral over [0,1)^3 of
                 [ G is constant on y_0, y_1, y_2, y_0 - 3y_1 + 3y_2 ]  dy,

giving  delta~_4 <= Theta_4(G)/2  by the same unordered / non-degenerate count as
in the block family.

Cell decomposition.  With y_i = (a_i + s_i)/m, a_i in Z_m and s_i in [0,1), the
fourth block index is (a_0 - 3a_1 + 3a_2 + e) mod m where e = floor(s_0 - 3s_1 +
3s_2) in {-3,...,3}.  The weight of a cell is the exact volume

    V(e) = vol{ s in [0,1]^3 : e <= s_0 - 3s_1 + 3s_2 < e+1 },

computed in closed form below, so the whole table is exact rational arithmetic.
"""

import sys
from fractions import Fraction
from functools import lru_cache
from math import gcd
from pathlib import Path


def cube_cdf(weights, t):
    """vol{ s in [0,1]^n : sum_i w_i s_i <= t } exactly, for arbitrary real w_i.

    Negative weights are handled by the substitution s_i -> 1 - s_i, which maps
    the cube to itself.  For positive weights the volume is the standard
    inclusion-exclusion formula
        (1/(n! prod w_i)) * sum_{S} (-1)^{|S|} max(0, t - sum_{i in S} w_i)^n.
    """
    n = len(weights)
    t = Fraction(t)
    pos = []
    for w in weights:
        w = Fraction(w)
        if w < 0:
            t -= w          # s_i -> 1 - s_i turns w_i s_i into w_i - w_i s_i
            pos.append(-w)
        elif w > 0:
            pos.append(w)
        else:
            raise ValueError("zero weight")
    if t <= 0:
        return Fraction(0)
    fac = 1
    for i in range(2, n + 1):
        fac *= i
    denom = Fraction(fac)
    for w in pos:
        denom *= w
    total = Fraction(0)
    for mask in range(1 << n):
        s = t
        bits = 0
        for i in range(n):
            if mask >> i & 1:
                s -= pos[i]
                bits += 1
        if s > 0:
            total += (-1) ** bits * s ** n
    v = total / denom
    assert 0 <= v <= 1, (weights, t, v)
    return v


@lru_cache(maxsize=None)
def carry_volumes_k4():
    """Exact ((e, volume)) for e = floor(s_0 - 3 s_1 + 3 s_2) on [0,1)^3."""
    w = (1, -3, 3)
    out = []
    for e in range(-3, 4):
        v = cube_cdf(w, e + 1) - cube_cdf(w, e)
        if v > 0:
            out.append((e, v))
    assert sum(v for _, v in out) == 1, "carry volumes must sum to 1"
    return tuple(out)


@lru_cache(maxsize=None)
def quad_weights_int(m):
    """Integer-scaled exact weights for the quadratic family at k = 4.

    Returns (den, items) with items a tuple of (sorted index tuple, int weight)
    and sum(weights) == den, so Theta_4(G) = (matched sum)/den.
    """
    vols = carry_volumes_k4()
    cell_den = 1
    for _, v in vols:
        cell_den = cell_den * v.denominator // gcd(cell_den, v.denominator)
    cells = [(e, int(v * cell_den)) for e, v in vols]
    acc = {}
    for a0 in range(m):
        for a1 in range(m):
            base = a0 - 3 * a1
            for a2 in range(m):
                b = base + 3 * a2
                for e, w in cells:
                    U = tuple(sorted({a0, a1, a2, (b + e) % m}))
                    acc[U] = acc.get(U, 0) + w
    den = cell_den * m * m * m
    items = tuple(sorted(acc.items()))
    assert sum(w for _, w in items) == den, "quadratic weights must sum to 1"
    return den, items


def theta(c, items=None, den=None):
    """Exact Theta_4 of the quadratic colouring G given by the word c."""
    if items is None:
        den, items = quad_weights_int(len(c))
    tot = 0
    for U, w in items:
        first = c[U[0]]
        if all(c[i] == first for i in U[1:]):
            tot += w
    return Fraction(tot, den)


def export(m, path):
    den, items = quad_weights_int(m)
    with open(path, "w") as f:
        f.write(f"{m} 4 {len(items)} {den}\n")
        for idx, w in items:
            f.write(f"{len(idx)} {w} " + " ".join(map(str, idx)) + "\n")
    return len(items), den


if __name__ == "__main__":
    if len(sys.argv) >= 3 and sys.argv[1] == "export":
        m, out = int(sys.argv[2]), sys.argv[3]
        n, den = export(m, Path(out))
        print(f"m={m} quadratic sets={n} den={den} -> {out}")
    else:
        print("carry volumes:", [(e, str(v)) for e, v in carry_volumes_k4()])
