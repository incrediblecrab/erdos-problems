"""Continuous block-boundary evaluator for Psi_k.

`exact.py` computes Psi_k for colourings that are constant on the m equal
blocks [i/m, (i+1)/m).  That forces every block boundary onto a multiple of
1/m, which is the binding constraint once the grid search gets close to the
Parrilo-Robertson-Saracino optimum: their 12 block lengths over 548 have gcd 1,
so a coarse grid cannot even represent them, and a fine grid costs O(m^2)
rational areas to tabulate.

Here the colouring is given directly by its boundaries

    0 = a_0 < a_1 < ... < a_B = 1,

block i being [a_i, a_{i+1}) with colour i mod 2 (consecutive blocks must
differ in colour, otherwise merge them).  Psi_k is then

    sum over monochromatic index tuples (i_0,...,i_{k-1}) of
        area { (x,t) in T : a_{i_j} <= x + j t < a_{i_j + 1} for all j }

with T = {x >= 0, t >= 0, x + (k-1) t <= 1}, exactly the region and
normalisation used by exact.triangle_weights_int.  Each summand is a convex
polygon, so the whole thing is computed by clipping.

This is a from-scratch second implementation: it shares no code with
exact.py -- different decomposition (continuous boundaries instead of carry
regions over a grid), its own clipper and its own area routine.  Agreement
between the two on a colouring both can express is therefore a real check.

Works in Fractions (exact, for certification) or floats (fast, for search).
"""

from fractions import Fraction
from math import gcd


# --------------------------------------------------------------------------
# convex polygon geometry (self-contained; nothing imported from exact.py)
# --------------------------------------------------------------------------

def clip(poly, n0, n1, c):
    """Intersect the convex polygon with the half-plane n0*x + n1*y <= c.

    `poly` is a list of (x, y) in counter-clockwise order.  Returns the new
    vertex list, possibly empty.  Exact when the inputs are Fractions.
    """
    if not poly:
        return poly
    out = []
    npts = len(poly)
    prev = poly[-1]
    dprev = n0 * prev[0] + n1 * prev[1] - c
    for idx in range(npts):
        cur = poly[idx]
        dcur = n0 * cur[0] + n1 * cur[1] - c
        if dcur <= 0:
            if dprev > 0:
                s = dprev / (dprev - dcur)
                out.append((prev[0] + s * (cur[0] - prev[0]),
                            prev[1] + s * (cur[1] - prev[1])))
            out.append(cur)
        elif dprev <= 0:
            s = dprev / (dprev - dcur)
            out.append((prev[0] + s * (cur[0] - prev[0]),
                        prev[1] + s * (cur[1] - prev[1])))
        prev, dprev = cur, dcur
    return out


def area(poly):
    """Shoelace area of a convex polygon given counter-clockwise."""
    if len(poly) < 3:
        return 0
    s = 0
    prev = poly[-1]
    for cur in poly:
        s += prev[0] * cur[1] - cur[0] * prev[1]
        prev = cur
    return s / 2 if s >= 0 else -s / 2


def _span(poly, j):
    """Range of x + j*t over the vertices of a convex polygon."""
    vals = [p[0] + j * p[1] for p in poly]
    return min(vals), max(vals)


# --------------------------------------------------------------------------
# Psi_k for a block colouring given by its boundaries
# --------------------------------------------------------------------------

def psi_blocks(bounds, k):
    """Psi_k for the colouring with the given boundaries.

    bounds: a sequence 0 = a_0 < a_1 < ... < a_B = 1 of Fractions or floats.
    Block i is [a_i, a_{i+1}), coloured i mod 2.

    Returns the exact area (a Fraction if the bounds are Fractions).
    """
    nb = len(bounds) - 1
    assert nb >= 1
    one = bounds[-1]
    zero = one - one                                  # 0 in the same type
    tri = [(zero, zero), (one, zero), (zero, one / (k - 1))]

    total = zero

    def rec(j, poly, parity):
        nonlocal total
        if j == k:
            total += area(poly)
            return
        lo, hi = _span(poly, j)
        for i in range(nb):
            if i % 2 != parity:
                continue
            if bounds[i + 1] <= lo or bounds[i] >= hi:
                continue
            q = poly
            if bounds[i] > lo:
                q = clip(q, -1, -j, -bounds[i])       # x + j t >= a_i
                if len(q) < 3:
                    continue
            if bounds[i + 1] < hi:
                q = clip(q, 1, j, bounds[i + 1])      # x + j t <= a_{i+1}
                if len(q) < 3:
                    continue
            rec(j + 1, q, parity)

    for i0 in range(nb):
        p = tri
        lo, hi = _span(tri, 0)
        if bounds[i0] > lo:
            p = clip(p, -1, 0, -bounds[i0])
        if bounds[i0 + 1] < hi:
            p = clip(p, 1, 0, bounds[i0 + 1])
        if len(p) >= 3:
            rec(1, p, i0 % 2)
    return total


# --------------------------------------------------------------------------
# conversions
# --------------------------------------------------------------------------

def bounds_of_word(word):
    """Boundaries (Fractions over m) of the grid colouring given by `word`.

    `word` is a string/sequence of 0/1 of length m.  Runs of equal colour are
    merged.  The leading colour is forced to 0 by complementing if needed;
    Psi_k is invariant under swapping the two colours.
    """
    c = [int(x) for x in word]
    m = len(c)
    if c[0] == 1:
        c = [1 - x for x in c]
    cuts = [0]
    for i in range(1, m):
        if c[i] != c[i - 1]:
            cuts.append(i)
    cuts.append(m)
    return [Fraction(x, m) for x in cuts]


def bounds_of_lengths(lengths):
    """Boundaries from consecutive block lengths (integers), normalised."""
    tot = sum(lengths)
    acc, out = 0, [Fraction(0)]
    for L in lengths:
        acc += L
        out.append(Fraction(acc, tot))
    return out


def lengths_of_bounds(bounds):
    """Integer block lengths and their common denominator."""
    den = 1
    for b in bounds:
        den = den * b.denominator // gcd(den, b.denominator)
    cuts = [int(b * den) for b in bounds]
    return [cuts[i + 1] - cuts[i] for i in range(len(cuts) - 1)], den


def word_of_bounds(bounds, m):
    """Round the boundaries onto a grid of m equal blocks, as a 0/1 string."""
    cuts = [round(float(b) * m) for b in bounds]
    cuts[0], cuts[-1] = 0, m
    out = []
    for i in range(len(cuts) - 1):
        out.append(str(i % 2) * max(0, cuts[i + 1] - cuts[i]))
    s = "".join(out)
    assert len(s) == m, (len(s), m)
    return s
