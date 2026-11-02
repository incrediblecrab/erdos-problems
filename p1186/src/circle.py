"""Phi_k for a 2-colouring of the circle R/Z given by arbitrary breakpoints.

`exact.py` scores the block colouring x -> c[floor(m*x/p)] of F_p, whose limit
as p -> infinity is

    Phi_k(S) = mu{(x,d) in (R/Z)^2 : x+jd in S   for all j<k}
             + mu{(x,d) in (R/Z)^2 : x+jd in S^c for all j<k},

with S a union of the equal blocks [i/m,(i+1)/m), and the bound is
delta~_k <= Phi_k(S)/2.  Here S is any finite union of arcs: the colouring is
given by its breakpoints

    0 = a_0 < a_1 < ... < a_B = 1,   B even,

arc i = [a_i, a_{i+1}) having colour i mod 2 (B must be even so that the
colours alternate around the circle).  A step function with rational breakpoints
of common denominator M is the equal-M-block word word_of_bounds(bounds, M), so
this family is the union over M of exact.py's; what it adds is that the
breakpoints can be moved continuously, which a grid search cannot do.

Method.  Take (x,d) in the unit square [0,1)^2.  The j-th term x+jd lies in
[0, 1+j), and x+jd mod 1 is in arc i exactly when x+jd lies in
[a_i + w, a_{i+1} + w) for some integer w in {0,...,j}.  Recursively clip the
square by these strips, keeping only strips of one colour, and sum the areas
of the surviving convex polygons.  This is the same decomposition as
blocks.psi_blocks (which handles the triangle of {1,...,n}), with the square
and the wrap offsets w in place of the triangle; the clipper and area routine
are blocks.py's, not exact.py's, so agreement with exact.phi on a word both can
express is a check between two different decompositions.

Works in Fractions (exact) or floats (for search).  This is a search and
cross-check tool: a bound is certified by exact.phi (or exact.phi_bits, above
the table cap) on the equal-block word, and cross-checked by this module
(check F of src/final_check.py), the exact brute-force count of
src/verify_zpx.c and the floating-point one of src/verify_zp.c.
"""

from bisect import bisect_right
from fractions import Fraction
from itertools import groupby

from blocks import area, clip


def _span(poly, j):
    vals = [p[0] + j * p[1] for p in poly]
    return min(vals), max(vals)


def phi_circle(bounds, k):
    """Phi_k of the circle colouring with the given breakpoints (see module doc)."""
    nb = len(bounds) - 1
    assert nb >= 2 and nb % 2 == 0, "need an even number of arcs"
    one = bounds[-1]
    zero = one - one
    square = [(zero, zero), (one, zero), (one, one), (zero, one)]
    # extended breakpoints over [0, k): arc index i + nb*w is [a_i + w, a_{i+1} + w)
    ext = [b + w for w in range(k) for b in bounds[:-1]] + [one * k]
    total = zero

    def rec(j, poly, parity):
        nonlocal total
        if j == k:
            total += area(poly)
            return
        lo, hi = _span(poly, j)
        first = max(0, bisect_right(ext, lo) - 1)
        last = min(len(ext) - 1, bisect_right(ext, hi))
        for g in range(first, last):
            if g % 2 != parity:          # nb even, so arc colour is g mod 2
                continue
            left, right = ext[g], ext[g + 1]
            if right <= lo or left >= hi:
                continue
            q = poly
            if left > lo:
                q = clip(q, -1, -j, -left)       # x + j d >= left
                if len(q) < 3:
                    continue
            if right < hi:
                q = clip(q, 1, j, right)         # x + j d <= right
                if len(q) < 3:
                    continue
            rec(j + 1, q, parity)

    for i0 in range(nb):
        p = clip(square, -1, 0, -bounds[i0]) if bounds[i0] > zero else square
        if bounds[i0 + 1] < one:
            p = clip(p, 1, 0, bounds[i0 + 1])
        if len(p) >= 3:
            rec(1, p, i0 % 2)
    return total


def _arcs(bounds, colour, shift):
    """Disjoint sorted intervals in [0,1) of colour `colour`, translated by -shift mod 1."""
    one = bounds[-1]
    out = []
    for i in range(len(bounds) - 1):
        if i % 2 != colour:
            continue
        lo = (bounds[i] - shift) % one
        hi = lo + (bounds[i + 1] - bounds[i])
        if hi <= one:
            out.append((lo, hi))
        else:
            out.append((lo, one))
            out.append((lo - lo, hi - one))
    out.sort()
    return out


def _meet(A, B):
    out, i, j = [], 0, 0
    while i < len(A) and j < len(B):
        lo, hi = max(A[i][0], B[j][0]), min(A[i][1], B[j][1])
        if lo < hi:
            out.append((lo, hi))
        if A[i][1] < B[j][1]:
            i += 1
        else:
            j += 1
    return out


def _g(bounds, k, d):
    """Measure of x in R/Z with x + jd (j < k) all in one colour class."""
    tot = 0
    for colour in (0, 1):
        cur = _arcs(bounds, colour, 0 * d)
        for j in range(1, k):
            cur = _meet(cur, _arcs(bounds, colour, j * d))
            if not cur:
                break
        tot += sum(hi - lo for lo, hi in cur)
    return tot


def phi_sweep(bounds, k):
    """Phi_k by a one-dimensional sweep over d, a second exact method.

    Phi_k = integral over d in [0,1) of g(d), with g(d) the measure of the x
    whose k terms x + jd share a colour.  The x-intervals that make up g(d)
    have endpoints a_i - jd mod 1, and their cyclic order changes only where
    two endpoints meet, a_i - jd = a_i' - j'd mod 1, i.e. at
    d = (a_i' - a_i + w)/(j - j') mod 1.  Between consecutive such d the order
    is fixed and g is linear; g is continuous; so the trapezoid rule over the
    critical points is exact.  This shares no geometry with phi_circle (no
    polygons, no clipping) and accepts any rational breakpoints, not only
    those on a common grid.
    """
    nb = len(bounds) - 1
    assert nb >= 2 and nb % 2 == 0
    one = bounds[-1]
    zero = one - one
    crit = {zero, one}
    pts = bounds[:-1]
    for q in range(1, k):
        for ai in pts:
            for bi in pts:
                diff = (ai - bi) % one
                for w in range(q):
                    crit.add((diff + w) / q)
    crit = sorted(crit)
    vals = [_g(bounds, k, d) for d in crit]
    return sum((crit[i + 1] - crit[i]) * (vals[i] + vals[i + 1])
               for i in range(len(crit) - 1)) / 2


def bounds_of_cyclic_word(word):
    """Breakpoints (Fractions) of an equal-block word, rotated to start a run.

    Rotating the circle preserves Phi_k, and complementing the colours does too,
    so the word is rotated to begin at a colour change and complemented to start
    with colour 0.  Returns (bounds, shift) with shift the rotation in blocks.
    """
    c = [int(ch) for ch in word]
    m = len(c)
    shift = next(i for i in range(m) if c[i] != c[i - 1])
    c = c[shift:] + c[:shift]
    if c[0] == 1:
        c = [1 - x for x in c]
    runs = [len(list(g)) for _, g in groupby(c)]
    assert len(runs) % 2 == 0
    cuts, acc = [Fraction(0)], 0
    for r in runs:
        acc += r
        cuts.append(Fraction(acc, m))
    return cuts, shift


def word_of_lengths(lengths):
    """Equal-block word of integer arc lengths (colour 0 first)."""
    assert len(lengths) % 2 == 0 and all(L > 0 for L in lengths)
    return "".join(str(i % 2) * L for i, L in enumerate(lengths))
