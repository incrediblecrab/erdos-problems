"""Exact rational evaluation of the monochromatic k-AP density of a block colouring.

Setting.  Erdos #1186 asks for delta_k, the least constant such that every
2-colouring of {1,...,n} contains at least (delta_k+o(1)) n^2 monochromatic
k-term APs, and for the finite-field analogue delta~_k in F_p.

Upper bounds on delta~_k come from explicit colourings.  The family used here is
the *block* (equivalently, interval / one-frequency Bohr) family: fix m and a
binary word c in {0,1}^m, and colour x in Z_p by c[floor(m*x/p)].  Equivalently,
colour by the indicator of S = union of [i/m,(i+1)/m) over {i : c_i = 1},
pulled back along x -> x/p in R/Z.

Because x -> x/p carries a k-AP in Z_p to a k-AP in R/Z exactly, and because the
p^2 pairs (x/p, d/p) equidistribute in (R/Z)^2 with discrepancy O(1/p) against
regions bounded by finitely many rational lines, the density of ordered
monochromatic pairs (x,d) converges to

    Phi_k(S) = mu{(x,d) : x+jd in S   for all j<k}
             + mu{(x,d) : x+jd in S^c for all j<k}.

Counting each non-degenerate AP once (the pairs (x,d) and (x+(k-1)d,-d) give the
same AP) and discarding the p degenerate d=0 pairs gives

    #mono k-APs = Phi_k(S)*p^2/2 + O(p),      hence   delta~_k <= Phi_k(S)/2.

This module computes Phi_k(c) in exact rational arithmetic.

Method.  Substitute X = m*x, D = m*d, then X = a+s and D = b+t with integers
a,b in [0,m) and s,t in [0,1).  Then

    floor(m*(x+jd) mod m) = (a + j*b + e_j) mod m,      e_j := floor(s + j*t).

e_0 = 0 and e_j in {0,...,j}.  The vector e is constant on each cell of the
arrangement of the lines s + j*t = r inside the unit square, and every cell is a
convex polygon with rational vertices, so its area is an exact rational.  Hence

    Phi_k(c) = (1/m^2) * sum_{a,b in Z_m} sum_e area(e) * mono(c, a, b, e),

a finite sum of rationals.  Whether an AP is monochromatic depends only on the
*set* of block indices it visits, so the weights are aggregated by that set.

No floating point is used anywhere in this file.
"""

from fractions import Fraction
from itertools import product
from functools import lru_cache
from math import gcd

ZERO = Fraction(0)
ONE = Fraction(1)
UNIT_SQUARE = [(ZERO, ZERO), (ONE, ZERO), (ONE, ONE), (ZERO, ONE)]


def clip_halfplane(poly, A, B, C):
    """Clip a convex polygon by the closed half-plane A*s + B*t + C <= 0.

    poly is a list of (Fraction, Fraction) vertices in order.  Sutherland-Hodgman,
    run entirely in exact rational arithmetic.
    """
    if not poly:
        return []
    out = []
    n = len(poly)
    for i in range(n):
        px, py = poly[i]
        qx, qy = poly[(i + 1) % n]
        fp = A * px + B * py + C
        fq = A * qx + B * qy + C
        if fp <= 0:
            out.append((px, py))
        if (fp < 0 < fq) or (fq < 0 < fp):
            lam = fp / (fp - fq)
            out.append((px + lam * (qx - px), py + lam * (qy - py)))
    return out


def polygon_area(poly):
    """Exact area of a simple polygon by the shoelace formula."""
    if len(poly) < 3:
        return ZERO
    acc = ZERO
    n = len(poly)
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        acc += x1 * y2 - x2 * y1
    return abs(acc) / 2


@lru_cache(maxsize=None)
def carry_regions(k):
    """Cells of the unit square on which e_j = floor(s + j*t) is constant.

    Returns a tuple of (e_vector, exact_area) pairs, only those of positive area.
    The areas sum to exactly 1.
    """
    out = []
    for e in product(*[range(j + 1) for j in range(k)]):
        poly = UNIT_SQUARE
        for j in range(1, k):
            # e_j <= s + j*t
            poly = clip_halfplane(poly, Fraction(-1), Fraction(-j), Fraction(e[j]))
            # s + j*t <= e_j + 1
            poly = clip_halfplane(poly, ONE, Fraction(j), Fraction(-(e[j] + 1)))
            if len(poly) < 3:
                poly = []
                break
        if poly:
            a = polygon_area(poly)
            if a > 0:
                out.append((e, a))
    assert sum(a for _, a in out) == 1, "carry-region areas must sum to 1"
    return tuple(out)


@lru_cache(maxsize=None)
def block_weights_int(m, k, mult=1):
    """Integer-scaled exact weights.  Returns (den, items) with items a tuple of
    (sorted index tuple, integer weight) and sum(weights) == den, so that

        Phi_k(c) = (sum of weights over sets on which c is constant) / den.

    Scaling by the common denominator keeps everything in machine integers
    without leaving exact arithmetic.

    `mult` selects which of two families of colourings of Z_p is being scored.
    In both, the k-AP index of the j-th term is

        i_j = (a + j*b + mult*e_j) mod m,       e_j = floor(s + j*t),

    with (a,b) uniform on Z_m^2 and (s,t) uniform on [0,1)^2, independent.

      mult = 1   is the block colouring  x -> c[floor(m*x/p)].  Here (a,s) and
                 (b,t) are the integer and fractional parts of m*x/p and m*d/p,
                 the carry e_j is the wrap of the j-th term past a block edge,
                 and the colouring is defined for *every* prime p.

      mult = -r  is the periodic colouring  x -> c[x mod m] with p = r (mod m).
                 Reducing x+jd mod p and then mod m subtracts p once per wrap of
                 (x+jd)/p, so the index picks up -r*e_j.  This is Lu-Peng's
                 family.  It needs gcd(r,m) = 1 for infinitely many such primes
                 to exist (Dirichlet), so gcd(mult,m) = 1; the resulting bound
                 is on liminf_p, i.e. it holds along that residue class, not for
                 every large p.

    mult = 1 and mult = m-1 are therefore the same family up to e -> -e, which
    is why the m = 20 optimum here coincides with Lu-Peng's B_20.
    """
    regs = carry_regions(k)
    cell_den = 1
    for _, a in regs:
        cell_den = cell_den * a.denominator // gcd(cell_den, a.denominator)
    cells = [(e, int(a * cell_den)) for e, a in regs]
    acc = {}
    for a in range(m):
        for b in range(m):
            for e, w in cells:
                U = tuple(sorted({(a + j * b + mult * e[j]) % m for j in range(k)}))
                acc[U] = acc.get(U, 0) + w
    den = cell_den * m * m
    items = tuple(sorted(acc.items()))
    assert sum(w for _, w in items) == den, "block weights must sum to 1"
    return den, items


@lru_cache(maxsize=None)
def block_weights(m, k, mult=1):
    """Aggregated exact weights for the block family.

    Returns a tuple of (index_set, weight) with

        Phi_k(c) = sum over pairs of  weight * [ c is constant on index_set ].

    index_set is a frozenset of block indices in Z_m; weight is a Fraction.
    The weights sum to exactly 1.
    """
    den, items = block_weights_int(m, k, mult)
    out = tuple((frozenset(U), Fraction(w, den)) for U, w in items)
    assert sum(w for _, w in out) == 1, "block weights must sum to 1"
    return out


@lru_cache(maxsize=None)
def zm_weights_int(m, k):
    """Integer-scaled exact weights for the periodic family over Z_m.

    Colour i in {1,...,n} by c[i mod m].  For an arithmetic progression
    (a, a+d, ..., a+(k-1)d) inside {1,...,n} the colours depend only on
    (a mod m, d mod m), and each residue pair occurs equally often up to
    O(n), so the number of monochromatic k-APs is Psi_k(c) n^2 + O(n) with

        Psi_k(c) = Phi^{Z_m}_k(c) / (2(k-1)),
        Phi^{Z_m}_k(c) = (1/m^2) #{(a,b) in Z_m^2 : c constant on a+jb, j<k}.

    The 1/(2(k-1)) is the area of the triangle {x,t >= 0, x+(k-1)t <= 1}, i.e.
    the density of k-APs in {1,...,n}; it is the same normalisation used by
    triangle_weights_int, so the two families are directly comparable.

    Note Phi^{Z_m}_k >= 1/m always, because every b = 0 pair is monochromatic;
    equality says exactly that Z_m has a 2-colouring with no monochromatic
    k-AP of nonzero common difference.

    Returns (den, items) with den = m^2; Phi = (matched sum)/den.
    """
    acc = {}
    for a in range(m):
        for b in range(m):
            U = tuple(sorted({(a + j * b) % m for j in range(k)}))
            acc[U] = acc.get(U, 0) + 1
    items = tuple(sorted(acc.items()))
    assert sum(w for _, w in items) == m * m, "Z_m weights must sum to m^2"
    return m * m, items


def phi_zm(c, k):
    """Exact Phi^{Z_m}_k of the periodic colouring c (length m, entries 0/1)."""
    m = len(c)
    den, items = zm_weights_int(m, k)
    tot = 0
    for U, w in items:
        first = c[U[0]]
        if all(c[i] == first for i in U[1:]):
            tot += w
    return Fraction(tot, den)


def psi_zm(c, k):
    """Exact Psi_k (the delta_k upper bound) of x -> c[x mod m]."""
    return phi_zm(c, k) / (2 * (k - 1))


@lru_cache(maxsize=None)
def carry_polys(k):
    """As carry_regions, but returning the clipped polygon as well."""
    out = []
    for e in product(*[range(j + 1) for j in range(k)]):
        poly = UNIT_SQUARE
        for j in range(1, k):
            poly = clip_halfplane(poly, Fraction(-1), Fraction(-j), Fraction(e[j]))
            poly = clip_halfplane(poly, ONE, Fraction(j), Fraction(-(e[j] + 1)))
            if len(poly) < 3:
                poly = []
                break
        if poly:
            a = polygon_area(poly)
            if a > 0:
                out.append((e, poly, a))
    return tuple(out)


@lru_cache(maxsize=None)
def triangle_weights_int(m, k):
    """Integer-scaled exact weights for the {1,...,n} setting.

    Colour i in {1,...,n} by c[floor(m*i/n)].  Writing x = a/n and t = d/n, the
    k-APs (a, a+d, ..., a+(k-1)d) inside {1,...,n} correspond to the triangle

        T = { (x,t) : x >= 0, t >= 0, x + (k-1)t <= 1 },   area 1/(2(k-1)),

    and the number of monochromatic k-APs is Psi_k(c)*n^2 + O(n) where

        Psi_k(c) = integral over T of [ c constant on blocks of x+jt, j<k ].

    There is no wrap-around here: the block index of x+jt is the honest integer
    floor(m(x+jt)) = a + j*b + e_j with e_j = floor(s+j*u) as before, but the
    cell (a,b) must additionally be clipped to the triangle.  Since the clipping
    line s + (k-1)u = m - a - (k-1)b has integer coefficients, all areas stay
    rational and are computed exactly.

    Returns (den, items) with items a tuple of (sorted index tuple, int weight);
    Psi_k(c) = (matched sum)/den.  The weights sum to den/(2(k-1)), not den.
    """
    polys = carry_polys(k)
    cell_den = 1
    for _, _, a in polys:
        cell_den = cell_den * a.denominator // gcd(cell_den, a.denominator)
    acc = {}
    extra_den = 1
    raw = []
    for a in range(m):
        bmax = (m - a) // (k - 1)
        for b in range(bmax + 1):
            R = m - a - (k - 1) * b          # s + (k-1)u <= R
            if R <= 0:
                continue
            for e, poly, full in polys:
                if R >= k:                    # clip is inactive
                    ar = full
                else:
                    q = clip_halfplane(poly, ONE, Fraction(k - 1), Fraction(-R))
                    if len(q) < 3:
                        continue
                    ar = polygon_area(q)
                    if ar == 0:
                        continue
                U = tuple(sorted({a + j * b + e[j] for j in range(k)}))
                raw.append((U, ar))
                extra_den = extra_den * ar.denominator // gcd(extra_den, ar.denominator)
    scale = extra_den
    for U, ar in raw:
        v = ar * scale
        assert v.denominator == 1
        acc[U] = acc.get(U, 0) + int(v)
    den = scale * m * m
    items = tuple(sorted(acc.items()))
    assert Fraction(sum(w for _, w in items), den) == Fraction(1, 2 * (k - 1)), \
        "triangle weights must sum to the triangle area 1/(2(k-1))"
    assert max(max(U) for U, _ in items) < m, "block index out of range"
    return den, items


@lru_cache(maxsize=None)
def triangle_cells_int(B, k):
    """As triangle_weights_int, but keeping the ORDERED tuple of block indices.

    triangle_weights_int deduplicates the k block indices into a set, which is
    all that is needed when the colour depends on the block alone.  For the
    two-scale family the colour of the j-th term of the progression depends on
    the pair (fine residue, coarse block), so the indices must stay in order.

    Returns (den, items) with items a tuple of (ordered k-tuple, int weight);
    the weights sum to den/(2(k-1)), the area of the triangle.
    """
    polys = carry_polys(k)
    acc = {}
    raw = []
    extra_den = 1
    for a in range(B):
        bmax = (B - a) // (k - 1)
        for b in range(bmax + 1):
            R = B - a - (k - 1) * b
            if R <= 0:
                continue
            for e, poly, full in polys:
                if R >= k:
                    ar = full
                else:
                    q = clip_halfplane(poly, ONE, Fraction(k - 1), Fraction(-R))
                    if len(q) < 3:
                        continue
                    ar = polygon_area(q)
                    if ar == 0:
                        continue
                idx = tuple(a + j * b + e[j] for j in range(k))
                raw.append((idx, ar))
                extra_den = extra_den * ar.denominator // gcd(extra_den, ar.denominator)
    for idx, ar in raw:
        v = ar * extra_den
        assert v.denominator == 1
        acc[idx] = acc.get(idx, 0) + int(v)
    den = extra_den * B * B
    items = tuple(sorted(acc.items()))
    assert Fraction(sum(w for _, w in items), den) == Fraction(1, 2 * (k - 1)), \
        "triangle cells must sum to the triangle area 1/(2(k-1))"
    assert max(max(t) for t, _ in items) < B, "block index out of range"
    return den, items


@lru_cache(maxsize=None)
def two_scale_weights_int(m, B, k):
    """Weights for the two-scale family on {1,...,n}.

    Colour i by C[i mod m][floor(B*i/n)], i.e. by a fine residue mod m and a
    coarse position in one of B equal blocks of {1,...,n}.  The colour table C
    is flattened to a single index  r * B + q  with r in Z_m, q in [B].

    For a progression (a, a+d, ..., a+(k-1)d) the residues (a mod m, d mod m)
    are asymptotically uniform on Z_m^2 and asymptotically independent of the
    scaled position (a/n, d/n) in the triangle, because for each fixed residue
    pair the admissible (a,d) form a sublattice of density 1/m^2 that
    equidistributes.  So the weights are the product of the Z_m enumeration
    and the triangle-cell enumeration:

        Psi_k(C) = sum over (alpha,beta) in Z_m^2 and triangle cells (i_0..i_k-1)
                   of  (area/m^2) * [ C constant on { (alpha+j*beta, i_j) } ].

    This family contains the block colourings (m = 1), the periodic colourings
    (B = 1), and products of the two.

    Returns (den, items) with items a tuple of (sorted index tuple, int weight).
    """
    cden, cells = triangle_cells_int(B, k)
    acc = {}
    for alpha in range(m):
        for beta in range(m):
            res = [(alpha + j * beta) % m for j in range(k)]
            for idx, w in cells:
                U = tuple(sorted({res[j] * B + idx[j] for j in range(k)}))
                acc[U] = acc.get(U, 0) + w
    den = cden * m * m
    items = tuple(sorted(acc.items()))
    assert Fraction(sum(w for _, w in items), den) == Fraction(1, 2 * (k - 1)), \
        "two-scale weights must sum to the triangle area 1/(2(k-1))"
    assert max(max(U) for U, _ in items) < m * B, "flattened index out of range"
    return den, items


def psi_two_scale(C, m, B, k):
    """Exact Psi_k for the two-scale colouring; C is a flat 0/1 list of length m*B."""
    assert len(C) == m * B
    den, items = two_scale_weights_int(m, B, k)
    tot = 0
    for U, w in items:
        first = C[U[0]]
        if all(C[i] == first for i in U[1:]):
            tot += w
    return Fraction(tot, den)


def psi(c, k, items=None, den=None):
    """Exact Psi_k of the block colouring c in the {1,...,n} setting."""
    if items is None:
        den, items = triangle_weights_int(len(c), k)
    tot = 0
    for U, w in items:
        first = c[U[0]]
        if all(c[i] == first for i in U[1:]):
            tot += w
    return Fraction(tot, den)


def phi(c, k, weights=None, mult=1):
    """Exact Phi_k of the block colouring c (a sequence of 0/1 of length m).

    See block_weights_int for the meaning of `mult`; mult=1 is the block
    colouring x -> c[floor(m*x/p)], valid for every prime p.
    """
    m = len(c)
    if weights is None:
        weights = block_weights(m, k, mult)
    tot = ZERO
    for U, w in weights:
        it = iter(U)
        first = c[next(it)]
        if all(c[i] == first for i in it):
            tot += w
    return tot


def delta_tilde_bound(c, k, weights=None):
    """The upper bound on delta~_k certified by the block colouring c: Phi_k(c)/2."""
    return phi(c, k, weights) / 2


def lambda_form(c, k, weights=None):
    """Exact Lambda_k(f) = E_{x,d} prod_{j<k} f(x+jd) for f = 2*1_S - 1.

    Derived from Phi via  Phi_k = 2^{1-k} * sum over even subsets T of Lambda_T,
    but computed here directly so that the two routes can be cross-checked.
    """
    m = len(c)
    regs = carry_regions(k)
    f = [1 if ci else -1 for ci in c]
    tot = ZERO
    for a in range(m):
        for b in range(m):
            for e, ar in regs:
                prod = 1
                for j in range(k):
                    prod *= f[(a + j * b + e[j]) % m]
                if prod > 0:
                    tot += ar
                else:
                    tot -= ar
    return tot / Fraction(m * m)


def mean(c):
    """Exact mean of f = 2*1_S - 1."""
    m = len(c)
    return Fraction(2 * sum(c) - m, m)
