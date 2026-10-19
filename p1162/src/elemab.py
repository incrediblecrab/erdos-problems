"""Exact counts of elementary abelian 2-subgroups of S_n.

a(n, k) = #{ H <= S_n : H is isomorphic to (Z/2)^k }.

Everything here is exact integer arithmetic.  The chain is

  (1) #Hom(F_2^k, S_n), from the exponential formula for group actions;
  (2) #Inj(F_2^k, S_n),  by Moebius inversion over the subspace lattice of F_2^k
      (equivalently q-binomial inversion at q = 2);
  (3) a(n, k) = #Inj(F_2^k, S_n) / |GL_k(2)|, because the map
      (injective hom) -> (its image) is |Aut(F_2^k)|-to-one.

Step (1).  For a finite group G the exponential generating function of the
number of G-actions on an n-set is

      exp( sum_{H <= G} x^[G:H] / [G:H] ),

because a transitive G-action on a labelled m-set is an index-m subgroup
(the point stabiliser) together with (m-1)! labellings of its cosets.  For
G = F_2^k the subgroups of index 2^j are the codimension-j subspaces, of
which there are the Gaussian binomial [k choose j]_2.  Writing the EGF as
exp(sum_m a_m x^m/m!) gives a_m = 0 unless m = 2^j, and

      a_{2^j} = [k choose j]_2 * (2^j - 1)! ,

and exp of an EGF obeys  h(n) = sum_m C(n-1, m-1) a_m h(n-m),  h(0) = 1.
Only j <= log2(n) contribute, so each h(n) costs O(log n) big-integer
multiplications: n can be pushed into the thousands.

Nothing in this module knows about permutations; `verify_elemab.py` checks it
against a brute-force enumeration that shares no code with it.
"""

from __future__ import annotations

import functools
from math import comb, factorial

# ---------------------------------------------------------------- q-binomials


@functools.lru_cache(maxsize=None)
def gaussian_binomial(k: int, j: int, q: int = 2) -> int:
    """[k choose j]_q -- the number of j-dimensional subspaces of F_q^k.

    Also the number of codimension-j subspaces, since [k,j]_q = [k,k-j]_q.
    """
    if j < 0 or j > k:
        return 0
    num = 1
    den = 1
    for i in range(j):
        num *= q ** (k - i) - 1
        den *= q ** (i + 1) - 1
    assert num % den == 0
    return num // den


@functools.lru_cache(maxsize=None)
def gl_order(k: int, q: int = 2) -> int:
    """|GL_k(F_q)| = prod_{t=0}^{k-1} (q^k - q^t).  |GL_0| = 1."""
    out = 1
    for t in range(k):
        out *= q**k - q**t
    return out


@functools.lru_cache(maxsize=None)
def galois_number(r: int, q: int = 2) -> int:
    """G_r = number of subspaces (= subgroups) of F_q^r, all dimensions.

    Uses the one-step recurrence [r,j] = [r,j-1] * (q^{r-j+1}-1)/(q^j-1), so a
    single G_r costs O(r) big-integer operations rather than O(r^2).
    """
    total = 1
    term = 1
    for j in range(1, r + 1):
        term = term * (q ** (r - j + 1) - 1)
        d = q**j - 1
        assert term % d == 0
        term //= d
        total += term
    return total


# ------------------------------------------------------- homomorphism counting


def hom_counts(k: int, n_max: int) -> list[int]:
    """[ #Hom(F_2^k, S_n) for n in 0..n_max ].

    Exponential formula; exact integers.
    """
    # a_m for m = 2^j, j = 0..min(k, log2 n_max)
    weights: dict[int, int] = {}
    j = 0
    while 2**j <= n_max and j <= k:
        m = 2**j
        weights[m] = gaussian_binomial(k, j) * factorial(m - 1)
        j += 1

    h = [0] * (n_max + 1)
    h[0] = 1
    for n in range(1, n_max + 1):
        total = 0
        for m, a_m in weights.items():
            if m <= n:
                total += comb(n - 1, m - 1) * a_m * h[n - m]
        h[n] = total
    return h


def hom_table(k_max: int, n_max: int) -> list[list[int]]:
    """table[k][n] = #Hom(F_2^k, S_n) for 0 <= k <= k_max, 0 <= n <= n_max."""
    return [hom_counts(k, n_max) for k in range(k_max + 1)]


def inj_from_hom(hom: list[list[int]], k: int, n: int) -> int:
    """#Inj(F_2^k, S_n) by q-binomial inversion of #Hom = sum [k,m]_2 #Inj."""
    total = 0
    for m in range(k + 1):
        d = k - m
        term = gaussian_binomial(k, m) * (2 ** comb(d, 2)) * hom[m][n]
        total += -term if d % 2 else term
    return total


def rank_counts(n: int, hom: list[list[int]] | None = None) -> list[int]:
    """[ a(n, k) for k in 0..floor(n/2) ] -- exact.

    a(n, 0) = 1 (the trivial subgroup).  a(n, k) = 0 for k > floor(n/2),
    which the code asserts rather than assumes: it is a sharp structural
    check on the whole inversion.
    """
    k_max = n // 2
    if hom is None:
        hom = hom_table(k_max + 1, n)
    out = []
    for k in range(k_max + 1):
        i_k = inj_from_hom(hom, k, n)
        g = gl_order(k)
        assert i_k >= 0, f"negative Inj count at n={n}, k={k}"
        assert i_k % g == 0, f"|GL_{k}(2)| does not divide Inj at n={n}"
        out.append(i_k // g)
    # structural check: rank floor(n/2)+1 is impossible inside S_n
    if k_max + 1 <= len(hom) - 1:
        over = inj_from_hom(hom, k_max + 1, n)
        assert over == 0, f"nonzero rank-{k_max + 1} count in S_{n}: {over}"
    return out


def elementary_abelian_total(n: int, hom: list[list[int]] | None = None) -> int:
    """e(n) = total number of elementary abelian 2-subgroups of S_n (incl. 1)."""
    return sum(rank_counts(n, hom))


# --------------------------------------------------------------------- driver

if __name__ == "__main__":
    import argparse
    import json
    import pathlib
    import time

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-max", type=int, default=64)
    ap.add_argument("--out", type=pathlib.Path, default=None)
    args = ap.parse_args()

    n_max = args.n_max
    t0 = time.time()
    hom = hom_table(n_max // 2 + 2, n_max)
    t1 = time.time()
    print(f"hom table up to k={n_max // 2 + 1}, n={n_max}: {t1 - t0:.1f}s")

    rows = {}
    for n in range(0, n_max + 1):
        rc = rank_counts(n, hom)
        rows[n] = rc
        tot = sum(rc)
        peak = max(range(len(rc)), key=lambda k: rc[k])
        print(
            f"n={n:4d}  e(n)=2^{tot.bit_length() - 1:<7d} "
            f"argmax k={peak:<4d} (k/n={peak / n if n else 0:.4f})  "
            f"maxrank={len(rc) - 1}"
        )
    print(f"total {time.time() - t0:.1f}s")

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(
            json.dumps({str(n): [str(v) for v in rc] for n, rc in rows.items()})
        )
        print(f"wrote {args.out}")
