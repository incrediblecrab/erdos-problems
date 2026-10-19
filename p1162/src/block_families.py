"""Exact counts of Frattini-closed 2-subgroups of S_n with every orbit of size <= 8.

Let H <= S_n be a 2-group with orbits O_1..O_t and restrictions P_i = H|O_i,
each a transitive 2-group.  Call H Frattini-closed if H >= N := prod_i Phi(P_i).
Then H/N is a subspace V of prod_i P_i/Phi(P_i) = F_2^D, D = sum_i d(P_i), whose
projection onto every block F_2^{d(P_i)} is onto, because a subgroup of P_i that
maps onto P_i/Phi(P_i) is P_i (Phi(P_i) is the set of non-generators).
Conversely every such V is H/N for exactly one such H, and H determines its
orbits, its P_i, and N.  Hence

  #{Frattini-closed H with orbit data (O_i, P_i)} = #{V <= F_2^D, onto every block}
      = sum_{U_i <= F_2^{d_i}} prod_i mu(U_i, F_2^{d_i}) G_{sum_i dim U_i}

by Moebius inversion on the subspace lattice, G_r being the Galois number.
By the q-binomial theorem the block weight sum_U mu(U, F_2^d) z^{dim U} is
w_d(z) = prod_{i<d} (z - 2^i), and the exponential formula gives

  B(n) = n! [x^n] L( exp( sum_b W_b(z) x^b/b! ) ),   W_b = sum_{P transitive on [b]} w_{d(P)},

with L the linear map z^D -> G_D.  |H| = 2^{dim V} prod_i |Phi(P_i)|, and
sum_k k [D,k]_2 = (D/2) G_D by the symmetry [D,k] = [D,D-k], so the mean of
log2|H| comes from the same recurrence carried with dual numbers.

The only input is the census of transitive 2-groups of degree 1, 2, 4, 8
written by src/census8.py, which counts the same family by brute force on
permutations for n <= 8.  Gate 1 below checks that the two agree.

Every Frattini-closed family counted here is a set of distinct subgroups of
S_n, so each count is a rigorous lower bound for f(n) = |Sub(S_n)|.
Everything about which family *dominates* f(n) is heuristic and is reported
as such.

Exit status 0 only if every gate passes.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import sys
import time
from collections import defaultdict
from fractions import Fraction
from math import comb, factorial

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from elemab import elementary_abelian_total, galois_number, gaussian_binomial, rank_counts  # noqa: E402

BLOCKS = (1, 2, 4, 8)
LOG2E = math.log2(math.e)

# A005432, number of subgroups of S_n, n = 0..18 (refs/oeis_A005432.txt)
A005432 = [
    1, 1, 2, 6, 30, 156, 1455, 11300, 151221, 1694723, 29594446, 404126228,
    10594925360, 175238308453, 5651774693595, 117053117995400,
    5320744503742316, 125889331236297288, 7598016157515302757,
]


# ------------------------------------------------------------ polynomials


def w_poly(d: int) -> list[int]:
    """prod_{i<d} (z - 2^i), coefficients from z^0 up."""
    poly = [1]
    for i in range(d):
        a = 2**i
        new = [0] * (len(poly) + 1)
        for k, c in enumerate(poly):
            new[k + 1] += c
            new[k] -= a * c
        poly = new
    return poly


def pmul(p: list[int], q: list[int]) -> list[int]:
    out = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a:
            for j, b in enumerate(q):
                out[i + j] += a * b
    return out


def padd_into(acc: list[int], p: list[int], scale: int) -> None:
    for i, c in enumerate(p):
        acc[i] += scale * c


def L(poly: list[int]) -> int:
    """z^D -> G_D: the number of subspaces of F_2^D."""
    return sum(c * galois_number(D) for D, c in enumerate(poly) if c)


def L_D(poly: list[int]) -> int:
    """z^D -> D G_D = 2 sum_{V <= F_2^D} dim V (by the symmetry [D,k] = [D,D-k])."""
    return sum(D * c * galois_number(D) for D, c in enumerate(poly) if c)


_SQ: dict[int, int] = {}


def L_sq(poly: list[int]) -> int:
    """z^D -> sum_{V <= F_2^D} (dim V)^2."""
    total = 0
    for D, c in enumerate(poly):
        if not c:
            continue
        if D not in _SQ:
            term, acc = 1, 0
            for k in range(1, D + 1):
                term = term * (2 ** (D - k + 1) - 1) // (2**k - 1)
                acc += k * k * term
            _SQ[D] = acc
        total += c * _SQ[D]
    return total


# ------------------------------------------------------------ block types


def load_types(census: dict, keep=None) -> dict[int, list[tuple[int, int, int]]]:
    """{b: [(d, log2|Phi|, count)]} from census8.json, optionally filtered."""
    out = {}
    for b in BLOCKS:
        rows = []
        for cls in census["transitive"][str(b)]["classes"]:
            if keep is None or keep(b, cls):
                rows.append((cls["d"], cls["log2_phi"], cls["class_size"]))
        if rows:
            out[b] = rows
    return out


def weights(types) -> dict[int, list[list[int]]]:
    """[W_b, W_b', W_b''] in z, the u-derivatives at u = 1 of sum_P u^{log2|Phi(P)|} w_{d(P)}(z)."""
    W = {}
    for b, rows in types.items():
        dmax = max(d for d, _, _ in rows)
        ders = [[0] * (dmax + 1) for _ in range(3)]
        for d, ph, cnt in rows:
            for r, factor in enumerate((1, ph, ph * (ph - 1))):
                padd_into(ders[r], w_poly(d), cnt * factor)
        W[b] = ders
    return W


def egf(W, n_max: int) -> list[list[list[int]]]:
    """F[n] = [F_n, F_n', F_n''] with F_n = n! [x^n] exp(sum_b W_b x^b/b!), a polynomial in z.

    Exponential formula F_n = sum_b C(n-1, b-1) W_b F_{n-b}, differentiated
    twice in u by Leibniz.
    """
    F: list[list[list[int]]] = [[[1], [0], [0]]]
    for n in range(1, n_max + 1):
        size = n // 2 + 1
        out = [[0] * size for _ in range(3)]
        for b, wd in W.items():
            if b > n:
                continue
            c = comb(n - 1, b - 1)
            prev = F[n - b]
            for r in range(3):
                for a in range(r + 1):
                    wa, fb = wd[a], prev[r - a]
                    k = c * comb(r, a)
                    for i, x in enumerate(wa):
                        if not x:
                            continue
                        kx = k * x
                        for j, y in enumerate(fb):
                            if y:
                                out[r][i + j] += kx * y
        while len(out[0]) > 1 and not any(o[-1] for o in out):
            for o in out:
                o.pop()
        F.append(out)
    return F


def order_distribution(types, n_max: int) -> list[dict[int, int]]:
    """Exact distribution of log2|H| over the family, for n <= n_max."""
    Wd: dict[int, dict[tuple[int, int], int]] = {}
    for b, rows in types.items():
        acc: dict[tuple[int, int], int] = defaultdict(int)
        for d, ph, cnt in rows:
            for D, c in enumerate(w_poly(d)):
                acc[(D, ph)] += cnt * c
        Wd[b] = {k: v for k, v in acc.items() if v}
    F: list[dict[tuple[int, int], int]] = [{(0, 0): 1}]
    for n in range(1, n_max + 1):
        acc = defaultdict(int)
        for b, wb in Wd.items():
            if b > n:
                continue
            c = comb(n - 1, b - 1)
            for (D1, p1), x in wb.items():
                for (D2, p2), y in F[n - b].items():
                    acc[(D1 + D2, p1 + p2)] += c * x * y
        F.append({k: v for k, v in acc.items() if v})
    out = []
    for n in range(n_max + 1):
        dist: dict[int, int] = defaultdict(int)
        for (D, ph), c in F[n].items():
            for k in range(D + 1):
                dist[k + ph] += c * gaussian_binomial(D, k)
        out.append({t: v for t, v in sorted(dist.items()) if v})
    return out


def dense_structures(c: dict[int, int], n_max: int) -> list[int]:
    """n! [x^n] exp(sum_b c_b x^b / b!): labelled set partitions with c_b colours per b-block."""
    h = [1] + [0] * n_max
    for n in range(1, n_max + 1):
        h[n] = sum(comb(n - 1, b - 1) * cb * h[n - b] for b, cb in c.items() if b <= n)
    return h


def log2i(x: int) -> float:
    return math.log2(x)


def log2_ratio(a: int, b: int) -> float:
    """log2(a/b) for positive integers with a/b near 1, to full double precision.

    log2i(a) - log2i(b) loses everything below 1e-16 of log2 a, which is
    1e-12 at the sizes here, so it reports rounding noise as a value.
    """
    return math.log1p(float(Fraction(a - b, b))) / math.log(2)


# ------------------------------------------------------------------- driver


def moments(Fn: list[list[int]]) -> tuple[int, int, int]:
    """(count, 2 * sum log2|H|, sum (log2|H|)^2) over the family at one n, exact.

    log2|H| = dim V + phi with phi = sum_i log2|Phi(P_i)|, so
    sum X^2 = sum dimV^2 + 2 sum phi dimV + sum phi^2.
    """
    F0, F1, F2 = Fn
    count = L(F0)
    s1_twice = 2 * L(F1) + L_D(F0)
    s2 = L_sq(F0) + L_D(F1) + L(F2) + L(F1)
    return count, s1_twice, s2


def mean_sd(count: int, s1_twice: int, s2: int) -> tuple[float, float]:
    mean = s1_twice / (2 * count)
    var = (4 * count * s2 - s1_twice * s1_twice) / (4 * count * count)
    return mean, math.sqrt(var) if var >= 0 else math.nan


def hist_moments(hist: dict[int, int]) -> tuple[int, int, int]:
    """Same triple from an exact histogram {log2|H|: count}."""
    return (
        sum(hist.values()),
        2 * sum(t * v for t, v in hist.items()),
        sum(t * t * v for t, v in hist.items()),
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--census", type=pathlib.Path, default=ROOT / "data" / "census8.json")
    ap.add_argument("--n-max", type=int, default=512)
    ap.add_argument("--dist-max", type=int, default=64)
    ap.add_argument("--out", type=pathlib.Path, default=ROOT / "data" / "block_families.json")
    args = ap.parse_args()
    n_max = args.n_max
    t0 = time.time()

    failures: list[str] = []

    def gate(ok: bool, label: str) -> None:
        print(f"  [{'pass' if ok else 'FAIL'}] {label}")
        if not ok:
            failures.append(label)

    census = json.loads(args.census.read_text())
    types_all = load_types(census)
    types_elab = load_types(census, lambda b, cls: cls["exponent"] <= 2)
    types_dense = load_types(census, lambda b, cls: b >= 2 and 2 * cls["d"] == b)
    c_dense = {b: sum(cnt for _, _, cnt in rows) for b, rows in types_dense.items()}
    c8 = c_dense[8]
    phibar8 = sum(ph * cnt for _, ph, cnt in types_dense[8]) / c8
    print(f"density-1/2 block types (d = b/2), conjugates per b-set: {c_dense}; their mean log2|Phi| at b=8: {phibar8}")

    W_all = weights(types_all)
    F = egf(W_all, n_max)
    print(f"EGF of the Frattini-closed family, n <= {n_max}: {time.time() - t0:.1f}s")

    print("\n== gate 0: block weights and subspace identities ==")
    w4 = w_poly(4)
    gate(w4 == [64, -120, 70, -15, 1], "w_4(z) = (z-1)(z-2)(z-4)(z-8) = z^4 - 15z^3 + 70z^2 - 120z + 64")
    gate(L(w4) == 1, "L(w_4) = 1: the only subspace of F_2^4 onto F_2^4 is F_2^4")
    gate(L(pmul(w4, w4)) == 2**16, "L(w_4^2) = 2^16: subspaces of F_2^4 x F_2^4 onto both factors = 4x4 matrices (Goursat)")
    gate(L(pmul(w_poly(2), w_poly(2))) == 16, "L(w_2^2) = 16 = 6 + 9 + 1 (Goursat on F_2^2 x F_2^2)")

    print("\n== gate 1: formula vs brute force on permutations (census8.json), n <= 8 ==")
    dist = order_distribution(types_all, max(8, args.dist_max))
    for n in range(1, 9):
        bf = census["frattini_closed"][str(n)]
        bf_hist = {int(k).bit_length() - 1: v for k, v in bf["order_histogram"].items()}
        mom = moments(F[n])
        gate(
            dist[n] == bf_hist and mom == hist_moments(bf_hist) and mom[1] == 2 * bf["sum_log2_order"],
            f"n={n}: {mom[0]} subgroups; order histogram and first two moments of log2|H| equal the brute force",
        )

    print(f"\n== gate 2: moments by u-derivatives vs the exact order distribution, n <= {args.dist_max} ==")
    bad = [n for n in range(args.dist_max + 1) if moments(F[n]) != hist_moments(dist[n])]
    gate(not bad, f"count, sum and sum of squares of log2|H| agree exactly for all n <= {args.dist_max}" + (f" (first failure n={bad[0]})" if bad else ""))

    print("\n== gate 3: elementary abelian block types reproduce elemab.py ==")
    Fe = egf(weights(types_elab), 20)
    for n in range(0, 21):
        got, want = L(Fe[n][0]), elementary_abelian_total(n)
        if n <= 15:
            gate(got == want, f"n={n}: {got} = e(n) (every orbit of an elementary abelian subgroup of S_{n} has size <= 8)")
        else:
            gate(0 < want - got, f"n={n}: e(n) exceeds the orbits-<=8 count by {want - got} (regular (Z/2)^4 on 16 points)")

    print("\n== gate 4: explicit disjoint families never exceed f(n) = A005432(n) ==")
    s3_count = [0] * (n_max + 1)
    for n in range(3, n_max + 1):
        s3_count[n] = comb(n, 3) * L(pmul([0, 1], F[n - 3][0]))
    for n in range(0, 19):
        b_n = L(F[n][0])
        tot = b_n + s3_count[n]
        gate(tot <= A005432[n], f"n={n:2d}: |B(n)| + |S3(n)| = {tot} <= f(n)   (share of f: B {b_n / A005432[n]:.3f}, S3 {s3_count[n] / A005432[n]:.3f})")

    # ------------------------------------------------ e(n): the saddle point
    print("\n== Finding A: e(n) against G_{floor(n/2)} M(n), M = set partitions into blocks of size 2 and 4 ==")
    stats = json.loads((ROOT / "data" / "elemab_stats.json").read_text())["rows"]
    log2_e = {r["n"]: r["log2_e"] for r in stats}
    M_FAR = 32768
    M = dense_structures({2: 1, 4: 1}, max(n_max, M_FAR + 8))
    mism = [n for n in range(2, 41) if n in log2_e and abs(log2i(elementary_abelian_total(n)) - log2_e[n]) > 1e-9]
    gate(not mism, "data/elemab_stats.json log2 e(n) equals elemab.py to 1e-9 for n = 2..40")

    def model_e(n: int) -> float:
        if n % 2 == 0:
            return log2i(galois_number(n // 2)) + log2i(M[n])
        return log2i(galois_number((n - 1) // 2)) + log2i(n * M[n - 1])

    resid = {n: log2_e[n] - model_e(n) for n in sorted(log2_e) if n >= 8 and n <= n_max}
    for n in (8, 9, 16, 17, 32, 33, 64, 128, 256, 512):
        if n in resid:
            print(f"  n={n:4d}  log2 e(n) - log2 G - log2 M = {resid[n]: .3e}")
    big = {n: r for n, r in resid.items() if n >= 128}
    if big:
        gate(max(abs(r) for r in big.values()) < 1e-4, f"|log2 e(n) - log2 G - log2 M| < 1e-4 for all {len(big)} stored n >= 128 (max {max(abs(r) for r in big.values()):.1e})")

    om = json.loads((ROOT / "data" / "orbit_model.json").read_text())
    pred_excess = {n: model_e(n) - n * n / 16 for n in log2_e if n <= n_max}

    def a_local(ex, n, h=8):
        f = lambda t: t * math.log2(t)  # noqa: E731
        return (ex[n + h] - 2 * ex[n] + ex[n - h]) / (f(n + h) - 2 * f(n) + f(n - h))

    diffs = [(r["n"], r["A_local"], a_local(pred_excess, r["n"])) for r in om["local_coefficients"]["8"] if r["n"] >= 64 and r["n"] + 8 in pred_excess and r["n"] - 8 in pred_excess]
    worst = None
    if diffs:
        tail = [abs(a - b) for n, a, b in diffs if n >= 128]
        if tail:
            worst = max(tail)
            gate(worst < 1e-4, f"A_local(n) of data/orbit_model.json equals the G*M prediction to {worst:.1e} for n = 128..{diffs[-1][0]}")
        for n, a, b in diffs:
            if n in (64, 128, 256, 504):
                two_term = 0.75 - math.sqrt(6) / 8 / math.sqrt(n)
                print(f"     n={n:4d}  stored A_local {a:.6f}   from G*M {b:.6f}   3/4 - (sqrt6/8)/sqrt(n) = {two_term:.6f}   (stored - that) * n^1.5 = {(a - two_term) * n**1.5:.3f}")

    # Beyond the stored range, from exact M(n) alone: the second difference of log2 G_{n/2} - n^2/16 is O(2^(-n/4)).
    # Ratios of exact integers keep the float error near 1e-16 where naive differences of log2 M would lose ~1e-3 at n = 8192.
    def a_local_M(n: int, h: int = 8) -> float:
        num = math.log2(M[n + h] * M[n - h] / (M[n] * M[n]))
        den = ((n + h) * math.log1p(h / n) + (n - h) * math.log1p(-h / n)) / math.log(2)
        return num / den

    c32 = 9 * math.sqrt(6) / 32
    far = [2**j for j in range(9, M_FAR.bit_length())]
    rem = {n: (a_local_M(n) - 0.75 + math.sqrt(6) / 8 / math.sqrt(n)) * n**1.5 for n in far}
    print("  from exact M(n): (A_local - 3/4 + (sqrt6/8) n^-1/2) n^1.5 = " + "  ".join(f"{n}: {v:.5f}" for n, v in rem.items()))
    rich = (math.sqrt(2) * rem[far[-1]] - rem[far[-2]]) / (math.sqrt(2) - 1)
    gate(abs(rich - c32) < 1e-3, f"its Richardson limit {rich:.5f} equals the saddle-point coefficient 9 sqrt6/32 = {c32:.5f}: A_local = 3/4 - (sqrt6/8) n^-1/2 + (9 sqrt6/32) n^-3/2 + O(n^-2), no n^-1 term")

    print("\n== Finding B: the rank of a random elementary abelian subgroup is the dimension of a random subspace of F_2^{floor(n/2)} ==")
    tv_rows = {}
    for n in (16, 32, 64, 128, 256):
        if n > n_max:
            continue
        a = rank_counts(n)
        e_n, m = sum(a), n // 2
        g_m = galois_number(m)
        tv = sum(abs(a[k] * g_m - gaussian_binomial(m, k) * e_n) for k in range(m + 1)) / (2 * e_n * g_m)
        tv_rows[n] = tv
        print(f"  n={n:4d}  TV(a(n,.)/e(n), [m,.]_2/G_m) = {tv:.3e}   log2(TV 2^(n/4) / n^2) = {math.log2(tv) + n / 4 - 2 * math.log2(n):6.2f}")
    tvn = sorted(tv_rows)
    gate(all(tv_rows[x] > tv_rows[y] for x, y in zip(tvn, tvn[1:])), "the total variation distance decreases strictly along n = " + ", ".join(map(str, tvn)))
    band = {n: math.log2(tv_rows[n]) + n / 4 - 2 * math.log2(n) for n in tvn if n >= 64}
    gate(all(-12 <= v <= -5 for v in band.values()), "TV = n^2 2^(-n/4) 2^c with c in [-12, -5] for n = 64..256: the rate of the 8-block correction s0^8 2^(-n/4)/1344")
    sdp = json.loads((ROOT / "data" / "sd_precision.json").read_text())["rows"]
    print("  data/sd_precision.json, |sd(n) - theta-law sd| at the same scale: " + "  ".join(f"n={r['n']}: {r['log10_error'] * math.log2(10) + r['n'] / 4 - 2 * math.log2(r['n']):.2f}" for r in sdp))

    # ----------------------------------------------- the families, n <= n_max
    print("\n== second-order constant alpha(X) = (log2 X(n) - n^2/16) / (n log2 n) ==")
    pure8 = {}
    Lw4 = {}
    wpow = [1]
    for m in range(1, n_max // 8 + 1):
        wpow = pmul(wpow, w4)
        n = 8 * m
        Lw4[m] = L(wpow)
        pure8[n] = factorial(n) // (factorial(m) * factorial(8) ** m) * c8**m * Lw4[m]
    rows = []
    for n in (16, 32, 64, 128, 256, 512):
        if n > n_max:
            continue
        lb = log2i(L(F[n][0]))
        den = n * math.log2(n)
        row = {
            "n": n,
            "alpha_e": (log2_e[n] - n * n / 16) / den,
            "alpha_B": (lb - n * n / 16) / den,
            "alpha_pure8": (log2i(pure8[n]) - n * n / 16) / den,
            "log2_B_over_e": lb - log2_e[n],
            "log2_pure8_over_e": log2i(pure8[n]) - log2_e[n],
        }
        rows.append(row)
        print(
            f"  n={n:4d}  e(n): {row['alpha_e']:.4f}   B(n): {row['alpha_B']:.4f}   pure 8-blocks: {row['alpha_pure8']:.4f}"
            f"   log2(B/e) = {row['log2_B_over_e']:8.1f}   log2(pure8/e) = {row['log2_pure8_over_e']:8.1f}"
        )
    cross = next((n for n in sorted(pure8) if n in log2_e and log2i(pure8[n]) > log2_e[n]), None)
    print(f"  smallest stored n where the pure 8-block family alone exceeds e(n): {cross}")

    print("\n== pure 8-block family: (log2 F8(n) - n^2/16 - (7/8) n log2 n) / n ==")
    lin8 = {n: (log2i(pure8[n]) - n * n / 16 - 7 / 8 * n * math.log2(n)) / n for n in (64, 128, 256, 512) if n in pure8}
    for n, v in lin8.items():
        print(f"  n={n:4d}  {v:.5f}")
    pred_lin8 = -7 / 8 * LOG2E + 3 / 8 - math.log2(factorial(8) / c8) / 8
    print(f"  Stirling: -(7/8) log2 e + 3/8 - (1/8) log2(8!/{c8}) = {pred_lin8:.5f}")
    if 512 in lin8:
        gate(abs(lin8[512] - pred_lin8) < 0.02, f"n=512 value is within 0.02 of the Stirling limit (O(1)/n remainder): {lin8[512] - pred_lin8:+.4f}")
    # Stirling to O(1/n) plus log2 G_{4m} - 4m^2 -> log2(theta/eta); what is left should be the next Stirling term.
    theta = sum(2.0 ** (-j * j) for j in range(-40, 41))
    eta = math.prod(1.0 - 2.0**-i for i in range(1, 200))
    stirling_rows = {}
    for n in (64, 128, 256, 512):
        if n in pure8:
            m = n // 8
            full = n * n / 16 + 7 / 8 * n * math.log2(n) + pred_lin8 * n + 1.5 + math.log2(theta / eta) - 7 * LOG2E / (12 * n)
            stirling_rows[n] = {"remainder": log2i(pure8[n]) - full, "next_stirling_term": LOG2E * (1 - 8.0**-3) / (360 * m**3), "log2_L_over_G": log2_ratio(Lw4[m], galois_number(4 * m))}
            print(f"  n={n:4d}  log2 F8(n) - [n^2/16 + (7/8) n log2 n + cn + 3/2 + log2(theta/eta) - 7 log2(e)/(12n)] = {stirling_rows[n]['remainder']: .3e}"
                  f"   next Stirling term log2(e)(1 - 8^-3)/(360 m^3) = {stirling_rows[n]['next_stirling_term']:.3e}   log2(L(w_4^m)/G_4m) = {stirling_rows[n]['log2_L_over_G']: .3e}")

    print("\n== B(n) against G_{n/2} Bstar(n), Bstar = set partitions into blocks of density 1/2 ==")
    Bstar = dense_structures(c_dense, n_max)
    ratio_rows = {}
    for n in (16, 32, 64, 128, 256, 512):
        if n <= n_max:
            ratio_rows[n] = log2_ratio(L(F[n][0]), galois_number(n // 2) * Bstar[n])
            print(f"  n={n:4d}  log2 B(n) - log2 G_(n/2) - log2 Bstar(n) = {ratio_rows[n]: .3e}  (from the exact integers)")
    tail_r = [abs(v) for n, v in ratio_rows.items() if n >= 128]
    if tail_r:
        gate(max(tail_r) < 1e-3, f"|log2 B - log2 G - log2 Bstar| < 1e-3 for n >= 128 (max {max(tail_r):.1e})")

    print("\n== order statistics inside B(n) ==")
    # If every orbit type were dense and dim V averaged D/2 = n/4, E log2|H| / n would be
    # 1/4 + sum_b f_b phibar_b / b, with phibar_b the mean log2|Phi| of the dense b-types.
    phibar = {b: Fraction(sum(ph * cnt for _, ph, cnt in rows), sum(cnt for _, _, cnt in rows)) for b, rows in types_dense.items()}
    print("  mean log2|Phi| of the dense types: " + ", ".join(f"b={b}: {v}" for b, v in phibar.items()))
    blocks_rows = {}
    for n in (8, 16, 32, 64, 128, 256, 512):
        if n > n_max:
            continue
        cnt, s1t, s2 = moments(F[n])
        mu, sd = mean_sd(cnt, s1t, s2)
        fx = {b: Fraction(b * comb(n, b) * L(pmul(W_all[b][0], F[n - b][0])), cnt * n) for b in W_all if b <= n}
        frac = {b: float(v) for b, v in fx.items()}
        gap = Fraction(s1t, 2 * cnt * n) - Fraction(1, 4) - sum(phibar[b] / b * fx.get(b, 0) for b in phibar)
        blocks_rows[n] = {"mean_log2_order": mu, "sd_log2_order": sd, "point_fraction": frac, "dense_identity_gap": float(gap)}
        print(
            f"  n={n:4d}  E[log2|H|]/n = {mu / n:.4f}  sd(log2|H|) = {sd:6.3f}  sd/n^(1/4) = {sd / n**0.25:.3f}   points in blocks of size "
            + "  ".join(f"{b}: {frac[b]:.3f}" for b in sorted(frac))
            + f"   E[log2|H|]/n - 1/4 - sum_b f_b phibar_b/b = {float(gap):.2e} (exact)"
        )
    print(f"  all points in 8-blocks: 1/4 + (1/8)(mean log2|Phi| = {phibar8}) = {0.25 + phibar8 / 8:.4f}")
    gate(all(v["sd_log2_order"] >= 0 for v in blocks_rows.values()), "every variance of log2|H| is nonnegative")

    print("\n== odd n: S3(n) = {H : A_3 on some 3-set is normal in H, H/A_3 in the Frattini-closed family} ==")
    s3_rows = {}
    for n in (9, 17, 33, 65, 129, 257, 511):
        if n <= n_max:
            s3_rows[n] = log2i(s3_count[n]) - log2i(L(F[n][0]))
            print(f"  n={n:4d}  log2(|S3(n)| / |B(n)|) = {s3_rows[n]: .3f}")
    s3_even = {n: log2i(s3_count[n]) - log2i(L(F[n][0])) for n in (8, 16, 32, 64, 128, 256, 512) if n <= n_max}
    print("  even n, where S3 needs a fixed point: " + "  ".join(f"n={n}: {v:.3f} ({v / (n / 4):.3f} per n/4)" for n, v in s3_even.items()))

    print("\n== conjectured formula vs A005432, n <= 18 (pre-asymptotic; information only) ==")
    conj = {}
    for n in range(2, 19):
        if n % 2 == 0:
            val = galois_number(n // 2) * Bstar[n]
        else:
            val = galois_number((n - 1) // 2) * (n * Bstar[n - 1] + comb(n, 3) * Bstar[n - 3])
        conj[n] = A005432[n] / val
    print("  f(n) / conjecture: " + "  ".join(f"{n}:{v:.2f}" for n, v in conj.items()))

    out = {
        "description": "Exact Frattini-closed block-family counts (src/block_families.py)",
        "dense_types": {str(b): v for b, v in c_dense.items()},
        "mean_log2_phi_degree8_dense": phibar8,
        "e_minus_GM_log2": {str(n): v for n, v in resid.items()},
        "A_local_max_abs_diff_n_ge_128": worst,
        "A_local_remainder_times_n15_from_M": {str(n): v for n, v in rem.items()},
        "A_local_n_minus_1_5_coefficient": {"richardson": rich, "saddle_point_9sqrt6_over_32": c32},
        "rank_law_tv_vs_random_subspace": {str(n): v for n, v in tv_rows.items()},
        "alpha": rows,
        "pure8_first_exceeds_e_at": cross,
        "pure8_linear_coefficient": {str(n): v for n, v in lin8.items()},
        "pure8_linear_prediction": pred_lin8,
        "pure8_stirling_remainder": {str(n): v for n, v in stirling_rows.items()},
        "B_minus_GBstar_log2": {str(n): v for n, v in ratio_rows.items()},
        "order_stats": {str(n): {"mean_log2_order": v["mean_log2_order"], "sd_log2_order": v["sd_log2_order"], "point_fraction": {str(b): f for b, f in v["point_fraction"].items()}, "dense_identity_gap": v["dense_identity_gap"]} for n, v in blocks_rows.items()},
        "odd_n_log2_S3_over_B": {str(n): v for n, v in s3_rows.items()},
        "even_n_log2_S3_over_B": {str(n): v for n, v in s3_even.items()},
        "f_over_conjecture_small_n": {str(n): v for n, v in conj.items()},
        "log2_B": {str(n): log2i(L(F[n][0])) for n in range(1, n_max + 1)},
        "gates_failed": failures,
        "seconds": round(time.time() - t0, 1),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=1) + "\n")
    shown = args.out.resolve().relative_to(ROOT) if args.out.resolve().is_relative_to(ROOT) else args.out
    print(f"\nwrote {shown}  ({time.time() - t0:.1f}s)")
    print(f"{len(failures)} gate(s) failed" if failures else "all gates passed")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
