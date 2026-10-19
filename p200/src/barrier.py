"""The Hardy-Littlewood count C(k,N) in closed form, and B(k,N) = 2^k k! C(k,N).

Correction, September 22, 2026. This module was titled "the sieve barrier" and argued that B(k,N) >= 1 everywhere shows no sieve upper bound yields L(N) <= c log N. The inference fails in both directions. The constant 2^k k! + o(1) is a theorem for fixed k as x -> infinity, with an o(1) that depends on k, so B < 1 at k = c log N would have proved nothing. And B >= 1 shows nothing either, because B overstates what a sieve gives: the Selberg main term is at most the number of candidates T0, and B > T0 at 69 of the 82 grid points with log10 N <= 10^6 and at all 19 records. The parenthetical that closed this docstring, calling B > 1 an a fortiori demonstration because granting the bound for every d only inflates B, had the direction backwards: inflating an upper bound cannot show that the true bound is large. The conclusion itself is true, and src/sieve_limits.py proves it for the actual Selberg main term; NOTES.md section 6 is the write-up. Every number below is computed as before, and results/barrier.json is unchanged.

Two things are computed and checked here.

1.  An asymptotic closed form for the Hardy-Littlewood count

        C(k,N) ~ A_k T_k E[g] * N^2 / (2 (k-1) W_k (log N)^k),
        A_k = prod_{p<=k} (1-1/p)^{1-k},

    validated against the exact heuristic.predicted_count in the range where the
    latter is computable.  The shape N^2/(2(k-1) log^k N) times a singular series
    is the same normalisation Tao-Teravainen state in Example 1.7 of
    arXiv:2107.02158, which is an independent check on the constant.

2.  B(k,N): the classical k-dimensional (Selberg / Halberstam-Richert) upper bound for a fixed admissible k-tuple,

        #{a <= x : a+h_1, ..., a+h_k all prime} <= (2^k k! + o(1)) S(H) x/(log x)^k,

    summed over the admissible d as though it held uniformly in k, gives B(k,N) = 2^k k! * C(k,N). We evaluate log B on a grid and report its minimum; src/sieve_limits.py compares B with T0 at these points.
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import heuristic as H
from sympy import primerange

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GAMMA = H.GAMMA


KEXACT = 2 * 10 ** 6   # above this, prime loops are replaced by asymptotics


def log_theta(k):
    """log W_k = theta(k).  Exact for k <= KEXACT, else PNT with the RS-size error."""
    if k <= KEXACT:
        return sum(math.log(p) for p in primerange(2, int(k) + 1))
    return float(k)   # PNT; measured error 0.07% at k = 2e6 (see validate_large_k)


def log_A(k):
    """log prod_{p<=k} (1-1/p)^{1-k}.  Exact for k <= KEXACT, else Mertens:
    prod_{p<=k} (1-1/p)^{-1} = e^gamma log k (1 + O(1/log^2 k))."""
    if k <= KEXACT:
        return (1 - k) * sum(math.log1p(-1.0 / p) for p in primerange(2, int(k) + 1))
    return (k - 1) * (GAMMA + math.log(math.log(k)))


def log_T_asymp(k):
    """log prod_{p>k} (1-k/p)(1-1/p)^{-k} ~ -(k^2-k)/2 * sum_{p>k} 1/p^2 ~ -k/(2 log k)."""
    if k <= 20000:
        _, lT, _ = H.constants(int(k))
        return lT
    return -(k * k - k) / 2.0 * (1.0 / (k * math.log(k)))


def log_Eg_asymp(k):
    """log E[g] = sum_{p>k} log[(1-1/p)(1-(k-1)/p)/(1-k/p)] ~ (k-1) sum_{p>k} 1/p^2."""
    if k <= 20000:
        return H.log_Eg(int(k))
    return (k - 1) / (k * math.log(k))


def log_C_asymp(k, L):
    """log of the asymptotic Hardy-Littlewood count, L = log N."""
    return (log_A(k) + log_T_asymp(k) + log_Eg_asymp(k) + 2 * L
            - math.log(2 * (k - 1)) - log_theta(k) - k * math.log(L))


def log_sieve_bound(k, L):
    """log of 2^k k! times the asymptotic count."""
    return k * math.log(2.0) + math.lgamma(k + 1) + log_C_asymp(k, L)


def validate_asymptotic():
    print("1.  asymptotic closed form vs the exact heuristic (same quantity)")
    print("%3s %10s %16s %16s %10s" % ("k", "log10 N", "exact log C", "asymp log C", "ratio"))
    rows = []
    for k, e in [(8, 8), (9, 8), (10, 8), (11, 9), (12, 9), (13, 9.5), (14, 9.7),
                 (15, 9.9), (17, 10.3), (19, 11.2), (22, 13.6), (26, 16.1)]:
        N = 10.0 ** e
        ex = H.predicted_count(k, N)
        if ex <= 0:
            continue
        le, la = math.log(ex), log_C_asymp(k, math.log(N))
        rows.append(math.exp(la - le))
        print("%3d %10.2f %16.6f %16.6f %10.4f" % (k, e, le, la, math.exp(la - le)))
    gm = math.exp(sum(math.log(r) for r in rows) / len(rows))
    print("  geometric mean asymp/exact = %.4f  (range %.4f .. %.4f)"
          % (gm, min(rows), max(rows)))
    print("  the asymptotic drops the t<N gain in the density integral and the")
    print("  Sum_m -> integral discretisation, so a few tens of percent is expected.")
    return gm


def validate_large_k():
    """The k>KEXACT formulas must agree with the exact ones where both are defined."""
    print("\n1b. large-k asymptotics vs exact prime sums (overlap region)")
    print("%8s %12s %12s %10s %12s %12s %10s"
          % ("k", "theta exact", "theta asym", "ratio", "logA exact", "logA asym", "ratio"))
    for k in (1000, 5000, 20000, 200000, 2 * 10 ** 6):
        te = sum(math.log(p) for p in primerange(2, k + 1))
        ta = float(k)
        ae = (1 - k) * sum(math.log1p(-1.0 / p) for p in primerange(2, k + 1))
        aa = (k - 1) * (GAMMA + math.log(math.log(k)))
        print("%8d %12.6g %12.6g %10.5f %12.6g %12.6g %10.5f"
              % (k, te, ta, ta / te, ae, aa, aa / ae))
    print("  theta(k)/k -> 1, so theta(k)=k is used above KEXACT; the residual is o(k)")
    print("  while the margin in log B grows like k log log k, so it cannot flip the sign.")
    print("  The asymptotic C understates the exact heuristic by ~2x (section 1), so B is understated by the same factor.")
    print("  B is not a sieve bound at k = c log N (see the module docstring), so its sign proves nothing either way.")


def main():
    validate_asymptotic()
    validate_large_k()

    print("\n2.  B(k,N) = 2^k k! C(k,N),  k = floor(c log N): the fixed-k sieve constant, extrapolated")
    print("    B is not a valid sieve bound at k = c log N, so neither sign in the log10 column proves anything; see src/sieve_limits.py.")
    print("%12s %6s %8s %14s" % ("log10 N", "c", "k", "log10 B(k,N)"))
    worst = None
    grid = []
    for e10 in [2, 3, 4, 6, 8, 12, 16, 20, 30, 50, 100, 300, 1000, 10 ** 4,
                10 ** 6, 10 ** 9, 10 ** 12, 10 ** 30, 10 ** 100]:
        L = e10 * math.log(10)
        for c in (0.2, 0.3, 0.5, 0.7, 0.9, 1.0):
            k = int(c * L)
            if k < 3:
                continue
            lb = log_sieve_bound(k, L) / math.log(10)
            grid.append((e10, c, k, lb))
            if worst is None or lb < worst[3]:
                worst = (e10, c, k, lb)
            if e10 in (4, 16, 100, 10 ** 4, 10 ** 9, 10 ** 30):
                print("%12g %6.1f %8d %14.4g" % (e10, c, k, lb))
    print("\n  minimum of log10 B over the whole grid: %.4g  at log10 N=%g, c=%.1f, k=%d"
          % (worst[3], worst[0], worst[1], worst[2]))
    print("  B >= 1 everywhere: %s" % all(g[3] >= 0 for g in grid))

    print("\n3.  the asymptotic condition for B < 1")
    print("    log B < 0  requires   log log k + log c + gamma + log 2 - 2 + 2/c < 0.")
    print("    min over c of (log c + 2/c) = 1 + log 2 at c = 2, so the requirement is")
    print("    log log k < 2 - gamma - 2 log 2 - 1 = %.4f, i.e. k < %.4f."
          % (2 - GAMMA - 2 * math.log(2) - 1, math.exp(math.exp(2 - GAMMA - 2 * math.log(2) - 1))))
    for c in (0.5, 1.0, 2.0):
        thr = 2 - GAMMA - math.log(2) - math.log(c) - 2 / c
        print("    c=%.1f: need log log k < %+.4f -> k < %.4f"
              % (c, thr, math.exp(math.exp(thr)) if thr < 5 else float("inf")))
    print("    No k >= 3 qualifies at any c, so B >= 1 for every c. That says nothing about the sieve itself, whose main term src/sieve_limits.py bounds directly.")

    print("\n4.  what the overshoot measures")
    print("    B = exp(k[log log k + log c + gamma + log 2 - 2 + 2/c] + o(k)), which grows like (log k)^k.")
    print("    This was once read as the factor by which the constant 2^k k! must be beaten. It is not: the (log k)^k comes from B/T0, about (2 e^(gamma-1) c log k)^k, the factor by which B exceeds the candidate count T0, which no sieve main term can exceed.")
    print("    src/sieve_limits.py measures the real shortfall: a proof needs about 2/c - 1 nats per term, and at any fixed level N^A the Selberg main term saves o(1) per term.")
    print("    An upper bound matching the Hardy-Littlewood count uniformly in k would give L(N) <= (2 + o(1)) log N / log log N and answer the problem.")

    json.dump({"grid": [{"log10N": g[0], "c": g[1], "k": g[2], "log10B": g[3]} for g in grid],
               "min_log10B": worst[3]},
              open(os.path.join(ROOT, "results", "barrier.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
