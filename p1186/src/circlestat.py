"""The exact stationary point behind a rounded circle word: is the continuum
optimum itself an equal-block word?

Phi_k of a circle colouring is piecewise quadratic in its breakpoints: the
polygons that make it up have vertices that are rational-linear in the
breakpoints, and their combinatorial type is constant on cells of a
polyhedral subdivision of breakpoint space.  On one cell Phi_k is therefore
an exact quadratic Q with rational coefficients, and the stationary point of
Q, the solution of a rational linear system, is rational.  If it lies in the
cell it is the local continuum optimum, and since it is rational it is itself
an equal-block word, of length D = the common denominator of its breakpoints.

This script, starting from a stored word (by default the best results/zp_k*_c
word for k):

  1. relaxes the word's breakpoints to the floating-point local minimum with
     the C evaluator (search only);
  2. snaps that point to a rational base point a0 of denominator 2^30 and
     takes the exact gradient g and Hessian H of Phi_k there by central
     differences with step h, using circle.phi_circle in Fractions.  On a
     quadratic piece these differences are exact; they are recomputed with
     step h/2 and must agree, a test that the whole stencil lies on one piece;
  3. solves H delta = -g in Fractions, a* = a0 + delta, and evaluates Phi_k(a*)
     exactly with both circle.phi_circle and circle.phi_sweep.  a* is on the
     same piece only if Phi_k(a*) equals the quadratic's prediction
     Phi(a0) + g.delta + delta.H.delta / 2; that equality is checked, and H is
     checked positive definite (exact Cholesky).

It prints D, Phi_k(a*), and the gap to the certified rounded word.  Whether
a* could be certified the repository's way (exact.phi or exact.phi_bits on the
length-D word, and src/verify_zpx) depends on D, which is measured, not
assumed.  Nothing here is stored in results/.

    python src/circlestat.py 4 [--word-file results/zp_k4_c.json] [--workers 4]
"""

import argparse
import glob
import json
import os
import sys
import time
from fractions import Fraction
from math import lcm
from multiprocessing import Pool

import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import circle      # noqa: E402
import circleopt   # noqa: E402

K = None


def phi_free(a):
    """Exact Phi_k at interior breakpoints a (a_0 = 0 and a_B = 1 implicit)."""
    return circle.phi_circle([Fraction(0)] + list(a) + [Fraction(1)], K)


def _init(k):
    global K
    K = k


def _eval(a):
    return phi_free(a)


def float_min(L, k):
    """Polish arc lengths to a floating-point local minimum (C evaluator)."""
    f = lambda v: circleopt.fast_phi_lengths(v, k)
    best = np.asarray(L, float)
    for _ in range(3):
        r = minimize(f, best, method="L-BFGS-B",
                     options={"eps": 1e-10, "ftol": 1e-16, "gtol": 1e-14, "maxiter": 20000})
        r = minimize(f, np.abs(r.x) / np.abs(r.x).sum(), method="Powell",
                     options={"xtol": 1e-12, "ftol": 1e-16, "maxfev": 200000})
        best = np.abs(r.x) / np.abs(r.x).sum()
    return best, f(best)


def derivatives(pool, a0, h):
    n = len(a0)
    pts = [tuple(a0)]
    E = lambda i, s: tuple(a0[t] + (s * h if t == i else 0) for t in range(n))
    for i in range(n):
        pts += [E(i, 1), E(i, -1)]
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for i, j in pairs:
        for si, sj in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
            pts.append(tuple(a0[t] + (si * h if t == i else 0) + (sj * h if t == j else 0)
                             for t in range(n)))
    vals = dict(zip(pts, pool.map(_eval, pts, chunksize=4)))
    f0 = vals[tuple(a0)]
    g = [(vals[E(i, 1)] - vals[E(i, -1)]) / (2 * h) for i in range(n)]
    H = [[Fraction(0)] * n for _ in range(n)]
    for i in range(n):
        H[i][i] = (vals[E(i, 1)] - 2 * f0 + vals[E(i, -1)]) / (h * h)
    for i, j in pairs:
        P = lambda si, sj: vals[tuple(a0[t] + (si * h if t == i else 0) + (sj * h if t == j else 0)
                                      for t in range(n))]
        H[i][j] = H[j][i] = (P(1, 1) - P(1, -1) - P(-1, 1) + P(-1, -1)) / (4 * h * h)
    return f0, g, H


def solve(H, b):
    n = len(b)
    A = [row[:] + [b[i]] for i, row in enumerate(H)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c] / A[c][c]
                A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [A[i][n] / A[i][i] for i in range(n)]


def positive_definite(H):
    n = len(H)
    A = [row[:] for row in H]
    for c in range(n):
        if A[c][c] <= 0:
            return False
        for r in range(c + 1, n):
            f = A[r][c] / A[c][c]
            for t in range(c, n):
                A[r][t] -= f * A[c][t]
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("k", type=int)
    ap.add_argument("--word-file", default=None)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--hexp", type=int, default=34, help="finite-difference step 2^-hexp")
    args = ap.parse_args()
    k = args.k
    t0 = time.time()
    path = args.word_file or os.path.join(ROOT, "results", f"zp_k{k}_c.json")
    recs = sorted((Fraction(r["value"]), key, r["word"]) for key, r in json.load(open(path)).items())
    val_word, key, word = recs[0]
    print(f"start: {os.path.relpath(path, ROOT)} m={key}, Phi = {val_word} = {float(val_word):.12f}",
          flush=True)
    bounds, _ = circle.bounds_of_cyclic_word(word)
    L, vf = float_min(np.diff([float(b) for b in bounds]), k)
    L = circleopt.clean(L)
    vf = circleopt.fast_phi_lengths(L, k)
    print(f"floating-point local minimum: Phi = {vf:.15f}, {L.size} arcs", flush=True)
    a_float = np.cumsum(L)[:-1]
    Q = 2 ** 30
    a0 = [Fraction(int(round(x * Q)), Q) for x in a_float]
    n = len(a0)
    print(f"{n} free breakpoints; exact derivatives need {2 * n * n + 1} exact evaluations per step",
          flush=True)
    with Pool(args.workers, initializer=_init, initargs=(k,)) as pool:
        h = Fraction(1, 2 ** args.hexp)
        f0, g, H = derivatives(pool, a0, h)
        print(f"  step 2^-{args.hexp}: done ({time.time() - t0:.0f}s)", flush=True)
        f0b, g2, H2 = derivatives(pool, a0, h / 2)
        same = g == g2 and H == H2
        print(f"  step 2^-{args.hexp + 1}: gradient and Hessian {'IDENTICAL' if same else 'DIFFER'}"
              f" ({time.time() - t0:.0f}s)", flush=True)
    _init(k)
    pd = positive_definite(H)
    delta = solve(H, [-x for x in g])
    a_star = [x + d for x, d in zip(a0, delta)]
    pred = f0 + sum(x * y for x, y in zip(g, delta)) / 1 + \
        sum(delta[i] * H[i][j] * delta[j] for i in range(n) for j in range(n)) / 2
    increasing = all(a_star[i] < a_star[i + 1] for i in range(n - 1)) and 0 < a_star[0] and a_star[-1] < 1
    exact_circle = phi_free(a_star) if increasing else None
    exact_sweep = circle.phi_sweep([Fraction(0)] + a_star + [Fraction(1)], k) if increasing else None
    D = 1
    for x in a_star:
        D = lcm(D, x.denominator)
    print(f"Hessian positive definite: {pd}; breakpoints increasing: {increasing}")
    print(f"common denominator D of a*: {D} ({len(str(D))} digits)")
    print(f"Phi(a*) by circle.phi_circle = {float(exact_circle) if exact_circle is not None else None!r}")
    print(f"Phi(a*) by circle.phi_sweep  equal: {exact_circle == exact_sweep}")
    print(f"quadratic prediction equals exact value: {pred == exact_circle}")
    if exact_circle is not None:
        print(f"Phi(a*) = {exact_circle}" if len(str(exact_circle)) < 200 else
              f"Phi(a*) has a {len(str(exact_circle.denominator))}-digit denominator")
        print(f"Phi(a*) = {float(exact_circle):.15f}; certified word {float(val_word):.15f};"
              f" gap {float(val_word - exact_circle):.3e}")
    print(f"total {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
