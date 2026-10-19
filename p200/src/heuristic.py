#!/usr/bin/env python3
"""heuristic.py -- Hardy-Littlewood prediction for the number of k-term
arithmetic progressions of primes with last term <= N.

Model.  For a fixed common difference d the k-tuple H_d = {0, d, ..., (k-1)d}
has singular series

    S(d) = prod_p  (1 - w_d(p)/p) / (1 - 1/p)^k,
    w_d(p) = #{ i*d mod p : 0 <= i < k } = 1 if p | d, else min(k, p),

and the prime k-tuple conjecture predicts

    #{ a : a+id prime for all i }  ~  S(d) * INT_{t0}^{X} prod_i dt/log(t+id).

S(d) vanishes unless every prime p <= k divides d, so write d = W_k * m with
W_k = prod_{p<=k} p.  Factoring S(W m) into a part independent of m and a
multiplicative correction,

    S(W m) = A_k * T_k * g(m),
    A_k = prod_{p<=k} (1-1/p)^(1-k),
    T_k = prod_{p>k}  (1-k/p) (1-1/p)^(-k),
    g(m) = prod_{p | m, p > k} (1-1/p) / (1-k/p).

Summing over m gives the predicted count.  Note E[g(m)] = 1 over m, because for
each p > k the factor averages to (1/p)(1-1/p)^(1-k) + (1-1/p)(1-k/p)(1-1/p)^(-k)
divided by the p-factor of T_k; A_k*T_k is therefore the mean singular series.

This file also contains the "expected count = 1" threshold estimate for a(k).
"""
import math

from sympy import primerange, factorint

GAMMA = 0.5772156649015328606


def primorial(k):
    r = 1
    for p in primerange(2, k + 1):
        r *= p
    return r


def log_A(k):
    return sum((1 - k) * math.log1p(-1.0 / p) for p in primerange(2, k + 1))


def log_T(k, plim=2 * 10 ** 6):
    """log prod_{k < p <= plim} (1-k/p)(1-1/p)^(-k), plus an estimate of the tail.

    The tail satisfies log factor = -(k^2-k)/(2p^2) + O(k^3/p^3), so it is
    approximated by -(k^2-k)/2 * sum_{p>plim} 1/p^2 ~ -(k^2-k)/(2 plim log plim).
    """
    s = 0.0
    for p in primerange(k + 1, plim + 1):
        s += math.log1p(-k / p) - k * math.log1p(-1.0 / p)
    tail = -(k * k - k) / (2.0 * plim * math.log(plim))
    return s + tail, tail


def g_of_m(m, k):
    r = 1.0
    for p in factorint(m):
        if p > k:
            r *= (1.0 - 1.0 / p) / (1.0 - k / p)
    return r


_LEG = {}


def _leg(n):
    if n not in _LEG:
        import numpy as np
        _LEG[n] = np.polynomial.legendre.leggauss(n)
    return _LEG[n]


def integrate_density(lo, hi, d, k, n=64):
    """Gauss-Legendre quadrature of prod_i dt/log(t+id) over [lo, hi]."""
    if hi <= lo:
        return 0.0
    import numpy as np
    x, w = _leg(n)
    mid, half = 0.5 * (lo + hi), 0.5 * (hi - lo)
    t = mid + half * x
    s = np.zeros_like(t)
    for i in range(k):
        s += np.log(np.log(t + i * d))
    return float(half * np.sum(w * np.exp(-s)))


def log_Eg(k, plim=2 * 10 ** 6):
    """log of the mean of g(m) over m:  prod_{p>k} (1-1/p)(1-(k-1)/p)/(1-k/p).

    The p-factor expands as 1 + (k-1)/p^2 + O(k^2/p^3), so the tail beyond plim
    is approximately (k-1) * sum_{p>plim} 1/p^2 ~ (k-1)/(plim log plim).
    """
    s = 0.0
    for p in primerange(k + 1, plim + 1):
        s += math.log1p(-1.0 / p) + math.log1p(-(k - 1) / p) - math.log1p(-k / p)
    return s + (k - 1) / (plim * math.log(plim))


_gsieve = {}


def g_table(k, M):
    """Array g[1..M] of the multiplicative weight, via a least-prime-factor sieve."""
    import numpy as np
    key = (k, M)
    if key in _gsieve:
        return _gsieve[key]
    g = np.ones(M + 1)
    comp = np.zeros(M + 1, dtype=bool)
    for p in range(2, M + 1):
        if comp[p]:
            continue
        comp[p * p::p] = True
        if p > k:
            g[p::p] *= (1.0 - 1.0 / p) / (1.0 - k / p)
    _gsieve.clear()
    _gsieve[key] = g
    return g


def _smooth_integral(k, N, t0, nd=128, nt=128):
    """(1/W) * INT_0^D INT_{t0}^{N-(k-1)d} prod_i dt/log(t+id) dd,  D=(N-t0)/(k-1).

    Used in place of the exact sum over d = W*m when the number of terms is large;
    the two are cross-checked against each other in the overlap region.
    """
    import numpy as np
    D = (N - t0) / (k - 1)
    if D <= 0:
        return 0.0
    xd, wd = _leg(nd)
    xt, wt = _leg(nt)
    dmid, dhalf = 0.5 * D, 0.5 * D
    d = dmid + dhalf * xd                      # (nd,)
    hi = N - (k - 1) * d                       # (nd,)
    tmid = 0.5 * (t0 + hi)
    thalf = 0.5 * (hi - t0)
    t = tmid[:, None] + thalf[:, None] * xt[None, :]   # (nd, nt)
    s = np.zeros_like(t)
    for i in range(k):
        s += np.log(np.log(t + i * d[:, None]))
    inner = thalf * np.sum(wt[None, :] * np.exp(-s), axis=1)
    return float(dhalf * np.sum(wd * inner))


_cache = {}


def constants(k):
    """(log A_k, log T_k, tail estimate in log T_k)."""
    if k not in _cache:
        lT, tail = log_T(k)
        _cache[k] = (log_A(k), lT, tail)
    return _cache[k]


MEXACT = 20000


def predicted_count(k, N, t0=None, force=None):
    """Predicted number of k-APs of primes with last term <= N and W_k | d."""
    W = primorial(k)
    lA, lT, _ = constants(k)
    if t0 is None:
        t0 = max(2.0, float(k))
    M = int((N - t0) // ((k - 1) * W))
    if M < 1:
        return 0.0
    mode = force or ("exact" if M <= MEXACT else "smooth")
    if mode == "exact":
        import numpy as np
        g = g_table(k, M)[1:M + 1]
        d = W * np.arange(1, M + 1, dtype=float)
        hi = N - (k - 1) * d
        x, w = _leg(64)
        mid = 0.5 * (t0 + hi)
        half = 0.5 * (hi - t0)
        t = mid[:, None] + half[:, None] * x[None, :]
        lg = np.zeros_like(t)
        for i in range(k):
            lg += np.log(np.log(t + i * d[:, None]))
        inner = half * np.sum(w[None, :] * np.exp(-lg), axis=1)
        return math.exp(lA + lT) * float(np.sum(g * inner))
    return math.exp(lA + lT + log_Eg(k)) * _smooth_integral(k, N, t0) / W


def threshold_ak(k, lo=None, hi=None):
    """Smallest N with predicted_count(k, N) >= 1 -- the heuristic estimate of a(k)."""
    lo = lo or float((k - 1) * primorial(k) + k)
    hi = hi or lo * 10.0
    while predicted_count(k, hi) < 1.0:
        lo, hi = hi, hi * 8.0
        if hi > 1e300:
            return float("inf")
    for _ in range(60):
        if hi / lo < 1.001:
            break
        mid = math.sqrt(lo * hi)
        if predicted_count(k, mid) >= 1.0:
            hi = mid
        else:
            lo = mid
    return math.sqrt(lo * hi)


def granville_ak(k):
    """Granville 2008, formula (2.1): last term of the smallest k-AP ~ (e^(1-gamma) k/2)^(k/2)."""
    return (math.exp(1 - GAMMA) * k / 2.0) ** (k / 2.0)


if __name__ == "__main__":
    # direct check of the singular series against a brute-force count
    import numpy as np
    X = 10 ** 7
    sieve = np.ones(X + 200, dtype=bool)
    sieve[:2] = False
    for i in range(2, int((X + 200) ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    print("singular-series spot checks: exact count of a<=X with a+i*d prime, i<k")
    print(f"{'k':>2} {'d':>8} {'X':>10} {'actual':>10} {'predicted':>12} {'ratio':>7}")
    for k, d in [(3, 6), (4, 30), (5, 30), (6, 30), (7, 210), (5, 60), (6, 90), (8, 210)]:
        W = primorial(k)
        assert d % W == 0, (k, d, W)
        L = X - (k - 1) * d
        m = sieve[0:L].copy()
        for i in range(1, k):
            m &= sieve[i * d:i * d + L]
        actual = int(m.sum())
        lA, lT, _ = constants(k)
        pred = math.exp(lA + lT) * g_of_m(d // W, k) * integrate_density(2.0, L, d, k)
        print(f"{k:>2} {d:>8} {X:>10} {actual:>10} {pred:>12.1f} {pred/max(actual,1):>7.3f}")
