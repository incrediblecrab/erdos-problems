#!/usr/bin/env python3
"""analysis.py -- L(N) = length of the longest AP of primes in {1,...,N}.

Produces:
  * the exact values of L(N) on its 26 known jump points,
  * the extremal ratio L(N)/log N,
  * an explicit rigorous upper bound on L(N) (see NOTES.md, Prop. 1-2),
  * comparison with Granville's heuristic a(k) ~ (e^(1-gamma) k/2)^(k/2)
    and with the Hardy-Littlewood "expected count = 1" threshold,
  * extrapolation of the conjectured decay rate of L(N)/log N.
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

from sympy import isprime, nextprime, primerange

import heuristic as H

GAMMA = H.GAMMA


def primorial(k):
    r = 1
    for p in primerange(2, k + 1):
        r *= p
    return r


def min_last_lower_bound(k):
    """Rigorous lower bound for the last term of any k-AP of primes (k >= 3).

    Case A (first term > k):  W_k | d, first term >= nextprime(k)
                              => last >= nextprime(k) + (k-1) W_k
    Case B (k prime, first term = k): (W_k/k) | d
                              => last >= k + (k-1) W_k / k
    """
    W = primorial(k)
    caseA = nextprime(k) + (k - 1) * W
    if isprime(k):
        caseB = k + (k - 1) * (W // k)
        return min(caseA, caseB)
    return caseA


_bound_table = [None, None, None]   # _bound_table[k] = min_last_lower_bound(k)


def _extend(kmax):
    """Build min_last_lower_bound(k) for k <= kmax with a running primorial."""
    k = len(_bound_table) - 1
    W = primorial(k) if k >= 2 else 1
    while k < kmax:
        k += 1
        if isprime(k):
            W *= k
            b = min(nextprime(k) + (k - 1) * W, k + (k - 1) * (W // k))
        else:
            b = nextprime(k) + (k - 1) * W
        _bound_table.append(b)


def rigorous_L_bound(N):
    """Largest k for which a k-AP of primes with last term <= N is not excluded.

    Monotone because min_last_lower_bound(k) grows at least like (k-1)W_k/k.
    """
    k = 3
    while True:
        _extend(k + 64)
        while k + 1 < len(_bound_table) and _bound_table[k + 1] <= N:
            k += 1
        if k + 1 < len(_bound_table) and _bound_table[k + 1] > N:
            return k


_TH = {}
try:
    _TH = json.load(open(os.path.join(ROOT, "results", "thresholds.json")))
except OSError:
    pass


def main():
    recs = json.load(open(os.path.join(ROOT, "data", "records.json")))
    a = {r["k"]: r["a_k"] for r in recs}
    kmax = max(a)

    print("=" * 108)
    print("TABLE 1.  L(N) = longest AP of primes in {1..N}.  L jumps to k at N = a(k).")
    print("=" * 108)
    print(f"{'k':>3} {'a(k)':>19} {'log a(k)':>9} {'k/log a(k)':>11} "
          f"{'rig.bound':>9} {'bound/k':>8} {'Granville':>11} {'a(k)/Gran':>10} "
          f"{'HL a(k)':>11} {'a(k)/HL':>8}")
    rows = []
    for k in range(3, kmax + 1):
        ak = a[k]
        lg = math.log(ak)
        rb = rigorous_L_bound(ak)
        gv = H.granville_ak(k)
        hl = _TH.get(str(k), {}).get("hl") or H.threshold_ak(k)
        rows.append((k, ak, lg, k / lg, rb, gv, hl))
        print(f"{k:>3} {ak:>19} {lg:>9.3f} {k/lg:>11.4f} {rb:>9} {rb/k:>8.3f} "
              f"{gv:>11.4g} {ak/gv:>10.3f} {hl:>11.4g} {ak/hl:>8.3f}")

    best = max(rows, key=lambda r: r[3])
    print(f"\nmax over 3<=k<={kmax} of L/log N  =  {best[3]:.4f}  at k={best[0]}, N={best[1]}")
    lo = [r for r in rows if r[0] >= 11]
    print(f"max over k>=11              =  {max(r[3] for r in lo):.4f}  "
          f"(k={max(lo, key=lambda r: r[3])[0]})")
    print(f"value at k={kmax}                =  {rows[-1][3]:.4f}")

    print()
    print("=" * 108)
    print("TABLE 2.  Rigorous upper bound vs truth, and the conjectured decay rate.")
    print("  rigorous:  L(N) <= max{k : min_last_lower_bound(k) <= N}   (Prop. 2, NOTES.md)")
    print("  conjecture (Granville):  L(N)/log N ~ 2 / (log(L/2) + 1 - gamma)")
    print("=" * 108)
    print(f"{'log10 N':>9} {'log N':>10} {'rigorous L<=':>13} {'bound/logN':>11} "
          f"{'conj. L':>9} {'conj. L/logN':>13}")
    for e in [2, 4, 8, 12, 16, 17, 20, 30, 50, 100, 300, 1000, 10 ** 4, 10 ** 6, 10 ** 9]:
        N = 10 ** e if e <= 20000 else None
        lgN = e * math.log(10)
        rb = rigorous_L_bound(N) if N is not None else None
        # solve k*(1-gamma+log(k/2))/2 = log N  for the Granville threshold
        k = 4.0
        for _ in range(300):
            k = 2 * lgN / (1 - GAMMA + math.log(k / 2))
        s = f"{e:>9} {lgN:>10.4g} "
        s += f"{rb:>13} {rb/lgN:>11.4f} " if rb is not None else f"{'-':>13} {'-':>11} "
        s += f"{k:>9.2f} {k/lgN:>13.4f}"
        print(s)

    print()
    print("Reading: the conjectured ratio falls off like 2/log log N.  To see it drop")
    print("below 1/2 needs N ~ 10^100; below 1/10 needs log log N ~ 20, i.e. N ~ exp(exp(20)).")

    json.dump({"rows": [{"k": r[0], "a_k": r[1], "log_ak": r[2], "ratio": r[3],
                         "rigorous_bound": r[4], "granville": r[5], "hl_threshold": r[6]}
                        for r in rows]}, open(os.path.join(ROOT, "results", "analysis.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
