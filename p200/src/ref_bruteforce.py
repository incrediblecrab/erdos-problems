#!/usr/bin/env python3
"""ref_bruteforce.py -- independent reference enumerator for k-term APs of primes.

Shares no code and no mathematical shortcut with src/apsearch.c: it loops over
every common difference d >= 1 and every start a >= 2 with a+(k-1)d <= B, using
only a plain sieve of Eratosthenes.  Its purpose is to validate apsearch.c,
whose speed depends on a structure theorem (W_k | d) that this program does not
assume.

Output: results/reference.json  mapping k -> {B, min_last, total}.
"""
import os
import json
import sys
import time

import numpy as np


def sieve(n):
    p = np.ones(n + 1, dtype=bool)
    p[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if p[i]:
            p[i * i::i] = False
    return p


def enumerate_aps(k, B, Pb):
    """All (a,d), d>=1, a+(k-1)d <= B, every term prime.  Returns (count, min_last)."""
    total = 0
    best = None
    dmax = (B - 2) // (k - 1)
    for d in range(1, dmax + 1):
        L = B - (k - 1) * d + 1          # a ranges over [0, L)
        m = Pb[0:L] & Pb[d:d + L]
        if k >= 3:
            m &= Pb[2 * d:2 * d + L]
        idx = np.flatnonzero(m)
        if idx.size == 0:
            continue
        for i in range(3, k):
            idx = idx[Pb[idx + i * d]]
            if idx.size == 0:
                break
        if idx.size:
            total += int(idx.size)
            last = int(idx[0]) + (k - 1) * d
            if best is None or last < best[0]:
                best = (last, int(idx[0]), d)
    return total, best


def main():
    plan = {3: 20000, 4: 50000, 5: 50000, 6: 200000, 7: 200000,
            8: 300000, 9: 300000, 10: 300000,
            11: 1000000, 12: 1000000, 13: 1000000}
    if len(sys.argv) > 1:
        plan = {int(sys.argv[1]): int(sys.argv[2])}
    Pb = sieve(max(plan.values()))
    out = {}
    for k in sorted(plan):
        B = plan[k]
        t = time.time()
        total, best = enumerate_aps(k, B, Pb)
        rec = {"B": B, "total": total, "seconds": round(time.time() - t, 2)}
        if best:
            rec["min_last"], rec["min_first"], rec["min_diff"] = best
        out[k] = rec
        print(f"k={k:2d} B={B:>8d} total={total:>12d} min_last={rec.get('min_last')} "
              f"({rec['seconds']}s)", flush=True)
    with open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "reference.json"), "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
