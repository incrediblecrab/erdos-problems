"""Measured vs Hardy-Littlewood predicted counts of k-APs of primes.

Reads results/count_k*.json produced by apsearch, converts the 2^(1/16) histogram
of last terms into cumulative counts C_k(N), and compares against heuristic.py.

The model covers Case A only (a > k, W_k | d), so Case B APs are subtracted.
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import heuristic as H

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def cumulative(hist, B):
    """[(N, exact count of APs with last term <= N)] from the bucket dict.

    apsearch prints each bucket keyed by its lower edge lo = (1+j/16)*2^e, holding
    last terms in [lo, next_lo).  So the running total through bucket i is the exact
    count of APs with last term <= next_lo - 1; through the final bucket it is <= B.
    """
    items = sorted((float(b), c) for b, c in hist.items())
    out, run = [], 0
    for i, (lo, c) in enumerate(items):
        run += c
        hi = (items[i + 1][0] - 1) if i + 1 < len(items) else B
        out.append((hi, run))
    return out


def main():
    rows = []
    for k in range(3, 30):
        path = os.path.join(ROOT, "results", "count_k%d.json" % k)
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            continue
        d = json.load(open(path))
        B, tot, cb = d["B"], d["total"], d["caseB_count"]
        cum = cumulative(d["hist"], B)
        assert cum[-1][1] == tot, (k, cum[-1][1], tot)
        # Case B APs all have last term k + (k-1)d with (W_k/k) | d; they are rare.
        # Compare at four points spread over the top two decades of the range.
        for frac in (0.01, 0.1, 0.4, 1.0):
            target = B * frac
            cand = [(e, c) for e, c in cum if e <= target and c >= 20]
            if not cand:
                continue
            Nedge, meas = cand[-1]
            pred = H.predicted_count(k, Nedge)
            rows.append((k, Nedge, meas, cb, pred, pred / meas))
    print("%3s %14s %12s %6s %14s %8s" % ("k", "N", "measured", "caseB", "HL predicted", "pred/obs"))
    for k, N, meas, cb, pred, r in rows:
        print("%3d %14.5g %12d %6d %14.6g %8.4f" % (k, N, meas, cb, pred, r))
    if rows:
        lr = [math.log(r[5]) for r in rows]
        gm = math.exp(sum(lr) / len(lr))
        sd = math.exp(math.sqrt(sum((x - math.log(gm)) ** 2 for x in lr) / (len(lr) - 1)))
        print("\n%d comparisons; geometric mean pred/obs = %.4f, geometric sd = %.4f"
              % (len(rows), gm, sd))
        print("range %.4f .. %.4f" % (min(r[5] for r in rows), max(r[5] for r in rows)))
    json.dump([{"k": r[0], "N": r[1], "measured": r[2], "caseB": r[3],
                "predicted": r[4], "ratio": r[5]} for r in rows],
              open(os.path.join(ROOT, "results", "count_comparison.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
