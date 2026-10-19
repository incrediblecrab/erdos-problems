"""Calibration test: is the spread of a(k) around the Hardy-Littlewood threshold
consistent with the model, or does the model drift?

a(k) is the *minimum* last term over a sparse, approximately Poisson set of k-APs,
so a(k)/threshold is expected to fluctuate by an O(1) multiplicative factor.  The
right test is not "is the ratio near 1" but "is the implied probability integral
transform uniform".

Model: let C_k(N) be the expected number of k-APs with last term <= N.  If those
APs behave like a Poisson process in N, then
    P(a(k) > N) = exp(-C_k(N)),   so  U_k := exp(-C_k(a(k)))  ~ Uniform(0,1).
U_k is computed directly from the heuristic -- no local power-law approximation.

A one-sample Kolmogorov-Smirnov test against Uniform(0,1) then asks whether the 24
observed records are jointly consistent with the model.  A systematic drift in the
model (e.g. a wrong constant, or a k-dependent bias) shows up as U clustering.
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import heuristic as H

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    recs = json.load(open(os.path.join(ROOT, "data", "records.json")))
    a = {r["k"]: r["a_k"] for r in recs}
    first = {r["k"]: r["first"] for r in recs}
    # The heuristic models Case A only (first term > k, W_k | d).  For k = 3, 5, 7
    # the record AP is a Case B progression (first term = k), which the model does
    # not cover, so those three k are excluded rather than scored as C_k = 0.
    skipped = [k for k in range(3, 27) if first[k] == k]
    rows = []
    for k in range(3, 27):
        if k in skipped:
            continue
        C = H.predicted_count(k, float(a[k]))
        U = math.exp(-C)
        rows.append((k, a[k], C, U))
    print("excluded (Case B record, outside the model): k = %s\n"
          % ", ".join(map(str, skipped)))
    print("%3s %20s %12s %9s" % ("k", "a(k)", "C_k(a(k))", "U_k"))
    for k, ak, C, U in rows:
        print("%3d %20d %12.5g %9.5f" % (k, ak, C, U))

    U = sorted(r[3] for r in rows)
    n = len(U)
    dplus = max((i + 1) / n - U[i] for i in range(n))
    dminus = max(U[i] - i / n for i in range(n))
    D = max(dplus, dminus)
    # asymptotic KS p-value (Kolmogorov distribution), adequate for n=24
    lam = (math.sqrt(n) + 0.12 + 0.11 / math.sqrt(n)) * D
    p = 2 * sum((-1) ** (j - 1) * math.exp(-2 * j * j * lam * lam) for j in range(1, 101))
    p = min(1.0, max(0.0, p))
    print("\nKS test of U_k against Uniform(0,1):  n=%d  D=%.4f  p=%.3f" % (n, D, p))
    print("mean U = %.4f (uniform expects 0.500)" % (sum(U) / n))
    print("count U<0.1: %d, U>0.9: %d (uniform expects %.1f each)"
          % (sum(1 for u in U if u < 0.1), sum(1 for u in U if u > 0.9), n * 0.1))

    # drift check: regress log C_k(a(k)) on k.  A correct model has no trend.
    ks = [r[0] for r in rows]
    y = [math.log(r[2]) for r in rows]
    mx, my = sum(ks) / n, sum(y) / n
    sxx = sum((x - mx) ** 2 for x in ks)
    slope = sum((ks[i] - mx) * (y[i] - my) for i in range(n)) / sxx
    resid = [y[i] - (my + slope * (ks[i] - mx)) for i in range(n)]
    s2 = sum(r * r for r in resid) / (n - 2)
    se = math.sqrt(s2 / sxx)
    print("\ndrift: slope of log C_k(a(k)) vs k = %+.4f +- %.4f  (t = %+.2f)"
          % (slope, se, slope / se))
    print("a nonzero slope would mean the model is systematically wrong in k.")

    json.dump([{"k": r[0], "a_k": r[1], "C": r[2], "U": r[3]} for r in rows],
              open(os.path.join(ROOT, "results", "calibration.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
