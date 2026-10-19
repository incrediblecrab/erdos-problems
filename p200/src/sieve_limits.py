"""sieve_limits.py -- what Selberg's sieve, the large sieve and Gallagher's larger sieve can say about Erdős #200.

Run:  ~/.venvs/main/bin/python src/sieve_limits.py [--plant larger|dusart|collapse]

Exit 0 means every check passed. Each --plant option introduces a known defect, and the run must then exit 1; nothing is written when a defect is planted.

Reads data/records.json, results/barrier.json and src/barrier.py (imported for B(k,N) and the Hardy-Littlewood count), plus the tracker database ../tracker-and-solver/data/erdos.db, opened read-only, if it is present. Writes results/sieve_limits.json, which src/final_check.py re-checks. NOTES.md sections 6 and 8.1 are the write-up.

Sections
  1. corpus    neighbours of #200 in the tracker, and the techniques it records for solved problems tagged 'primes'
  2. records   at N = a(k)-1, where no k-AP of primes exists, a lower bound for the Selberg main term: exact denominator at level N, closed form at level max(T0, N)
  3. grid      the same lower bound in closed form, for k = floor(c log N) up to N = 10^(10^6)
  4. dusart    the explicit bound for sum 1/p that section 3 uses beyond 2e8, tested at every prime below 2e8
  5. barrier   barrier.py's B(k,N) = 2^k k! C(k,N) against the number of Case A candidates
  6. larger    Gallagher's larger sieve applied to the progression itself
"""
import argparse
import bisect
import json
import math
import os
import sqlite3
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DB = os.path.join(os.path.dirname(ROOT), "tracker-and-solver", "data", "erdos.db")
OUT = os.path.join(ROOT, "results", "sieve_limits.json")

sys.path.insert(0, HERE)
import barrier  # noqa: E402

# Meissel-Mertens constant, and Dusart, arXiv:1002.0442v1, Theorem 6.10 (read from the rendered page):
#   sum_{p<=x} 1/p - ln ln x - B <= 1/(10 ln^2 x) + 4/(15 ln^3 x)   for x >= 10372.
MERTENS_B = 0.26149721284764278
P0 = 2 * 10 ** 8      # primes are enumerated exactly up to here
EXACT_E10 = 10 ** 4   # at N = 10^e with e <= this, T0 is computed in exact integer arithmetic
XI_SMALL = 20000      # the DFS is cross-checked against brute force at this xi

ap = argparse.ArgumentParser(description="Limits of sieve methods on Erdős #200; see NOTES.md section 6.")
ap.add_argument("--plant", choices=["larger", "dusart", "collapse"], default=None,
                help="introduce a known defect; the run must then exit 1")
PLANT = ap.parse_args().plant

FAILS = []


def check(name, ok, detail=""):
    print(("  PASS  " if ok else "  FAIL  ") + name + ("  (" + detail + ")" if detail else ""))
    if not ok:
        FAILS.append(name)


def dusart_eps(lnx):
    return 0.0 if PLANT == "dusart" else 1.0 / (10 * lnx ** 2) + 4.0 / (15 * lnx ** 3)


# ---------------------------------------------------------------- primes

def primes_upto(n):
    s = np.ones(n // 2 + 1, dtype=bool)
    s[0] = False
    for i in range(1, math.isqrt(n) // 2 + 1):
        if s[i]:
            p = 2 * i + 1
            s[p * p // 2::p] = False
    odd = 2 * np.nonzero(s)[0] + 1
    return np.concatenate(([2], odd[odd <= n])).astype(np.int64)


t_start = time.time()
P = primes_upto(P0)
PF = P.astype(np.float64)
THETA = np.cumsum(np.log(PF))            # THETA[i] = theta(P[i])
LOGPHI = np.cumsum(np.log1p(-1.0 / PF))  # LOGPHI[i] = log prod_{p<=P[i]} (1-1/p)
INV = np.cumsum(1.0 / PF)                # INV[i] = sum_{p<=P[i]} 1/p
print("sieved %d primes up to %.0e in %.1fs" % (len(P), P0, time.time() - t_start))


def idx(x):
    """Index of the largest prime <= x."""
    return int(np.searchsorted(P, x, side="right")) - 1


def theta(k):
    return float(THETA[idx(k)])


def log_phi_ratio(k):
    """log phi(W_k)/W_k."""
    return float(LOGPHI[idx(k)])


def next_prime(k):
    return int(P[idx(k) + 1])


def primorial(k):
    return math.prod(P[: idx(k) + 1].tolist())


def case_a(N, k):
    """(M, T0), T0 = #{(a, m) : a > k, m >= 1, a + (k-1) m W_k <= N}, the Case A candidates of Prop. 1, exactly."""
    W = primorial(k)
    M = (N - k - 1) // ((k - 1) * W)
    if M <= 0:
        return 0, 0
    return M, M * (N - k) - (k - 1) * W * M * (M + 1) // 2


def log_T0_asymp(L, k):
    """log T0 = 2 log N - log(2(k-1)) - theta(k) + O(1/M) at N = e^L. None where M < e^40, since the O(1/M) could then matter."""
    if L - theta(k) - math.log(k - 1) < 40:
        return None
    return 2 * L - math.log(2 * (k - 1)) - theta(k)


def log_T0_at(e10, k):
    """log T0 at N = 10^e10: exact integers up to EXACT_E10, the asymptotic beyond. None if there are no candidates."""
    if e10 != int(e10):
        raise ValueError("N = 10^%r is not an integer" % e10)
    e10 = int(e10)
    if e10 <= EXACT_E10:
        _, T0 = case_a(10 ** e10, k)
        return math.log(T0) if T0 > 0 else None
    lt = log_T0_asymp(e10 * math.log(10), k)
    if lt is None:
        raise ValueError("asymptotic log T0 is not accurate at log10 N = %d, k = %d" % (e10, k))
    return lt


_S_ALL = {}


def S_bound(k, logxi, xi=None):
    """Upper bound for S = sum_{k<p<=xi} k/(p-k): exact below P0, Dusart's Theorem 6.10 beyond. Pass the integer xi when it is known."""
    i0 = idx(k) + 1
    if xi is not None:
        if xi >= P0:
            raise ValueError("integer xi beyond the prime table")
        return float(np.sum(k / (PF[i0: idx(xi) + 1] - k)))
    if logxi <= math.log(P0):
        # a prime within rounding of xi is included, which can only raise the bound
        return float(np.sum(k / (PF[i0: idx(math.exp(logxi) * (1 + 1e-12)) + 1] - k)))
    if k not in _S_ALL:
        _S_ALL[k] = float(np.sum(k / (PF[i0:] - k)))
    tail_recip = math.log(logxi) + MERTENS_B + dusart_eps(logxi) - float(INV[-1])
    return _S_ALL[k] + k * P0 / (P0 - k) * tail_recip


def log_EJ(J, S):
    """log E_J(S) = log sum_{j<=J} S^j / j!."""
    if J == 0 or S <= 0:
        return 0.0
    j = np.arange(J + 1, dtype=np.float64)
    logfact = np.concatenate(([0.0], np.cumsum(np.log(np.arange(1, J + 1, dtype=np.float64)))))
    t = j * math.log(S) - logfact
    m = t.max()
    return float(m + math.log(np.exp(t - m).sum()))


def J_exact(k, xi):
    """Largest j with (p_k^+)^j <= xi, for integer xi: a squarefree q <= xi has at most this many prime factors above k."""
    p = next_prime(k)
    j, v = 0, p
    while v <= xi:
        j, v = j + 1, v * p
    return j


def J_from_log(k, logxi):
    """The same from log xi in floating point. The 1e-9 can only raise J, which only weakens the bound."""
    return int(logxi / math.log(next_prime(k)) + 1e-9)


# ---------------------------------------------------------------- exact G* by DFS

def gstar_table(k, xi):
    """Primes up to xi with g(p) = 1/(p-1) for p <= k and k/(p-k) for p > k, plus prefix sums."""
    n = idx(xi) + 2                            # one prime beyond xi for the P[j+1] lookahead
    ps = P[:n].tolist()
    pf = PF[:n]
    g = np.where(pf <= k, 1.0 / (pf - 1.0), k / np.maximum(pf - k, 1e-300))
    pref = np.concatenate(([0.0], np.cumsum(g)))
    return ps, g.tolist(), pref.tolist()


def gstar_dfs(k, xi):
    """G*(xi) = sum over squarefree q <= xi of prod_{p|q} g(p). The last prime is summed by prefix sums."""
    ps, g, pref = gstar_table(k, xi)
    off = 0 if PLANT == "collapse" else 1
    sys.setrecursionlimit(10000)

    def rec(q, w, i):
        lim = xi // q
        jmax = bisect.bisect_right(ps, lim) - 1
        tot = w
        if jmax < i:
            return tot
        j = i
        while j <= jmax and ps[j] * ps[j + 1] <= lim:
            tot += rec(q * ps[j], w * g[j], j + 1)
            j += 1
        return tot + w * (pref[jmax + off] - pref[j])

    return rec(1, 1.0, 0)


def gstar_brute(k, xi):
    spf = list(range(xi + 1))
    for p in range(2, math.isqrt(xi) + 1):
        if spf[p] == p:
            for m in range(p * p, xi + 1, p):
                if spf[m] == m:
                    spf[m] = p
    tot = 1.0
    for q in range(2, xi + 1):
        n, w, ok = q, 1.0, True
        while n > 1:
            p = spf[n]
            n //= p
            if n % p == 0:
                ok = False
                break
            w *= 1.0 / (p - 1) if p <= k else k / (p - k)
        if ok:
            tot += w
    return tot


RESULTS = {}

# ================================================================ 1. corpus
print("\n1. corpus: neighbours of #200 in the tracker, and how the solved ones were settled")
NEIGHBOURS = [
    # id, question (paraphrased from the statement), how many values must be prime at once (annotation), solved?
    (219, "arbitrarily long APs of primes?  yes, Green-Tao", "k fixed", 1),
    (1187, "monochromatic prime k-APs in every finite colouring?  yes, from Green-Tao", "k fixed", 1),
    (457, "all p <= (2+eps) log n divide prod_{i<=log n} (n+i), infinitely often?  yes, by a CRT construction", "-", 1),
    (1202, "sifting (p-1)/2 classes modulo k primes always leaves <= eps n?  no, by a construction", "-", 1),
    (1140, "infinitely many n with n-2x^2 prime for all 2x^2 < n?  no", "~sqrt(n/2)", 1),
    (1141, "infinitely many n with n-k^2 prime for all (k,n)=1, k^2 < n?  no", "~sqrt(n)", 1),
    (1142, "any n > 105 with n-2^k prime for all 1 < 2^k < n?", "log2 n", 0),
    (236, "#{k : n-2^k prime} = o(log n)?", "log2 n", 0),
    (852, "longest run of distinct consecutive prime gaps below x: h(x) = o(log x)?", "~log x", 0),
    (1181, "q(n, log n) < (1-c) log^2 n for all large n?", "-", 0),
    (141, "k consecutive primes in AP, for every k?", "k fixed", 0),
    (200, "longest AP of primes in [1,N] has length o(log N)?", "~log N", 0),
]
if not os.path.exists(DB):
    print("  SKIP  tracker database not found at %s: nothing checked, corpus recorded as null" % DB)
    RESULTS["corpus"] = None
else:
    con = sqlite3.connect("file:%s?mode=ro" % DB, uri=True)
    meta = dict(con.execute("SELECT * FROM meta").fetchall())
    rows = []
    for pid, question, plen, want in NEIGHBOURS:
        st, solved = con.execute("SELECT status, is_solved FROM problems WHERE id=?", (pid,)).fetchone()
        rows.append({"id": pid, "status": st, "solved": solved, "question": question, "pattern_length": plen})
        print("  #%-5d %-16s pattern %-11s %s" % (pid, st, plen, question))
    check("tracker statuses are as described", all(r["solved"] == w for r, (_, _, _, w) in zip(rows, NEIGHBOURS)))
    tech, models, nsolved = {}, set(), 0
    for tj, model in con.execute(
            "SELECT s.techniques, s.model FROM problems p JOIN problem_tags t ON t.problem_id=p.id AND t.tag='primes' "
            "JOIN solutions s ON s.problem_id=p.id WHERE p.is_solved=1"):
        nsolved += 1
        models.add(model)
        for t in json.loads(tj or "[]"):
            tech[t] = tech.get(t, 0) + 1
    con.close()
    print("  technique labels on the %d solved problems tagged 'primes' (written by %s, not by the site):"
          % (nsolved, ", ".join(sorted(models))))
    for t, n in sorted(tech.items(), key=lambda x: (-x[1], x[0])):
        print("    %2d  %s" % (n, t))
    RESULTS["corpus"] = {"tracker_last_build": meta.get("last_build"), "neighbours": rows,
                         "solved_primes_problems": nsolved, "techniques": tech,
                         "technique_label_models": sorted(models)}

# ================================================================ 2. records
print("\n2. records: can Selberg's sieve prove the known fact that no k-AP of primes ends below a(k)?")
print("   T0 = Case A candidates (a,d); G* = the Selberg denominator at xi = sqrt(N), computed exactly, an upper bound for every d.")
print("   log(T0/G*) > 0 means the main term alone exceeds 1, so the sieve cannot prove it.")
print("   bound = log((W/phi) E_J(S)), the closed form used in section 3.")
print("   joint = log(T0) - bound at level max(T0, N), for the sieve over pairs (a, m), whose natural level is at most T0.")

for kk in (5, 12, 26):
    a, b = gstar_dfs(kk, XI_SMALL), gstar_brute(kk, XI_SMALL)
    check("DFS G* equals brute force at xi=%d, k=%d" % (XI_SMALL, kk), abs(a - b) <= 1e-9 * b,
          "%.10g vs %.10g" % (a, b))

records = {r["k"]: r["a_k"] for r in json.load(open(os.path.join(ROOT, "data", "records.json")))}
print("%4s %18s %8s %3s %8s %8s %8s %8s %9s %8s %8s %6s" %
      ("k", "N=a(k)-1", "xi", "J", "log T0", "log G*", "bound", "T0/G*", "T0/bound", "joint", "log B", "secs"))
rec_rows = []
for k in range(8, 27):
    N = records[k] - 1
    M, T0 = case_a(N, k)
    xi = math.isqrt(N)
    ts = time.time()
    G = gstar_dfs(k, xi)
    secs = time.time() - ts
    J = J_exact(k, xi)
    S = S_bound(k, math.log(xi), xi)
    bound = -log_phi_ratio(k) + log_EJ(J, S)
    xj = math.isqrt(max(T0, N))
    Jj = J_exact(k, xj)
    Sj = S_bound(k, math.log(xj), xj if xj < P0 else None)
    joint = math.log(T0) + log_phi_ratio(k) - log_EJ(Jj, Sj)
    lB = barrier.log_sieve_bound(k, math.log(N))
    row = {"k": k, "N": N, "xi": xi, "M": M, "T0": T0, "log_T0": math.log(T0), "J": J, "S": S,
           "log_Gstar": math.log(G), "log_Gstar_bound": bound, "log_main_lower": math.log(T0) - math.log(G),
           "log_main_closed": math.log(T0) - bound, "xi_joint": xj, "J_joint": Jj, "S_joint": Sj,
           "log_main_joint": joint, "log_B": lB}
    rec_rows.append(row)
    print("%4d %18d %8.2g %3d %8.3f %8.3f %8.3f %8.3f %9.3f %8.3f %8.3f %6.1f"
          % (k, N, xi, J, row["log_T0"], row["log_Gstar"], bound, row["log_main_lower"],
             row["log_main_closed"], joint, lB, secs))
check("exact G* <= (W/phi) E_J(S) at every record", all(r["log_Gstar"] <= r["log_Gstar_bound"] + 1e-9 for r in rec_rows))
check("Selberg main term > 1 at every record (the sieve cannot prove L(a(k)-1) < k)",
      all(r["log_main_lower"] > 0 for r in rec_rows),
      "min log(T0/G*) = %.3f" % min(r["log_main_lower"] for r in rec_rows))
check("the closed form alone shows it too (T0 > (W/phi) E_J(S) at every record)",
      all(r["log_main_closed"] > 0 for r in rec_rows),
      "min %.3f" % min(r["log_main_closed"] for r in rec_rows))
check("the joint sieve over (a, m) cannot either, at level max(T0, N), by the closed form",
      all(r["log_main_joint"] > 0 for r in rec_rows),
      "min %.3f, level up to N^%.3f" % (min(r["log_main_joint"] for r in rec_rows),
                                       max(2 * math.log(r["xi_joint"]) / math.log(r["N"]) for r in rec_rows)))
check("barrier.py's B exceeds the candidate count T0 at every record",
      all(r["log_B"] > r["log_T0"] for r in rec_rows))
RESULTS["records"] = rec_rows

# ================================================================ 3. grid
print("\n3. grid: k = floor(c log N); Selberg or large sieve at level D = N^A, xi = sqrt(D)")
print("   need = log T / k, the nats per term a proof must save; sieve = log E_J(S) / k, the most it can save;")
print("   HL = (log T - log C_HL) / k, the per-term deficit the Hardy-Littlewood model predicts (heuristic);")
print("   the model's expected count is exp(k (need - HL)), so HL > need is where it expects no k-AP.")
print("%9s %5s %9s %3s %11s %11s %7s %7s %7s %8s %9s" %
      ("log10 N", "c", "k", "A", "log T", "log E_J", "need", "sieve", "HL", "margin", "A needed"))
grid, asymp_err = [], []
for e in (10, 16, 30, 100, 1000, 10 ** 4, 10 ** 5, 10 ** 6):
    L = e * math.log(10)
    for c in (0.3, 0.5, 0.7, 0.9, 1.0):
        k = int(c * L)
        lT0 = log_T0_at(e, k)
        if lT0 is None:
            continue
        if e <= EXACT_E10 and log_T0_asymp(L, k) is not None:
            asymp_err.append(abs(log_T0_asymp(L, k) - lT0))
        lT = lT0 + log_phi_ratio(k)
        hl = (lT - barrier.log_C_asymp(k, L)) / k if k <= 20000 else None

        def logE(A):
            logxi = A * L / 2
            return log_EJ(J_from_log(k, logxi), S_bound(k, logxi))

        lo, hi = 1.0, 2.0
        while logE(hi) < lT and hi < 1e4:
            lo, hi = hi, 2 * hi
        for _ in range(40):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if logE(mid) < lT else (lo, mid)
        for A in (1, 2):
            lE = logE(A)
            row = {"log10N": e, "c": c, "k": k, "A": A, "T0_exact": e <= EXACT_E10, "log_T": lT, "log_EJ": lE,
                   "J": J_from_log(k, A * L / 2), "need_per_term": lT / k, "sieve_per_term": lE / k,
                   "hl_per_term": hl, "margin": lT - lE, "A_needed": hi}
            grid.append(row)
            print("%9g %5.1f %9d %3d %11.4g %11.4g %7.3f %7.3f %7s %8.3g %9.3g"
                  % (e, c, k, A, lT, lE, lT / k, lE / k, "%.3f" % hl if hl else "-", lT - lE, hi))
check("the asymptotic log T0 used beyond 10^%d agrees with exact integers wherever both apply" % EXACT_E10,
      bool(asymp_err) and max(asymp_err) < 1e-9,
      "%d points, max |diff| %.2g" % (len(asymp_err), max(asymp_err) if asymp_err else float("nan")))
check("sieve vacuous (log T > log E_J) at every grid point, A = 1 and 2", all(r["margin"] > 0 for r in grid))
check("level needed exceeds N^2 at every grid point", all(r["A_needed"] > 2 for r in grid),
      "min A needed %.3g" % min(r["A_needed"] for r in grid))
big = [r for r in grid if r["log10N"] >= 1000 and r["A"] == 1]
check("per-term sieve supply falls as N grows (log10 N >= 10^3, A = 1, each c)",
      all(r1["sieve_per_term"] > r2["sieve_per_term"] for r1, r2 in zip(big, big[5:]) if r1["c"] == r2["c"]))
RESULTS["grid"] = grid

# ================================================================ 4. dusart
print("\n4. dusart: the bound for sum 1/p used above for xi > %.0e, tested wherever it can be" % P0)
ii = np.nonzero(P >= 10372)[0]
lnx = np.log(PF[ii])
lhs = INV[ii] - np.log(lnx) - MERTENS_B
rhs = 0.0 if PLANT == "dusart" else 1 / (10 * lnx ** 2) + 4 / (15 * lnx ** 3)
viol = int(np.sum(lhs > rhs))
check("Dusart Thm 6.10 upper bound holds at every prime 10372 <= p <= %.0e" % P0, viol == 0,
      "%d violations of %d" % (viol, len(ii)))
RESULTS["dusart"] = {"primes_tested": len(ii), "violations": viol, "range": [10372, P0]}

# ================================================================ 5. barrier
print("\n5. barrier.py: is B(k,N) = 2^k k! C(k,N) a possible Selberg output?  The main term is <= T0, since G >= 1.")
bj = json.load(open(os.path.join(ROOT, "results", "barrier.json")))["grid"]
cmp_rows, asymp_ratio = [], []
for r in bj:
    if r["log10N"] > 10 ** 6:
        continue                                 # k is beyond the prime table
    L = r["log10N"] * math.log(10)
    k = r["k"]
    lT0 = log_T0_at(r["log10N"], k)
    if lT0 is None:
        continue
    logxi = L / 2
    lE = log_EJ(J_from_log(k, logxi), S_bound(k, logxi))
    x = {"log10N": r["log10N"], "c": r["c"], "k": k, "log10B": r["log10B"], "log10T0": lT0 / math.log(10),
         "log_main_lower": lT0 + log_phi_ratio(k) - lE}
    cmp_rows.append(x)
    if k >= 10 ** 5:
        lk = math.log(k)
        pred = (k * math.log(2 * math.exp(barrier.GAMMA - 1) * r["c"] * lk) - k / (2 * lk)
                + 0.5 * math.log(2 * math.pi * k) - barrier.GAMMA - math.log(lk)) / math.log(10)
        asymp_ratio.append(pred / (x["log10B"] - x["log10T0"]))
over = [x for x in cmp_rows if x["log10B"] > x["log10T0"]]
gap = [x["log10B"] - x["log10T0"] for x in cmp_rows]
print("  %d of %d points in results/barrier.json with log10 N <= 10^6 have B > T0; log10(B/T0) ranges %.3g .. %.3g"
      % (len(over), len(cmp_rows), min(gap), max(gap)))
for x in cmp_rows:
    if x["log10N"] in (2, 16, 100, 10 ** 4, 10 ** 6) and x["c"] in (0.5, 0.7, 1.0):
        print("    log10 N=%-8g c=%.1f k=%-8d log10 B=%-12.5g log10 T0=%-12.5g" %
              (x["log10N"], x["c"], x["k"], x["log10B"], x["log10T0"]))
print("  points with B <= T0, where the comparison is inconclusive:")
for x in cmp_rows:
    if x["log10B"] <= x["log10T0"]:
        print("    log10 N=%-8g c=%.1f k=%-4d log10 B=%-9.4g log10 T0=%-9.4g" %
              (x["log10N"], x["c"], x["k"], x["log10B"], x["log10T0"]))
strong = [x for x in cmp_rows if x["c"] >= 0.7 or x["log10N"] >= 300]
check("B > T0 wherever c >= 0.7 or log10 N >= 300", all(x["log10B"] > x["log10T0"] for x in strong),
      "%d points" % len(strong))
check("log10(B/T0) matches (2 e^(gamma-1) c log k)^k e^(-k/(2 log k)) sqrt(2 pi k)/(e^gamma log k) to 1% for k >= 10^5",
      bool(asymp_ratio) and all(abs(q - 1) < 0.01 for q in asymp_ratio),
      "%d points, ratios %s" % (len(asymp_ratio), ", ".join("%.5f" % q for q in asymp_ratio)))
check("yet the sieve is vacuous at every such point with candidates (T/E_J > 1 at level N)",
      all(x["log_main_lower"] > 0 for x in cmp_rows),
      "min log(T/E_J) = %.3g" % min(x["log_main_lower"] for x in cmp_rows))
RESULTS["barrier_vs_T0"] = cmp_rows

# ================================================================ 6. larger sieve
print("\n6. Gallagher's larger sieve on the progression itself (Case A):")
print("   nu(p) = 1 for p <= k (p | d), nu(p) = k for k < p <= z.  Claim: bound < k  iff  theta(k) > log N, for every z.")
agree = total = 0
for k in range(3, 81):
    for e in range(2, 41):
        L = e * math.log(10)
        for z in (k, 2 * k, k * k, 10 ** 4, 10 ** 7):
            if z < k:
                continue
            th_k, th_z = theta(k), theta(z)
            num = th_z - L
            den = (th_k / k if PLANT == "larger" else th_k) + (th_z - th_k) / k - L
            proves = den > 0 and num / den < k
            total += 1
            agree += proves == (th_k > L)
print("  %d of %d (k, N, z) cases agree" % (agree, total))
check("larger sieve reproduces exactly the primorial condition", agree == total)
RESULTS["larger_sieve"] = {"cases": total, "agree": agree}

RESULTS["fails"] = FAILS
if PLANT:
    print("\nplanted defect %r: %s not written" % (PLANT, os.path.relpath(OUT, ROOT)))
else:
    with open(OUT, "w") as f:
        json.dump(RESULTS, f, indent=1)
        f.write("\n")
    print("\nwrote %s" % os.path.relpath(OUT, ROOT))
print("%.0fs; " % (time.time() - t_start) + ("%d checks failed: %s" % (len(FAILS), FAILS) if FAILS else "all checks passed"))
sys.exit(1 if FAILS else 0)
