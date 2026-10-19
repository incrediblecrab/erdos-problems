"""final_check.py -- re-derive every claim in NOTES.md from the stored artifacts.

This is the trusted validation path.  It shares no code with apsearch.c and does not
use heuristic.py; it re-tests primality independently, re-checks the structure theorem
against all 26 known records, re-checks the rigorous upper bound against them, and
re-derives the range of N over which L(N) has been independently determined here.

Section G re-derives the sieve-limit claims of NOTES.md sections 6 and 8.1 from data/records.json, results/barrier.json and results/sieve_limits.json, without importing src/sieve_limits.py or numpy: it has its own prime sieve up to sqrt(a(26)), which makes it the slow part of this script.

Exit code is nonzero if any check fails.
"""
import array
import bisect
import itertools
import json
import math
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FAIL = []
# Dusart, arXiv:1002.0442v1, Theorem 6.10: sum_{p<=x} 1/p <= ln ln x + B + 1/(10 ln^2 x) + 4/(15 ln^3 x) for x >= 10372
MERTENS_B = 0.26149721284764278
EULER_GAMMA = 0.57721566490153286


def check(name, cond, detail=""):
    print("  [%s] %s%s" % ("ok" if cond else "FAIL", name, ("  " + detail) if detail else ""))
    if not cond:
        FAIL.append(name)


def is_prime(n):
    """Deterministic Miller-Rabin for n < 3.3e24 (first 13 primes as bases)."""
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41):
        if n % p == 0:
            return n == p
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41):
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def primorial(k):
    r = 1
    for p in primes_upto(k):
        r *= p
    return r


def primes_array(n):
    """All primes <= n: the sieve of primes_upto, collected with itertools.compress so that n ~ 1e8 is practical."""
    s = bytearray([1]) * (n + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, math.isqrt(n) + 1):
        if s[i]:
            s[i * i::i] = bytes(len(range(i * i, n + 1, i)))
    return array.array("q", itertools.compress(range(n + 1), s))


def log_EJ(J, S):
    """log sum_{j<=J} S^j/j!, the closed-form bound of NOTES.md section 6.1."""
    if J == 0:
        return 0.0
    t = [j * math.log(S) - math.lgamma(j + 1) for j in range(J + 1)]
    m = max(t)
    return m + math.log(math.fsum(math.exp(x - m) for x in t))


def section_g():
    print("\nG. limits of sieve methods (NOTES.md sections 6 and 8.1), re-derived without src/sieve_limits.py")
    path = os.path.join(ROOT, "results", "sieve_limits.json")
    if not os.path.exists(path):
        check("results/sieve_limits.json present", False, "missing")
        return
    R = json.load(open(path))
    check("sieve_limits.py recorded no failed checks", R["fails"] == [], str(R["fails"]))
    a = {r["k"]: r["a_k"] for r in json.load(open(os.path.join(ROOT, "data", "records.json")))}
    rec = {r["k"]: r for r in R["records"]}

    t0 = time.time()
    P1 = math.isqrt(a[26] - 1)                   # the largest xi = sqrt(N) at a record
    P = primes_array(P1)
    print("     %d primes up to %d in %.1fs" % (len(P), P1, time.time() - t0))
    INV = math.fsum(1.0 / p for p in P)
    memo = {}

    def upto(x):
        return bisect.bisect_right(P, x)

    def cached(key, f):
        if key not in memo:
            memo[key] = f()
        return memo[key]

    def theta(x):
        return cached(("theta", x), lambda: math.fsum(math.log(p) for p in P[:upto(x)]))

    def log_W_phi(k):
        """log W_k / phi(W_k)."""
        return cached(("wphi", k), lambda: -math.fsum(math.log1p(-1.0 / p) for p in P[:upto(k)]))

    def S_to(k, x):
        """An upper bound for sum_{k<p<=x} k/(p-k): exact up to P1, Dusart's Theorem 6.10 beyond, using p/(p-k) < P1/(P1-k) there."""
        if x <= P1:
            return math.fsum(k / (p - k) for p in P[upto(k):upto(x)])
        L = math.log(x)
        tail = math.log(L) + MERTENS_B + 1 / (10 * L * L) + 4 / (15 * L ** 3) - INV
        return cached(("S", k), lambda: math.fsum(k / (p - k) for p in P[upto(k):])) + k * P1 / (P1 - k) * tail

    def J_of(k, x):
        """The most prime factors above k that a squarefree q <= x can have."""
        p = P[upto(k)]
        j, v = 0, p
        while v <= x:
            j, v = j + 1, v * p
        return j

    def case_a(N, k):
        """(M, T0): T0 = #{(a, m) : a > k, m >= 1, a + (k-1) m W_k <= N}, in exact integers."""
        W = math.prod(P[:upto(k)])
        M = (N - k - 1) // ((k - 1) * W)
        return (M, M * (N - k) - (k - 1) * W * M * (M + 1) // 2) if M > 0 else (0, 0)

    def log_T0_at(e10, k):
        """log T0 at N = 10^e10: exact integers up to e10 = 10^4, 2 log N - log(2(k-1)) - theta(k) beyond (NaN where that is inaccurate)."""
        if e10 <= 10 ** 4:
            T0 = case_a(10 ** e10, k)[1]
            return math.log(T0) if T0 else None
        L = e10 * math.log(10)
        if L - theta(k) - math.log(k - 1) < 40:
            return float("nan")
        return 2 * L - math.log(2 * (k - 1)) - theta(k)

    print("   G1. the 19 records, N = a(k)-1 for k = 8..26, where no k-AP of primes exists")
    bad, level_n, joint, lt0 = [], [], [], {}
    for k in range(8, 27):
        r, N = rec[k], a[k] - 1
        M, T0 = case_a(N, k)
        lt0[k] = math.log(T0)
        xi = math.isqrt(N)
        J, S = J_of(k, xi), S_to(k, xi)
        bound = log_W_phi(k) + log_EJ(J, S)
        xj = math.isqrt(max(T0, N))
        lj = lt0[k] - log_W_phi(k) - log_EJ(J_of(k, xj), S_to(k, xj))
        if ((N, M, T0, J) != (r["N"], r["M"], r["T0"], r["J"]) or abs(S - r["S"]) > 1e-9 * S
                or abs(bound - r["log_Gstar_bound"]) > 1e-9 * bound):
            bad.append(k)
        level_n.append((lt0[k] - bound, k, r["log_Gstar"] <= bound + 1e-9))
        joint.append((lj, k, lj <= r["log_main_joint"] + 1e-9 * abs(lj)))
    check("T0, M, J, S and log (W/phi) E_J(S) recomputed at every record match the stored values", not bad,
          "mismatch at k = %s" % bad if bad else "")
    check("NOTES: T0 > (W/phi) E_J(S) at every record, so the Selberg main term at level N exceeds 1; least e^0.78 at k = 10",
          min(level_n)[0] > 0 and math.floor(100 * min(level_n)[0]) == 78 and min(level_n)[1] == 10,
          "min log %.5f at k=%d" % min(level_n)[:2])
    small = [k for k in rec if case_a(a[k] - 1, k)[1] < a[k] - 1]
    lvl26 = 2 * math.log(math.isqrt(max(case_a(a[26] - 1, 26)[1], a[26] - 1))) / math.log(a[26] - 1)
    check("NOTES: M = 1 and T0 = 190, 189, 188 at k = 8, 9, 10; T0 < N exactly at k = 8, 9, 10, 13; the joint level reaches N^1.376 at k = 26",
          [case_a(a[k] - 1, k) for k in (8, 9, 10)] == [(1, 190), (1, 189), (1, 188)] and sorted(small) == [8, 9, 10, 13]
          and round(lvl26, 3) == 1.376, "T0 < N at k = %s; N^%.4f" % (sorted(small), lvl26))
    check("the stored exact G* never exceeds that closed form", all(x[2] for x in level_n))
    check("the joint sieve's closed form at level max(T0, N) exceeds 1 at every record, and is at most the stored value",
          min(joint)[0] > 0 and all(x[2] for x in joint), "min log %.3f at k=%d" % min(joint)[:2])
    ml = min((r["log_main_lower"], k) for k, r in rec.items())
    check("NOTES: log(T0/G*) is least at k = 10, 1.419, and is 38.9 at k = 26",
          (round(ml[0], 3), ml[1]) == (1.419, 10) and round(rec[26]["log_main_lower"], 1) == 38.9,
          "%.3f at k=%d; %.3f at k=26" % (ml[0], ml[1], rec[26]["log_main_lower"]))
    rb = [(rec[k]["log_B"] - lt0[k]) / math.log(10) for k in rec]
    check("NOTES: barrier.py's B exceeds T0 at all 19 records, by 10^4.80 to 10^11.27",
          min(rb) > 0 and math.floor(100 * min(rb)) == 480 and round(max(rb), 2) == 11.27,
          "10^%.4f .. 10^%.4f" % (min(rb), max(rb)))

    print("   G2. the grid k = floor(c log N), log10 N = 10 .. 10^6, level N^A with A = 1, 2")
    grid = R["grid"]
    badT = [g for g in grid if not abs(log_T0_at(g["log10N"], g["k"]) - log_W_phi(g["k"]) - g["log_T"])
            <= 1e-9 * abs(g["log_T"])]
    check("log T recomputed at all %d grid rows matches (exact integers to 10^(10^4), asymptotic beyond)" % len(grid),
          not badT, "%d mismatched" % len(badT))
    ax = []
    for g in grid:
        if g["log10N"] <= 10 ** 4 and g["A"] == 1:
            L = g["log10N"] * math.log(10)
            if L - theta(g["k"]) - math.log(g["k"] - 1) >= 40:
                ax.append(abs(2 * L - math.log(2 * (g["k"] - 1)) - theta(g["k"]) - log_T0_at(g["log10N"], g["k"])))
    check("NOTES: exact log T0 at 30 of the 40 (c, N) points; the asymptotic form agrees to within 5e-11, rounding error, at the 14 where both apply",
          sum(1 for g in grid if g["A"] == 1 and g["T0_exact"]) == 30 and len(ax) == 14 and max(ax) < 5e-11,
          "%d points, max |diff| %.3g" % (len(ax), max(ax) if ax else float("nan")))
    mm = min(grid, key=lambda g: g["log_T"] - g["log_EJ"])
    check("NOTES: the main term exceeds 1 at every grid row, by at least e^3.05 (log10 N = 10, c = 1, A = 2)",
          all(g["log_T"] - g["log_EJ"] > 0 for g in grid) and math.floor(100 * (mm["log_T"] - mm["log_EJ"])) == 305
          and (mm["log10N"], mm["c"], mm["A"]) == (10, 1.0, 2),
          "min %.4f at log10 N=%g, c=%.1f, A=%d" % (mm["log_T"] - mm["log_EJ"], mm["log10N"], mm["c"], mm["A"]))
    spot = [g for g in grid if g["log10N"] == 10]
    sbad = []
    for g in spot:
        xi = 10 ** (5 * g["A"])
        lE = log_EJ(J_of(g["k"], xi), S_to(g["k"], xi))
        if not (lE >= g["log_EJ"] - 1e-9 * abs(lE) and g["log_T"] - lE > 0):
            sbad.append((g["c"], g["A"]))
    check("log E_J re-derived at the %d rows with log10 N = 10, where margins are smallest: at least the stored value, and still below log T"
          % len(spot), not sbad, str(sbad) if sbad else "")
    mA = min(grid, key=lambda g: g["A_needed"])
    check("NOTES: the level needed exceeds N^2 everywhere, least N^2.3398 at log10 N = 10, c = 1",
          mA["A_needed"] > 2 and round(mA["A_needed"], 4) == 2.3398 and (mA["log10N"], mA["c"]) == (10, 1.0),
          "N^%.4f at log10 N=%g, c=%.1f" % (mA["A_needed"], mA["log10N"], mA["c"]))
    brk = []
    for g in spot:
        if g["A"] == 1:
            lo, hi = (10 ** (5 * A) for A in (2.33, 2.35))
            brk.append((g["c"], g["log_T"] - log_EJ(J_of(g["k"], lo), S_to(g["k"], lo)) > 0,
                        g["log_T"] - log_EJ(J_of(g["k"], hi), S_to(g["k"], hi)) > 0))
    check("NOTES: re-derived here, the closed form still exceeds 1 at level N^2.33 at every log10 N = 10 row, and at c = 1 no longer does at N^2.35",
          all(b[1] for b in brk) and [b[2] for b in brk if b[0] == 1.0] == [False], str(brk))
    at6 = [round(g["log_T"] / g["k"], 3) for g in sorted(grid, key=lambda g: g["c"]) if g["log10N"] == 10 ** 6 and g["A"] == 1]
    sup = {A: {c: [g["log_EJ"] / g["k"] for g in sorted(grid, key=lambda g: g["log10N"]) if g["A"] == A and g["c"] == c]
               for c in (0.3, 0.5, 0.7, 0.9, 1.0)} for A in (1, 2)}
    flat = {A: [x for v in sup[A].values() for x in v] for A in (1, 2)}
    falls = all(s[i] > s[i + 1] for A in (1, 2) for s in sup[A].values() for i in range(1, len(s) - 1))
    hl = next(g for g in grid if (g["log10N"], g["c"], g["A"]) == (16, 0.7, 1))
    check("NOTES: per term a proof needs 5.668, 3.001, 1.858, 1.223, 1.001 nats at log10 N = 10^6; the sieve supplies 0.23-1.62 at level N and 0.42-2.79 at N^2, falling from log10 N = 16 on; 1.951 needed and 0.62 supplied at k = 25, log10 N = 16",
          at6 == [5.668, 3.001, 1.858, 1.223, 1.001] and all(0.23 <= x <= 1.62 for x in flat[1])
          and all(0.42 <= x <= 2.79 for x in flat[2]) and falls and hl["k"] == 25
          and round(hl["log_T"] / 25, 3) == 1.951 and round(hl["log_EJ"] / 25, 2) == 0.62,
          "need %s; supply %.4f-%.4f, %.4f-%.4f; falls %s" % (at6, min(flat[1]), max(flat[1]), min(flat[2]), max(flat[2]), falls))

    print("   G3. barrier.py's B(k,N) = 2^k k! C(k,N) against the candidate count T0")
    pts = []
    for r in json.load(open(os.path.join(ROOT, "results", "barrier.json")))["grid"]:
        if r["log10N"] <= 10 ** 6:
            lt = log_T0_at(int(r["log10N"]), r["k"])
            if lt is not None:
                pts.append((r, lt / math.log(10)))
    over = [r for r, t in pts if r["log10B"] > t]
    under = [r for r, t in pts if not r["log10B"] > t]
    strong = [(r, t) for r, t in pts if r["c"] >= 0.7 or r["log10N"] >= 300]
    check("NOTES: B > T0 at 69 of the 82 points with log10 N <= 10^6", (len(over), len(pts)) == (69, 82),
          "%d of %d" % (len(over), len(pts)))
    check("NOTES: including all 57 with c >= 0.7 or log10 N >= 300",
          len(strong) == 57 and all(r["log10B"] > t for r, t in strong), "%d points" % len(strong))
    check("NOTES: the 13 others have c <= 0.5 and log10 N <= 100",
          len(under) == 13 and all(r["c"] <= 0.5 and r["log10N"] <= 100 for r in under))
    gmax = max(r["log10B"] - t for r, t in pts)
    check("NOTES: log10(B/T0) reaches 2.92e6", round(gmax / 1e6, 2) == 2.92, "%.6g" % gmax)
    ratios = []
    for r, t in pts:
        if r["k"] >= 10 ** 5:
            k, lk = r["k"], math.log(r["k"])
            pred = (k * math.log(2 * math.exp(EULER_GAMMA - 1) * r["c"] * lk) - k / (2 * lk)
                    + 0.5 * math.log(2 * math.pi * k) - EULER_GAMMA - math.log(lk))
            ratios.append((r["log10B"] - t) * math.log(10) / pred)
    check("NOTES: at the 6 points with k >= 10^5, log(B/T0) matches (2e^(gamma-1) c log k)^k e^(-k/(2 log k)) sqrt(2 pi k)/(e^gamma log k) to within 0.021%",
          len(ratios) == 6 and max(abs(q - 1) for q in ratios) < 0.00021, ", ".join("%.6f" % q for q in ratios))
    stored = {(x["log10N"], x["c"]): x["log_main_lower"] for x in R["barrier_vs_T0"]}
    redo = []
    for r, t in pts:
        if r["log10N"] <= 16:
            xi = math.isqrt(10 ** int(r["log10N"]))
            lm = t * math.log(10) - log_W_phi(r["k"]) - log_EJ(J_of(r["k"], xi), S_to(r["k"], xi))
            redo.append((lm, abs(lm - stored[(r["log10N"], r["c"])]) <= 1e-9 * max(1.0, abs(lm))))
    mb = min(stored.items(), key=lambda x: x[1])
    check("main term at level N re-derived at the %d points with log10 N <= 16 matches the stored values" % len(redo),
          all(x[1] for x in redo))
    check("NOTES: at level N it exceeds 1 at all 82 points, least e^2.40 at log10 N = 2, c = 0.9",
          mb[1] > 0 and round(mb[1], 2) == 2.40 and mb[0] == (2, 0.9) and abs(min(x[0] for x in redo) - mb[1]) < 1e-9,
          "%.4f at %s" % (mb[1], mb[0]))

    print("   G4. Gallagher's larger sieve on the progression")
    agree = total = 0
    for k in range(3, 81):
        for e in range(2, 41):
            L = e * math.log(10)
            for z in (k, 2 * k, k * k, 10 ** 4, 10 ** 7):
                den = theta(k) + (theta(z) - theta(k)) / k - L
                total += 1
                agree += (den > 0 and (theta(z) - L) / den < k) == (theta(k) > L)
    ls = R["larger_sieve"]
    check("bound < k exactly when theta(k) > log N, in all %d (k, N, z) cases" % total,
          agree == total == ls["cases"] == ls["agree"], "%d of %d here; %d of %d stored" % (agree, total, ls["agree"], ls["cases"]))

    print("   G5. Dusart's bound for sum 1/p, which the grid relies on beyond 2e8")
    acc, n, viol = 0.0, 0, 0
    for p in P:
        acc += 1.0 / p
        if p >= 10372:
            L = math.log(p)
            n += 1
            viol += acc - math.log(L) - MERTENS_B > 1 / (10 * L * L) + 4 / (15 * L ** 3)
    check("holds at every prime 10372 <= p <= %d, re-tested here" % P1, viol == 0, "%d violations of %d" % (viol, n))
    d = R["dusart"]
    check("and sieve_limits.py found none among the %d primes up to 2e8" % d["primes_tested"],
          d["violations"] == 0 and d["primes_tested"] == 11077665 and d["range"] == [10372, 200000000])

    print("   G6. the tracker survey")
    c = R["corpus"]
    if c is None:
        print("     skipped: the tracker database was absent when sieve_limits.py ran")
    else:
        want = {219: "PROVED (LEAN)", 1187: "SOLVED", 457: "PROVED (LEAN)", 1202: "SOLVED", 1140: "DISPROVED",
                1141: "DISPROVED (LEAN)", 1142: "OPEN", 236: "OPEN", 852: "OPEN", 1181: "OPEN", 141: "OPEN", 200: "OPEN"}
        print("     tracker build %s, as recorded by sieve_limits.py" % c["tracker_last_build"])
        check("NOTES: neighbour statuses as tabulated", {x["id"]: x["status"] for x in c["neighbours"]} == want)
        tech = {"sieve method": 9, "analytic estimate": 6, "explicit construction": 6, "elementary number theory": 5,
                "reduction to known theorem": 3, "covering / density argument": 1, "fourier analysis": 1,
                "greedy algorithm": 1}
        check("NOTES: technique labels on the 16 solved 'primes' problems as quoted",
              c["solved_primes_problems"] == 16 and c["techniques"] == tech
              and c["technique_label_models"] == ["claude-sonnet-5", "gemini-3.8-flash"])
    print("     section G took %.0fs" % (time.time() - t0))


def main():
    recs = json.load(open(os.path.join(ROOT, "data", "records.json")))
    a = {r["k"]: r["a_k"] for r in recs}
    first = {r["k"]: r["first"] for r in recs}
    diff = {r["k"]: r["diff"] for r in recs}

    print("A. the 26 records are genuine k-APs of distinct primes")
    allok = True
    for r in recs:
        k, f, d, last = r["k"], r["first"], r["diff"], r["a_k"]
        terms = [f + i * d for i in range(k)]
        ok = (terms[-1] == last and len(set(terms)) == k
              and all(is_prime(t) for t in terms))
        allok &= ok
        if not ok:
            print("     k=%d FAILED" % k)
    check("all 26 records verified by independent Miller-Rabin", allok)
    check("a(k) strictly increasing", all(a[k] < a[k + 1] for k in range(1, 26)))

    print("\nB. structure theorem (Prop. 1) holds for every record with k >= 3")
    allok = True
    for k in range(3, 27):
        f, d, W = first[k], diff[k], primorial(k)
        caseA = (f > k) and (d % W == 0)
        caseB = (is_prime(k) and f == k and d % (W // k) == 0)
        if not (caseA ^ caseB):
            allok = False
            print("     k=%d: caseA=%s caseB=%s  (f=%d d=%d W=%d)" % (k, caseA, caseB, f, d, W))
    check("every record is exactly one of Case A / Case B", allok)
    cb = [k for k in range(3, 27) if first[k] == k]
    check("Case B records identified", cb == [3, 5, 7], "k = %s" % cb)

    print("\nC. rigorous upper bound (Prop. 2) is satisfied by every record")
    allok = True
    for k in range(3, 27):
        W = primorial(k)
        nxt = next(p for p in range(k + 1, 4 * k) if is_prime(p))
        lb = nxt + (k - 1) * W
        if is_prime(k):
            lb = min(lb, k + (k - 1) * (W // k))
        if a[k] < lb:
            allok = False
            print("     k=%d: a(k)=%d < lower bound %d" % (k, a[k], lb))
    check("a(k) >= min_last_lower_bound(k) for k=3..26", allok)

    # the weaker but cleaner form: W_k <= a(k) * k/(k-1)
    allok = all(primorial(k) <= a[k] * k / (k - 1) for k in range(3, 27))
    check("W_k <= a(k)*k/(k-1)  (i.e. theta(L) <= log N + 1/(L-1))", allok)

    print("\nD. range over which L(N) is independently determined by this work")
    runs = {}
    for fn in os.listdir(os.path.join(ROOT, "results")):
        if not fn.endswith(".json"):
            continue
        try:
            d = json.load(open(os.path.join(ROOT, "results", fn)))
        except Exception:
            continue
        if isinstance(d, dict) and "B" in d and "min_last" in d:
            k = d["k"]
            if k not in runs or d["B"] > runs[k]["B"]:
                runs[k] = d
    for k in sorted(runs):
        r = runs[k]
        check("exhaustive k=%d to B=%d reproduces a(%d)" % (k, r["B"], k),
              r["min_last"] == a[k], "found %d" % r["min_last"])

    # ref_bruteforce.py is assumption-free (no structure theorem, no shared code)
    ref = json.load(open(os.path.join(ROOT, "results", "reference.json")))
    for ks, r in sorted(ref.items(), key=lambda t: int(t[0])):
        k = int(ks)
        check("brute force k=%d to B=%d reproduces a(%d)" % (k, r["B"], k),
              r["min_last"] == a[k], "found %d" % r["min_last"])

    kmax = max(runs)
    enum_total = sum(r["total"] for k, r in runs.items() if k >= 8)
    enum_caseB = sum(r["caseB_count"] for k, r in runs.items() if k >= 8)
    check("no Case B progression among the %d enumerated for k>=8" % enum_total,
          enum_caseB == 0, "caseB total %d" % enum_caseB)
    # A (k+1)-AP with last term <= B contains two k-APs with last term <= B
    # (terms 1..k and terms 2..k+1), so total == 1 rules out any (k+1)-AP <= B.
    top = runs[kmax]
    check("k=%d run has total==1, so no %d-AP exists below B=%d"
          % (kmax, kmax + 1, top["B"]), top["total"] == 1)
    verified_to = top["B"]
    covered = set(runs) | {int(x) for x in ref}
    missing = [k for k in range(3, kmax + 1) if k not in covered]
    check("a(k) reproduced by some run for every k in 3..%d" % kmax,
          not missing, "missing: %s" % missing if missing else "no gaps")
    print("\n  => L(N) independently re-derived for all N <= %d  (%.3g)"
          % (verified_to, verified_to))
    print("     L(N) = max{k : a(k) <= N}; beyond this range a(%d)..a(26) are cited,"
          % (kmax + 1))
    print("     not recomputed here, and L(N) is unknown for N >= a(27).")

    print("\nE. extremal ratio")
    rows = [(k, a[k], k / math.log(a[k])) for k in range(3, 27)]
    best = max(rows, key=lambda r: r[2])
    print("     max k/log a(k) over 3<=k<=26 = %.4f at k=%d" % (best[2], best[0]))
    print("     value at k=26                = %.4f" % rows[-1][2])
    check("ratio at k=26 is below 0.71", rows[-1][2] < 0.71, "%.4f" % rows[-1][2])

    print("\nF. Hardy-Littlewood model vs exhaustively measured counts")
    cpath = os.path.join(ROOT, "results", "count_comparison.json")
    if not os.path.exists(cpath):
        check("count_comparison.json present", False, "missing")
    else:
        comp = json.load(open(cpath))
        # Re-derive the measured column straight from the raw histograms, with a
        # cumulative sum written independently of src/compare_counts.py.  The
        # histogram keys are bucket LOWER edges, so every bucket with lo <= N
        # lies entirely at or below N.
        raw, bad = {}, []
        for k in sorted({int(r["k"]) for r in comp}):
            f = os.path.join(ROOT, "results", "count_k%d.json" % k)
            raw[k] = json.load(open(f))
        for r in comp:
            k, N = int(r["k"]), float(r["N"])
            tot = sum(v for lo, v in raw[k]["hist"].items() if float(lo) <= N)
            if tot != int(r["measured"]):
                bad.append((k, N, tot, r["measured"]))
        check("measured counts re-derived from raw histograms",
              not bad, "%d rows, %d mismatched" % (len(comp), len(bad)))
        # each k's largest tabulated N must recover that run's own total
        tot_bad = [k for k in raw
                   if max((int(r["measured"]) for r in comp if int(r["k"]) == k),
                          default=-1) != raw[k]["total"]]
        check("largest N per k recovers the run total", not tot_bad, str(tot_bad))

        ratios = [r["predicted"] / r["measured"] for r in comp]
        lg = [math.log(x) for x in ratios]
        gm = math.exp(sum(lg) / len(lg))
        gsd = math.exp(math.sqrt(sum((x - math.log(gm)) ** 2 for x in lg) / (len(lg) - 1)))
        ks = sorted({int(r["k"]) for r in comp})
        print("     %d points, k=%d..%d, N=%.3g..%.3g, measured %d..%d"
              % (len(comp), ks[0], ks[-1],
                 min(r["N"] for r in comp), max(r["N"] for r in comp),
                 min(int(r["measured"]) for r in comp),
                 max(int(r["measured"]) for r in comp)))
        print("     geometric mean pred/obs = %.4f, geometric sd = %.4f, range %.4f..%.4f"
              % (gm, gsd, min(ratios), max(ratios)))
        # the claim made in NOTES.md section 4
        check("geometric mean within 3% of 1", abs(gm - 1.0) < 0.03, "%.4f" % gm)
        check("every point within 15% of prediction",
              max(max(ratios), 1 / min(ratios)) < 1.15,
              "worst %.4f" % max(max(ratios), 1 / min(ratios)))
        check("no Case B progression found in any counting run",
              all(int(r["caseB"]) == 0 for r in comp)
              and all(v["caseB_count"] == 0 for v in raw.values()), "")

    section_g()

    print("\n%s" % ("ALL CHECKS PASSED" if not FAIL else "FAILURES: %s" % FAIL))
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
