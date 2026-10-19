#!/usr/bin/env python3
"""verify_records.py -- independent verification of the AP-record data.

Builds the table of a(k) = least possible last term of a k-term arithmetic
progression of primes from OEIS A005115 / A113827 / A093364 (parsed from the
raw OEIS text in refs/, not hard-coded), then checks it three ways:

 1. arithmetic:  first + (k-1)*diff == a(k)
 2. primality:   every one of the k terms is prime, tested twice with
                 independent implementations (sympy.isprime, and a
                 deterministic Miller-Rabin with the first 12 prime bases,
                 which is proved correct for n < 3.317e24)
 3. cross-source: the progression agrees with Norman Luhn's independently
                 maintained table (pzktupel.de), which states each record in
                 the unrelated form  (x + y*n) * p# + c.

Writes data/records.json.
"""
import os
import json
import re
import sys

from sympy import isprime, primerange

REFS = "refs/"


def oeis_terms(seqid):
    txt = open(REFS + "oeis_%s.txt" % seqid).read()
    body = "".join(re.findall(r"^%[STU] " + seqid + r" (.*)$", txt, re.M))
    return [int(x) for x in body.replace(" ", "").strip(",").split(",") if x]


def mr_is_prime(n):
    """Deterministic Miller-Rabin; the first 12 prime bases decide n < 3.317e24."""
    if n < 2:
        return False
    bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    for p in bases:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in bases:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def primorial(n):
    r = 1
    for p in primerange(2, n + 1):
        r *= p
    return r


def parse_luhn():
    """Parse the 'AP-k with minimal end' table:  (x + y*n) * p# + c  ->  (a, d)."""
    import html
    s = open(REFS + "luhn.html", encoding="utf-8", errors="replace").read()
    s = re.sub(r"<script.*?</script>|<style.*?</style>", "", s, flags=re.S)
    t = html.unescape(re.sub(r"<[^>]+>", " ", s))
    t = re.sub(r"\s+", " ", t)
    # anchor on the table header row, not the table-of-contents link of the same name
    i = t.find("with minimal end k Primes Ends with")
    if i < 0:
        raise RuntimeError("Luhn minimal-end table not found")
    j = t.find("Record history for smallest AP", i)
    seg = t[i:j if j > 0 else len(t)]
    out = {}
    # forms:  "k (x + y * n) * p# + c  end" / "k m * n * p# + c  end" / "k n * p# + c  end"
    pat = re.compile(
        r"(?<![\d.])(\d{1,2}) \(? *(?:(\d+) *\+ *)?(?:(\d+) *\u2022 *)?n *\)? *"
        r"\u2022 *(\d+)# *\+ *(\d+) (\d+) ")
    for m in pat.finditer(seg):
        k = int(m.group(1))
        x = int(m.group(2) or 0)
        y = int(m.group(3) or 1)
        prim = primorial(int(m.group(4)))
        c = int(m.group(5))
        end = int(m.group(6))
        out[k] = (x * prim + c, y * prim, end)
    if len(out) < 20:
        raise RuntimeError("Luhn parse matched only %d rows; regex is stale" % len(out))
    return out


def main():
    a = oeis_terms("A005115")
    first = oeis_terms("A113827")
    diff = oeis_terms("A093364")
    assert len(a) == len(first) == len(diff), (len(a), len(first), len(diff))
    luhn = parse_luhn()

    recs, fails = [], 0
    for i, ak in enumerate(a):
        k = i + 1
        f, d = first[i], diff[i]
        chk = {}
        chk["arith"] = (f + (k - 1) * d == ak)
        terms = [f + j * d for j in range(k)]
        chk["sympy_prime"] = all(isprime(t) for t in terms)
        chk["millerrabin"] = all(mr_is_prime(t) for t in terms)
        chk["distinct"] = (len(set(terms)) == k)
        if k in luhn:
            chk["luhn_agrees"] = (luhn[k] == (f, d, ak))
        else:
            chk["luhn_agrees"] = None
        ok = all(v for v in chk.values() if v is not None)
        fails += not ok
        recs.append({"k": k, "a_k": ak, "first": f, "diff": d,
                     "checks": chk, "ok": ok})
        flag = "OK " if ok else "FAIL"
        lu = {True: "yes", False: "NO", None: "-"}[chk["luhn_agrees"]]
        print(f"k={k:2d} a(k)={ak:>18d} first={f:>16d} d={d:>15d} "
              f"luhn={lu:>3} {flag}")

    print(f"\nrecords={len(recs)} failures={fails}")
    json.dump(recs, open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "records.json"), "w"), indent=1)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
