#!/usr/bin/env python3
"""explicit_n0.py -- the explicit pieces of the conditional threshold N_0(k) of NOTES.md section 7.1.

NOTES.md section 7.1 shows that no published source states the constants needed to turn
Tao-Teraevaeinen's Theorem 1.6 (arXiv:2107.02158v4) into an explicit N_0(k).  It writes the
threshold as

    N_0(k) = max{ N_1(k), N_2(k), exp exp((16 (k-1) C_k / S_k)^(1/c_k)) }

where C_k, c_k and N_1(k) are the unstated constants of that theorem, and computes here the two
pieces that ARE explicit:

  S_k     the singular series prod_p beta_p of TT Example 1.7, as a rigorous enclosure
          [lo, hi]: the product over p <= P exactly, times a tail factor in
          [exp(-(k-1)^2 / (2 (1 - (k-1)/P) P)), 1].  For p >= k,
          log beta_p = log(1-x) - (k-1) log(1-1/p) with x = (k-1)/p lies in [-x^2/(2(1-x)), 0],
          and sum_{p>P} 1/p^2 < sum_{n>P} 1/n^2 < 1/P.
  N_2(k)  the threshold beyond which the tuples with a proper prime power p^j (j >= 2) among
          their terms cost less than S_k/(16(k-1)) of the main term:
          log N_2(k) = max(30 log 2, log 2k, L*), where L* >= 2k solves
          L/2 - k log L = log(32 k (k-1) / S_k_lo).

It also checks numerically, on the range computed, the elementary bound
log(1/S_k) <= 4 k log k (all k >= 3) used in the tower-height lemma of section 7.1.

Nothing here evaluates C_k, c_k or N_1(k): no source states them (NOTES.md section 7.1).
Writes results/explicit_n0.json.  src/final_check.py section H recomputes every value
independently, in floating point without mpmath, and does not import this file.
"""
import json
import os

import mpmath as mp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
P = 10 ** 6
KS = list(range(3, 27))
mp.mp.dps = 30


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


PR = primes_upto(P)


def log_beta(p, k):
    """log beta_p for the k-AP forms x + i y, i = 0..k-1 (TT Example 1.7)."""
    p = mp.mpf(p)
    if p < k:
        return -mp.log(p) - (k - 1) * mp.log1p(-1 / p)
    return mp.log1p(-(k - 1) / p) - (k - 1) * mp.log1p(-1 / p)


def singular_series(k):
    """(lo, hi): a rigorous enclosure of S_k = prod_p beta_p."""
    hi = mp.fsum(log_beta(p, k) for p in PR)
    tail = -mp.mpf((k - 1) ** 2) / (2 * (1 - mp.mpf(k - 1) / P) * P)
    return mp.exp(hi + tail), mp.exp(hi)


def log_n2(k, s_lo):
    """log N_2(k), natural log."""
    rhs = mp.log(32 * k * (k - 1) / s_lo)
    h = lambda L: L / 2 - k * mp.log(L) - rhs
    lo = mp.mpf(2 * k)
    if h(lo) >= 0:
        Lstar = lo
    else:
        hi = lo
        while h(hi) < 0:
            hi *= 2
        Lstar = mp.findroot(h, (lo, hi), solver="bisect")
    return max(30 * mp.log(2), mp.log(2 * k), Lstar)


def main():
    rows = []
    for k in KS:
        lo, hi = singular_series(k)
        ln2 = log_n2(k, lo)
        rows.append({"k": k, "S_lo": float(lo), "S_hi": float(hi), "log_N2": float(ln2),
                     "log_inv_S_over_klogk": float(-mp.log(lo) / (k * mp.log(k)))})
        print("k=%2d  S_k in [%.7f, %.7f]  log N_2(k) = %.4f" % (k, lo, hi, ln2))
    worst = max(r["log_inv_S_over_klogk"] for r in rows)
    print("max over k of log(1/S_k)/(k log k) = %.4f (the analytic bound used in NOTES 7.1 is 4)" % worst)
    out = {"P": P, "ks": KS, "rows": rows, "max_log_inv_S_over_klogk": worst,
           "sources": {"tao_teravainen": "arXiv:2107.02158v4", "manners": "arXiv:1811.00718v2"}}
    with open(os.path.join(ROOT, "results", "explicit_n0.json"), "w") as f:
        json.dump(out, f, indent=1)
        f.write("\n")
    return 0 if worst < 4 else 1


if __name__ == "__main__":
    raise SystemExit(main())
