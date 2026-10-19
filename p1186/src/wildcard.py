"""The wildcard construction: a Z_44 base with free cosets, for k = 5.

Lu and Peng [LuPe12, section 4] found their k = 4 bound by noticing that the
two halves of their good colouring B22 of Z_22 differ in exactly one bit.  That
gives B11 = (1,1,1,0,1,*,0,1,0,0,0), a colouring of Z_11 with one wildcard, and
a recursion B11 |x| B_t on Z_{11t} whose limit is m_4(Z_n) -> 1/12.

The same thing happens one level up for k = 5, with a different base.  The
optimal periodic colourings found at m = 88 and m = 132 both reduce to the SAME
44-bit pattern with exactly four free positions,

    100?0111010001?1101111011?1000101110?0010000
       ^         ^          ^          ^
       3        14         25         36          (an AP of difference 11)

so the base is Z_44 with a wildcard coset structure, and the measured optima obey

    m_5(Z_{44t}) = (10t+1)/(484t),

which is what the recursion  m_5(Z_{44t}) <= (40 + 4*m_5(Z_t))/44^2  gives when
the inner colouring is 5-AP-free.  This script tests that law: it fixes the 40
determined bits and searches only the 4t free ones, then certifies the result in
exact rational arithmetic.

    python src/wildcard.py --t 2,3,4,5,6 --tabu 20000 --restarts 12
"""

import argparse
import glob
import json
import os
import sys
import time
from fractions import Fraction

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exact                                          # noqa: E402
from relax import Objective                           # noqa: E402

BASE = "10010111010001011011110110100010111010010000"
FREE = (3, 14, 25, 36)


def inner_min(ft):
    """Best known m_5(Z_ft), the inner term of the recursion.

    Exhaustive for small ft (where it is the true minimum), otherwise the best
    stored periodic optimum for that modulus.  Returns None if neither is
    available, so callers can tell "unknown" from "large".
    """
    import itertools
    if ft <= 18:
        best = None
        for bits in itertools.product((0, 1), repeat=ft):
            phi = exact.phi_zm(list(bits), 5)
            if best is None or phi < best:
                best = phi
        return best
    best = None
    for fn in glob.glob(os.path.join("results", "zm_k5*.json")):
        try:
            data = json.load(open(fn))
        except Exception:
            continue
        recs = list(data.values()) if isinstance(data, dict) else data
        for rec in recs:
            if not isinstance(rec, dict) or rec.get("m") != ft:
                continue
            w = rec.get("word")
            if w and len(w) == ft:
                phi = exact.phi_zm([int(ch) for ch in w], 5)
                if best is None or phi < best:
                    best = phi
    return best


def predicted(t):
    """The wildcard recursion for the Z_44 base with a free coset of size 4:

        m_5(Z_{44t}) <= ((44-4) + 4^2 * m_5(Z_{4t})) / 44^2
                      = (5 + 2 * m_5(Z_{4t})) / 242,

    whose fixed point is 1/48.  This supersedes the naive law (10t+1)/(484t),
    which assumes the inner colouring is 5-AP-free and is therefore wrong
    whenever it is not -- e.g. t=4 (m_5(Z_16)=5/64) and t=6 (m_5(Z_24)=5/72).
    Returns None when the inner value is unknown.
    """
    x = inner_min(4 * t)
    if x is None:
        return None
    return (Fraction(44 - 4) + 16 * x) / (44 * 44)


def derive_free_positions():
    """Re-derive BASE and FREE from the stored optima instead of trusting them.

    Returns (base, free) or None if there is not enough stored data.
    """
    import glob
    best = {}
    for fn in glob.glob(os.path.join("results", "zm_k5*.json")):
        try:
            data = json.load(open(fn))
        except Exception:
            continue
        for v in data.values():
            m = int(v["m"])
            p = Fraction(v["phi"])
            if m % 44 == 0 and m > 44 and (m not in best or p < best[m][0]):
                best[m] = (p, v["word"])
    if not best:
        return None
    free = set()
    for m, (_, w) in best.items():
        t = m // 44
        blocks = [w[i * 44:(i + 1) * 44] for i in range(t)]
        free |= {j for j in range(44) if len({b[j] for b in blocks}) > 1}
    fixed = {}
    for m, (_, w) in best.items():
        t = m // 44
        for j in range(44):
            if j in free:
                continue
            col = w[j]
            if fixed.setdefault(j, col) != col:
                return None                          # the bases disagree
    return fixed, tuple(sorted(free))


def restricted_tabu(obj, free_idx, base_bits, iters, rng, tabu_len):
    """Tabu search that is only allowed to flip the free positions."""
    m = obj.m
    v = np.array([1 - 2 * b for b in base_bits], dtype=np.float64)
    mask = np.zeros(m, dtype=bool)
    mask[free_idx] = True
    cur = obj.value(v)
    best_v, best_val = v.copy(), cur
    until = np.zeros(m, dtype=np.int64)
    for it in range(iters):
        gains = obj.flip_gains(v)
        allowed = mask & (until <= it)
        allowed |= mask & ((cur + gains) < best_val - 1e-15)
        if not allowed.any():
            allowed = mask.copy()
        idx = np.flatnonzero(allowed)
        j = idx[np.argmin(gains[idx])]
        cur += gains[j]
        v[j] = -v[j]
        until[j] = it + tabu_len
        if cur < best_val - 1e-15:
            best_val, best_v = cur, v.copy()
    return best_v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--t", default="2,3,4,5,6")
    ap.add_argument("--tabu", type=int, default=20000)
    ap.add_argument("--restarts", type=int, default=12)
    ap.add_argument("--seed", type=int, default=20260202)
    ap.add_argument("--store", default="results/wildcard_k5.json")
    args = ap.parse_args()

    k = 5
    derived = derive_free_positions()
    if derived is None:
        base_map, free = {j: BASE[j] for j in range(44)}, FREE
        print("could not re-derive from results/; using the hard-coded base")
    else:
        base_map, free = derived
        print("re-derived from the stored optima:")
    base_str = "".join(base_map.get(j, "?") for j in range(44))
    print("  base  %s" % base_str)
    print("  free  %s   (differences %s)"
          % (list(free), [free[i + 1] - free[i] for i in range(len(free) - 1)]))
    if derived is not None:
        agree = (tuple(free) == FREE and
                 all(base_map[j] == BASE[j] for j in base_map))
        print("  agrees with the hard-coded base and free set: %s" % agree)

    store = {}
    if args.store and os.path.exists(args.store):
        store = json.load(open(args.store))

    rng = np.random.default_rng(args.seed)
    ok = True
    for t in [int(x) for x in args.t.split(",")]:
        m = 44 * t
        t0 = time.time()
        obj = Objective(m, k, setting="zm")
        free_idx = [i * 44 + j for i in range(t) for j in free]
        base_bits = [int(base_map.get(j, "0")) for i in range(t) for j in range(44)]
        best_phi, best_c = None, None
        for r in range(args.restarts):
            bb = list(base_bits)
            for i in free_idx:
                bb[i] = int(rng.integers(2))
            v = restricted_tabu(obj, free_idx, bb, args.tabu, rng, max(3, len(free_idx) // 3))
            c = [(1 if x < 0 else 0) for x in v]
            for j in range(m):                       # the fixed bits must be untouched
                if (j % 44) not in free:
                    assert c[j] == int(base_map[j % 44]), "restricted search moved a fixed bit"
            phi = exact.phi_zm(c, k)                 # exact certification
            if best_phi is None or phi < best_phi:
                best_phi, best_c = phi, c
        pred = predicted(t)
        verdict = ("MATCHES the law" if best_phi == pred
                   else ("BELOW the law" if best_phi < pred else "above the law"))
        print("  t=%-3d m=%-5d  found %-14s = %.10f   law %-14s = %.10f   %s   [%.0fs]"
              % (t, m, best_phi, float(best_phi), pred, float(pred), verdict, time.time() - t0))
        if best_phi > pred:
            ok = False
        store[str(m)] = {
            "k": k, "m": m, "t": t,
            "phi": str(best_phi), "psi": str(best_phi / (2 * (k - 1))),
            "psi_float": float(best_phi / (2 * (k - 1))),
            "predicted": str(pred),
            "matches_law": best_phi == pred,
            "ap_free": best_phi == Fraction(1, m),
            "word": "".join(str(b) for b in best_c),
            "family": "wildcard",
        }
        if args.store:
            json.dump(store, open(args.store, "w"), indent=1, sort_keys=True)
    print("\nlaw reproduced at every t tested: %s" % ok)
    return 0


if __name__ == "__main__":
    sys.exit(main())
