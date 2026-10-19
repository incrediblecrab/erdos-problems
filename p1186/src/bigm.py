"""Seeded local search in the periodic family at large moduli.

Random restarts are hopeless once m is in the hundreds, but the periodic family
has a free lunch: tiling a colouring of Z_m to Z_{tm} preserves Phi exactly,
because the colour pattern of a progression depends only on (a mod m, b mod m)
and that pair is uniform on Z_m^2 either way.  So the tiled copy of the best
known Z_m word is a starting point at Z_{tm} whose value is already the record,
and tabu search can only improve on it.

    python src/bigm.py 5 --base 44 --mult 2,3,4 --tabu 40000 --restarts 8

Writes to results/zm_k<k>_big.json in the same schema as periodic.py.
"""

import argparse
import json
import os
import sys
import time
from fractions import Fraction

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exact                                          # noqa: E402
from relax import Objective                           # noqa: E402


def best_word(k, m):
    """Best stored periodic word at modulus m, merging every shard."""
    best = None
    import glob
    for fn in glob.glob(os.path.join("results", "zm_k%d*.json" % k)):
        try:
            with open(fn) as fh:
                data = json.load(fh)
        except Exception:
            continue
        for v in data.values():
            if int(v["m"]) != m:
                continue
            phi = Fraction(v["phi"])
            if best is None or phi < best[0]:
                best = (phi, v["word"])
    return best


def tabu(obj, v0, iters, rng, tabu_len):
    """Plain tabu over sign vectors; returns the best vector seen."""
    v = v0.copy()
    cur = obj.value(v)
    best_v, best_val = v.copy(), cur
    until = np.zeros(obj.m, dtype=np.int64)
    for it in range(iters):
        gains = obj.flip_gains(v)
        allowed = until <= it
        # aspiration: a move that beats the incumbent is always allowed
        allowed |= (cur + gains) < best_val - 1e-15
        if not allowed.any():
            allowed[:] = True
        idx = np.flatnonzero(allowed)
        j = idx[np.argmin(gains[idx])]
        cur += gains[j]
        v[j] = -v[j]
        until[j] = it + tabu_len
        if cur < best_val - 1e-15:
            best_val, best_v = cur, v.copy()
    return best_v, best_val


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("k", type=int)
    ap.add_argument("--base", type=int, required=True,
                    help="modulus whose best stored word is tiled as the seed")
    ap.add_argument("--mult", default="2,3,4",
                    help="comma separated tiling multipliers")
    ap.add_argument("--tabu", type=int, default=40000)
    ap.add_argument("--restarts", type=int, default=6)
    ap.add_argument("--seed", type=int, default=20260101)
    ap.add_argument("--store", default=None)
    args = ap.parse_args()

    k = args.k
    seed = best_word(k, args.base)
    if seed is None:
        print("no stored word at m=%d for k=%d" % (args.base, k))
        return 1
    seed_phi, seed_str = seed
    print("seed  m=%d  Phi=%s = %.12f" % (args.base, seed_phi, float(seed_phi)))
    base = [int(ch) for ch in seed_str]
    assert exact.phi_zm(base, k) == seed_phi, "stored seed does not re-certify"
    print("      seed re-certified exactly")

    store = {}
    if args.store and os.path.exists(args.store):
        with open(args.store) as fh:
            store = json.load(fh)

    rng = np.random.default_rng(args.seed)
    for t in [int(x) for x in args.mult.split(",")]:
        m = args.base * t
        t0 = time.time()
        obj = Objective(m, k, setting="zm")
        tiled = base * t
        phi_tiled = exact.phi_zm(tiled, k)
        print("\nm=%d  tiled seed Phi=%s = %.12f  (tiling must preserve Phi: %s)"
              % (m, phi_tiled, float(phi_tiled),
                 "yes" if phi_tiled == seed_phi else "NO -- tiling argument is wrong"))

        best_c, best_phi = tiled, phi_tiled
        starts = [np.array([1 - 2 * b for b in tiled], dtype=np.float64)]
        for _ in range(args.restarts - 1):
            v = np.array([1 - 2 * b for b in tiled], dtype=np.float64)
            flip = rng.random(m) < 0.15               # perturb 15% of the seed
            v[flip] *= -1
            starts.append(v)
        for si, v0 in enumerate(starts):
            v, _ = tabu(obj, v0, args.tabu, rng, max(3, m // 8))
            c = [(1 if x < 0 else 0) for x in v]
            phi = exact.phi_zm(c, k)                  # exact certification
            if phi < best_phi:
                best_phi, best_c = phi, c
                print("    start %d: Phi=%s = %.12f  *" % (si, phi, float(phi)))
        dt = time.time() - t0
        psi = best_phi / (2 * (k - 1))
        print("  m=%-5d BEST Phi=%-16s = %.12f   delta_%d <= %s = %.12f   [%.0fs]"
              % (m, best_phi, float(best_phi), k, psi, float(psi), dt))
        if best_phi < seed_phi:
            print("    IMPROVES on the m=%d seed" % args.base)
        store["%d" % m] = {
            "k": k, "m": m,
            "phi": str(best_phi), "psi": str(psi),
            "psi_float": float(psi),
            "ap_free": best_phi == Fraction(1, m),
            "word": "".join(str(b) for b in best_c),
            "seeded_from": args.base,
        }
        if args.store:
            with open(args.store, "w") as fh:
                json.dump(store, fh, indent=1, sort_keys=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
