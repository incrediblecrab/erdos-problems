"""Continuous breakpoint search for F_p colourings, with exact certification.

The block searches in relax.py fix every colour change of the circle R/Z to a
multiple of 1/m.  This script frees the breakpoints: it starts from stored
block words, minimises circle.phi_circle over the arc lengths in floating
point (L-BFGS-B and Powell, then random perturbations of the incumbent), rounds the result
onto grids of M equal blocks for every M in a range, polishes the integer
cuts of the best few grids by single +-1 moves,
and re-derives the value of the resulting word in exact rational arithmetic
with exact.phi.  Only that exact value is recorded; the float optimum is a
proposal.

    python src/circleopt.py 4 --mmax 600 --store results/zp_k4_c.json

The stored record has the same shape as relax.py's (word, value, bound), so
src/final_check.py re-certifies it like any other zp_* word.
"""

import argparse
import json
import os
import sys
import time
from fractions import Fraction
from glob import glob

import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import circle  # noqa: E402
import exact   # noqa: E402


def phi_of_lengths(L, k):
    L = np.abs(np.asarray(L, dtype=float))
    L = L / L.sum()
    b = np.concatenate(([0.0], np.cumsum(L)))
    b[-1] = 1.0
    return circle.phi_circle(b.tolist(), k)


def optimise(L0, k, rng, hops, scale, maxfev=20000):
    """L-BFGS-B (finite differences) then a bounded Powell polish from L0, then
    `hops` perturb-and-reoptimise steps (basin hopping)."""
    def run(L):
        r = minimize(phi_of_lengths, L, args=(k,), method="L-BFGS-B",
                     options={"eps": 1e-9, "ftol": 1e-15, "gtol": 1e-12, "maxiter": 5000})
        r = minimize(phi_of_lengths, np.abs(r.x) / np.abs(r.x).sum(), args=(k,),
                     method="Powell",
                     options={"xtol": 1e-8, "ftol": 1e-13, "maxfev": maxfev})
        L = np.abs(r.x) / np.abs(r.x).sum()
        return L, phi_of_lengths(L, k)
    best_L, best_v = run(L0)
    for _ in range(hops):
        trial = best_L * np.exp(rng.normal(0.0, scale, size=best_L.size))
        L, v = run(trial / trial.sum())
        if v < best_v - 1e-13:
            best_L, best_v = L, v
    return best_L, best_v


def cuts_on_grid(L, M):
    cum = np.concatenate(([0.0], np.cumsum(L)))
    cuts = [int(round(c * M)) for c in cum]
    cuts[0], cuts[-1] = 0, M
    return cuts


def phi_of_cuts(cuts, k):
    M = cuts[-1]
    return circle.phi_circle([c / M for c in cuts], k)


def polish(cuts, k):
    """Greedy single-cut +-1 moves on the M-grid (float); keeps cuts strictly increasing."""
    cuts = list(cuts)
    val = phi_of_cuts(cuts, k)
    improved = True
    while improved:
        improved = False
        for i in range(1, len(cuts) - 1):
            for step in (-1, 1):
                c = cuts[i] + step
                if not (cuts[i - 1] < c < cuts[i + 1]):
                    continue
                trial = cuts[:i] + [c] + cuts[i + 1:]
                v = phi_of_cuts(trial, k)
                if v < val - 1e-15:
                    cuts, val, improved = trial, v, True
    return cuts, val


def stored_words(k, top):
    """The `top` best stored zp words for k, over all zp_k{k}*.json shards."""
    recs = []
    for path in glob(os.path.join(ROOT, "results", f"zp_k{k}*.json")):
        for key, rec in json.load(open(path)).items():
            if "word" in rec and "value" in rec:
                recs.append((Fraction(rec["value"]), key, rec["word"]))
    recs.sort()
    seen, out = set(), []
    for v, key, w in recs:
        if w in seen:
            continue
        seen.add(w)
        out.append((v, key, w))
    return out[:top]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("k", type=int)
    ap.add_argument("--top", type=int, default=5, help="number of stored words to start from")
    ap.add_argument("--hops", type=int, default=20, help="basin-hopping steps per start")
    ap.add_argument("--scale", type=float, default=0.3, help="log-normal perturbation size")
    ap.add_argument("--mmin", type=int, default=100)
    ap.add_argument("--mmax", type=int, default=600)
    ap.add_argument("--polish", type=int, default=15, help="polish the best N rounded grids")
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--store", default=None)
    args = ap.parse_args()
    k = args.k
    rng = np.random.default_rng(args.seed)

    starts = stored_words(k, args.top)
    best = None
    for v0, key, w in starts:
        t = time.time()
        bounds, _ = circle.bounds_of_cyclic_word(w)
        L0 = np.diff([float(b) for b in bounds])
        L, v = optimise(L0, k, rng, args.hops, args.scale)
        print(f"start m={key:>4} Phi={float(v0):.10f} -> continuous {v:.10f}"
              f"  arcs={L.size}  {time.time() - t:.0f}s", flush=True)
        if best is None or v < best[1]:
            best = (L, v)
    L, vcont = best
    print(f"best continuous Phi = {vcont:.12f}; arc lengths:", np.round(L, 6).tolist(), flush=True)

    rounded = []
    for M in range(args.mmin, args.mmax + 1):
        cuts = cuts_on_grid(L, M)
        if any(cuts[i + 1] <= cuts[i] for i in range(len(cuts) - 1)):
            continue
        rounded.append((phi_of_cuts(cuts, k), M, cuts))
    rounded.sort()
    grid = None
    for v_round, M, cuts in rounded[:args.polish]:
        cuts, v = polish(cuts, k)
        print(f"  M={M:4d} rounded Phi={v_round:.12f} polished {v:.12f}", flush=True)
        if grid is None or v < grid[1] - 1e-15:
            grid = (cuts, v)
    cuts, vgrid = grid
    word = circle.word_of_lengths([cuts[i + 1] - cuts[i] for i in range(len(cuts) - 1)])
    M = len(word)
    val = exact.phi([int(ch) for ch in word], k)          # the certificate
    assert abs(float(val) - vgrid) < 1e-9, (float(val), vgrid)
    assert circle.phi_circle([Fraction(c, M) for c in cuts], k) == val
    print(f"certified m={M}: Phi = {val} = {float(val):.12f}, delta~_{k} <= {val / 2}", flush=True)

    best_stored = starts[0][0]
    print(f"best stored block word: Phi = {best_stored} = {float(best_stored):.12f};"
          f" {'IMPROVES' if val < best_stored else 'does not improve'}", flush=True)
    if args.store and val < best_stored:
        db = json.load(open(args.store)) if os.path.exists(args.store) else {}
        db[str(M)] = {
            "word": word,
            "value": str(val),
            "value_float": float(val),
            "bound": str(val / 2),
            "bound_float": float(val / 2),
            "found_by": "circleopt",
            "exhaustive": False,
        }
        with open(args.store, "w") as fh:
            json.dump(db, fh, indent=1, sort_keys=True)
        print(f"stored in {args.store}")


if __name__ == "__main__":
    main()
