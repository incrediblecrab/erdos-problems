"""Continuous optimisation of block boundaries for Psi_k, with exact certification.

The grid searches in relax.py pin every boundary to a multiple of 1/m.  That is
the binding constraint near the Parrilo-Robertson-Saracino value: their twelve
block lengths over 548 have gcd 1, so no coarse grid represents them, and the
O(m^2) exact weight table makes a fine grid expensive.

Here the free variables are the block lengths themselves.  Search runs in
floating point (blocks.psi_blocks evaluates the same piecewise-polynomial in
either arithmetic); every reported value is then re-derived in exact rational
arithmetic from integer block lengths.

    python src/optblocks.py 3 --blocks 13-30 --starts 200 --store results/blocks_k3.json

Reference values (Psi_k normalisation, = delta_k):
    k=3  117/2192 = 0.05337591  [PRS08]     random 1/16   = 0.0625
    k=4    1/72   = 0.01388889  [LuPe12]    random 1/48   = 0.02083333
    k=5    1/304  = 0.00328947  [LuPe12]    random 1/128  = 0.0078125
"""

import argparse
import json
import os
import sys
import time
from fractions import Fraction

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import blocks                                         # noqa: E402

PUBLISHED = {3: Fraction(117, 2192), 4: Fraction(1, 72), 5: Fraction(1, 304)}


def random_area(k):
    """Psi_k of a random colouring: 2^(1-k) * area(T)."""
    return Fraction(1, (k - 1) * 2 ** k)


# --------------------------------------------------------------------------
# float objective over an unconstrained parametrisation
# --------------------------------------------------------------------------

def bounds_from_z(z):
    """Softmax: unconstrained z -> positive block lengths summing to 1 -> bounds."""
    w = np.exp(z - z.max())
    u = w / w.sum()
    b = np.concatenate(([0.0], np.cumsum(u)))
    b[-1] = 1.0
    return b


def f_of_z(z, k):
    return float(blocks.psi_blocks(bounds_from_z(z).tolist(), k))


def refine(z0, k, maxiter=400):
    """L-BFGS-B on the softmax parametrisation (scipy's finite differences)."""
    res = minimize(f_of_z, z0, args=(k,), method="L-BFGS-B",
                   options={"maxiter": maxiter, "eps": 1e-7, "maxfun": 100000})
    return res.x, res.fun


# --------------------------------------------------------------------------
# integer snap + exact certification
# --------------------------------------------------------------------------

def snap(bounds, den):
    """Round boundaries to integer block lengths over `den`; drop empty blocks.

    Blocks that round to zero length are deleted.  Deleting one block would
    merge its two neighbours, which have the same colour, so blocks are removed
    in pairs: an empty block and one neighbour.  This keeps the alternation.
    """
    cuts = [int(round(b * den)) for b in bounds]
    cuts[0], cuts[-1] = 0, den
    cuts = sorted(cuts)
    lens = [cuts[i + 1] - cuts[i] for i in range(len(cuts) - 1)]
    while 0 in lens:
        i = lens.index(0)
        j = i + 1 if i + 1 < len(lens) else i - 1
        if len(lens) <= 2:
            break
        merged = lens[i] + lens[j]
        lo, hi = min(i, j), max(i, j)
        lens = lens[:lo] + lens[hi + 1:]
        if lens:
            lens[max(0, lo - 1)] += merged
        else:
            lens = [merged]
    assert sum(lens) == den, (sum(lens), den)
    return lens


def hillclimb(lens, k, rounds=60):
    """Move one unit of length across a boundary while it helps (float arithmetic)."""
    den = sum(lens)
    lens = list(lens)
    best = float(blocks.psi_blocks([c / den for c in np.cumsum([0] + lens)], k))
    for _ in range(rounds):
        improved = False
        for i in range(len(lens)):
            for j in (i - 1, i + 1):
                if j < 0 or j >= len(lens) or lens[i] <= 1:
                    continue
                trial = list(lens)
                trial[i] -= 1
                trial[j] += 1
                v = float(blocks.psi_blocks([c / den for c in np.cumsum([0] + trial)], k))
                if v < best - 1e-15:
                    best, lens, improved = v, trial, True
        if not improved:
            break
    return lens, best


def certify(lens, k):
    """Exact Psi_k from integer block lengths (Fractions throughout)."""
    return blocks.psi_blocks(blocks.bounds_of_lengths(lens), k)


def best_exact(bounds, k, dens=(548, 1000, 2192, 5000, 10000), climb=True):
    """Snap to several denominators, hill-climb, certify exactly; return the best."""
    best = None
    for den in dens:
        lens = snap(bounds, den)
        if len(lens) < 2:
            continue
        if climb:
            lens, _ = hillclimb(lens, k)
        g = 0
        for L in lens:
            g = np.gcd(g, L)
        if g > 1:
            lens = [L // int(g) for L in lens]
        val = certify(lens, k)
        if best is None or val < best[0]:
            best = (val, lens)
    return best


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------

def load(path):
    if path and os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    return {}


def save(path, data):
    if not path:
        return
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(data, fh, indent=1, sort_keys=True)
    os.replace(tmp, path)


def parse_range(s):
    out = []
    for part in s.split(","):
        if "-" in part:
            a, b = part.split("-")
            out.extend(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("k", type=int)
    ap.add_argument("--blocks", default="13-28")
    ap.add_argument("--starts", type=int, default=60)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--maxiter", type=int, default=400)
    ap.add_argument("--store", default="")
    ap.add_argument("--seed-lengths", default="",
                    help="comma-separated block lengths to use as one start")
    args = ap.parse_args()

    k = args.k
    rng = np.random.default_rng(args.seed)
    data = load(args.store)
    pub = PUBLISHED.get(k)
    rnd = random_area(k)

    seeds = []
    if args.seed_lengths:
        seeds.append([int(x) for x in args.seed_lengths.split(",")])

    for B in parse_range(args.blocks):
        t0 = time.time()
        best = None
        starts = []
        # perturbations of any supplied seed, padded/truncated to B blocks
        for s in seeds:
            if len(s) <= B:
                pad = s + [max(1, sum(s) // (10 * B))] * (B - len(s))
            else:
                pad = s[:B]
            z = np.log(np.array(pad, dtype=float) / sum(pad))
            starts.append(z)
            for _ in range(min(10, args.starts // 4)):
                starts.append(z + rng.normal(0, 0.35, B))
        while len(starts) < args.starts:
            starts.append(rng.normal(0, 1.0, B))

        for z0 in starts:
            z, v = refine(z0, k, args.maxiter)
            if best is None or v < best[0]:
                best = (v, z)

        got = best_exact(bounds_from_z(best[1]).tolist(), k)
        if got is None:
            continue
        val, lens = got
        key = str(B)
        old = data.get(key)
        if old is None or Fraction(*map(int, old["value"].split("/"))) > val:
            data[key] = {"blocks": len(lens), "lengths": lens,
                         "den": sum(lens), "value": "%d/%d" % (val.numerator, val.denominator),
                         "float": float(val)}
            save(args.store, data)
            mark = "NEW"
        else:
            mark = "   "
        rel_p = float(val / pub) if pub else float("nan")
        print("%s B=%3d -> %3d blocks/%-6d  %-22s = %.10f  %.4fx random  %.4fx published  [%.0fs]"
              % (mark, B, len(lens), sum(lens),
                 "%d/%d" % (val.numerator, val.denominator), float(val),
                 float(val / rnd), rel_p, time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
