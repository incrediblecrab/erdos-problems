"""Search the two-scale family on {1,...,n}: colour i by C[i mod m][floor(B i/n)].

This family contains both constructions that hold the published records:

    m = 1        the block colourings of Parrilo-Robertson-Saracino (delta_3)
    B = 1        the periodic colourings searched by periodic.py

so seeding the search with the best member of each parent family guarantees the
result is never worse than either.  See exact.two_scale_weights_int for the
derivation; the evaluator is checked against direct counting over {1,...,n}.

    python src/twoscale.py 4 --fine 1-12 --coarse 2-40 --store results/two_k4.json

Reference values (delta_k normalisation):
    k=3  117/2192 = 0.05337591  [PRS08]    random 1/16  = 0.0625
    k=4    1/72   = 0.01388889  [LuPe12]   random 1/48  = 0.02083333
    k=5    1/304  = 0.00328947  [LuPe12]   random 1/128 = 0.0078125
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

PUBLISHED = {3: Fraction(117, 2192), 4: Fraction(1, 72), 5: Fraction(1, 304)}


def parse_range(s):
    out = []
    for part in s.split(","):
        if "-" in part:
            a, b = part.split("-")
            out.extend(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return out


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


def best_word(path, key_field="word"):
    """Best stored word from a results file keyed by modulus."""
    data = load(path)
    best = None
    for v in data.values():
        w = v.get(key_field)
        if not w:
            continue
        val = v.get("psi_float", v.get("float"))
        if val is None:
            continue
        if best is None or val < best[0]:
            best = (val, w)
    return best[1] if best else None


def best_from_shards(pattern, length):
    """Best word of the given length across all result shards matching pattern.

    Results are sharded across several json files written by concurrent
    sweeps, and a sweep rewrites its whole file on every save, so a single
    file is not guaranteed to hold the best word for a given modulus.
    """
    best = None
    for path in sorted(glob.glob(pattern)):
        for v in load(path).values():
            w = v.get("word")
            if not w or len(w) != length:
                continue
            val = v.get("psi_float", v.get("float"))
            if val is not None and (best is None or val < best[0]):
                best = (val, w)
    return best[1] if best else None


def seeds_for(m, B, k, rng, nrandom, per_override=None):
    """Seed vectors in {-1,+1}^(m*B): parent families first, then random."""
    out = []

    def flat(fn):
        return np.array([1.0 if fn(r, q) else -1.0
                         for r in range(m) for q in range(B)])

    blk = best_from_shards("results/n_k%d*.json" % k, B)
    per = per_override or best_from_shards("results/zm_k%d*.json" % k, m)
    if m == 1:
        per = "0"

    if blk:
        out.append(flat(lambda r, q: blk[q] == "1"))
    if per:
        out.append(flat(lambda r, q: per[r] == "1"))
    if blk and per:
        out.append(flat(lambda r, q: (blk[q] == "1") ^ (per[r] == "1")))
        # independent coarse colouring per residue, started from the block word
        out.append(flat(lambda r, q: (blk[q] == "1") ^ (per[r] == "1") ^ (r % 2 == 1)))
    base = list(out)
    for v in base:
        for frac in (0.03, 0.10):
            w = v.copy()
            w[rng.random(w.size) < frac] *= -1
            out.append(w)
    while len(out) < nrandom + len(base) * 3:
        out.append(np.where(rng.random(m * B) < 0.5, -1.0, 1.0))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("k", type=int)
    ap.add_argument("--fine", default="1-12")
    ap.add_argument("--coarse", default="2-30")
    ap.add_argument("--restarts", type=int, default=12)
    ap.add_argument("--tabu", type=int, default=8000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--maxvars", type=int, default=900)
    ap.add_argument("--store", default="")
    ap.add_argument("--periodic", default="",
                    help="explicit fine-scale seed word of length m")
    args = ap.parse_args()

    k = args.k
    rng = np.random.default_rng(args.seed)
    data = load(args.store)
    pub = PUBLISHED.get(k)
    rnd = Fraction(1, (k - 1) * 2 ** k)
    overall = None

    for m in parse_range(args.fine):
        for B in parse_range(args.coarse):
            if m * B > args.maxvars:
                continue
            t0 = time.time()
            try:
                obj = Objective((m, B), k, "two")
            except (MemoryError, AssertionError) as exc:
                print("   m=%2d B=%3d skipped: %s" % (m, B, exc), flush=True)
                continue
            best = None
            for v in seeds_for(m, B, k, rng, args.restarts,
                               args.periodic or None):
                v, _ = obj.descend(v)
                v, val = obj.tabu(v, args.tabu, rng)
                if best is None or val < best[0]:
                    best = (val, v.copy())
            C = [1 if x > 0 else 0 for x in best[1]]
            psi = exact.psi_two_scale(C, m, B, k)      # exact certification
            assert abs(float(psi) - best[0]) < 1e-9, (psi, best[0])
            key = "%d_%d" % (m, B)
            old = data.get(key)
            oldv = Fraction(*map(int, old["psi"].split("/"))) if old else None
            if oldv is None or psi < oldv:
                data[key] = {"m": m, "B": B, "k": k,
                             "colour": "".join(map(str, C)),
                             "psi": "%d/%d" % (psi.numerator, psi.denominator),
                             "psi_float": float(psi)}
                save(args.store, data)
                mark = "NEW"
            else:
                mark = "   "
            if overall is None or psi < overall[0]:
                overall = (psi, m, B)
            rel = (" %.4fx published" % float(psi / pub)) if pub else ""
            print("%s m=%2d B=%3d vars=%4d terms=%8d  %-18s %.9f  %.4fx random%s  [%.0fs]"
                  % (mark, m, B, m * B, obj.nterms,
                     "%d/%d" % (psi.numerator, psi.denominator), float(psi),
                     float(psi / rnd), rel, time.time() - t0), flush=True)

    if overall:
        print("best: Psi = %s at m=%d B=%d" % overall, flush=True)


if __name__ == "__main__":
    main()
