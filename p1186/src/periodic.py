"""Search the periodic family: colour i in {1,...,n} by c[i mod m].

The bound this family certifies is

    delta_k <= Psi_k = Phi^{Z_m}_k(c) / (2(k-1)),
    Phi^{Z_m}_k(c) = (1/m^2) #{(a,b) in Z_m^2 : c constant on a, a+b, ..., a+(k-1)b}.

Every b = 0 pair is monochromatic, so Phi^{Z_m}_k >= 1/m, with equality exactly
when c is a 2-colouring of Z_m with no monochromatic k-AP of nonzero common
difference.  Such a colouring therefore certifies

    delta_k <= 1 / (2(k-1)m),

and the whole problem for this family is to make m as large as possible.

This is a different family from the block colourings searched by relax.py
--setting n, and for k >= 4 it is much stronger: the block family bottoms out
near 0.0179 for k=4 while Z_11 gives 1/66 = 0.01515.

    python src/periodic.py 5 2-80 --store results/zm_k5.json

Reference values (delta_k normalisation):
    k=3  117/2192 = 0.05337591  [PRS08]    random 1/16  = 0.0625
    k=4    1/72   = 0.01388889  [LuPe12]   random 1/48  = 0.02083333
    k=5    1/304  = 0.00328947  [LuPe12]   random 1/128 = 0.0078125
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


def search_m(m, k, restarts, tabu_iters, rng):
    """Best vertex found for Phi^{Z_m}_k; returns (float value, word)."""
    obj = Objective(m, k, "zm")
    floor = 1.0 / m                     # b = 0 pairs, unavoidable
    best = None
    for _ in range(restarts):
        v = np.where(rng.random(m) < 0.5, -1.0, 1.0)
        v, _ = obj.descend(v)
        v, val = obj.tabu(v, tabu_iters, rng)
        if best is None or val < best[0]:
            best = (val, v.copy())
        if best[0] <= floor + 1e-12:     # AP-free: cannot do better
            break
    word = "".join("1" if x > 0 else "0" for x in best[1])
    return best[0], word


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("k", type=int)
    ap.add_argument("ms")
    ap.add_argument("--restarts", type=int, default=40)
    ap.add_argument("--tabu", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--store", default="")
    args = ap.parse_args()

    k = args.k
    rng = np.random.default_rng(args.seed)
    data = load(args.store)
    pub = PUBLISHED.get(k)
    rnd = Fraction(1, (k - 1) * 2 ** k)          # random colouring
    best_overall = None

    for m in parse_range(args.ms):
        t0 = time.time()
        _, word = search_m(m, k, args.restarts, args.tabu, rng)
        c = [int(x) for x in word]
        phi = exact.phi_zm(c, k)                  # exact certification
        psi = phi / (2 * (k - 1))
        apfree = phi == Fraction(1, m)
        key = str(m)
        old = data.get(key)
        oldv = Fraction(*map(int, old["phi"].split("/"))) if old else None
        if oldv is None or phi < oldv:
            data[key] = {"m": m, "k": k, "word": word,
                         "phi": "%d/%d" % (phi.numerator, phi.denominator),
                         "psi": "%d/%d" % (psi.numerator, psi.denominator),
                         "psi_float": float(psi), "ap_free": apfree}
            save(args.store, data)
            mark = "NEW"
        else:
            mark = "   "
        if best_overall is None or psi < best_overall[0]:
            best_overall = (psi, m)
        note = " AP-FREE" if apfree else ""
        rel = (" %.4fx published" % float(psi / pub)) if pub else ""
        print("%s m=%3d Phi=%-14s Psi=%-14s %.9f  %.4fx random%s%s  [%.0fs]"
              % (mark, m, "%d/%d" % (phi.numerator, phi.denominator),
                 "%d/%d" % (psi.numerator, psi.denominator), float(psi),
                 float(psi / rnd), rel, note, time.time() - t0), flush=True)

    if best_overall:
        print("best over the sweep: Psi = %s at m = %d" % best_overall, flush=True)


if __name__ == "__main__":
    main()
