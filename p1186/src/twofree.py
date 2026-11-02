"""Free coarse breakpoints for colourings of {1,...,n}: a search for delta_3, delta_4.

The F_p bounds of NOTES.md section 11 improved when the colour changes of the
circle were freed from the equal-block grid.  This asks whether the same idea
helps on {1,...,n}, in two families:

  m = 1   block colourings x -> colour of the block of x/n containing it, with
          free block boundaries (the family of src/blocks.py and
          src/optblocks.py, in which [PRS08]'s 117/2192 lives);
  m > 1   two-scale colourings, i -> C[i mod m][coarse block of i/n], with the
          coarse boundaries free instead of equal (src/twoscale.py fixes them
          on an equal grid; its best k = 4 value is 2917/209088 at m = 66,
          B = 4).

Each start alternates a continuous optimisation of the boundaries (L-BFGS-B
then Powell, colour table fixed) with first-improvement single-bit flips of
the colour table (boundaries fixed), until neither improves, then tries random
multi-bit kicks.  The evaluator is the C port src/twof.c in floating point,
checked here against exact.psi_two_scale and blocks.psi_blocks before use.
Nothing this script prints is a certificate; it reports floating-point values
against the target (117/2192 for k = 3, 1/72 for k = 4), and a result would
need exact certification before it could be claimed.

    python src/twofree.py 4 --m 66 --nb 2-8 --starts 12
    python src/twofree.py 3 --m 1 --nb 2-30 --starts 20
"""

import argparse
import ctypes
import glob
import json
import os
import sys
import time
from fractions import Fraction
from multiprocessing import Pool

import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import blocks  # noqa: E402
import exact   # noqa: E402

TARGET = {3: Fraction(117, 2192), 4: Fraction(1, 72), 5: Fraction(1, 384)}
_LIB = None


def lib():
    global _LIB
    if _LIB is None:
        so, src = os.path.join(HERE, "libtwof.so"), os.path.join(HERE, "twof.c")
        if not os.path.exists(so) or os.path.getmtime(so) < os.path.getmtime(src):
            assert os.system(f"cc -O3 -shared -fPIC -o '{so}' '{src}'") == 0
        L = ctypes.CDLL(so)
        L.psi_two.restype = ctypes.c_double
        L.psi_two.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_int,
                              ctypes.POINTER(ctypes.c_double), ctypes.c_int]
        _LIB = L
    return _LIB


def psi(C, m, bounds, k):
    b = np.ascontiguousarray(bounds, dtype=np.float64)
    return lib().psi_two(bytes(C), m, b.size - 1,
                         b.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), k)


def bounds_of(L):
    L = np.abs(np.asarray(L, dtype=float))
    L = L / L.sum()
    b = np.concatenate(([0.0], np.cumsum(L)))
    b[-1] = 1.0
    return b


def self_check():
    worst = 0.0
    n = 0
    for path in glob.glob(os.path.join(ROOT, "results", "two_k*.json")):
        k = int(os.path.basename(path).split("two_k")[1][0])
        for _, r in list(json.load(open(path)).items())[:6]:
            C = [int(c) for c in r["colour"]]
            v = psi(C, r["m"], np.linspace(0, 1, r["B"] + 1), k)
            worst = max(worst, abs(v - float(Fraction(r["psi"]))))
            n += 1
    rng = np.random.default_rng(0)
    for _ in range(30):
        k, nb = int(rng.integers(3, 7)), int(rng.integers(1, 13))
        b = bounds_of(rng.dirichlet(np.ones(nb)))
        worst = max(worst, abs(psi([i % 2 for i in range(nb)], 1, b, k)
                               - blocks.psi_blocks(b.tolist(), k)))
        n += 1
    assert worst < 1e-12, worst
    return n, worst


def opt_bounds(C, m, L, k):
    f = lambda v: psi(C, m, bounds_of(v), k)
    r = minimize(f, L, method="L-BFGS-B", options={"eps": 1e-9, "ftol": 1e-15, "maxiter": 2000})
    r = minimize(f, np.abs(r.x) / np.abs(r.x).sum(), method="Powell",
                 options={"xtol": 1e-9, "ftol": 1e-14, "maxfev": 6000})
    L = np.abs(r.x) / np.abs(r.x).sum()
    return L, f(L)


def flips(C, m, L, k, val):
    b = bounds_of(L)
    improved = True
    while improved:
        improved = False
        for i in range(len(C)):
            C[i] ^= 1
            v = psi(C, m, b, k)
            if v < val - 1e-15:
                val, improved = v, True
            else:
                C[i] ^= 1
    return C, val


def run_start(task):
    k, m, nb, seed, C0, kicks = task
    t = time.time()
    rng = np.random.default_rng(seed)
    C = list(C0) if C0 is not None else [int(x) for x in rng.integers(0, 2, m * nb)]
    if m == 1:
        C = [i % 2 for i in range(nb)]
    L = rng.dirichlet(np.ones(nb) * 4) if seed else np.ones(nb) / nb
    val = psi(C, m, bounds_of(L), k)

    def descend(C, L, val):
        while True:
            L, v1 = opt_bounds(C, m, L, k)
            if m == 1:
                return C, L, v1
            C, v2 = flips(C, m, L, k, v1)
            if v2 >= val - 1e-14 and v1 >= val - 1e-14:
                return C, L, min(v1, v2)
            val = v2

    C, L, val = descend(C, L, val)
    for _ in range(kicks):
        C2, L2 = list(C), L * np.exp(rng.normal(0, 0.2, L.size))
        if m > 1:
            for i in rng.choice(len(C2), size=max(1, len(C2) // 20), replace=False):
                C2[i] ^= 1
        C2, L2, v2 = descend(C2, L2 / L2.sum(), psi(C2, m, bounds_of(L2), k))
        if v2 < val - 1e-14:
            C, L, val = C2, L2, v2
    return val, m, nb, seed, "".join(map(str, C)), L.tolist(), time.time() - t


def stored_two(k, m, nb):
    """Stored two-scale colour tables with this m and at most nb coarse blocks,
    refined onto nb equal blocks (q -> q*B // nb), as starting points."""
    out = []
    for path in glob.glob(os.path.join(ROOT, "results", f"two_k{k}*.json")):
        for key, r in json.load(open(path)).items():
            if r["m"] != m or r["B"] > nb:
                continue
            B = r["B"]
            C = [int(r["colour"][res * B + (q * B) // nb]) for res in range(m) for q in range(nb)]
            out.append((Fraction(r["psi"]), key, C))
    out.sort()
    return out[:3]


def parse_range(s):
    if "-" in s:
        a, b = s.split("-")
        return list(range(int(a), int(b) + 1))
    return [int(x) for x in s.split(",")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("k", type=int)
    ap.add_argument("--m", default="1", help="fine moduli, e.g. 1 or 11,22,66")
    ap.add_argument("--nb", default="2-8", help="numbers of coarse blocks")
    ap.add_argument("--starts", type=int, default=8, help="random starts per (m, nb)")
    ap.add_argument("--kicks", type=int, default=6)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    k = args.k
    n, worst = self_check()
    print(f"C evaluator vs exact.psi_two_scale / blocks.psi_blocks: {n} cases, max |diff| = {worst:.1e}",
          flush=True)
    rng = np.random.default_rng(args.seed)
    tasks = []
    for m in parse_range(args.m):
        for nb in parse_range(args.nb):
            if m == 1 and nb % 1:
                continue
            for v, key, C in (stored_two(k, m, nb) if m > 1 else []):
                tasks.append((k, m, nb, 0, C, args.kicks))
            for _ in range(args.starts):
                tasks.append((k, m, nb, int(rng.integers(1, 2 ** 31)), None, args.kicks))
    target = TARGET[k]
    print(f"{len(tasks)} starts; target {target} = {float(target):.10f}", flush=True)
    res = []
    with Pool(args.workers) as pool:
        for r in pool.imap_unordered(run_start, tasks):
            res.append(r)
            v, m, nb, seed, C, L, dt = r
            print(f"  m={m:>3} nb={nb:>2} seed={seed:>10}: Psi={v:.10f}"
                  f" ({100 * (v / float(target) - 1):+.3f}% vs target) {dt:.0f}s", flush=True)
    res.sort()
    print("best 10:")
    for v, m, nb, seed, C, L, _ in res[:10]:
        print(f"  {v:.12f}  m={m} nb={nb} seed={seed}  {100 * (v / float(target) - 1):+.4f}%"
              f"  lengths={np.round(L, 6).tolist()}" + (f"  C={C}" if m > 1 else ""), flush=True)
    v = res[0][0]
    print(f"best floating-point Psi_{k} = {v:.12f} vs target {float(target):.12f}:"
          f" {'BELOW target (needs exact certification)' if v < float(target) - 1e-12 else 'not below target'}",
          flush=True)


if __name__ == "__main__":
    main()
