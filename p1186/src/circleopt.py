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

--wide (added November 2, 2026) runs a broader search with the C port of
circle.phi_circle (src/circlef.c, about 30 times faster, checked against the
Python evaluator before use):

  stage A  many short local optimisations: the stored words, and random
           Dirichlet arc lengths for every even number of arcs in
           --bmin..--bmax, unconstrained and in the symmetric subspaces
           (antipodal c(x+1/2) = 1-c(x), mirror x -> -x, or both);
  stage B  the best --keep distinct stage-A optima are basin-hopped without
           symmetry constraints: log-normal perturbations, splitting an arc
           in three and merging three arcs into one, each followed by a
           full local optimisation, keeping only improvements;
  round    the best --round final optima are rounded onto every grid
           --mmin..--mmax and the best --polish grids are polished.

    python src/circleopt.py 5 --wide --bmin 2 --bmax 30 --random 6 --workers 4

The rounded word is certified by exact.phi when its weight table fits
(exact.phi_table_fits) and by exact.phi_bits (the same formula, streamed) otherwise, and must agree with
circle.phi_circle in Fractions.  Floating-point values in the transcript are
search output.
"""

import argparse
import ctypes
import json
import os
import sys
import time
from fractions import Fraction
from glob import glob
from multiprocessing import Pool

import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import circle  # noqa: E402
import exact   # noqa: E402

_LIB = None


def _lib():
    """The C port of circle.phi_circle (search only), built on first use."""
    global _LIB
    if _LIB is None:
        so = os.path.join(HERE, "libcirclef.so")
        src = os.path.join(HERE, "circlef.c")
        if not os.path.exists(so) or os.path.getmtime(so) < os.path.getmtime(src):
            rc = os.system(f"cc -O3 -shared -fPIC -o '{so}' '{src}'")
            assert rc == 0, "could not build src/libcirclef.so"
        lib = ctypes.CDLL(so)
        lib.phi_lengths.restype = ctypes.c_double
        lib.phi_lengths.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_int]
        lib.phi_circle.restype = ctypes.c_double
        lib.phi_circle.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_int]
        _LIB = lib
    return _LIB


def fast_phi_lengths(L, k):
    L = np.ascontiguousarray(np.asarray(L, dtype=np.float64))
    return _lib().phi_lengths(L.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), L.size, k)


def fast_phi_cuts(cuts, k):
    M = cuts[-1]
    b = np.ascontiguousarray(np.asarray(cuts, dtype=np.float64) / M)
    b[-1] = 1.0
    return _lib().phi_circle(b.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), b.size - 1, k)


def check_fast(k, trials=200, seed=0):
    """The C port must agree with circle.phi_circle to 1e-12 on random colourings."""
    rng = np.random.default_rng(seed)
    worst = 0.0
    for _ in range(trials):
        nb = 2 * int(rng.integers(1, 21))
        L = rng.dirichlet(np.ones(nb))
        b = np.concatenate(([0.0], np.cumsum(L)))
        b[-1] = 1.0
        worst = max(worst, abs(fast_phi_lengths(L, k) - circle.phi_circle(b.tolist(), k)))
    assert worst < 1e-12, f"C evaluator disagrees with circle.py by {worst}"
    return worst


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


# ------------------------------------------------------------------ --wide

def expand(v, sym):
    """Free parameters -> arc lengths (colour 0 first) for a symmetry class.

    none        the lengths themselves
    anti        c(x + 1/2) = 1 - c(x): lengths l + l with len(l) odd
    mirror      c(-x) = c(x) about the centre of arc 0: L_i = L_{B-i}
    antimirror  both: the mirror pattern on an odd half, repeated twice
    """
    v = list(np.abs(np.asarray(v, dtype=float)))
    if sym == "none":
        return np.array(v)
    if sym == "anti":
        return np.array(v + v)
    if sym == "mirror":
        h = len(v) - 1
        return np.array([v[0]] + v[1:h] + [v[h]] + v[1:h][::-1])
    if sym == "antimirror":
        half = [v[0]] + v[1:] + v[1:][::-1]
        return np.array(half + half)
    raise ValueError(sym)


def nfree(B, sym):
    """Number of free parameters for B arcs in a symmetry class, or None."""
    if sym == "none":
        return B
    if sym == "anti":
        return B // 2 if (B // 2) % 2 == 1 else None
    if sym == "mirror":
        return B // 2 + 1
    if sym == "antimirror":
        return (B // 2 + 1) // 2 if (B // 2) % 2 == 1 else None
    return None


def local_opt(v0, k, sym, maxfev, full=True):
    """L-BFGS-B then (if full) Powell on phi(expand(v)), with the C evaluator."""
    f = lambda v: fast_phi_lengths(expand(v, sym), k)
    r = minimize(f, v0, method="L-BFGS-B",
                 options={"eps": 1e-9, "ftol": 1e-15, "gtol": 1e-12, "maxiter": 3000})
    v = np.abs(r.x) / np.abs(r.x).sum()
    if full:
        r = minimize(f, v, method="Powell",
                     options={"xtol": 1e-9, "ftol": 1e-14, "maxfev": maxfev})
        v = np.abs(r.x) / np.abs(r.x).sum()
    L = expand(v, sym)
    L = L / L.sum()
    return L, fast_phi_lengths(L, k)


def clean(L, tol=1e-7):
    """Drop arcs that collapsed to (near) zero by merging them with their neighbours."""
    L = list(L)
    changed = True
    while changed and len(L) > 2:
        changed = False
        for i, x in enumerate(L):
            if x < tol:
                n = len(L)
                a, c = (i - 1) % n, (i + 1) % n
                merged = L[a] + x + L[c]
                # rotation and colour swap preserve Phi, so start at the merged
                # arc and continue in cyclic order after arc c
                L = [merged] + [L[(c + 1 + t) % n] for t in range(n - 3)]
                changed = True
                break
    L = np.array(L)
    return L / L.sum()


def hop_move(L, rng, bmax):
    """One basin-hopping proposal: perturb, split an arc in three, or merge three."""
    L = list(L)
    B = len(L)
    u = rng.random()
    if u < 0.2 and B + 2 <= bmax:
        i = int(rng.integers(B))
        eps = rng.uniform(0.05, 0.3)
        s = rng.uniform(0.1, 0.9) * (1 - eps)
        L[i:i + 1] = [L[i] * s, L[i] * eps, L[i] * (1 - eps - s)]
        move = "split"
    elif u < 0.4 and B >= 4:
        w = 1.0 / np.asarray(L)
        i = int(rng.choice(B, p=w / w.sum()))
        L = L[i:] + L[:i]                       # rotate: arc i first
        L = [L[-1] + L[0] + L[1]] + L[2:-1]
        move = "merge"
    else:
        scale = rng.choice([0.05, 0.15, 0.4])
        L = list(np.asarray(L) * np.exp(rng.normal(0.0, scale, size=B)))
        move = f"perturb{scale}"
    L = np.asarray(L)
    return L / L.sum(), move


def stage_a(task):
    k, sym, B, seed, v0, maxfev = task
    t = time.time()
    try:
        L, v = local_opt(np.asarray(v0), k, sym, maxfev)
        L = clean(L)
        v = fast_phi_lengths(L, k)
    except Exception as exc:              # a failed start is reported, not fatal
        return (float("inf"), sym, B, seed, None, f"error {exc}", 0.0)
    return (v, sym, B, seed, L.tolist(), "", time.time() - t)


def stage_b(task):
    k, L0, seed, hops, bmax, maxfev = task
    rng = np.random.default_rng(seed)
    L, v = local_opt(np.asarray(L0), k, "none", maxfev)
    L = clean(L)
    v = fast_phi_lengths(L, k)
    hist = []
    for h in range(hops):
        trial, move = hop_move(L, rng, bmax)
        try:
            T, tv = local_opt(trial, k, "none", maxfev)
            T = clean(T)
            tv = fast_phi_lengths(T, k)
        except Exception:
            continue
        if tv < v - 1e-13:
            hist.append((h, move, len(L), len(T), v, tv))
            L, v = T, tv
    return (v, L.tolist(), seed, hist)


def round_and_polish(L, k, mmin, mmax, npolish):
    rounded = []
    for M in range(mmin, mmax + 1):
        cuts = cuts_on_grid(np.asarray(L), M)
        if any(cuts[i + 1] <= cuts[i] for i in range(len(cuts) - 1)):
            continue
        rounded.append((fast_phi_cuts(cuts, k), M, cuts))
    rounded.sort()
    out = []
    for v_round, M, cuts in rounded[:npolish]:
        cuts, v = polish_fast(cuts, k)
        out.append((v, M, cuts, v_round))
    out.sort()
    return out


def polish_fast(cuts, k):
    """polish() with the C evaluator, plus paired moves of two adjacent cuts."""
    cuts = list(cuts)
    val = fast_phi_cuts(cuts, k)
    improved = True
    while improved:
        improved = False
        n = len(cuts)
        for i in range(1, n - 1):
            for step in (-1, 1, -2, 2):
                c = cuts[i] + step
                if not (cuts[i - 1] < c < cuts[i + 1]):
                    continue
                trial = cuts[:i] + [c] + cuts[i + 1:]
                v = fast_phi_cuts(trial, k)
                if v < val - 1e-15:
                    cuts, val, improved = trial, v, True
        for i in range(1, n - 2):                     # shift an arc by one block
            for step in (-1, 1):
                a, b = cuts[i] + step, cuts[i + 1] + step
                if not (cuts[i - 1] < a < b < cuts[i + 2]):
                    continue
                trial = cuts[:i] + [a, b] + cuts[i + 2:]
                v = fast_phi_cuts(trial, k)
                if v < val - 1e-15:
                    cuts, val, improved = trial, v, True
    return cuts, val


def certify_word(word, k):
    """Exact Phi of an equal-block word: exact.phi (with phi_bits agreeing) when
    the weight table fits in memory, exact.phi_bits beyond; circle.py must agree."""
    c = [int(ch) for ch in word]
    M = len(c)
    vb = exact.phi_bits(c, k)
    if exact.phi_table_fits(M, k):
        val = exact.phi(c, k)
        assert val == vb, (val, vb)
        how = "exact.phi (= exact.phi_bits)"
    else:
        val, how = vb, "exact.phi_bits"
    bounds, _ = circle.bounds_of_cyclic_word(word)
    assert circle.phi_circle(bounds, k) == val, "circle.py disagrees"
    return val, how


def wide(args):
    k = args.k
    t0 = time.time()
    print(f"C evaluator vs circle.py: max |diff| = {check_fast(k):.1e} on 200 random colourings",
          flush=True)
    rng = np.random.default_rng(args.seed)
    starts = stored_words(k, args.top)
    tasks = []
    for v0, key, w in starts:
        bounds, _ = circle.bounds_of_cyclic_word(w)
        L0 = np.diff([float(b) for b in bounds])
        tasks.append((k, "none", L0.size, f"stored m={key}", L0.tolist(), args.maxfev_a))
    syms = args.sym.split(",")
    for B in range(max(2, args.bmin + args.bmin % 2), args.bmax + 1, 2):
        for sym in syms:
            n = nfree(B, sym)
            if n is None or n < 1:
                continue
            for r in range(args.random):
                seed = int(rng.integers(2 ** 31))
                v0 = np.random.default_rng(seed).dirichlet(np.ones(n))
                tasks.append((k, sym, B, seed, v0.tolist(), args.maxfev_a))
    print(f"stage A: {len(tasks)} starts, B = {args.bmin}..{args.bmax}, symmetries {syms},"
          f" {args.random} random per (B, symmetry), {args.workers} workers", flush=True)
    with Pool(args.workers) as pool:
        resA = []
        for res in pool.imap_unordered(stage_a, tasks):
            resA.append(res)
            v, sym, B, seed, L, err, dt = res
            arcs = len(L) if L else 0
            print(f"  A {sym:<10} B={B:>2} seed={seed}: Phi={v:.10f} arcs={arcs} {dt:.0f}s {err}",
                  flush=True)
    resA = [r for r in resA if r[4] is not None]
    resA.sort(key=lambda r: r[0])
    print("stage A best 15:", flush=True)
    for v, sym, B, seed, L, _, _ in resA[:15]:
        print(f"  {v:.12f}  {sym:<10} B={B:>2} arcs={len(L)} seed={seed}", flush=True)
    keep, seen = [], []
    for r in resA:
        if any(abs(r[0] - s) < 1e-9 and len(r[4]) == ls for s, ls in seen):
            continue
        seen.append((r[0], len(r[4])))
        keep.append(r)
        if len(keep) >= args.keep:
            break
    tasksB = [(k, r[4], int(rng.integers(2 ** 31)), args.hops, args.bmax_hop, args.maxfev_b)
              for r in keep]
    print(f"stage B: {len(tasksB)} basins x {args.hops} hops", flush=True)
    with Pool(args.workers) as pool:
        resB = pool.map(stage_b, tasksB)
    resB.sort(key=lambda r: r[0])
    for (v, L, seed, hist), r in zip(resB, [None] * len(resB)):
        print(f"  B Phi={v:.12f} arcs={len(L)} seed={seed} improving hops={len(hist)}", flush=True)
        for h, move, b0, b1, va, vb in hist:
            print(f"      hop {h:>3} {move:<12} arcs {b0}->{b1}  {va:.12f} -> {vb:.12f}", flush=True)
    print(f"best continuous Phi = {resB[0][0]:.12f} (floating point, search output); arc lengths:",
          np.round(resB[0][1], 6).tolist(), flush=True)

    cands = []
    for v, L, seed, _ in resB[:args.round]:
        out = round_and_polish(L, k, args.mmin, args.mmax, args.polish)
        for vv, M, cuts, vr in out[:3]:
            print(f"  basin {v:.12f}: M={M:5d} rounded {vr:.12f} polished {vv:.12f}", flush=True)
        cands += out[:3]
    cands.sort()
    vgrid, M, cuts, _ = cands[0]
    word = circle.word_of_lengths([cuts[i + 1] - cuts[i] for i in range(len(cuts) - 1)])
    val, how = certify_word(word, k)
    assert abs(float(val) - vgrid) < 1e-9, (float(val), vgrid)
    print(f"certified m={M}: Phi = {val} = {float(val):.12f}, delta~_{k} <= {val / 2}"
          f"  [{how}; circle.py agrees]", flush=True)
    best_stored = stored_words(k, 1)[0][0]
    print(f"best stored word: Phi = {best_stored} = {float(best_stored):.12f};"
          f" {'IMPROVES' if val < best_stored else 'does not improve'}"
          f"  ({time.time() - t0:.0f}s total)", flush=True)
    if args.store and val < best_stored:
        store(args.store, word, val, "circleopt --wide")


def store(path, word, val, found_by):
    db = json.load(open(path)) if os.path.exists(path) else {}
    db[str(len(word))] = {
        "word": word,
        "value": str(val),
        "value_float": float(val),
        "bound": str(val / 2),
        "bound_float": float(val / 2),
        "found_by": found_by,
        "exhaustive": False,
    }
    with open(path, "w") as fh:
        json.dump(db, fh, indent=1, sort_keys=True)
    print(f"stored in {path}", flush=True)


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
    ap.add_argument("--wide", action="store_true", help="many-start search with the C evaluator")
    ap.add_argument("--bmin", type=int, default=2, help="--wide: fewest arcs in random starts")
    ap.add_argument("--bmax", type=int, default=30, help="--wide: most arcs in random starts")
    ap.add_argument("--bmax-hop", type=int, default=40, help="--wide: most arcs a split may create")
    ap.add_argument("--random", type=int, default=4, help="--wide: random starts per (B, symmetry)")
    ap.add_argument("--sym", default="none,anti,mirror,antimirror", help="--wide: symmetry classes")
    ap.add_argument("--keep", type=int, default=8, help="--wide: stage-A basins to hop from")
    ap.add_argument("--round", type=int, default=3, help="--wide: final basins to round")
    ap.add_argument("--maxfev-a", type=int, default=4000)
    ap.add_argument("--maxfev-b", type=int, default=20000)
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    if args.wide:
        return wide(args)
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
