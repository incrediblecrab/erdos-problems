"""Continuous (multilinear) relaxation of the block-colouring objective.

Both objectives used in this project have the shape

    F(c) = (1/den) * sum_U  w_U * [ c is constant on the index set U ]

where U ranges over *sets* of distinct block indices (see exact.block_weights_int
and exact.triangle_weights_int).  Writing the colouring as v_i = +-1, the
indicator of "constant on U" is

    [const on U] = ( prod_{i in U} (1+v_i) + prod_{i in U} (1-v_i) ) / 2^{|U|},

which is exact at every vertex v in {-1,+1}^m and is *multilinear*: every
variable occurs with degree at most one in every term.  Consequently

  * F extends smoothly to the box [-1,1]^m,
  * the minimum over the box is attained at a vertex (multilinear functions are
    linear along each coordinate, so a minimiser can always be pushed to a
    corner without increasing the value), and
  * along coordinate l the function is F = alpha_l + beta_l * v_l with
    beta_l = dF/dv_l independent of v_l, so a single gradient evaluation gives
    the *exact* gain of every one of the m possible single flips:

        F(flip l) - F(v) = -2 * v_l * beta_l.

That last identity turns steepest-descent local search into one gradient per
step instead of m re-evaluations, and lets L-BFGS-B explore the interior of the
box before being rounded back to a vertex at no loss.

Nothing here is a certificate.  Every word this module proposes is re-certified
in exact rational arithmetic by src/certify.py through src/exact.py.
"""

import argparse
import json
import os
import sys
from fractions import Fraction

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import exact  # noqa: E402


class Objective:
    """Multilinear extension of Phi_k (setting 'zp') or Psi_k (setting 'n')."""

    def __init__(self, m, k, setting="zp", mult=1):
        self.m = m
        self.k = k
        self.setting = setting
        self.mult = mult
        if setting == "zp":
            den, items = exact.block_weights_int(m, k, mult)
        elif setting == "n":
            den, items = exact.triangle_weights_int(m, k)
        elif setting == "zm":
            den, items = exact.zm_weights_int(m, k)
        elif setting == "two":
            # m is the pair (fine modulus, coarse block count)
            fine, coarse = m
            den, items = exact.two_scale_weights_int(fine, coarse, k)
            self.m = fine * coarse
        else:
            raise ValueError("setting must be 'zp', 'n', 'zm' or 'two'")
        self.den = den
        # group the terms by |U| so each group is a dense integer array
        groups = {}
        for U, w in items:
            groups.setdefault(len(U), ([], []))
            groups[len(U)][0].append(U)
            groups[len(U)][1].append(w)
        self.groups = []
        for s in sorted(groups):
            idx, wts = groups[s]
            self.groups.append((
                s,
                np.asarray(idx, dtype=np.int64),            # (T, s)
                np.asarray(wts, dtype=np.float64) / den / (2.0 ** s),
            ))
        self.nterms = sum(g[1].shape[0] for g in self.groups)

    def value_grad(self, v):
        """Return (F(v), grad F(v)) for v in [-1,1]^m."""
        val = 0.0
        grad = np.zeros(self.m, dtype=np.float64)
        for s, idx, w in self.groups:
            p = 1.0 + v[idx]                      # (T, s)
            q = 1.0 - v[idx]
            # leave-one-out products via prefix/suffix scans (no division)
            pre_p = np.ones_like(p)
            pre_q = np.ones_like(q)
            np.cumprod(p[:, :-1], axis=1, out=pre_p[:, 1:])
            np.cumprod(q[:, :-1], axis=1, out=pre_q[:, 1:])
            suf_p = np.ones_like(p)
            suf_q = np.ones_like(q)
            np.cumprod(p[:, :0:-1], axis=1, out=suf_p[:, -2::-1])
            np.cumprod(q[:, :0:-1], axis=1, out=suf_q[:, -2::-1])
            full_p = pre_p[:, -1] * p[:, -1]
            full_q = pre_q[:, -1] * q[:, -1]
            val += float(np.dot(w, full_p + full_q))
            loo = (pre_p * suf_p - pre_q * suf_q) * w[:, None]
            # bincount is far faster than np.add.at for this scatter-add
            grad += np.bincount(idx.ravel(), weights=loo.ravel(),
                                minlength=self.m)
        return val, grad

    def value(self, v):
        return self.value_grad(v)[0]

    def flip_gains(self, v):
        """Exact change in F from flipping each coordinate, valid at vertices."""
        _, g = self.value_grad(v)
        return -2.0 * v * g

    def descend(self, v, max_steps=100000):
        """Steepest-descent 1-flip local search from a vertex.  Returns (v, F)."""
        val, _ = self.value_grad(v)
        for _ in range(max_steps):
            gains = self.flip_gains(v)
            l = int(np.argmin(gains))
            if gains[l] >= -1e-13:
                break
            v = v.copy()
            v[l] = -v[l]
            val += gains[l]
        return v, self.value(v)


    def tabu(self, v, iters, rng, tenure=None, best_val=None):
        """Tabu search over vertices.

        Each iteration costs exactly one gradient evaluation, because
        flip_gains() returns the exact gain of all m candidate moves at once.
        The best non-tabu move is taken even when it worsens the value, which is
        what lets the walk leave the deep, narrow basins that defeat plain
        descent; the aspiration criterion overrides the tabu list whenever a
        move would set a new global best.
        """
        m = self.m
        if tenure is None:
            tenure = max(3, m // 8)
        v = v.copy()
        val, _ = self.value_grad(v)
        best_v, bv = v.copy(), val
        if best_val is not None and best_val < bv:
            bv = best_val
        tabu_until = np.zeros(m, dtype=np.int64)
        for it in range(1, iters + 1):
            gains = self.flip_gains(v)
            allowed = tabu_until < it
            aspire = (val + gains) < bv - 1e-13
            ok = allowed | aspire
            if not ok.any():
                ok = allowed
                if not ok.any():
                    tabu_until[:] = 0
                    ok = np.ones(m, dtype=bool)
            cand = np.where(ok, gains, np.inf)
            l = int(np.argmin(cand))
            v[l] = -v[l]
            val += gains[l]
            tabu_until[l] = it + tenure + int(rng.integers(0, max(2, tenure)))
            if val < bv - 1e-13:
                bv, best_v = val, v.copy()
        return best_v, bv


def word_of(v):
    return "".join("1" if x > 0 else "0" for x in v)


def stretch(word, r):
    """Subdivide every block r ways.  This denotes the *same* subset of [0,1),
    so both Phi and Psi are preserved exactly -- it is the correct way to move
    a record from resolution m to resolution r*m."""
    return "".join(ch * r for ch in word)


def repeat(word, r):
    """Tile the word r times.  This is the pullback of the colouring under
    x -> r*x on the torus, which preserves Phi exactly (that map is a measure-
    preserving r^2-to-1 cover of the (x,d) torus) but *not* Psi, because x -> rx
    does not preserve the triangle x + (k-1)t <= 1.  Useful as a zp seed only."""
    return word * r


def vec_of(word):
    return np.array([1.0 if ch == "1" else -1.0 for ch in word], dtype=np.float64)


def certify(word, k, setting):
    """Re-derive the value in exact rational arithmetic (the actual certificate)."""
    c = [1 if ch == "1" else 0 for ch in word]
    if setting == "zp":
        val = exact.phi(c, k)
        return val, val / 2
    val = exact.psi(c, k)
    return val, val


def polish_pairs(obj, v, rounds=4):
    """Greedy 2-flip descent; 1-flip is already exhausted by descend()."""
    v, val = obj.descend(v)
    m = obj.m
    for _ in range(rounds):
        improved = False
        base = val
        for i in range(m):
            v2 = v.copy()
            v2[i] = -v2[i]
            cand, cval = obj.descend(v2)
            if cval < base - 1e-13:
                v, val, base = cand, cval, cval
                improved = True
        if not improved:
            break
    return v, val


def search(m, k, setting, restarts, seeds=(), rng=None, pairs=False,
           tabu_iters=0, rounds=None, obj=None):
    """Multi-start L-BFGS-B in the box plus iterated tabu search over vertices.

    L-BFGS-B supplies diverse starting vertices (randomised rounding of interior
    points); tabu search does the actual optimisation.  Every candidate returned
    here is a float-precision proposal only -- it is re-certified exactly by the
    caller before being recorded.
    """
    from scipy.optimize import minimize

    obj = obj or Objective(m, k, setting)
    rng = rng or np.random.default_rng(12345)
    best = (float("inf"), None)

    def consider(v):
        nonlocal best
        v, val = obj.descend(np.sign(v) + (v == 0))
        if pairs:
            v, val = polish_pairs(obj, v)
        if val < best[0]:
            best = (val, v.copy())
        return v

    for w in seeds:                      # warm starts (e.g. blown-up records)
        consider(vec_of(w))

    fun = lambda x: obj.value_grad(x)    # noqa: E731
    for _ in range(restarts):
        x0 = rng.uniform(-1.0, 1.0, size=m)
        r = minimize(fun, x0, jac=True, method="L-BFGS-B",
                     bounds=[(-1.0, 1.0)] * m,
                     options={"maxiter": 2000, "ftol": 1e-15, "gtol": 1e-12})
        consider(r.x)

    if tabu_iters:
        start = best[1].copy()
        nrounds = rounds if rounds is not None else max(1, restarts // 4)
        for _ in range(nrounds):
            v, val = obj.tabu(start, tabu_iters, rng, best_val=best[0])
            if val < best[0]:
                best = (val, v.copy())
            # perturb the incumbent to reseed the next round
            start = best[1].copy()
            flip = rng.choice(m, size=max(1, m // 6), replace=False)
            start[flip] = -start[flip]
    return obj, best[0], word_of(best[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("k", type=int)
    ap.add_argument("ms", help="comma list and/or a-b ranges of m")
    ap.add_argument("--setting", default="zp", choices=("zp", "n"))
    ap.add_argument("--restarts", type=int, default=200)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--pairs", action="store_true")
    ap.add_argument("--tabu", type=int, default=0,
                    help="tabu iterations per round (0 disables)")
    ap.add_argument("--rounds", type=int, default=None,
                    help="tabu rounds (default: restarts//4)")
    ap.add_argument("--store", default=None, help="results json to update")
    args = ap.parse_args()

    ms = []
    for part in args.ms.split(","):
        if "-" in part:
            a, b = part.split("-")
            ms.extend(range(int(a), int(b) + 1))
        else:
            ms.append(int(part))

    store = args.store or f"results/{args.setting}_k{args.k}.json"
    db = {}
    if os.path.exists(store):
        db = json.load(open(store))

    rng = np.random.default_rng(args.seed)
    for m in ms:
        seeds = []
        for d in range(1, m):
            if m % d == 0 and str(d) in db:
                w = db[str(d)]["word"]
                seeds.append(stretch(w, m // d))       # same subset of [0,1)
                if args.setting == "zp":
                    seeds.append(repeat(w, m // d))    # Phi-preserving only
        if str(m) in db:
            seeds.append(db[str(m)]["word"])
        obj, val, word = search(m, args.k, args.setting, args.restarts,
                                seeds=seeds, rng=rng, pairs=args.pairs,
                                tabu_iters=args.tabu, rounds=args.rounds)
        exact_val, bound = certify(word, args.k, args.setting)
        assert abs(float(exact_val) - val) < 1e-9, (m, float(exact_val), val)
        prev = db.get(str(m))
        better = prev is None or Fraction(exact_val) < Fraction(prev["value"])
        flag = "NEW " if better else "    "
        print(f"{flag}m={m:4d} k={args.k} {args.setting} terms={obj.nterms:8d} "
              f"val={str(exact_val):>24s} = {float(exact_val):.10f} "
              f"bound={float(bound):.10f}", flush=True)
        if better:
            db[str(m)] = {
                "word": word,
                "value": str(exact_val),
                "value_float": float(exact_val),
                "bound": str(bound),
                "bound_float": float(bound),
                "found_by": "relax",
                "exhaustive": False,
            }
            with open(store, "w") as fh:
                json.dump(db, fh, indent=1, sort_keys=True)


if __name__ == "__main__":
    main()
