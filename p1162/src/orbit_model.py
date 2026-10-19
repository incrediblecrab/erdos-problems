"""Where the second-order term of log_2 e(n) comes from.

`asymptotics.py` measures  excess(n) = log_2 e(n) - n^2/16  and fits it to
A n log_2 n + B n + C.  This module explains the A, and gives a falsifiable
prediction that the measurement can be tested against.

The model
---------
Let H = F_2^k act on [n].  Split [n] by orbit size; all orbits have size
2^j for some j >= 0.  The number of homomorphisms F_2^k -> S_n whose orbits
*all* have size m = 2^j is exactly

      n!/t! * ( [k choose j]_2 / m )^t ,            t = n/m,

because a partition of [n] into t unordered blocks of size m can be made in
n!/(m!^t t!) ways and each block carries [k,j]_2 * (m-1)! transitive
actions.  Dividing by |GL_k(2)| = 2^{k^2}(1+o(1)) turns homomorphisms into
(an upper bound for) subgroups.  Taking logs with t = n/2^j,

  log_2 N_j(n,k)/|GL_k(2)|
      = (1 - 2^-j) n log_2 n  +  j k n / 2^j  -  k^2  +  O_j(n).

Maximising the n^2 part over k gives k = j n / 2^{j+1} and the value
j^2 n^2 / 2^{2j+2}:

      j = 1:  k = n/4,    n^2/16,     n log n coefficient 1/2
      j = 2:  k = n/4,    n^2/16,     n log n coefficient 3/4   <-- wins
      j = 3:  k = 3n/16,  9n^2/256    (loses already at order n^2)

So j = 1 and j = 2 tie at the leading order n^2/16 -- the first-order
constant of Pyber and of RoTr25 Theorem 1 -- and j = 2 wins the tie-break at
order n log n.  Allowing mixtures does not help: if c_j is the fraction of
points in orbits of size 2^j then the n^2 coefficient is
kappa * sum_j c_j j/2^j - kappa^2 with kappa = k/n, and sum_j c_j j/2^j <= 1/2
with equality iff c_j is supported on {1,2}; on that face the n log n
coefficient is 1/2 + c_2/4, maximised at c_2 = 1.

  PREDICTION A:  log_2 e(n) = n^2/16 + (3/4) n log_2 n + O(n),
                 i.e. the fitted A converges to 0.75.

  PREDICTION B:  the approach is slow.  Moving a fraction eps of the points
                 from 4-orbits into 2-orbits costs only (eps/4) n log_2 n but
                 buys entropy, so A_local(n) approaches 3/4 from below and the
                 limit has to be extrapolated rather than read off.

                 [The first version of this file guessed eps ~ c/log_2 n.
                 That is FALSE: measured eps*log_2 n falls 1.589 -> 0.924 over
                 n = 64..512 while eps*sqrt(n) only drifts 2.119 -> 2.323, and
                 a held-out scan over A_inf + c n^-p selects p = 0.48.  The
                 guessed exponent is kept here as a rejected hypothesis; the
                 exponent is selected by held-out error, not asserted.  The
                 first version also said "off 4-orbits"; eps is the 2-orbit
                 fraction.  NOTES.md section 9.2 derives the value 1/2
                 heuristically from the saddle point of the count M(n) of set
                 partitions into blocks of size 2 and 4: eps ~ sqrt(6/n), and
                 A_local = 3/4 - (sqrt6/8) n^-1/2 + O(n^-3/2).]

  PREDICTION C:  at the dominant rank k = n/4 the exact expected fraction of
                 points lying in orbits of size 4 tends to 1, slowly, and the
                 residual mass sits on orbits of size 2 (not 8 or more).

Predictions B and C are what make this testable at n <= 512: A itself
converges like a power of n and cannot be read off directly.

All logarithms are base 2.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import sys
from math import comb, factorial

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from asymptotics import log2_int, ratio  # noqa: E402
from elemab import gaussian_binomial, gl_order, hom_counts  # noqa: E402


def orbit_profile(n: int, k: int, row: list[int] | None = None) -> dict:
    """Exact orbit-size profile of a uniformly random hom F_2^k -> S_n.

    For a species with EGF exp(A(x)), marking orbits of size m gives
    (a_m x^m/m!) exp(A(x)), so summed over all homomorphisms the number of
    orbits of size m is C(n,m) * a_m * h(n-m).  Hence

        E[#orbits of size m] = C(n,m) a_m h(n-m) / h(n).

    Returns expected orbit counts and the expected fraction of the n points
    that lie in orbits of each size.  Exact rationals, reported as floats.
    """
    if row is None:
        row = hom_counts(k, n)
    h_n = row[n]
    prof = {}
    total_points = 0.0
    j = 0
    while 2**j <= n:
        m = 2**j
        a_m = gaussian_binomial(k, j) * factorial(m - 1)
        cnt = ratio(comb(n, m) * a_m * row[n - m], h_n)
        if cnt > 0:
            prof[m] = {"exp_orbits": cnt, "exp_point_fraction": m * cnt / n}
            total_points += m * cnt
        j += 1
    return {
        "n": n,
        "k": k,
        "log2_hom": log2_int(h_n),
        "profile": prof,
        # exact identity sum_m m * E[#orbits_m] = n; a self-check on the code
        "point_total_check": total_points,
        "point_total_error": abs(total_points - n),
    }


def pure_type_log2(n: int, j: int, k: int) -> float | None:
    """log_2 ( N_j(n,k) / |GL_k(2)| ), exactly, or None if 2^j does not divide n."""
    m = 2**j
    if n % m:
        return None
    t = n // m
    qb = gaussian_binomial(k, j)
    if qb == 0:
        return None
    num = factorial(n) // factorial(t) * (qb**t)
    den = (m**t) * gl_order(k)
    return log2_int(num) - log2_int(den)


def pure_type_scan(n: int, j_max: int = 4) -> list[dict]:
    """For each pure orbit type 2^j, the best rank and the excess over n^2/16."""
    out = []
    for j in range(1, j_max + 1):
        if n % (2**j):
            continue
        best = None
        for k in range(1, n // 2 + 1):
            v = pure_type_log2(n, j, k)
            if v is None:
                continue
            if best is None or v > best[1]:
                best = (k, v)
        if best is None:
            continue
        k, v = best
        out.append(
            {
                "j": j,
                "orbit_size": 2**j,
                "best_k": k,
                "best_k_over_n": k / n,
                "log2_value": v,
                "excess_over_n2_16": v - n * n / 16.0,
                "excess_per_n_log2n": (v - n * n / 16.0) / (n * math.log2(n)),
                "predicted_k_over_n": j / 2 ** (j + 1),
                "predicted_n2_coef": j * j / 2 ** (2 * j + 2),
                "predicted_nlogn_coef": 1 - 2.0**-j,
            }
        )
    return out


def local_coefficients(rows: list[dict], h: int) -> list[dict]:
    """Second-difference estimate of A in  excess = A n log2 n + B n + C.

    The second difference over an arithmetic triple (n-h, n, n+h) annihilates
    both B n and C exactly, so

        A_local(n) = D2[excess] / D2[n log2 n]

    with no contamination from the linear or constant term.
    """
    by_n = {r["n"]: r["excess"] for r in rows}
    out = []
    for n in sorted(by_n):
        if n - h not in by_n or n + h not in by_n:
            continue
        d2y = by_n[n + h] - 2 * by_n[n] + by_n[n - h]
        f = lambda t: t * math.log2(t)  # noqa: E731
        d2x = f(n + h) - 2 * f(n) + f(n - h)
        out.append(
            {
                "n": n,
                "h": h,
                "A_local": d2y / d2x,
                "inv_log2n": 1.0 / math.log2(n),
            }
        )
    return out


def fit_local(local: list[dict], n_min: int) -> dict:
    """Least squares of A_local against 1/log2 n.  Intercept should be 3/4."""
    import numpy as np

    sel = [r for r in local if r["n"] >= n_min]
    if len(sel) < 3:
        return {"n_min": n_min, "points": len(sel), "error": "too few points"}
    X = np.array([[r["inv_log2n"], 1.0] for r in sel])
    y = np.array([r["A_local"] for r in sel])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = float(np.max(np.abs(X @ coef - y)))
    return {
        "n_min": n_min,
        "n_max": max(r["n"] for r in sel),
        "points": len(sel),
        "slope_inv_log2n": float(coef[0]),
        "intercept": float(coef[1]),
        "predicted_intercept": 0.75,
        "intercept_error": abs(float(coef[1]) - 0.75),
        "max_abs_residual": resid,
    }


# --------------------------------------------------------- held-out selection

def _np():
    import numpy as np

    return np


BASES = {
    "1/log2 n": lambda n: 1.0 / _np().log2(n),
    "1/log2^2 n": lambda n: 1.0 / _np().log2(n) ** 2,
    "n^-1/4": lambda n: n**-0.25,
    "n^-1/2": lambda n: n**-0.5,
    "n^-1": lambda n: 1.0 / n,
    "n^-1/2 log2 n": lambda n: n**-0.5 * _np().log2(n),
    "log2 n / n": lambda n: _np().log2(n) / n,
}


def holdout_extrapolation(
    local: list[dict], train_hi: int = 300, test_lo: int = 400, n_min: int = 96
) -> dict:
    """Extrapolate A_local to n = infinity, scoring on data not used to fit.

    A_local(n) still drifts at n = 512, so the limit has to be extrapolated,
    and an extrapolation basis chosen by in-sample fit proves nothing: with
    seven candidate shapes one of them will land near any target.  So fit each
    candidate on n <= train_hi only and score it on n >= test_lo, which it
    never saw.  Also scan a free exponent p in A_inf + c n^-p the same way, so
    that the decay rate is selected by held-out error rather than by hand.

    Predicted A_inf = 3/4 from the orbit model at the top of this file.
    """
    try:
        np = _np()
    except ImportError:
        return {"error": "numpy not installed", "points": 0}
    sel = [r for r in local if r["n"] >= n_min]
    N = np.array([r["n"] for r in sel], dtype=float)
    Y = np.array([r["A_local"] for r in sel])
    tr, te = N <= train_hi, N >= test_lo
    if tr.sum() < 4 or te.sum() < 3:
        return {"error": "not enough points to hold out", "points": len(sel)}

    def fit(basis, mask_fit, mask_score):
        X = lambda v: np.stack([basis(v), np.ones_like(v)], 1)  # noqa: E731
        c, *_ = np.linalg.lstsq(X(N[mask_fit]), Y[mask_fit], rcond=None)
        return (
            float(c[1]),
            float(c[0]),
            float(np.max(np.abs(X(N[mask_fit]) @ c - Y[mask_fit]))),
            float(np.max(np.abs(X(N[mask_score]) @ c - Y[mask_score]))),
        )

    named = []
    for nm, b in BASES.items():
        a_inf, coef, r_tr, r_te = fit(b, tr, te)
        named.append(
            {
                "basis": nm,
                "A_inf": a_inf,
                "coef": coef,
                "train_maxres": r_tr,
                "holdout_maxres": r_te,
                "A_inf_error_vs_0.75": abs(a_inf - 0.75),
            }
        )
    named.sort(key=lambda d: d["holdout_maxres"])

    scan = []
    for p in [round(x, 2) for x in np.arange(0.05, 1.51, 0.01)]:
        a_inf, coef, r_tr, r_te = fit(lambda v, p=p: v**-p, tr, te)
        scan.append({"p": p, "A_inf": a_inf, "coef": coef, "holdout_maxres": r_te})
    best = min(scan, key=lambda d: d["holdout_maxres"])

    return {
        "n_min": n_min,
        "train_range": [float(N[tr].min()), float(N[tr].max())],
        "test_range": [float(N[te].min()), float(N[te].max())],
        "n_train": int(tr.sum()),
        "n_holdout": int(te.sum()),
        "ranked_bases": named,
        "best_free_exponent": best,
        "best_free_exponent_error_vs_0.75": abs(best["A_inf"] - 0.75),
        "predicted_A_inf": 0.75,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description="structural model for the n log n term")
    ap.add_argument("--stats", type=pathlib.Path, default=pathlib.Path("data/elemab_stats.json"))
    ap.add_argument("--outdir", type=pathlib.Path, default=pathlib.Path("data"))
    ap.add_argument("--profile-n", type=int, nargs="*", default=[64, 128, 256, 512])
    ap.add_argument("--pure-n", type=int, nargs="*", default=[64, 128, 256])
    args = ap.parse_args()

    result: dict = {}

    print("== Prediction C: exact orbit profile at the dominant rank k = n/4 ==")
    profiles = []
    for n in args.profile_n:
        k = n // 4
        p = orbit_profile(n, k)
        profiles.append(p)
        frac = {m: v["exp_point_fraction"] for m, v in p["profile"].items()}
        print(
            f"n={n:4d} k={k:4d}  "
            + "  ".join(f"|orb|={m}: {f:.6f}" for m, f in sorted(frac.items()) if f > 1e-9)
            + f"   [check sum m*E = {p['point_total_check']:.6f}, err {p['point_total_error']:.2e}]",
            flush=True,
        )
    result["orbit_profiles"] = profiles

    print("\n== Model check: pure orbit types, best rank and excess ==")
    pures = {}
    for n in args.pure_n:
        rows = pure_type_scan(n)
        pures[str(n)] = rows
        for r in rows:
            print(
                f"n={n:4d} orbits of size {r['orbit_size']:2d}: "
                f"best k/n={r['best_k_over_n']:.5f} (pred {r['predicted_k_over_n']:.5f})  "
                f"excess/(n log2 n)={r['excess_per_n_log2n']:+.5f} "
                f"(pred {r['predicted_nlogn_coef']:.5f} if n^2 coef is 1/16; "
                f"pred n^2 coef {r['predicted_n2_coef']:.6f})",
                flush=True,
            )
    result["pure_types"] = pures

    if args.stats.exists():
        stats = json.loads(args.stats.read_text())
        rows = stats["rows"]
        print("\n== Prediction B: local A from second differences ==")
        locals_all = {}
        for h in (8, 16, 32):
            loc = local_coefficients(rows, h)
            if not loc:
                continue
            locals_all[str(h)] = loc
            print(f"  h={h:3d}  A_local: n={loc[0]['n']} -> {loc[0]['A_local']:.5f}   "
                  f"n={loc[-1]['n']} -> {loc[-1]['A_local']:.5f}")
        result["local_coefficients"] = locals_all

        print("\n== Extrapolating A_local to n -> infinity, scored out of sample ==")
        loc16 = locals_all.get("16", [])
        hv = holdout_extrapolation(loc16)
        result["holdout"] = hv
        if "error" in hv:
            print("  " + hv["error"])
        else:
            print(
                f"  fit on n in [{hv['train_range'][0]:.0f},{hv['train_range'][1]:.0f}] "
                f"({hv['n_train']} pts); scored on n in "
                f"[{hv['test_range'][0]:.0f},{hv['test_range'][1]:.0f}] "
                f"({hv['n_holdout']} pts never fitted)"
            )
            print(f"    {'basis':16s} {'A_inf':>9s} {'train':>10s} {'HELD-OUT':>10s}")
            for r in hv["ranked_bases"]:
                print(
                    f"    {r['basis']:16s} {r['A_inf']:9.5f} "
                    f"{r['train_maxres']:10.2e} {r['holdout_maxres']:10.2e}"
                )
            b = hv["best_free_exponent"]
            print(
                f"  free exponent selected by held-out error: p={b['p']:.2f}  "
                f"A_inf={b['A_inf']:.5f}  (predicted 0.75, error "
                f"{hv['best_free_exponent_error_vs_0.75']:.5f}, "
                f"held-out maxres {b['holdout_maxres']:.2e})"
            )
    else:
        print(f"\n(no {args.stats}; run asymptotics.py first for predictions B)")

    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "orbit_model.json").write_text(json.dumps(result, indent=1))
    print(f"\nwrote {args.outdir / 'orbit_model.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
