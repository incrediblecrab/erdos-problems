"""Asymptotics of the exact elementary abelian counts.

Three questions, all answered from exact integers produced by `elemab.py`:

  1. Where does the rank distribution of a random elementary abelian
     2-subgroup of S_n sit?  (This is the candidate "statistical theorem on
     their order" of Erdos #1162: rank = log_2 of the order.)

  2. How big is log_2 e(n) - n^2/16, where e(n) is the total number of
     elementary abelian 2-subgroups?  RoTr25 Theorem 1 brackets the true
     log_2 |Sub(S_n)| - n^2/16 between alpha*n*log n and beta*n^{3/2}.

  3. How much of that excess is already present in a *single* maximal
     elementary abelian subgroup E = (Z/2)^{floor(n/2)} -- i.e. in the Galois
     number G_{floor(n/2)} = |Sub(E)|, the classical Pyber lower bound?

All logarithms are to base 2, matching RoTr25 ("our logarithms are to the
base 2", p. 1).
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from elemab import galois_number, hom_table, rank_counts  # noqa: E402


def log2_int(x: int) -> float:
    """log2 of an arbitrarily large positive integer, to ~15 significant digits."""
    if x <= 0:
        raise ValueError("log2 of non-positive integer")
    b = x.bit_length()
    if b <= 900:
        return math.log2(x)
    shift = b - 900
    return shift + math.log2(x >> shift)


def ratio(num: int, den: int) -> float:
    """num/den as a float, safe for integers far beyond the float range."""
    if num == 0:
        return 0.0
    sign = -1.0 if (num < 0) != (den < 0) else 1.0
    num, den = abs(num), abs(den)
    return sign * math.exp2(log2_int(num) - log2_int(den))


def rank_stats(n: int, counts: list[int]) -> dict:
    """Summary of the distribution k -> a(n,k)/e(n)."""
    e = sum(counts)
    mode = max(range(len(counts)), key=lambda k: counts[k])
    s1 = sum(k * c for k, c in enumerate(counts))
    s2 = sum(k * k * c for k, c in enumerate(counts))
    mean = ratio(s1, e)
    var = ratio(s2, e) - mean * mean
    # smallest symmetric window around the mode holding >= 1 - 1e-6 of the mass
    width = 0
    while True:
        lo, hi = mode - width, mode + width
        mass = sum(counts[max(0, lo) : hi + 1])
        if ratio(e - mass, e) < 1e-6 or width > len(counts):
            break
        width += 1
    return {
        "n": n,
        "max_rank": len(counts) - 1,
        "log2_e": log2_int(e),
        "excess": log2_int(e) - n * n / 16.0,
        "mode": mode,
        "mode_over_n": mode / n if n else 0.0,
        "mean": mean,
        "mean_over_n": mean / n if n else 0.0,
        "sd": math.sqrt(max(var, 0.0)),
        "window_1e6": width,
        "e_bits": e.bit_length(),
    }


def galois_stats(r_values: list[int]) -> list[dict]:
    out = []
    for r in r_values:
        g = galois_number(r)
        out.append(
            {
                "r": r,
                "log2_G": log2_int(g),
                "excess_over_r2_4": log2_int(g) - r * r / 4.0,
                "parity": "even" if r % 2 == 0 else "odd",
            }
        )
    return out


def fit_excess(rows: list[dict], n_min: int) -> dict:
    """Least squares of excess(n) against [n log2 n, n, 1]."""
    import numpy as np

    sel = [r for r in rows if r["n"] >= n_min]
    A = np.array([[r["n"] * math.log2(r["n"]), r["n"], 1.0] for r in sel])
    y = np.array([r["excess"] for r in sel])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = float(np.max(np.abs(A @ coef - y)))
    return {
        "n_min": n_min,
        "n_max": max(r["n"] for r in sel),
        "points": len(sel),
        "coef_n_log2n": float(coef[0]),
        "coef_n": float(coef[1]),
        "coef_const": float(coef[2]),
        "max_abs_residual": resid,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-max", type=int, default=512)
    ap.add_argument("--dense-to", type=int, default=64)
    ap.add_argument("--step", type=int, default=16)
    ap.add_argument("--outdir", type=pathlib.Path, default=pathlib.Path("data"))
    args = ap.parse_args()

    schedule = sorted(
        set(range(2, args.dense_to + 1))
        | set(range(args.dense_to, args.n_max + 1, args.step))
        | {args.n_max}
    )

    t0 = time.time()
    hom = hom_table(args.n_max // 2 + 2, args.n_max)
    print(f"hom table to n={args.n_max}: {time.time() - t0:.1f}s", flush=True)

    rows = []
    small_exact: dict[str, list[str]] = {}
    for n in schedule:
        t = time.time()
        counts = rank_counts(n, hom)
        st = rank_stats(n, counts)
        st["seconds"] = round(time.time() - t, 2)
        rows.append(st)
        if n <= 40:
            small_exact[str(n)] = [str(c) for c in counts]
        print(
            f"n={n:4d} log2 e={st['log2_e']:12.3f} excess={st['excess']:+10.3f} "
            f"mode={st['mode']:4d} mode/n={st['mode_over_n']:.5f} "
            f"mean/n={st['mean_over_n']:.5f} sd={st['sd']:.3f} "
            f"win={st['window_1e6']:3d} ({st['seconds']}s)",
            flush=True,
        )

    gal = galois_stats(
        sorted({r for p in range(2, 10) for r in (2**p - 1, 2**p, 2**p + 1)}
               | set(range(2, 65)))
    )

    fits = [fit_excess(rows, m) for m in (32, 64, 128, 256) if m < args.n_max]

    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "elemab_stats.json").write_text(
        json.dumps({"rows": rows, "fits": fits}, indent=1)
    )
    (args.outdir / "rank_counts_small.json").write_text(json.dumps(small_exact))
    (args.outdir / "galois.json").write_text(json.dumps(gal, indent=1))

    print("\nfits of  log2 e(n) - n^2/16  ~  A n log2 n + B n + C")
    for f in fits:
        print(
            f"  n in [{f['n_min']},{f['n_max']}]  A={f['coef_n_log2n']:.5f} "
            f"B={f['coef_n']:+.4f} C={f['coef_const']:+.3f} "
            f"maxres={f['max_abs_residual']:.4f}"
        )

    print("\nGalois numbers: log2 G_r - r^2/4")
    for r in (64, 128, 256, 512):
        row = next(g for g in gal if g["r"] == r)
        print(f"  r={r:5d} ({row['parity']:4s})  {row['excess_over_r2_4']:+.9f}")
    for r in (63, 127, 255, 511):
        row = next(g for g in gal if g["r"] == r)
        print(f"  r={r:5d} ({row['parity']:4s})  {row['excess_over_r2_4']:+.9f}")
    print(f"\ntotal {time.time() - t0:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
