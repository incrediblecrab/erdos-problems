#!/usr/bin/env python3
"""Numeric checks of the constants in the Theorem B reconstruction (NOTES.md, section 4).

The proofs in NOTES.md are short analytic arguments; this script only guards against a slip in the constants they state. It checks, in floating point with a reported margin:

1. the sum bound S(n, α) = Σ_{a=1}^{n-1} (a(n-a))^{-α} ≤ 2^{2α} n^{1-2α} / (1-α), for 1/2 < α < 1 and every 2 ≤ n ≤ N_MAX;
2. the mean bound μ_n = Σ_{1 ≤ a < n/2} (a(n-a))^{-1/2-ε} ≤ 2^{2ε} n^{-2ε} / (1/2 - ε), for every 3 ≤ n ≤ N_MAX;
3. the table of k(ε) = ⌊1/(2ε)⌋ + 1 and the resulting eventual bound 2k - 1 on the ordered representation count, with the exponent check 2εk > 1 done in exact rational arithmetic.

It also prints the ratio 2^{2α} / ((1-α) B(1-α, 1-α)) that bound 1 approaches as n → ∞, where B is the Beta function, to show how much slack the constant has. Nothing here proves an asymptotic statement.

Usage: python theorem_b_check.py [--nmax N] [--out DIR] [--plant NAME]. Writes theorem_b_check.json to DIR (default ../results). Exit 0 means every checked inequality held with the float64 margin below.

--plant swaps in a deliberately wrong constant and writes nothing, to show the check can fail: "denominators" drops the 1/(1-α) and 1/(1/2-ε) factors from bounds 1 and 2, and "k" uses k = ⌊1/(2ε)⌋, which makes 2εk = 1 at ε = 1/200. Each plant must exit 1.
"""
import argparse
import hashlib
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ALPHAS = [0.505, 0.51, 0.55, 0.6, 0.75, 0.9, 0.99]
EPSILONS = [Fraction(1, 200), Fraction(1, 100), Fraction(1, 20), Fraction(1, 10), Fraction(1, 5), Fraction(1, 4), Fraction(2, 5), Fraction(49, 100)]
# The relative float64 error of a sum of at most N_MAX positive terms is below N_MAX * 2**-52 < 5e-12 for N_MAX = 20000. An inequality counts as verified only if it holds with a margin of 1e-10, twenty times that.
MARGIN = 1e-10


PLANTS = ("denominators", "k")


def sum_bound_check(alpha: float, nmax: int, plant: str | None = None) -> dict:
    worst_ratio, worst_n = math.inf, None
    for n in range(2, nmax + 1):
        a = np.arange(1, n, dtype=np.float64)
        lhs = float(np.sum((a * (n - a)) ** (-alpha)))
        rhs = 2 ** (2 * alpha) * n ** (1 - 2 * alpha) / (1 if plant == "denominators" else 1 - alpha)
        ratio = rhs / lhs
        if ratio < worst_ratio:
            worst_ratio, worst_n = ratio, n
    limit = 2 ** (2 * alpha) / ((1 - alpha) * math.exp(2 * math.lgamma(1 - alpha) - math.lgamma(2 - 2 * alpha)))
    return {"alpha": alpha, "min_rhs_over_lhs": worst_ratio, "at_n": worst_n, "asymptotic_ratio": limit, "holds": worst_ratio > 1 + MARGIN}


def mean_bound_check(eps: Fraction, nmax: int, plant: str | None = None) -> dict:
    e = float(eps)
    worst_ratio, worst_n = math.inf, None
    for n in range(3, nmax + 1):
        a = np.arange(1, (n + 1) // 2, dtype=np.float64)  # 1 <= a < n/2
        mu = float(np.sum((a * (n - a)) ** (-0.5 - e)))
        rhs = 2 ** (2 * e) * n ** (-2 * e) / (1 if plant == "denominators" else 0.5 - e)
        ratio = rhs / mu
        if ratio < worst_ratio:
            worst_ratio, worst_n = ratio, n
    return {"epsilon": str(eps), "min_rhs_over_mu": worst_ratio, "at_n": worst_n, "holds": worst_ratio > 1 + MARGIN}


def k_table(plant: str | None = None) -> list:
    rows = []
    for eps in EPSILONS:
        k = math.floor(1 / (2 * eps)) + (0 if plant == "k" else 1)
        rows.append({"epsilon": str(eps), "k": k, "two_eps_k": str(2 * eps * k), "summable": 2 * eps * k > 1, "eventual_bound_on_sumRep": 2 * k - 1, "density_exponent": str(Fraction(1, 2) - eps)})
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--nmax", type=int, default=20000)
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent.parent / "results")
    ap.add_argument("--plant", choices=PLANTS, help="use a deliberately wrong constant; nothing is written")
    args = ap.parse_args()
    if args.nmax > 20000:
        sys.exit("--nmax above 20000 would exceed the float64 margin argued in the header")
    sums = [sum_bound_check(a, args.nmax, args.plant) for a in ALPHAS]
    means = [mean_bound_check(e, args.nmax, args.plant) for e in EPSILONS]
    table = k_table(args.plant)
    for r in sums:
        print(f"sum bound   alpha={r['alpha']:<6} min RHS/LHS = {r['min_rhs_over_lhs']:.6f} at n={r['at_n']:<6} asymptotic ratio {r['asymptotic_ratio']:.6f}  {'ok' if r['holds'] else 'FAILS'}")
    for r in means:
        print(f"mean bound  eps={r['epsilon']:<6} min RHS/mu  = {r['min_rhs_over_mu']:.6f} at n={r['at_n']:<6}  {'ok' if r['holds'] else 'FAILS'}")
    for r in table:
        print(f"k table     eps={r['epsilon']:<6} k={r['k']:<4} 2*eps*k={r['two_eps_k']:<8} sumRep eventually <= {r['eventual_bound_on_sumRep']:<4} density exponent {r['density_exponent']}  {'ok' if r['summable'] else 'FAILS'}")
    ok = all(r["holds"] for r in sums) and all(r["holds"] for r in means) and all(r["summable"] for r in table)
    if args.plant:
        print(f"plant {args.plant}: " + ("NOT CAUGHT, every check still holds" if ok else "caught"))
        return 0 if ok else 1
    args.out.mkdir(parents=True, exist_ok=True)
    script_sha256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (args.out / "theorem_b_check.json").write_text(json.dumps({"nmax": args.nmax, "float64_margin": MARGIN, "sum_bound": sums, "mean_bound": means, "k_table": table, "all_hold": ok, "script_sha256": script_sha256}, indent=2) + "\n", encoding="utf-8")
    print("all checks hold" if ok else "SOME CHECK FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
