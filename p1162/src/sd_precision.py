"""High-precision check of the rank standard deviation against its closed form.

`asymptotics.py` stores the standard deviation of the elementary abelian rank distribution as a double in `data/elemab_stats.json`, computed through floating-point logarithms and a variance that cancels about four digits, so it is far noisier than the actual convergence from n=256 on (NOTES.md section 3): read off the doubles, the deviation from the closed form looks like 2.8e-10 at n=256 and 1.2e-09 at n=512, and non-monotone.  Both figures are rounding noise in the artifact, not convergence error.

This module recomputes the variance from the exact integer rank counts in
mpmath at 60 decimal places and compares it with

    sd = sqrt( sum_j j^2 2^{-j^2} / sum_j 2^{-j^2} )
       = 0.849305961022908280743697585694...

The true agreement is 16 decimal places at n=256 and 35 at n=512.  The error
roughly squares at each doubling of n, so it decays like 2^{-cn}, c ~ 0.23 --
not polynomially.

This is convergence of an exact finite quantity to a closed form at four
values of n.  It is strong numerical evidence and it is not a proof.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

import mpmath as mp

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from elemab import rank_counts  # noqa: E402

DPS = 60


def closed_form() -> mp.mpf:
    """sqrt( sum_{j in Z} j^2 2^{-j^2} / sum_{j in Z} 2^{-j^2} )."""
    theta = mp.nsum(lambda j: mp.mpf(2) ** (-j * j), [-mp.inf, mp.inf])
    second = mp.nsum(lambda j: j * j * mp.mpf(2) ** (-j * j), [-mp.inf, mp.inf])
    return mp.sqrt(second / theta)


def exact_sd(n: int) -> mp.mpf:
    """Standard deviation of the rank distribution, from exact integer counts."""
    counts = [mp.mpf(c) for c in rank_counts(n)]
    total = sum(counts)
    mean = sum(k * counts[k] for k in range(len(counts))) / total
    var = sum((k - mean) ** 2 * counts[k] for k in range(len(counts))) / total
    return mp.sqrt(var)


def decimal_places(diff: mp.mpf) -> int:
    """Largest d with |diff| < 0.5 * 10^-d."""
    if diff == 0:
        return DPS
    d = 0
    while diff < mp.mpf("0.5") * mp.mpf(10) ** (-(d + 1)) and d < DPS:
        d += 1
    return d


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, nargs="*", default=[64, 128, 256, 512])
    ap.add_argument("--out", type=pathlib.Path, default=None)
    args = ap.parse_args()

    mp.mp.dps = DPS
    target = closed_form()
    print(f"closed form (30 dp): {mp.nstr(target, 30)}")
    print()

    rows = []
    prev_log = None
    for n in sorted(args.n):
        sd = exact_sd(n)
        diff = abs(sd - target)
        dp = decimal_places(diff)
        log10 = float(mp.log10(diff)) if diff > 0 else float("-inf")
        ratio = log10 / prev_log if prev_log not in (None, 0) else None
        prev_log = log10
        rows.append(
            {
                "n": n,
                "sd_30dp": mp.nstr(sd, 30),
                "abs_error": mp.nstr(diff, 3),
                "decimal_places": dp,
                "log10_error": log10,
                "log10_ratio_vs_prev": ratio,
            }
        )
        r = f"{ratio:.2f}" if ratio is not None else "  --"
        print(
            f"  n={n:4d}  sd={mp.nstr(sd, 30):<34s}"
            f"  |err|={mp.nstr(diff, 3):>10s}  {dp:2d} dp  log10 ratio {r}"
        )

    print()
    print("The log10 ratio settling near 2 means the error squares at each")
    print("doubling of n, i.e. it decays like 2^(-cn), not polynomially.")

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(
            json.dumps(
                {"dps": DPS, "closed_form_30dp": mp.nstr(target, 30), "rows": rows},
                indent=2,
            )
            + "\n"
        )
        print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
