#!/usr/bin/env python
"""Export an exact integer weight table for (m, k) so that src/search can optimise it.

Two settings:

  zp   the finite-field problem.  Weights come from exact.block_weights_int and
       give Phi_k(c), with  delta~_k <= Phi_k(c)/2.
  n    the {1,...,n} problem.  Weights come from exact.triangle_weights_int and
       give Psi_k(c), with  delta_k <= Psi_k(c).

Every weight is an integer and the objective is (matched sum)/DEN, with no
rounding anywhere.  search.c only *searches* this table; whatever word it
returns is re-certified from scratch by src/certify.py.

Usage:  export_weights.py <m> <k> <out-file> [--setting zp|n]
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exact import block_weights_int, triangle_weights_int  # noqa: E402


def export(m, k, path, setting="zp"):
    if setting == "zp":
        den, items = block_weights_int(m, k)
    elif setting == "n":
        den, items = triangle_weights_int(m, k)
    else:
        raise ValueError(setting)
    with open(path, "w") as f:
        f.write(f"{m} {k} {len(items)} {den}\n")
        for idx, w in items:
            f.write(f"{len(idx)} {w} " + " ".join(map(str, idx)) + "\n")
    return len(items), den


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("m", type=int)
    ap.add_argument("k", type=int)
    ap.add_argument("out")
    ap.add_argument("--setting", choices=("zp", "n"), default="zp")
    a = ap.parse_args()
    n, den = export(a.m, a.k, a.out, a.setting)
    print(f"m={a.m} k={a.k} setting={a.setting} sets={n} den={den} -> {a.out}")
