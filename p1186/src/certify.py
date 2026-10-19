#!/usr/bin/env python
"""Certify a block colouring: recompute its objective and bound in exact arithmetic.

Takes a binary word (as produced by src/search) and recomputes everything from
scratch through exact.py, in rational arithmetic.  Nothing reported by the C
optimiser is trusted; this is the certification path.

  --setting zp   reports Phi_k and the bound delta~_k <= Phi_k/2   (finite field)
  --setting n    reports Psi_k and the bound delta_k  <= Psi_k     ({1,...,n})

Usage:  certify.py <word> <k> [--setting zp|n] [--json]
"""
import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exact import phi, psi, lambda_form, mean  # noqa: E402


def certify(word, k, setting="zp"):
    c = [int(ch) for ch in word]
    assert set(c) <= {0, 1}, "word must be binary"
    m = len(c)
    r = {"m": m, "k": k, "setting": setting, "word": word, "mu": str(mean(c))}
    if setting == "zp":
        v = phi(c, k)
        bound = v / 2
        rand = Fraction(1, 2 ** k)
        r["objective"] = "Phi"
        if k <= 6:
            lam = lambda_form(c, k)
            r["lambda"] = str(lam)
            r["lambda_float"] = float(lam)
            if k == 3:
                r["identity_ok"] = bool(v == (1 + 3 * mean(c) ** 2) / 4)
            elif k == 4:
                r["identity_ok"] = bool(v == (1 + 6 * mean(c) ** 2 + lam) / 8)
    else:
        v = psi(c, k)
        bound = v
        rand = Fraction(1, (k - 1) * 2 ** k)
        r["objective"] = "Psi"
    r.update({
        "value": str(v), "value_float": float(v),
        "bound": str(bound), "bound_float": float(bound),
        "random_bound": str(rand), "random_bound_float": float(rand),
        "beats_random": bool(bound < rand),
    })
    return r


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("word")
    ap.add_argument("k", type=int)
    ap.add_argument("--setting", choices=("zp", "n"), default="zp")
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args()
    r = certify(a.word, a.k, a.setting)
    if a.json:
        print(json.dumps(r))
    else:
        sym = "delta~_%d" % a.k if a.setting == "zp" else "delta_%d" % a.k
        print(f"m={r['m']} k={a.k} setting={a.setting} word={a.word}")
        print(f"  {r['objective']}_{a.k:<10s} = {r['value']} = {r['value_float']:.12f}")
        print(f"  {sym} <= {' ' * max(0, 8 - len(sym))}{r['bound']} = {r['bound_float']:.12f}")
        print(f"  random bound   = {r['random_bound']} = {r['random_bound_float']:.12f}")
        print(f"  beats random   = {r['beats_random']}")
        if "lambda" in r:
            print(f"  Lambda_{a.k}       = {r['lambda']} = {r['lambda_float']:.12f}")
        if "identity_ok" in r:
            print(f"  identity check = {r['identity_ok']}")
