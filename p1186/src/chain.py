"""Build and certify one level of the wildcard chain, for any (k, b, f, r).

Positions of Z_{bt} lying over the free coset F = {r + (b/f)y} of Z_b are
exactly {r + step*z : z in Z_{ft}} with step = b/f, so the construction is

    c(x) = BASE[x mod b]          if x % step != r
    c(r + step*z) = D[z]          for a colouring D of Z_{ft}

where D is supplied (here: the k-AP-free base itself, giving the deepest level
that is still cheap to certify).  Phi is then certified exactly and compared
against the recursion ((b-f) + f^2 * m_k(Z_{ft})) / b^2.

    python src/chain.py --k 6 --b 86 --f 2 --r 33
"""

import argparse
import json
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exact                                          # noqa: E402
from cosetprobe import best_word, safe_coset          # noqa: E402


def build(base, f, r, D):
    b = len(base)
    step = b // f
    t = len(D) // f
    m = b * t
    assert len(D) == f * t
    word = [base[x % b] for x in range(m)]
    nfree = 0
    for z in range(f * t):
        word[(r + step * z) % m] = D[z]
        nfree += 1
    assert nfree == f * t
    # every position not congruent to r mod step must still carry the base
    for x in range(m):
        if x % step != r % step:
            assert word[x] == base[x % b], "construction moved a fixed bit"
    return word


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--b", type=int, required=True)
    ap.add_argument("--f", type=int, required=True)
    ap.add_argument("--r", type=int, required=True)
    ap.add_argument("--store", default=None)
    a = ap.parse_args()

    got = best_word(a.k, a.b)
    assert got is not None, "no stored base word at m=%d" % a.b
    phi_base, base = got
    step = a.b // a.f
    F = {(a.r + step * y) % a.b for y in range(a.f)}
    ok, nd = safe_coset(base, a.k, F)
    print("base   m=%d  Phi=%s=%.10f   coset f=%d r=%d  safe=%s (%d dangerous APs)"
          % (a.b, phi_base, float(phi_base), a.f, a.r, ok, nd))
    assert ok, "refusing to build a chain on an unsafe coset"

    D = list(base)                                    # inner colouring = the base
    word = build(base, a.f, a.r, D)
    m = len(word)
    phi = exact.phi_zm(word, a.k)
    psi = phi / (2 * (a.k - 1))

    pred = Fraction((a.b - a.f) + a.f * a.f * phi_base, a.b * a.b)
    fixed = Fraction(1, a.b + a.f)
    print("chain  m=%d  Phi=%s=%.10f" % (m, phi, float(phi)))
    print("       recursion predicts %s=%.10f   %s"
          % (pred, float(pred), "MATCH" if phi == pred else "MISMATCH"))
    print("       delta_%d <= %s = %.10f" % (a.k, psi, float(psi)))
    print("       chain fixed point 1/(b+f) = %s -> delta_%d <= %s = %.10f (limit, NOTES.md section 3.1)"
          % (fixed, a.k, fixed / (2 * (a.k - 1)), float(fixed / (2 * (a.k - 1)))))
    assert phi == pred, "recursion and explicit construction disagree"

    if a.store:
        rec = {"k": a.k, "m": m, "b": a.b, "f": a.f, "r": a.r,
               "phi": str(phi), "psi": str(psi), "psi_float": float(psi),
               "word": "".join(str(v) for v in word), "ap_free": False}
        data = {}
        if os.path.exists(a.store):
            data = json.load(open(a.store))
        data[str(m)] = rec
        json.dump(data, open(a.store, "w"), indent=1, sort_keys=True)
        print("stored -> %s" % a.store)


if __name__ == "__main__":
    main()
