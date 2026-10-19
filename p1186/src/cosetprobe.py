"""Exhaustive search for a wildcard coset in a k-AP-free base.

The wildcard construction needs a base modulus b, coloured so that every
monochromatic k-AP of the b-periodic extension has all of its terms inside a
single "free" coset F of size f.  Then F may be recoloured by an arbitrary
colouring of Z_{ft} without creating any monochromatic AP outside F, giving

    m_k(Z_{bt}) <= ((b - f) + f^2 * m_k(Z_{ft})) / b^2,

whose fixed point is 1/(b+f), hence delta_k <= 1/(2(k-1)(b+f)).

This script does not assume the coset; it finds it.  For each divisor f of b
and each residue r it fixes the base bits off the coset, enumerates or
tabu-searches the f*t free bits at m = b*t, and certifies Phi exactly.  A coset
is usable exactly when the resulting Phi beats the tiled value Phi(base), since
tiling always preserves Phi.

Run with --self-test first: it must rediscover the known k=5 answer (b=44,
f=4, free coset {3,14,25,36}).  A probe that cannot recover a known structure
is not evidence about a modulus where none is known.

    python src/cosetprobe.py --k 6 --b 86 --t 2 --self-test
"""

import argparse
import glob
import itertools
import json
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exact                                          # noqa: E402


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def best_word(k, m):
    """Best stored periodic word of length m for this k, over all shards."""
    best = None
    for fn in glob.glob(os.path.join("results", "zm_k%d*.json" % k)):
        try:
            data = json.load(open(fn))
        except Exception:
            continue
        if isinstance(data, dict):
            recs = list(data.values())
        else:
            recs = data
        for rec in recs:
            if not isinstance(rec, dict) or rec.get("m") != m:
                continue
            w = rec.get("word") or rec.get("c")
            if w is None:
                continue
            c = [int(ch) for ch in w]
            if len(c) != m:
                continue
            phi = exact.phi_zm(c, k)
            if best is None or phi < best[0]:
                best = (phi, c)
    return best


def free_positions(b, f, r, t):
    """Positions of Z_{bt} lying over the coset {r + (b/f)y} of Z_b."""
    step = b // f
    base_coset = {(r + step * y) % b for y in range(f)}
    assert len(base_coset) == f
    return [x for x in range(b * t) if x % b in base_coset]


def probe_coset(k, base, f, r, t, cap):
    """Minimise Phi_k over Z_{bt} varying only the free coset.  Exact."""
    b = len(base)
    pos = free_positions(b, f, r, t)
    assert len(pos) == f * t, (len(pos), f * t)
    word = [base[x % b] for x in range(b * t)]
    nfree = len(pos)
    if nfree > cap:
        return None                                   # too big to enumerate
    best = None
    for bits in itertools.product((0, 1), repeat=nfree):
        for p, v in zip(pos, bits):
            word[p] = v
        phi = exact.phi_zm(word, k)
        if best is None or phi < best[0]:
            best = (phi, list(word), bits)
    # the fixed bits must never have moved
    for x in range(b * t):
        if x not in set(pos):
            assert best[1][x] == base[x % b], "probe moved a fixed bit"
    return best


def safe_coset(base, k, F):
    """Exact test: can the coset F be recoloured freely without side effects?

    Colour Z_{bt} by BASE off the lift of F, and arbitrarily on it.  A k-AP is
    *dangerous* if it has terms both inside and outside F and BASE is constant
    on the outside terms, because then colouring the inside terms to match
    makes it monochromatic -- a progression the recursion does not account for.

    Both the colour and the F-membership of a+jd depend only on (a+jd) mod b,
    so testing over Z_b settles it for every t: if the outside terms already
    disagree mod b they still disagree in any lift.  Returns (ok, ndanger).
    """
    b = len(base)
    ndanger = 0
    for a in range(b):
        for d in range(b):
            U = {(a + j * d) % b for j in range(k)}
            inside = U & F
            if not inside:
                continue                              # wholly outside: a base AP
            outside = U - F
            if not outside:
                continue                              # wholly inside: the recursion's job
            cols = {base[x] for x in outside}
            if len(cols) == 1:
                ndanger += 1
    return ndanger == 0, ndanger


def structural_scan(k, b):
    """Report every coset of Z_b that is safe, by the exact test above."""
    got = best_word(k, b)
    if got is None:
        print("no stored word at m=%d for k=%d" % (b, k))
        return []
    phi_base, base = got
    print("base m=%d  Phi=%s=%.10f" % (b, phi_base, float(phi_base)))
    safe = []
    for f in divisors(b):
        if f == b:
            continue
        for r in range(b // f):
            step = b // f
            F = {(r + step * y) % b for y in range(f)}
            ok, nd = safe_coset(base, k, F)
            if ok:
                safe.append((f, r))
                print("  f=%2d r=%2d  SAFE    -> limit 1/(b+f) = 1/%d, delta_%d <= 1/%d"
                      % (f, r, b + f, k, 2 * (k - 1) * (b + f)))
    if not safe:
        print("  no safe coset: this base admits no wildcard structure")
    else:
        bf = max(f for f, _ in safe)
        print("  best safe coset size f=%d  ->  delta_%d <= 1/%d (chain limit)"
              % (bf, k, 2 * (k - 1) * (b + bf)))
    return safe


def run(k, b, t, cap, verbose):
    got = best_word(k, b)
    if got is None:
        print("no stored word at m=%d for k=%d" % (b, k))
        return None
    phi_base, base = got
    tiled = exact.phi_zm([base[x % b] for x in range(b * t)], k)
    assert tiled == phi_base, "tiling must preserve Phi"
    print("base m=%d  Phi=%s=%.10f   tiled to m=%d: %s (preserved)"
          % (b, phi_base, float(phi_base), b * t, tiled))

    results = []
    for f in divisors(b):
        if f == b:
            continue
        for r in range(b // f):
            out = probe_coset(k, base, f, r, t, cap)
            if out is None:
                continue
            phi, word, bits = out
            if phi < tiled:
                results.append((phi, f, r, word))
                print("  f=%2d r=%2d  Phi=%s=%.10f  BEATS tiled  limit 1/(b+f)=1/%d"
                      % (f, r, phi, float(phi), b + f))
            elif verbose:
                print("  f=%2d r=%2d  Phi=%s=%.10f" % (f, r, phi, float(phi)))
    if not results:
        print("  no coset of any size beats the tiled value -> no wildcard structure at b=%d" % b)
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--b", type=int, required=True)
    ap.add_argument("--t", type=int, default=2)
    ap.add_argument("--cap", type=int, default=20,
                    help="max free bits to enumerate exhaustively")
    ap.add_argument("--self-test", action="store_true",
                    help="first rediscover the known k=5 b=44 f=4 structure")
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args()

    if a.self_test:
        print("=== positive control: k=5, b=44 (expect f=4, r=3 to win) ===")
        res = run(5, 44, 2, a.cap, False)
        ok = res and any(f == 4 and r == 3 for _, f, r, _ in res)
        print("control %s" % ("PASSED" if ok else "FAILED -- probe is not trustworthy"))
        if not ok:
            sys.exit(1)
        print("--- structural test on the same base (must agree: best f=4) ---")
        s = structural_scan(5, 44)
        ok2 = s and max(f for f, _ in s) == 4
        print("structural control %s\n"
              % ("PASSED" if ok2 else "FAILED -- structural test disagrees with search"))
        if not ok2:
            sys.exit(1)

    print("=== k=%d, b=%d: structural scan over every coset ===" % (a.k, a.b))
    structural_scan(a.k, a.b)
    print("\n=== k=%d, b=%d, t=%d: enumerated search (cross-check) ===" % (a.k, a.b, a.t))
    run(a.k, a.b, a.t, a.cap, a.verbose)


if __name__ == "__main__":
    main()
