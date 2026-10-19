"""Assemble the wildcard-recursion bases behind the claims into results/.

Each base is a colouring of Z_b with a safe coset F = r + (b/f)Z_b (NOTES.md
section 3).  The words come from two sources:

  * Rabung's quadratic-residue colourings (src/rabung.c), with F = {0} for a
    prime and F = {0, p} for a zipped Z_2p.  Lu and Peng's B_11 and B_37 are
    exactly these for p = 11 and p = 37.
  * the tabu-search optima already stored in results/zm_k*.json (B_44, B_86).

Every stored word is cyclic k-AP-free as a whole (so H1 holds for any F), and
F is re-checked with cosetprobe.safe_coset before anything is written.

    python src/bases.py            # writes results/safebase_k*.json
"""

import json
import os
import subprocess
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from apfree import check                              # noqa: E402
from cosetprobe import best_word, safe_coset          # noqa: E402

# (k, b, f, r, source); source is ("rabung", p, zipped) or ("zm", m)
BASES = [
    (4, 11, 1, 0, ("rabung", 11, False)),
    (5, 37, 1, 0, ("rabung", 37, False)),
    (5, 44, 4, 3, ("zm", 44)),
    (6, 86, 2, 33, ("zm", 86)),
    (6, 226, 2, 0, ("rabung", 113, True)),
    (7, 617, 1, 0, ("rabung", 617, False)),
    (8, 1642, 2, 0, ("rabung", 821, True)),
]


def word_of(k, src):
    if src[0] == "rabung":
        args = [os.path.join(HERE, "rabung"), "word", str(src[1])]
        if src[2]:
            args.append("zip")
        return subprocess.run(args, capture_output=True, text=True,
                              check=True).stdout.strip(), \
            "rabung word %d%s" % (src[1], " zip" if src[2] else "")
    got = best_word(k, src[1])
    assert got is not None, "no stored word at m=%d for k=%d" % (src[1], k)
    return "".join(str(x) for x in got[1]), "results/zm_k%d.json m=%d" % (k, src[1])


def main():
    out = {}
    for k, b, f, r, src in BASES:
        w, where = word_of(k, src)
        c = [int(x) for x in w]
        assert len(c) == b
        F = {(r + (b // f) * y) % b for y in range(f)}
        bad = check(w, k)
        ok, nd = safe_coset(c, k, F)
        assert bad == 0, "k=%d b=%d: %d monochromatic k-APs" % (k, b, bad)
        assert ok, "k=%d b=%d: coset unsafe (%d dangerous APs)" % (k, b, nd)
        bound = Fraction(1, 2 * (k - 1) * (b + f))
        out.setdefault(k, {})["%d/%d" % (b, f)] = {
            "k": k, "b": b, "f": f, "r": r, "word": w, "source": where,
            "bound": str(bound), "bound_float": float(bound)}
        print("k=%d b=%4d f=%d r=%2d  AP-free, coset safe  ->  delta_%d <= %s = %.10f  [%s]"
              % (k, b, f, r, k, bound, float(bound), where))
    for k, recs in sorted(out.items()):
        path = os.path.join(ROOT, "results", "safebase_k%d.json" % k)
        with open(path, "w") as fh:
            json.dump(recs, fh, indent=1, sort_keys=True)
        print("wrote", os.path.relpath(path, ROOT))


if __name__ == "__main__":
    main()
