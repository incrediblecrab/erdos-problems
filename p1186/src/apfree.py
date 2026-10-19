"""Exact decision procedure: does Z_m admit a 2-colouring with no monochromatic k-AP?

Such a colouring, used periodically on {1,...,n} (colour i by c[i mod m]),
certifies

    delta_k <= 1 / (2(k-1) m),

because the only monochromatic progressions left are the degenerate ones with
common difference divisible by m, and those have density 1/(2(k-1)m).
See exact.zm_weights_int for the derivation of the normalisation.

The question "is Z_m k-AP-free 2-colourable?" is a tiny SAT instance: one
Boolean per residue, and for every progression two clauses saying it is not
all-0 and not all-1.  A solver therefore settles it either way -- a model is an
explicit colouring, and UNSAT is a proof that no colouring exists.  That is
strictly stronger than the tabu search in periodic.py, which can only find.

    python src/apfree.py 5 2-150

Every SAT model returned here is re-checked by a plain independent loop before
being reported, so a solver bug cannot manufacture a bound.
"""

import argparse
import json
import os
import sys

from pysat.formula import CNF
from pysat.solvers import Cadical153

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def progressions(m, k):
    """Distinct index sets of k-APs in Z_m with nonzero common difference."""
    seen = set()
    for a in range(m):
        for b in range(1, m):
            U = frozenset((a + j * b) % m for j in range(k))
            if len(U) > 1:
                seen.add(U)
    return sorted(seen, key=lambda s: (len(s), sorted(s)))


def check(word, k):
    """Independent re-check: count monochromatic k-APs with b != 0 in Z_m."""
    c = [int(x) for x in word]
    m = len(c)
    bad = 0
    for a in range(m):
        for b in range(1, m):
            if len({c[(a + j * b) % m] for j in range(k)}) == 1:
                bad += 1
    return bad


def decide(m, k, conf_budget=0):
    """Return (status, word). status in {'SAT', 'UNSAT', 'UNKNOWN'}."""
    cnf = CNF()
    for U in progressions(m, k):
        lits = [i + 1 for i in sorted(U)]
        cnf.append(lits)                      # not all 0
        cnf.append([-x for x in lits])        # not all 1
    cnf.append([-1])                          # colour swap symmetry: c[0] = 0
    with Cadical153(bootstrap_with=cnf) as s:
        if conf_budget > 0:
            s.conf_budget(conf_budget)
            res = s.solve_limited()
        else:
            res = s.solve()
        if res is None:
            return "UNKNOWN", None
        if not res:
            return "UNSAT", None
        model = set(x for x in s.get_model() if x > 0)
        word = "".join("1" if (i + 1) in model else "0" for i in range(m))
    bad = check(word, k)
    assert bad == 0, "solver model is not AP-free: %d violations" % bad
    return "SAT", word


def parse_range(s):
    out = []
    for part in s.split(","):
        if "-" in part:
            a, b = part.split("-")
            out.extend(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("k", type=int)
    ap.add_argument("ms")
    ap.add_argument("--budget", type=int, default=0,
                    help="conflict budget per instance; 0 = unlimited")
    ap.add_argument("--store", default="")
    args = ap.parse_args()

    data = {}
    if args.store and os.path.exists(args.store):
        with open(args.store) as fh:
            data = json.load(fh)

    best = None
    for m in parse_range(args.ms):
        status, word = decide(m, args.k, args.budget)
        if status == "SAT":
            best = m
            data[str(m)] = {"m": m, "k": args.k, "word": word,
                            "bound_num": 1, "bound_den": 2 * (args.k - 1) * m,
                            "bound": 1.0 / (2 * (args.k - 1) * m)}
            print("m=%3d  SAT    delta_%d <= 1/%-6d = %.9f  %s"
                  % (m, args.k, 2 * (args.k - 1) * m,
                     1.0 / (2 * (args.k - 1) * m), word), flush=True)
        else:
            data[str(m)] = {"m": m, "k": args.k, "status": status}
            print("m=%3d  %s" % (m, status), flush=True)
        if args.store:
            tmp = args.store + ".tmp"
            with open(tmp, "w") as fh:
                json.dump(data, fh, indent=1, sort_keys=True)
            os.replace(tmp, args.store)

    if best:
        print("largest AP-free m in this range: %d  ->  delta_%d <= 1/%d"
              % (best, args.k, 2 * (args.k - 1) * best), flush=True)


if __name__ == "__main__":
    main()
