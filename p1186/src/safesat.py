"""Exact decision procedure for wildcard-recursion bases.

A base is a colouring of Z_b off the subgroup F = (b/f)Z_b (translating a
coset to F is free, so F is always the subgroup here).  It feeds the wildcard
recursion of NOTES.md section 3 when

    (H1) every k-AP of Z_b with nonzero difference and no term in F is
         non-monochromatic, and
    (H2) every k-AP with terms both in and outside F is non-constant on its
         terms outside F (cosetprobe.safe_coset).

Colours on F never enter either condition.  NOTES.md section 3 proves that
(H1)+(H2) give m_k(Z_bt) = ((b-f) + f^2 m_k(Z_ft)) / b^2 for every t and every
colouring of Z_ft, hence

    delta_k <= 1 / (2(k-1)(b+f)).

Both conditions are clauses over the colours outside F, so a SAT solver
settles each (k, b, f) either way: a model is an explicit base, UNSAT is a
proof that none exists.  apfree.py is the special case "F empty".

    python src/safesat.py 5 2-200 --above 48
    python src/safesat.py 3 2-60 --above 4      # must find nothing: PRS08 Thm 4

Every model is re-checked twice before being reported: by violations() below
and by cosetprobe.safe_coset, which was written independently of this file.

--symbreak L adds lex-leader constraints on the first L free positions for
every affine map x -> u*x + s (u a unit mod b, s in F).  Those maps fix F and
permute the clauses, so the constraints only discard symmetric copies; the
answer is unchanged, which the scans in NOTES.md check against the plain run.

    python src/safesat.py 6 221 --above 228 --symbreak 24 --solver kissat404
"""

import argparse
import json
import os
import sys
from math import gcd

from pysat.formula import CNF
from pysat.solvers import Solver

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from apfree import parse_range                        # noqa: E402
from cosetprobe import divisors, safe_coset           # noqa: E402


def subgroup(b, f):
    return set(range(0, b, b // f))


def outside_sets(b, f, k):
    """Distinct outside-term sets that must be non-constant, or None if some
    progression has a single outside term (then no colouring works)."""
    F = subgroup(b, f)
    seen = set()
    for a in range(b):
        for d in range(1, b):
            U = {(a + j * d) % b for j in range(k)}
            O = frozenset(U - F)
            if not O:
                continue
            if len(O) == 1:
                return None
            seen.add(O)
    return seen


def violations(word, k, F):
    """Independent re-check: (#H1 failures, #H2 failures), colours on F ignored."""
    b = len(word)
    h1 = h2 = 0
    for a in range(b):
        for d in range(1, b):
            U = {(a + j * d) % b for j in range(k)}
            O = U - F
            if O and len({word[x] for x in O}) == 1:
                if O == U:
                    h1 += 1
                else:
                    h2 += 1
    return h1, h2


def symmetry_maps(b, f):
    """Images x -> u*x + s (as lists) for u a unit mod b, s in F, not the
    identity.  Each fixes F setwise and sends progressions with nonzero
    difference to progressions with nonzero difference, so it permutes the
    clauses."""
    F = sorted(subgroup(b, f))
    units = [u for u in range(1, b) if gcd(u, b) == 1]
    return [[(u * x + s) % b for x in range(b)]
            for u in units for s in F if (u, s) != (1, 0)]


def lex_leader(cnf, var, free, perms, L):
    """For every permutation sigma: colours of free[:L] <=lex their images.

    The solution set is invariant under every sigma and under colour swap, so
    the lex-least member of any orbit satisfies all of these at once together
    with the colour-swap unit clause on free[0], the first position in the
    same order.  Auxiliary variable e_i means "equal on the first i pairs"."""
    top = len(free)
    for img in perms:
        pairs = [(var[x], var[img[x]]) for x in free[:L]]
        live = [i for i, (x, y) in enumerate(pairs) if x != y]
        prev = 0
        for i in live:
            x, y = pairs[i]
            guard = [-prev] if prev else []
            cnf.append(guard + [-x, y])
            if i == live[-1]:
                break
            top += 1
            cnf.append(guard + [-x, top])
            cnf.append(guard + [y, top])
            prev = top
    return top


def encode_cnf(k, b, f, symbreak=0, perms=None):
    """Return the SAT encoding for bases as (cnf, free).

    The clauses are exactly the ones used by decide(): for every outside-term
    set O, at least one colour 0 and at least one colour 1, plus the colour-swap
    unit clause.  If symbreak is nonzero, add the same lex-leader clauses that
    decide() uses.  Return (None, None) for FORCED instances."""
    sets = outside_sets(b, f, k)
    if sets is None:
        return None, None
    F = subgroup(b, f)
    free = [x for x in range(b) if x not in F]
    var = {x: i + 1 for i, x in enumerate(free)}
    cnf = CNF()
    for O in sets:
        lits = [var[x] for x in sorted(O)]
        cnf.append(lits)
        cnf.append([-v for v in lits])
    cnf.append([-1])                                  # colour-swap symmetry
    if symbreak:
        lex_leader(cnf, var, free,
                   symmetry_maps(b, f) if perms is None else perms, symbreak)
    return cnf, free


def decide(k, b, f, conf_budget=0, symbreak=0, solver="cadical153", perms=None):
    """Return (status, word); status in {'SAT', 'UNSAT', 'UNKNOWN', 'FORCED'}.

    FORCED means infeasible before solving: some progression has exactly one
    term outside F, and a single term is always constant.  perms overrides the
    symmetry group; only self_test() uses it, to plant a non-symmetry."""
    cnf, free = encode_cnf(k, b, f, symbreak, perms)
    if cnf is None:
        return "FORCED", None
    F = subgroup(b, f)
    var = {x: i + 1 for i, x in enumerate(free)}
    with Solver(name=solver, bootstrap_with=cnf) as s:
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
    # F positions are irrelevant; store them as 0
    word = [0 if x in F else (1 if var[x] in model else 0) for x in range(b)]
    h1, h2 = violations(word, k, F)
    ok, nd = safe_coset(word, k, F)
    assert (h1, h2) == (0, 0) and ok and nd == 0, \
        "solver model fails the conditions: H1=%d H2=%d safe_coset=%d" % (h1, h2, nd)
    return "SAT", "".join(map(str, word))


def candidates(b, above):
    """(b, f) with f | b, 1 <= f < b, b + f > above."""
    return [(b, f) for f in divisors(b) if f < b and b + f > above]


def self_test(L=24):
    """Symmetry breaking must not change any answer, and the comparison must
    be able to see it if it did.  Returns True on pass."""
    import random
    cases = [(k, b, f) for k, B in ((3, 30), (4, 60), (5, 90))
             for b in range(2, B) for _, f in candidates(b, 0)]
    plain = {c: decide(*c)[0] for c in cases}
    broken = {c: decide(*c, symbreak=L)[0] for c in cases}
    agree = sum(plain[c] == broken[c] for c in cases)
    sat = sum(v == "SAT" for v in plain.values())
    print("self-test: %d instances (%d SAT); symmetry breaking agrees on %d"
          % (len(cases), sat, agree))
    # Planted defect: random permutations of the free positions are not
    # symmetries, so lex-leader constraints for them may cut every model.
    rng = random.Random(1186)
    flipped = 0
    for c in cases:
        if plain[c] != "SAT":
            continue
        k, b, f = c
        F = subgroup(b, f)
        free = [x for x in range(b) if x not in F]
        perms = []
        for _ in range(20):
            sh = free[:]
            rng.shuffle(sh)
            img = list(range(b))
            for x, y in zip(free, sh):
                img[x] = y
            perms.append(img)
        if decide(k, b, f, symbreak=L, perms=perms)[0] == "UNSAT":
            flipped += 1
    print("self-test: planted non-symmetries turn %d of %d SAT instances UNSAT"
          % (flipped, sat))
    return agree == len(cases) and flipped > 0


def main():
    if "--self-test" in sys.argv:
        sys.exit(0 if self_test() else 1)
    ap = argparse.ArgumentParser()
    ap.add_argument("k", type=int)
    ap.add_argument("bs", help="range of b, e.g. 2-200 (or --self-test alone)")
    ap.add_argument("--above", type=int, default=0,
                    help="only test (b, f) with b + f > ABOVE")
    ap.add_argument("--budget", type=int, default=0,
                    help="conflict budget per instance; 0 = unlimited")
    ap.add_argument("--symbreak", type=int, default=0, metavar="L",
                    help="lex-leader symmetry breaking on the first L positions")
    ap.add_argument("--solver", default="cadical153",
                    help="any pysat solver name, e.g. kissat404, cadical195")
    ap.add_argument("--f", default="",
                    help="comma-separated f values to keep (default: all)")
    ap.add_argument("--store", default="")
    args = ap.parse_args()

    data = {}
    if args.store and os.path.exists(args.store):
        with open(args.store) as fh:
            data = json.load(fh)

    best = None
    counts = {}
    for b in parse_range(args.bs):
        for b_, f in candidates(b, args.above):
            if args.f and f not in {int(v) for v in args.f.split(",")}:
                continue
            status, word = decide(args.k, b, f, args.budget,
                                  args.symbreak, args.solver)
            counts[status] = counts.get(status, 0) + 1
            key = "%d/%d" % (b, f)
            rec = {"k": args.k, "b": b, "f": f, "status": status}
            if args.symbreak or args.solver != "cadical153":
                rec.update({"symbreak": args.symbreak, "solver": args.solver})
            if status == "SAT":
                den = 2 * (args.k - 1) * (b + f)
                rec.update({"word": word, "bound": "1/%d" % den})
                if best is None or b + f > best[0] + best[1]:
                    best = (b, f)
                print("b=%3d f=%3d  SAT      b+f=%d  delta_%d <= 1/%d  %s"
                      % (b, f, b + f, args.k, den, word), flush=True)
            elif status != "FORCED":
                print("b=%3d f=%3d  %s" % (b, f, status), flush=True)
            data[key] = rec
            if args.store:
                tmp = args.store + ".tmp"
                with open(tmp, "w") as fh:
                    json.dump(data, fh, indent=1, sort_keys=True)
                os.replace(tmp, args.store)

    print("k=%d  b in %s  b+f > %d:  %s" % (args.k, args.bs, args.above,
          "  ".join("%s=%d" % kv for kv in sorted(counts.items()))))
    if best:
        b, f = best
        print("best base: b=%d f=%d  ->  delta_%d <= 1/%d"
              % (b, f, args.k, 2 * (args.k - 1) * (b + f)), flush=True)
    else:
        print("no base found above b+f=%d" % args.above, flush=True)


if __name__ == "__main__":
    main()
