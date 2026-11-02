"""Determine g(10) exactly: classify C(8, 24) and C(9, 35), then show that no 10-point set has 49 triangles of one area.

    python src/g10_exact.py            prints the classification and the bound; exit 0 iff the run is complete
    python src/g10_exact.py --write    also writes results/g10_classes.json (--write-to PATH writes elsewhere)
    python src/g10_exact.py --full-lemma    only checks the lemma without the symmetry reduction (about 95 s; not
                                            run by final_check.py)

Needs kissat on PATH (Homebrew) and python-sat. Conventions are those of g7_exact.py and g9_exact.py: an "edge" is a
triangle of doubled area exactly 1 after one scaling map, C(n, t) is the set of n-point sets with at least t edges up
to affine maps of determinant +-1, and the written steps are in NOTES.md section 4.

Growth rules. A rule starts from a class C(m, t_m) that is already classified completely and adds points one at a
time; step k says that the new point has at least a_k edges with pairs of the earlier points and that the k-point
set has at least t_k edges. Before each step the script checks that no line {X : det(q - p, X - p) = +-1} of any base
serves as many pairs as the new point needs, max(a_k, t_k - edges(base)), so the new point is an intersection of two
distinct lines: finitely many candidates, each counted exactly. The rules used, with the classes they start from:

    a2  C(7, 16)                                 -> 8: t >= 24
    a3  C(6, 10)  -> 7: t >= 15                  -> 8: t >= 24
    b5  C(5, 6)   -> 6: a >= 3, t >= 9   -> 7: t >= 15 -> 8: t >= 24
    b4  C(4, 3)   -> 5: a >= 2, t >= 5   -> 6: a >= 3, t >= 8 -> 7: a >= 4, t >= 15 -> 8: t >= 24

C(4, 3) is the unit square and the triangle (0,0), (1,0), (0,3) with its centroid (the four-point lemma of g7_exact.py).

Eight points (written). An 8-point set with e >= 24 edges is in C(8, 25) if e >= 25. If e = 24, the degrees sum to
72; a vertex of degree at most 8 leaves a 7-subset with at least 16 edges (rule a2), and otherwise every degree is
exactly 9.

Lemma (machine-checked here, not written). Every 3-uniform hypergraph on 8 vertices with 24 edges and every degree 9
contains a chain of rule a3, b5 or b4: nested subsets whose edge counts and new-vertex degrees meet the rule's bounds.
The script encodes the negation as CNF and kissat must report it unsatisfiable. The lemma uses no geometry. Before
encoding, a written symmetry reduction: the codegrees sum to 3 * 24 = 72 > 2 * 28, so some pair lies in at least 3
edges, and after relabelling the pair is {0, 1} and (0,1,2), (0,1,3), (0,1,4) are edges. Three controls guard the
encoding: the degree constraints alone are satisfiable and the model is checked to be 9-regular with 24 edges; with
rule b4 left out the formula is satisfiable and a direct search over subsets, written separately from the CNF,
confirms that the model has no a3 or b5 chain and does have a b4 chain; and every 24-edge configuration found has a
chain of one of the rules a2, a3, b5, b4 under that search.

So C(8, 24) is C(8, 25) together with whatever rules a2, a3, b5 and b4 reach; the script checks that each rule reaches
the same 8-point classes.

Nine points (written). A 9-point set with e = 35 or 36 edges (37 is impossible, g(9) = 36) has a vertex of degree at
most floor(3e/9), and deleting it leaves at least e - floor(3e/9) = 24 edges, so C(9, 35) is C(8, 24) plus a point
adding at least 35 - edges(base) edges.

Ten points (written). A 10-point set with e >= 49 edges has a vertex of degree at most floor(3e/10), and
e - floor(3e/10) >= 35 for every e >= 49, so it is a member of C(9, 35) plus a point adding at least 49 - 36 = 13
edges. If no such extension exists, g(10) <= 48, and the stored 10-point configuration with 48 gives g(10) = 48.
"""
import itertools
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter
from fractions import Fraction as Fr
from pathlib import Path

from pysat.card import CardEnc, EncType, ITotalizer
from pysat.formula import CNF, IDPool

from g7_exact import canonical, edges, intersections, lines

ROOT = Path(__file__).resolve().parent.parent
TEN = 49  # the 10-point threshold this script rules out
NINE = 35
EIGHT = 24
RULES = {  # name: (start size, start threshold, [(a_k, t_k) for the following sizes up to 8])
    "a2": (7, 16, [(0, 24)]),
    "a3": (6, 10, [(0, 15), (0, 24)]),
    "b5": (5, 6, [(3, 9), (0, 15), (0, 24)]),
    "b4": (4, 3, [(2, 5), (3, 8), (4, 15), (0, 24)]),
}
LEMMA_RULES = ("a3", "b5", "b4")
CANDIDATES = {}  # label -> number of candidate points tried by grow


def grow(bases, a, t, label):
    """All base + {X} with X adding at least a edges and the new set having at least t, up to equivalence."""
    found, candidates = {}, 0
    for base in bases:
        eb = edges(base)
        need = max(a, t - eb)
        ld = lines(base)
        if max(ld.values()) >= need:
            raise RuntimeError(f"{label}: a line serves {max(ld.values())} pairs, the new point needs {need}")
        for pt in intersections(ld):
            if pt in base:
                continue
            candidates += 1
            new = base + (pt,)
            e = edges(new)
            if e - eb >= a and e >= t:
                found[canonical(new)] = e
    CANDIDATES[label] = candidates
    print(f"  {label}: {len(bases)} bases, {candidates} candidate points, {len(found)} classes, edge counts "
          f"{sorted(Counter(found.values()).items())}", flush=True)
    return found


def run_rule(name, start):
    m, _, steps = RULES[name]
    current = start
    for i, (a, t) in enumerate(steps):
        current = grow(list(current), a, t, f"rule {name}, {m + i + 1} points, a >= {a}, t >= {t}")
    return current


# ---- the hypergraph lemma -----------------------------------------------------------------------------------------

def chains(hyperedges, n, rule):
    """Direct search, independent of the CNF: does the hypergraph contain a chain of the rule?"""
    m, tm, steps = RULES[rule]
    E = set(hyperedges)
    count = lambda X: sum(1 for q in itertools.combinations(sorted(X), 3) if q in E)
    deg = lambda X, y: sum(1 for q in itertools.combinations(sorted(X), 3) if y in q and q in E)
    level = {frozenset(X) for X in itertools.combinations(range(n), m) if count(X) >= tm}
    for k, (a, t) in zip(range(m + 1, n + 1), steps):
        level = {X | {y} for X in level for y in range(n) if y not in X
                 if count(X | {y}) >= t and deg(X | {y}, y) >= a}
    return bool(level)


class Formula:
    """CNF over the 56 triples of {0..7}. geq(X, t, y) is a literal forced true whenever at least t triples of X
    (containing y, if given) are edges; pysat's totalizer only forces outputs upward, which is all the rule clauses
    need, since they only ever derive membership of a chain."""

    def __init__(self, n):
        self.n, self.pool, self.clauses, self.tot = n, IDPool(), [], {}
        self.x = {q: self.pool.id(("x",) + q) for q in itertools.combinations(range(n), 3)}

    def geq(self, X, t, y=None):
        lits = [self.x[q] for q in itertools.combinations(sorted(X), 3) if y is None or y in q]
        if t <= 0:
            return None
        if t > len(lits):
            return False
        key = (tuple(sorted(X)), y)
        if key not in self.tot:
            tot = ITotalizer(lits=lits, ubound=len(lits), top_id=self.pool.top)
            self.pool.occupy(self.pool.top + 1, tot.top_id)
            self.pool.top = tot.top_id
            self.clauses += tot.cnf.clauses
            self.tot[key] = tot
        return self.tot[key].rhs[t - 1]

    def exactly(self, lits, k):
        self.clauses += CardEnc.equals(lits, bound=k, vpool=self.pool, encoding=EncType.seqcounter).clauses

    def no_chain(self, rule):
        """Clauses saying that the hypergraph has no chain of the rule."""
        m, tm, steps = RULES[rule]
        n = self.n
        c = {X: self.pool.id((rule, X)) for k in range(m, n + 1) for X in itertools.combinations(range(n), k)}
        for X in itertools.combinations(range(n), m):
            lit = self.geq(X, tm)
            if lit is not False:
                self.clauses.append([c[X]] + ([-lit] if lit else []))
        for k, (a, t) in zip(range(m + 1, n + 1), steps):
            for X in itertools.combinations(range(n), k):
                le = self.geq(X, t)
                if le is False:
                    continue
                for y in X:
                    ld = self.geq(X, a, y)
                    if ld is False:
                        continue
                    rest = tuple(v for v in X if v != y)
                    self.clauses.append([c[X], -c[rest]] + ([-le] if le else []) + ([-ld] if ld else []))
        self.clauses.append([-c[tuple(range(n))]])


def regular_formula(rules, wlog=True):
    f = Formula(8)
    f.exactly(list(f.x.values()), 24)
    for v in range(8):
        f.exactly([f.x[q] for q in f.x if v in q], 9)
    if wlog:
        for q in ((0, 1, 2), (0, 1, 3), (0, 1, 4)):
            f.clauses.append([f.x[q]])
    for r in rules:
        f.no_chain(r)
    return f


def kissat(f):
    """Returns (satisfiable, list of edges of the model or None)."""
    exe = shutil.which("kissat") or "/opt/homebrew/bin/kissat"
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "f.cnf")
        CNF(from_clauses=f.clauses).to_file(path)
        r = subprocess.run([exe, "-q", path], capture_output=True, text=True)
    if r.returncode == 20:
        return False, None
    if r.returncode != 10:
        raise RuntimeError(f"kissat exited {r.returncode}")
    model = {int(v) for line in r.stdout.splitlines() if line.startswith("v") for v in line.split()[1:]}
    return True, sorted(q for q, v in f.x.items() if v in model)


def lemma():
    ok = True
    t0 = time.time()
    sat, model = kissat(regular_formula(()))
    regular = sat and len(model) == 24 and all(sum(v in q for q in model) == 9 for v in range(8))
    print(f"lemma control: the degree constraints alone are satisfiable with a 9-regular 24-edge model: {regular}")
    ok &= regular
    sat, model = kissat(regular_formula(("a3", "b5")))
    regular2 = sat and len(model) == 24 and all(sum(v in q for q in model) == 9 for v in range(8))
    found = {r: chains(model, 8, r) for r in RULES} if sat else {}
    good = sat and regular2 and not found["a3"] and not found["b5"] and found["b4"]
    print(f"lemma control: without rule b4 the formula is satisfiable, and the direct search finds in the model "
          f"chains {sorted(r for r, v in found.items() if v)} (no a3, no b5, a b4 chain), and that model is "
          f"9-regular with 24 edges: {good}")
    ok &= good
    sat, _ = kissat(regular_formula(LEMMA_RULES))
    print(f"lemma: every 9-regular 24-edge hypergraph on 8 vertices has a chain of rule a3, b5 or b4 "
          f"(kissat: {'SATISFIABLE' if sat else 'UNSATISFIABLE'}), {time.time() - t0:.1f} s")
    ok &= not sat
    ok &= 3 * 24 > 2 * 28  # the symmetry reduction: some pair has codegree >= 3
    return ok


# ---- driver -------------------------------------------------------------------------------------------------------

def as_json(classes):
    return [{"edges": e, "points": [[str(x), str(y)] for x, y in c]}
            for c, e in sorted(classes.items(), key=lambda kv: (-kv[1], kv[0]))]


def hypergraph(points):
    return [q for q in itertools.combinations(range(len(points)), 3)
            if abs((points[q[1]][0] - points[q[0]][0]) * (points[q[2]][1] - points[q[0]][1])
                   - (points[q[1]][1] - points[q[0]][1]) * (points[q[2]][0] - points[q[0]][0])) == 1]


def main(argv):
    start = time.time()
    ok = True
    if "--full-lemma" in argv:
        sat, _ = kissat(regular_formula(LEMMA_RULES, wlog=False))
        print(f"lemma without the symmetry reduction: kissat {'SATISFIABLE' if sat else 'UNSATISFIABLE'}, "
              f"{time.time() - start:.1f} s")
        return 1 if sat else 0
    parse = lambda e: tuple((Fr(x), Fr(y)) for x, y in e["points"])
    g7 = json.loads((ROOT / "results" / "g7_classes.json").read_text())
    g9 = json.loads((ROOT / "results" / "g9_classes.json").read_text())
    base = {
        "a2": [parse(e) for e in g9["seven_points_at_least_16"]],
        "a3": [parse(e) for e in g9["six_points_at_least_10"]],
        "b5": [parse(e) for e in g7["five_points_at_least_6"]],
        "b4": [((Fr(0), Fr(0)), (Fr(1), Fr(0)), (Fr(1), Fr(1)), (Fr(0), Fr(1))),
               ((Fr(0), Fr(0)), (Fr(1), Fr(0)), (Fr(0), Fr(3)), (Fr(1, 3), Fr(1)))],
    }
    sizes = {"a2": (19, 16), "a3": (15, 10), "b5": (4, 6), "b4": (2, 3)}
    for r, bs in base.items():
        good = len(bs) == sizes[r][0] and all(edges(b) >= sizes[r][1] for b in bs)
        ok &= good
        if not good:
            print(f"the starting classes of rule {r} are not the stored ones")
    c825 = {parse(e): e["edges"] for e in g9["eight_points_at_least_25"]}

    # averaging arithmetic for the written steps (e only goes up to 36 at 9 points, but check the whole range)
    ok &= all(e - 3 * e // 9 >= EIGHT for e in range(NINE, 85))
    ok &= all(e - 3 * e // 10 >= NINE for e in range(TEN, 121))
    ok &= 3 * 24 // 8 == 9 and 24 - 8 >= 16  # eight points: a degree <= 8 leaves >= 16, else all degrees are 9

    ok &= lemma()
    if not ok:  # the averaging arithmetic, the starting classes or the lemma failed: nothing below would be a proof
        print("a precondition failed; stopping")
        return 1

    print("eight points with at least 24 edges:")
    reached = {r: run_rule(r, base[r]) for r in RULES}
    eight = dict(c825)
    for r in RULES:
        eight.update(reached[r])
    same = all(set(v) == set(eight) for v in reached.values())
    print(f"C(8, 24): {len(eight)} classes, edge counts {sorted(Counter(eight.values()).items())}; every rule reaches "
          f"the same classes, which contain C(8, 25): {same and set(c825) <= set(eight)}")
    ok &= same and set(c825) <= set(eight)
    exact24 = [c for c, e in eight.items() if e == 24]
    covered = all(any(chains(hypergraph(c), 8, r) for r in RULES) for c in exact24)
    nine_regular = sum(1 for c in exact24 if all(sum(v in q for q in hypergraph(c)) == 9 for v in range(8)))
    print(f"control: each of the {len(exact24)} classes with exactly 24 has a chain under the direct search: "
          f"{covered}; {nine_regular} of them are 9-regular")
    ok &= covered

    nine = grow(list(eight), 0, NINE, f"C(9, {NINE})")
    ring9 = canonical(tuple((Fr(x, 2), Fr(y)) for x, y in [(x, y) for x in range(3) for y in range(3)
                                                           if (x, y) != (1, 1)] + [(3, 0)]))
    ok &= nine.get(ring9) == 36 and max(nine.values()) == 36 and sum(1 for e in nine.values() if e == 36) == 2
    mult = {k: max(max(lines(c).values()) for c in v) for k, v in (("eight", eight), ("nine", nine))}
    print(f"largest line multiplicity: {mult}; C(9, 35) has exactly two classes with 36: "
          f"{sum(1 for e in nine.values() if e == 36) == 2}")
    print(f"control: the stored 9-point configuration with 36 (3 x 3 grid minus its centre, plus (3, 0)) is in "
          f"C(9, 35): {nine.get(ring9) == 36}")

    ten = grow(list(nine), 0, TEN, f"C(10, {TEN})")
    best, top = 0, {}
    for b in nine:
        for pt in intersections(lines(b)):
            if pt not in b:
                e = edges(b + (pt,))
                if e > best:
                    best, top = e, {}
                if e == best:
                    top[canonical(b + (pt,))] = e
    print(f"largest count among all intersection-point extensions of C(9, {NINE}): {best}, attained by {len(top)} "
          f"class(es)")
    small = json.loads((ROOT / "results" / "small_n.json").read_text())
    stored10 = next(b for b in small["best"] if b["n"] == 10)
    pts10 = tuple((Fr(x, stored10["D"]), Fr(y)) for x, y in stored10["points"])
    in_top = canonical(pts10) in top and stored10["count"] == 48
    print(f"control: the stored 10-point configuration with 48 is one of them: {in_top}")
    ok &= not ten and best == 48 and in_top
    # the classification statistics NOTES.md states, so that the prose and the gate cannot drift apart
    stats = (len(eight), len(exact24), len(nine), nine_regular, mult, len(top))
    stated = (24, 18, 3, 0, {"eight": 4, "nine": 4}, 1)
    print(f"stated statistics (|C(8, 24)|, exactly 24, |C(9, 35)|, 9-regular, line multiplicities, 48-classes): "
          f"{stats == stated}")
    ok &= stats == stated

    target = None
    if "--write" in argv:
        target = ROOT / "results" / "g10_classes.json"
    elif "--write-to" in argv:
        target = Path(argv[argv.index("--write-to") + 1])
    if target is not None:
        out = {
            "doubled_area": 1,
            "eight_points_at_least_24": as_json(eight),
            "nine_points_at_least_35": as_json(nine),
            "ten_points_at_least_49": as_json(ten),
            "nine_regular_eight_point_sets": nine_regular,
            "largest_line_multiplicity": mult,
            "ten_point_candidates": CANDIDATES[f"C(10, {TEN})"],
            "largest_ten_point_extension": best,
            "ten_point_extensions_with_48": as_json(top),
        }
        target.write_text(json.dumps(out, indent=1) + "\n")
    print(f"elapsed {time.time() - start:.1f} s")
    if ok:
        print("no 10-point set has 49 triangles of one area, so g(10) <= 48; with the stored configuration, g(10) = 48")
        return 0
    print("the controls, the lemma or the enumeration did not come out as expected")
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
