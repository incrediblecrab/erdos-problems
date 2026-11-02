"""Cross-check src/g10_exact.py, sharing no code with it, with src/g7_exact.py or with src/g9_exact.py.

    python src/g10_check.py [results/g10_classes.json] [--box] [--skip-lemma]      exit 0 iff every check passes

Needs cadical on PATH (Homebrew). The integer helpers come from src/g9_box_check.py, itself independent of the
exact scripts. Three checks, and a fourth with --box:

1. Listed counts. Every configuration in the file, scaled to integer coordinates, has the stated number of triangles
   of the common doubled area and no larger class of another area.
2. Re-derived growth. The rules a2, a3, b5 and b4 of g10_exact.py are carried out again in integer arithmetic from
   the listed 5-, 6- and 7-point classes (g7_classes.json, g9_classes.json) and from the two 4-point bases written
   here; every rule must reach exactly the listed 8-point area profiles, no line may serve as many pairs as a new
   point needs, extending the listed 8-point classes must give exactly the listed 9-point profiles, and no extension
   of those may reach 49, the largest being 48.
3. The hypergraph lemma, with a CNF encoding written here (its own sequential counters and chain clauses) and a
   different solver, cadical. The 9-regular 24-edge hypergraphs on 8 vertices are split into 21 cases: some pair has
   codegree 4, 5 or 6 (relabelled to {0, 1}, with third points 2, 3, ...), or every codegree is at most 3 and the link of vertex 0 (the 9 pairs
   {u, v} with {0, u, v} an edge, a graph with maximum degree 3) is one of the 18 isomorphism classes listed here.
   Every case must be unsatisfiable together with "no chain of rule a3, b5 or b4". As a control, the 18 link classes
   are enumerated over all labelled graphs, and dropping rule b4 must make some case satisfiable.
4. (--box, not run by final_check.py) src/box_scan lists every 8-subset of the 6 x 6 box with at least 24 triangles of
   one doubled area; each must have the area profile of a listed configuration.

Checks 3 and 4 can refute the classification or the lemma; check 2 guards against implementation errors, not against
an error in the written argument.
"""
import itertools
import json
import math
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
from pathlib import Path

import networkx as nx

from g9_box_check import candidate_points, largest_class, profile, to_integers, unit_count

ROOT = Path(__file__).resolve().parent.parent
THREADS = 3
RULES = {
    "a2": (7, 16, [(0, 24)]),
    "a3": (6, 10, [(0, 15), (0, 24)]),
    "b5": (5, 6, [(3, 9), (0, 15), (0, 24)]),
    "b4": (4, 3, [(2, 5), (3, 8), (4, 15), (0, 24)]),
}
LEVELS = (("eight_points_at_least_24", 8, 24), ("nine_points_at_least_35", 9, 35), ("ten_points_at_least_49", 10, 49))


def step(bases, a, t):
    """Integer growth step: every base plus a candidate point adding >= a triangles of the unit, total >= t."""
    out, worst_ok = {}, True
    for points, unit in bases:
        have = unit_count(points, unit)
        cands, worst = candidate_points(points, unit)
        worst_ok &= worst < max(a, t - have)
        for x, y in cands:
            if (x, y) in points:
                continue
            den = math.lcm(x.denominator, y.denominator)
            scaled = [(px * den, py * den) for px, py in points] + [(int(x * den), int(y * den))]
            u = unit * den * den
            n = unit_count(scaled, u)
            if n - have >= a and n >= t:
                out.setdefault(profile(scaled, u), (scaled, u))
    return list(out.values()), worst_ok


def best_extension(bases):
    best = 0
    for points, unit in bases:
        cands, _ = candidate_points(points, unit)
        for x, y in cands:
            if (x, y) in points:
                continue
            den = math.lcm(x.denominator, y.denominator)
            scaled = [(px * den, py * den) for px, py in points] + [(int(x * den), int(y * den))]
            best = max(best, unit_count(scaled, unit * den * den))
    return best


# ---- check 3: the lemma, encoded from scratch ------------------------------------------------------------------

class Cnf:
    def __init__(self):
        self.nv, self.cls = 0, []

    def new(self):
        self.nv += 1
        return self.nv

    def counter(self, lits, both):
        """Sequential counter; returns s with s[j] (1 <= j <= len) implied by 'at least j of lits true'. With both,
        s[j] also implies it, so s can state exact counts."""
        n = len(lits)
        prev = [None] * (n + 1)
        for i, x in enumerate(lits, start=1):
            cur = [None] + [self.new() for _ in range(n)]
            for j in range(1, n + 1):
                if j <= i - 1 and prev[j]:
                    self.cls.append([-prev[j], cur[j]])
                if j == 1:
                    self.cls.append([-x, cur[1]])
                elif prev[j - 1]:
                    self.cls.append([-x, -prev[j - 1], cur[j]])
                if both:
                    if j > i:
                        self.cls.append([-cur[j]])
                    else:
                        keep = [prev[j]] if j <= i - 1 and prev[j] else []
                        self.cls.append([-cur[j], x] + keep)
                        if j >= 2:
                            self.cls.append([-cur[j]] + ([prev[j - 1]] if prev[j - 1] else []) + keep)
            prev = cur
        return prev

    def solve(self, exe):
        with tempfile.NamedTemporaryFile("w", suffix=".cnf", delete=False) as fh:
            fh.write(f"p cnf {self.nv} {len(self.cls)}\n")
            fh.writelines(" ".join(map(str, c)) + " 0\n" for c in self.cls)
        r = subprocess.run([exe, "-q", fh.name], capture_output=True, text=True)
        Path(fh.name).unlink()
        if r.returncode not in (10, 20):
            raise RuntimeError(f"{exe} exited {r.returncode}")
        return r.returncode == 10


def link_classes():
    """Graphs on {1..7} with 9 edges and maximum degree 3, one per isomorphism class (networkx isomorphism tests
    within buckets of equal degree sequence and triangle count)."""
    reps = {}
    for E in itertools.combinations(itertools.combinations(range(1, 8), 2), 9):
        deg = [0] * 8
        for a, b in E:
            deg[a] += 1
            deg[b] += 1
        if max(deg) > 3:
            continue
        g = nx.Graph(E)
        g.add_nodes_from(range(1, 8))
        bucket = reps.setdefault((tuple(sorted(deg[1:])), sum(nx.triangles(g).values())), [])
        if not any(nx.is_isomorphic(g, h) for h, _ in bucket):
            bucket.append((g, E))
    return [E for b in reps.values() for _, E in b]


def build_case(case, rules):
    """The CNF of one case: a 9-regular 24-edge hypergraph in the case with no chain of the given rules."""
    f = Cnf()
    trip = {q: f.new() for q in itertools.combinations(range(8), 3)}
    every = f.counter(list(trip.values()), both=True)
    f.cls += [[every[24]], [-every[25]]]
    for v in range(8):
        s = f.counter([trip[q] for q in trip if v in q], both=True)
        f.cls += [[s[9]], [-s[10]]]
    if isinstance(case, int):  # codegree of {0, 1} exactly case >= 4, its third points relabelled 2, ..., case + 1
        f.cls += [[trip[(0, 1, w)] if w <= case + 1 else -trip[(0, 1, w)]] for w in range(2, 8)]
    else:
        for p in itertools.combinations(range(8), 2):
            s = f.counter([trip[q] for q in trip if set(p) <= set(q)], both=False)
            f.cls.append([-s[4]])
        for u, w in itertools.combinations(range(1, 8), 2):
            f.cls.append([trip[(0, u, w)] if (u, w) in case else -trip[(0, u, w)]])
    memo = {}

    def at_least(X, t, y=None):
        lits = [trip[q] for q in itertools.combinations(X, 3) if y is None or y in q]
        if t <= 0:
            return None
        if t > len(lits):
            return False
        if (X, y) not in memo:
            memo[(X, y)] = f.counter(lits, both=False)
        return memo[(X, y)][t]

    for name in rules:
        m, tm, steps = RULES[name]
        reach = {X: f.new() for k in range(m, 9) for X in itertools.combinations(range(8), k)}
        for X in itertools.combinations(range(8), m):
            lit = at_least(X, tm)
            if lit is not False:
                f.cls.append([reach[X]] + ([-lit] if lit else []))
        for k, (a, t) in zip(range(m + 1, 9), steps):
            for X in itertools.combinations(range(8), k):
                le = at_least(X, t)
                if le is False:
                    continue
                for y in X:
                    ld = at_least(X, a, y)
                    if ld is False:
                        continue
                    sub = tuple(v for v in X if v != y)
                    f.cls.append([reach[X], -reach[sub]] + ([-le] if le else []) + ([-ld] if ld else []))
        f.cls.append([-reach[tuple(range(8))]])
    return f


def lemma_case(case, rules, exe):
    return build_case(case, rules).solve(exe)


def main(argv):
    box, skip_lemma = "--box" in argv, "--skip-lemma" in argv
    argv = [a for a in argv if a not in ("--box", "--skip-lemma")]
    path = Path(argv[0]) if argv else ROOT / "results" / "g10_classes.json"
    listed = json.loads(path.read_text())
    g7 = json.loads((ROOT / "results" / "g7_classes.json").read_text())
    g9 = json.loads((ROOT / "results" / "g9_classes.json").read_text())
    ok = True

    levels = {}
    for key, k, t in LEVELS:
        levels[key] = []
        for entry in listed[key]:
            pts, unit = to_integers(entry)
            n = unit_count(pts, unit)
            if len(set(pts)) != k or n != entry["edges"] or n < t or largest_class(pts) != n:
                print(f"{key}: stored count {entry['edges']} not reproduced (recount {n})")
                ok = False
            levels[key].append((pts, unit))
    print("listed counts: " + ", ".join(f"{key} {len(levels[key])}" for key, _, _ in LEVELS))
    ok &= not listed["ten_points_at_least_49"]

    start = {
        "a2": [to_integers(e) for e in g9["seven_points_at_least_16"]],
        "a3": [to_integers(e) for e in g9["six_points_at_least_10"]],
        "b5": [to_integers(e) for e in g7["five_points_at_least_6"]],
        "b4": [([(0, 0), (1, 0), (1, 1), (0, 1)], 1), ([(0, 0), (3, 0), (0, 9), (1, 3)], 9)],
    }
    ok &= [unit_count(*b) for b in start["b4"]] == [4, 3]
    mine8 = {profile(p, u) for p, u in levels["eight_points_at_least_24"]}
    mine8 |= {profile(*to_integers(e)) for e in g9["eight_points_at_least_25"]}
    ok &= len(mine8) == 23  # 24 listed configurations, two sharing a profile (NOTES.md section 4)
    if len(mine8) != len(levels["eight_points_at_least_24"]):
        print(f"  note: {len(levels['eight_points_at_least_24']) - len(mine8)} profile collision(s) among the listed "
              f"8-point configurations, which profile comparisons cannot tell apart")
    for name, (m, tm, steps) in RULES.items():
        current, worst_ok = start[name], True
        for a, t in steps:
            current, w = step(current, a, t)
            worst_ok &= w
        got = {profile(p, u) for p, u in current}
        same = got == mine8
        print(f"re-derived rule {name}: {len(got)} eight-point profiles, equal to the listed ones: {same}; every "
              f"line serves fewer pairs than needed: {worst_ok}")
        ok &= same and worst_ok
    nine, worst_ok = step(levels["eight_points_at_least_24"], 0, 35)
    same9 = {profile(p, u) for p, u in nine} == {profile(p, u) for p, u in levels["nine_points_at_least_35"]}
    ten, worst10 = step(levels["nine_points_at_least_35"], 0, 49)
    best = best_extension(levels["nine_points_at_least_35"])
    print(f"re-derived nine points: {len(nine)} profiles, equal to the listed ones: {same9}; ten points with >= 49: "
          f"{len(ten)}; largest extension {best}; lines fine: {worst_ok and worst10}")
    ok &= same9 and not ten and best == 48 == listed["largest_ten_point_extension"] and worst_ok and worst10

    if skip_lemma:  # used only by final_check.py's planted-defect run, which tests check 2
        print("lemma skipped")
        print("PASS" if ok else "FAIL")
        return 0 if ok else 1
    exe = shutil.which("cadical") or "/opt/homebrew/bin/cadical"
    links = link_classes()
    cases = [4, 5, 6] + [set(E) for E in links]
    with ThreadPoolExecutor(THREADS) as ex:
        res = list(ex.map(lambda c: lemma_case(c, ("a3", "b5", "b4"), exe), cases))
    print(f"lemma by cadical: {len(links)} link classes, {len(cases)} cases, satisfiable: {sum(res)}")
    ok &= len(links) == 18 and not any(res)
    with ThreadPoolExecutor(THREADS) as ex:
        ctrl = list(ex.map(lambda c: lemma_case(c, ("a3", "b5"), exe), cases))
    print(f"control: without rule b4, {sum(ctrl)} of the {len(cases)} cases are satisfiable")
    ok &= sum(ctrl) == 17  # the recorded control value (NOTES.md section 4); any change needs a look

    if box:
        r = subprocess.run([str(ROOT / "src" / "box_scan"), "6", "6", "8", "24", str(THREADS)], capture_output=True,
                           text=True, check=True)
        found = set()
        for line in r.stdout.splitlines():
            if line and not line.startswith("done"):
                f = line.split()
                pts = [tuple(map(int, s.split(","))) for s in f[2:]]
                found.add(profile(pts, int(f[0])))
        print(f"box 6x6, 8-subsets with >= 24: {len(found)} profiles, {len(found - mine8)} unlisted")
        ok &= not (found - mine8)

    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
