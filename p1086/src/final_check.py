"""Re-derive every claim in results/claims.json; exits 0 only if every check passes, the planted defects included.

    python src/final_check.py

Needs src/grid_count, src/brute_count, src/small_search and src/box_scan built (see README.md), results/raw/ as
written by src/lattice_table.py, and numpy, sympy and z3. Takes 3.5 to 7 minutes on a shared machine.
"""
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"
PY = sys.executable
BRUTE_SIZES = [(2, 2), (3, 3), (4, 6), (7, 7), (10, 10), (13, 4), (16, 16), (20, 20), (2, 40), (3, 100), (30, 30),
               (17, 53)]
REBUILD = [(25, 25), (50, 50), (100, 100), (200, 200), (3, 3333)]
outcomes = []


def check(name, ok, detail=""):
    outcomes.append(ok)
    print(f"{'PASS' if ok else 'FAIL'}  {name}{'  ' + detail if detail else ''}", flush=True)


def run(args, **kw):
    return subprocess.run([str(a) for a in args], capture_output=True, text=True, cwd=ROOT, **kw)


def last(r):
    lines = (r.stdout.strip() or r.stderr.strip() or "(no output)").splitlines()
    return lines[-1].strip()


def histogram(path):
    return {int(d): int(c) for d, c in (l.split() for l in Path(path).read_text().splitlines() if l and l[0] != "#")}


def main():
    claims = json.loads((ROOT / "results" / "claims.json").read_text())
    lattice = json.loads((ROOT / "results" / "lattice.json").read_text())
    tmp = Path(tempfile.mkdtemp(prefix="p1086check-"))
    try:
        # 1. grid_count against the triple-enumeration brute force, bin by bin
        for a, b in BRUTE_SIZES:
            out = tmp / f"gc_{a}x{b}.txt"
            out.write_text(run([SRC / "grid_count", a, b]).stdout)
            r = run([PY, SRC / "brute_check.py", "--compare", a, b, out])
            check(f"grid_count = brute force on the {a}x{b} grid", r.returncode == 0, last(r))
        planted = tmp / "gc_planted.txt"
        lines = (tmp / "gc_10x10.txt").read_text().splitlines()
        lines = [f"{l.split()[0]} {int(l.split()[1]) + 6}" if l.startswith("7 ") else l for l in lines]
        planted.write_text("\n".join(lines) + "\n")
        r = run([PY, SRC / "brute_check.py", "--compare", 10, 10, planted])
        check("planted: brute_check rejects a histogram with 6 extra triangles at D=7", r.returncode == 1)
        for a, b in [(25, 25), (50, 50), (100, 100)]:
            blob = run([SRC / "brute_count", a, b]).stdout
            stored = (ROOT / "results" / "raw" / f"grid_{a}x{b}.txt").read_text()
            check(f"brute_count (C triple enumeration) = results/raw/grid_{a}x{b}.txt byte for byte",
                  blob == stored and blob.startswith(f"# a={a} b={b} ") and len(blob.splitlines()) > 1)

        # 2. the stored histograms: rebuilt byte for byte, then rechecked against lattice.json
        by_grid = {(x["a"], x["b"]): x for x in lattice["runs"]}
        for a, b in REBUILD:
            blob = run([SRC / "grid_count", a, b, 0, 10]).stdout.encode()
            same = hashlib.sha256(blob).hexdigest() == by_grid[(a, b)]["sha256"]
            check(f"grid_count rebuilds the {a}x{b} histogram to the SHA-256 recorded in lattice.json", same)
        r = run([PY, SRC / "verify_lattice.py"])
        check("verify_lattice: lattice.json agrees with the raw histograms", r.returncode == 0, last(r))
        tree = tmp / "vl"
        (tree / "src").mkdir(parents=True)
        (tree / "results" / "raw").mkdir(parents=True)
        shutil.copy(SRC / "verify_lattice.py", tree / "src")
        for f in (ROOT / "results" / "raw").iterdir():
            (tree / "results" / "raw" / f.name).symlink_to(f)
        bad = tree / "results" / "raw" / "grid_25x25.txt"
        bad.unlink()
        lines = (ROOT / "results" / "raw" / "grid_25x25.txt").read_text().splitlines()
        bad.write_text("\n".join(f"12 {int(l.split()[1]) + 6}" if l.startswith("12 ") else l for l in lines) + "\n")
        doctored = json.loads(json.dumps(lattice))
        for x in doctored["runs"]:
            if (x["a"], x["b"]) == (25, 25):
                x["sha256"] = hashlib.sha256(bad.read_bytes()).hexdigest()
        (tree / "results" / "lattice.json").write_text(json.dumps(doctored))
        r = run([PY, tree / "src" / "verify_lattice.py"])
        check("planted: verify_lattice rejects 6 extra triangles at D=12 in the 25x25 grid, hash updated",
              r.returncode == 1)

        # 3. g(5) = 7, g(6) = 12 and g(7) = 18
        r = run([PY, SRC / "g5_certificate.py", 8])
        check("g5_certificate 8: no 5-point system with 8 equal areas is consistent over C",
              r.returncode == 0 and ", 0 consistent" in r.stdout, last(r))
        r = run([PY, SRC / "g5_certificate.py", 7])
        check("g5_certificate 7 (control): a real configuration with 7 is found", r.returncode == 0, last(r))
        r = run([PY, SRC / "g6_exact.py"], env={**os.environ, "PYTHONPATH": str(SRC)})
        check("g6_exact: classification of 7-configurations and the sixth point give g(6) = 12",
              r.returncode == 0 and "g(6) = max(12, 12) = 12" in r.stdout, last(r))
        r = run([PY, SRC / "g6_brute.py"])
        check("g6_brute: independent scan of sixth points, largest count 12", r.returncode == 0, last(r))
        wrong = tmp / "g6_wrong_R.py"
        wrong.write_text((SRC / "g6_exact.py").read_text().replace(
            "R = [(0, 0), (1, 0), (0, 1), (1, 1), (1, -1)]", "R = [(0, 0), (1, 0), (0, 1), (1, 1), (2, -1)]", 1))
        r = run([PY, wrong], env={**os.environ, "PYTHONPATH": str(SRC)})
        check("planted: g6_exact rejects a wrong reference configuration", r.returncode == 1)
        r = run([PY, SRC / "g7_no21.py"])
        check("g7_no21: no 7-point system has 21 triangles of one nonzero area",
              r.returncode == 0 and "g(7) <= 20" in r.stdout, last(r))
        r = run([PY, SRC / "g7_exact.py", "--write-to", tmp / "g7_classes.json"])
        g7_ok = r.returncode == 0 and "g(7) = 18" in r.stdout
        check("g7_exact: no 7-point set has 19 triangles of one nonzero area, so g(7) = 18", g7_ok, last(r))
        stored_classes = (ROOT / "results" / "g7_classes.json").read_text()
        check("g7_exact: results/g7_classes.json is what the enumeration writes",
              g7_ok and (tmp / "g7_classes.json").read_text() == stored_classes)
        r = run([PY, SRC / "g7_box_check.py"])
        check("g7_box_check: the 6x6 box scan finds only listed configurations, and the four-point lemma holds "
              "on the 7x7 box", r.returncode == 0, last(r))
        wrong = tmp / "g7_wrong_centroid.py"
        wrong.write_text((SRC / "g7_exact.py").read_text().replace(
            "(Fr(1, 3), Fr(1)))", "(Fr(1, 3), Fr(2)))", 1))
        r = run([PY, wrong])
        check("planted: g7_exact rejects a wrong centroid base", r.returncode == 1)
        doctored = json.loads(stored_classes)
        doctored["six_points_at_least_11"] = doctored["six_points_at_least_11"][1:]
        (tmp / "g7_classes_short.json").write_text(json.dumps(doctored))
        r = run([PY, SRC / "g7_box_check.py", tmp / "g7_classes_short.json"])
        check("planted: g7_box_check rejects a classification with one 6-point configuration dropped",
              r.returncode == 1)

        # 3b. g(9) = 36
        r = run([PY, SRC / "g9_exact.py", "--write-to", tmp / "g9_classes.json"])
        g9_ok = g7_ok and r.returncode == 0 and "g(9) = 36" in r.stdout
        check("g9_exact: no 9-point set has 37 triangles of one nonzero area, so g(9) = 36", g9_ok, last(r))
        stored_g9 = (ROOT / "results" / "g9_classes.json").read_text()
        check("g9_exact: results/g9_classes.json is what the enumeration writes",
              g9_ok and (tmp / "g9_classes.json").read_text() == stored_g9)
        r = run([PY, SRC / "g9_box_check.py"])
        g9_ok &= r.returncode == 0
        check("g9_box_check: 6x6 box scans, re-derived extensions and the z3 Type II check agree with g9_classes.json",
              r.returncode == 0, last(r))
        source = (SRC / "g9_exact.py").read_text()
        for name, old, new in (("an averaging step claimed for 36 instead of 37", "(9, 37, 25))", "(9, 36, 25))"),
                               ("Type II solutions dropped", "sols = sp.solve_poly_system(eqs, *gens)",
                                "sols = sp.solve_poly_system(eqs, *gens)[:0]"),
                               ("a Type II solution coordinate perturbed by 1/1000",
                                "sols = {tuple(sp.sympify(v) for v in sol) for sol in (sols or [])}",
                                "sols = {tuple(sp.sympify(v) + (sp.Rational(1, 1000) if i == 0 else 0) "
                                "for i, v in enumerate(sol)) for sol in (sols or [])}")):
            wrong = SRC / f"_planted_g9_{len(name)}.py"  # in src/ so that its imports resolve
            try:
                wrong.write_text(source.replace(old, new, 1))
                r = run([PY, wrong, "--write-to", tmp / "g9_planted.json"])
            finally:
                wrong.unlink()
            check(f"planted: g9_exact rejects {name}", source.count(old) == 1 and r.returncode == 1)
        doctored = json.loads(stored_g9)
        doctored["eight_points_at_least_25"] = doctored["eight_points_at_least_25"][:-1]
        (tmp / "g9_classes_short.json").write_text(json.dumps(doctored))
        r = run([PY, SRC / "g9_box_check.py", tmp / "g9_classes_short.json"])
        check("planted: g9_box_check rejects a classification with one 8-point configuration dropped",
              r.returncode == 1)

        # 4. small-n configurations
        small = json.loads((ROOT / "results" / "small_n.json").read_text())
        r = run([PY, SRC / "recount.py", ROOT / "results" / "small_n.json"])
        check("recount: every stored small-n configuration", r.returncode == 0, last(r))
        doctored = json.loads(json.dumps(small))
        doctored["best"][5]["count"] += 1
        (tmp / "small_bad.json").write_text(json.dumps(doctored))
        r = run([PY, SRC / "recount.py", tmp / "small_bad.json"])
        check("planted: recount rejects an overstated count", r.returncode == 1)
        r = run([PY, SRC / "recount.py", ROOT / "results" / "lower_search.json"])
        check("recount: every configuration kept by lower_search", r.returncode == 0, last(r))
        searched = json.loads((ROOT / "results" / "lower_search.json").read_text())["best"]
        doctored = json.loads((ROOT / "results" / "lower_search.json").read_text())
        doctored["best"][3]["count"] += 1
        (tmp / "lower_bad.json").write_text(json.dumps(doctored))
        r = run([PY, SRC / "recount.py", tmp / "lower_bad.json"])
        check("planted: recount rejects an overstated lower_search count", r.returncode == 1)

        # 5. claims.json against what was just recomputed
        from recount import classes
        best = {e["n"]: e for e in small["best"]}
        upper = {6: 12, 7: 18 if g7_ok else 20}
        for n in range(8, 17):
            upper[n] = upper[n - 1] * n // (n - 3)
            if n == 9 and g9_ok:
                upper[9] = min(upper[9], 36)
        for e in claims["exact"]:
            n, g = e["n"], e["g"]
            ok = best[n]["count"] == g if n in best else (n, g) == (3, 1)
            check(f"claim g({n}) = {g}: lower bound attained by a stored configuration", ok)
            if n >= 7:
                check(f"claim g({n}) = {g}: upper bound {upper[n]} from g(7) = 18, g(9) <= 36 and averaging",
                      upper[n] == g)
        for e in claims["small_n"]:
            n = e["n"]
            recounted = max(classes([tuple(p) for p in best[n]["points"]]).values())
            check(f"claim {e['lower']} <= g({n}) <= {e['upper']}",
                  recounted == e["lower"] and upper[n] == e["upper"] and e["lower"] <= e["upper"])
        for e in searched:
            n = e["n"]
            claimed = next((c["g"] for c in claims["exact"] if c["n"] == n), None)
            claimed = claimed if claimed is not None else next(c["lower"] for c in claims["small_n"] if c["n"] == n)
            check(f"lower_search at n = {n}: best {e['count']}, not above the claimed lower bound {claimed}",
                  e["count"] <= claimed and e["improved"] == (e["count"] > e["previous"]))
        for e in claims["lattice"]:
            h = histogram(ROOT / "results" / "raw" / f"grid_{e['a']}x{e['b']}.txt")
            ok = h.get(e["D"], 0) == e["count"] and max(h.values()) == e["count"] and e["n"] == e["a"] * e["b"]
            check(f"claim g({e['n']}) >= {e['count']} ({e['a']}x{e['b']} grid, doubled area {e['D']})", ok)
        r1 = claims["reproduced"]["DST09_T1"]
        h = histogram(ROOT / "results" / "raw" / f"grid_{r1['m']}x{r1['m']}.txt")
        ratio = h[1] / r1["m"] ** 4
        check(f"reproduced DST09: T_1/n^2 = {ratio:.6f} at m = {r1['m']}, 6/pi^2 = {6 / math.pi**2:.6f}",
              h[1] == r1["T1"] and abs(ratio - 6 / math.pi**2) < 1e-4)
        for s in claims["informal"]:
            print(f"NOT MACHINE-CHECKED  {s}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print(f"{sum(outcomes)} of {len(outcomes)} checks passed")
    return 0 if all(outcomes) else 1


if __name__ == "__main__":
    sys.path.insert(0, str(SRC))
    sys.exit(main())
