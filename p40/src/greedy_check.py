#!/usr/bin/env python3
"""Reproduce the greedy B2[2] computation posted in the erdosproblems.com thread for problem 158, extend it, and check it without trusting the program that produced it.

The thread comment (Zeraoulia Rafik, February 1, 2026) reports a_2000 = 7,445,662 for the greedy B2[2] set and |A ∩ [1, N]| / √N ≈ 0.733 at N = a_2000. The note it links (Zenodo, DOI 10.5281/zenodo.18452185, file v7) tabulates |A ∩ [1, N]| at seven values of N in its Table 1. By `ErdosProblem40.strong_iff_forall_B2` a B2[2] set with liminf |A ∩ [1, N]| / √N > 0 would make the answer set of problem 40 empty, so the claim bears on problem 40. No finite computation decides a liminf; everything here is illustrative.

Steps, each of which must pass:
1. Compile src/greedy_b2g.c into a temporary directory and run it under a CPU-time limit set with setrlimit, which the kernel enforces.
2. Known answer: g = 1 is the Mian-Chowla sequence, OEIS A005282; the first G1_TERMS terms must equal refs/b005282.txt.
3. A separate Python implementation of the greedy rule must reproduce the first PY_TERMS terms for g = 2.
4. The whole g = 2 output must be a B2[2] set, counted from all pairwise sums.
5. The thread's two numbers, and the seven counts and ratios of the note's Table 1, must be reproduced.
6. The CPU-time limit must stop a run when lowered to 1 s.

Writes results/greedy_b2g.json and results/greedy_b2g_g2.txt. With --plant, one defect is put in before the checks: "nongreedy" and "notb2" alter the g = 2 output, and "note" alters one count claimed in the note's table. The script must then print "caught" and exit 1 without writing anything.
"""
import argparse
import bisect
import hashlib
import json
import math
import resource
import signal
import subprocess
import sys
import tempfile
import time
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src" / "greedy_b2g.c"
BFILE = ROOT / "refs" / "b005282.txt"
OUT_JSON = ROOT / "results" / "greedy_b2g.json"
OUT_LIST = ROOT / "results" / "greedy_b2g_g2.txt"

K = 10000
VMAX = 400_000_000
G1_TERMS = 2000
G1_VMAX = 100_000_000
PY_TERMS = 500
CPU_LIMIT_S = 900
CHECKPOINTS = (100, 200, 500, 1000, 2000, 5000, 10000)
THREAD_K, THREAD_A, THREAD_RATIO = 2000, 7445662, 0.733
# Table 1 of the note: (N, |A ∩ [1, N]|, |A ∩ [1, N]| / √N as printed to four places).
NOTE_TABLE = ((10, 6, 1.8974), (100, 17, 1.7), (1000, 48, 1.5179), (10_000, 127, 1.27), (100_000, 332, 1.0499), (1_000_000, 870, 0.87), (7_445_662, 2000, 0.733))
PLANTS = ("nongreedy", "notb2", "note")


def run_greedy(binary, g, k, vmax, out, cpu_limit):
    def limit():
        resource.setrlimit(resource.RLIMIT_CPU, (cpu_limit, cpu_limit))

    t0 = time.monotonic()
    p = subprocess.run([str(binary), str(g), str(k), str(vmax), str(out)], preexec_fn=limit, capture_output=True, text=True)
    return p.returncode, time.monotonic() - t0, p.stderr.strip()


def read_terms(path):
    return [int(line) for line in path.read_text().split()]


def greedy_python(g, k):
    """The greedy rule written again from the definition: accept m when every sum it creates still has at most g representations a + b with a <= b."""
    terms, reps = [], Counter()
    m = 0
    while len(terms) < k:
        m += 1
        sums = [a + m for a in terms] + [2 * m]
        if all(reps[s] < g for s in sums):
            reps.update(sums)
            terms.append(m)
    return terms


def max_representations(terms):
    """Largest number of representations n = a + b with a <= b over all n: sort all pairwise sums and take the longest run of equal values."""
    a = np.asarray(terms, dtype=np.int64)
    k = len(a)
    sums = np.empty(k * (k + 1) // 2, dtype=np.int64)
    pos = 0
    for i in range(k):
        sums[pos:pos + k - i] = a[i] + a[i:]
        pos += k - i
    sums.sort()
    starts = np.flatnonzero(np.diff(sums)) + 1
    runs = np.diff(np.concatenate(([0], starts, [len(sums)])))
    return int(runs.max())


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--plant", choices=PLANTS, help="put one defect into the g = 2 output; the checks must catch it")
    args = ap.parse_args()

    problems = []
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        binary = tmp / "greedy_b2g"
        cc = subprocess.run(["cc", "-O2", "-Wall", "-Wextra", "-Werror", "-o", str(binary), str(SRC)], capture_output=True, text=True)
        if cc.returncode != 0:
            sys.exit(f"compile failed:\n{cc.stderr}")

        # Step 2: known answer for g = 1.
        rc, t1, err1 = run_greedy(binary, 1, G1_TERMS, G1_VMAX, tmp / "g1.txt", CPU_LIMIT_S)
        g1 = read_terms(tmp / "g1.txt") if rc == 0 else []
        if not BFILE.exists():
            sys.exit(f"{BFILE.relative_to(ROOT)} missing; run refs/fetch.sh")
        oeis = [int(line.split()[1]) for line in BFILE.read_text().splitlines() if line.strip() and not line.startswith("#")][:G1_TERMS]
        g1_ok = rc == 0 and g1 == oeis and len(oeis) == G1_TERMS
        if not g1_ok:
            problems.append(f"g = 1: exit {rc}, first {G1_TERMS} terms differ from OEIS A005282")

        # Step 1 for g = 2.
        rc, t2, err2 = run_greedy(binary, 2, K, VMAX, tmp / "g2.txt", CPU_LIMIT_S)
        child = resource.getrusage(resource.RUSAGE_CHILDREN)
        if rc != 0:
            sys.exit(f"g = 2 run failed with exit {rc}: {err2}")
        g2 = read_terms(tmp / "g2.txt")

        # Step 6: the CPU limit stops a run when lowered.
        rc_trip, t_trip, _ = run_greedy(binary, 2, 4 * K, 4_000_000_000 // 2, tmp / "trip.txt", 1)
        trip_ok = rc_trip == -signal.SIGXCPU
        if not trip_ok:
            problems.append(f"CPU limit of 1 s did not stop the run: exit {rc_trip} after {t_trip:.1f} s")

    note_table = list(NOTE_TABLE)
    if args.plant == "nongreedy":
        g2[PY_TERMS // 2] += 1
    elif args.plant == "notb2":
        g2 = sorted(g2 + [5])
    elif args.plant == "note":
        n, c, r = note_table[5]
        note_table[5] = (n, c + 1, r)

    # Step 3: independent reimplementation on a prefix.
    t0 = time.monotonic()
    py = greedy_python(2, PY_TERMS)
    t_py = time.monotonic() - t0
    if g2[:PY_TERMS] != py:
        problems.append(f"g = 2: first {PY_TERMS} terms differ from the Python implementation")

    # Step 4: the output is a B2[2] set.
    strictly_increasing = all(x < y for x, y in zip(g2, g2[1:])) and g2[0] >= 1
    t0 = time.monotonic()
    maxrep = max_representations(g2)
    t_b2 = time.monotonic() - t0
    if not strictly_increasing or maxrep > 2:
        problems.append(f"g = 2: not a B2[2] set (increasing {strictly_increasing}, max representations {maxrep})")

    # Step 5: the thread's numbers.
    a_thread = g2[THREAD_K - 1]
    ratio_thread = THREAD_K / math.sqrt(a_thread)
    if a_thread != THREAD_A or round(ratio_thread, 3) != THREAD_RATIO:
        problems.append(f"thread numbers not reproduced: a_{THREAD_K} = {a_thread}, ratio {ratio_thread:.6f}")
    note_rows = []
    for n, count_claimed, ratio_claimed in note_table:
        count = bisect.bisect_right(g2, n)
        note_rows.append({"N": n, "count_claimed": count_claimed, "count_computed": count, "ratio_claimed": ratio_claimed, "ratio_computed": round(count / math.sqrt(n), 6)})
        if count != count_claimed or round(count / math.sqrt(n), 4) != ratio_claimed:
            problems.append(f"note's Table 1 not reproduced at N = {n}: claimed {count_claimed} ({ratio_claimed}), computed {count} ({count / math.sqrt(n):.6f})")

    if problems:
        for p in problems:
            print("FAIL:", p)
        if args.plant:
            print(f"plant {args.plant}: caught")
        sys.exit(1)
    if args.plant:
        print(f"plant {args.plant}: NOT caught")
        sys.exit(2)

    # |A ∩ [1, N]| / √N is smallest just before each new term: at N = a_{j+1} - 1 the count is j.
    rows = []
    for k in CHECKPOINTS:
        lo = k // 2
        dips = [j / math.sqrt(g2[j] - 1) for j in range(lo, k)]
        row = {"k": k, "a_k": g2[k - 1], "ratio_at_a_k": round(k / math.sqrt(g2[k - 1]), 6), "window": [g2[lo - 1], g2[k - 1]], "min_ratio_in_window": round(min(dips), 6)}
        rows.append(row)
    for prev, row in zip(rows, rows[1:]):
        row["exponent_from_previous"] = round(math.log(row["a_k"] / prev["a_k"]) / math.log(row["k"] / prev["k"]), 4)

    OUT_LIST.write_text("".join(f"{x}\n" for x in g2))
    result = {
        "description": "Greedy B2[2] set a_1 = 1, a_{k+1} = least m > a_k keeping at most 2 representations a + b (a <= b) of every n. Illustrative only: no finite computation decides a liminf.",
        "thread_claim": {"source": "https://www.erdosproblems.com/forum/thread/158, comment of February 1, 2026", "k": THREAD_K, "a_k_claimed": THREAD_A, "a_k_computed": a_thread, "ratio_claimed": THREAD_RATIO, "ratio_computed": round(ratio_thread, 6), "reproduced": True},
        "note_table": {"source": "R. Zeraoulia, Zenodo, DOI 10.5281/zenodo.18452185, file erdos_problem_158_greedy_computation_note_v7.pdf, Table 1", "rows": note_rows, "reproduced": True},
        "g2": {"terms": K, "last": g2[-1], "vmax": VMAX, "seconds": round(t2, 2), "list": str(OUT_LIST.relative_to(ROOT)), "list_sha256": hashlib.sha256(OUT_LIST.read_bytes()).hexdigest(), "checkpoints": rows},
        "checks": {
            "g1_matches_oeis_A005282": {"terms": G1_TERMS, "a_last": g1[-1], "seconds": round(t1, 2), "ok": g1_ok},
            "python_reimplementation": {"terms": PY_TERMS, "seconds": round(t_py, 2), "ok": True},
            "b2_property": {"max_representations": maxrep, "seconds": round(t_b2, 2), "ok": True},
            "cpu_limit": {"limit_s": CPU_LIMIT_S, "trip_limit_s": 1, "trip_exit": rc_trip, "trip_seconds": round(t_trip, 2), "ok": trip_ok},
            "peak_rss_children_bytes": child.ru_maxrss,
        },
        "source_sha256": hashlib.sha256(SRC.read_bytes()).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    OUT_JSON.write_text(json.dumps(result, indent=2) + "\n")
    for row in rows:
        print(row)
    print(f"all checks passed; wrote {OUT_JSON.relative_to(ROOT)} and {OUT_LIST.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
