"""Exact equal-area counts for integer grids, via grid_count, summarised into results/lattice.json.

For each grid it records T_N = #triangles of doubled area N for the N that matter here, the best N, and a SHA-256 of
grid_count's full output so a rerun can be compared byte for byte. The raw histograms go to results/raw/ (not kept in git).

    python src/lattice_table.py            all runs (measured once: 12 min 41 s wall-clock on 10 threads, shared machine)
    python src/lattice_table.py --quick    skip the 800 x 800 grid
"""
import hashlib
import json
import math
import os
import subprocess
import sys

from sympy import divisor_sigma

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RAW = os.path.join(ROOT, "results", "raw")
FIXED_N = [1, 2, 6, 12, 60, 120, 360, 840, 2520]
C6 = 6 / math.pi**2
EULER_GAMMA = 0.57721566490153286061


def run(a, b, threads=10):
    path = os.path.join(RAW, f"grid_{a}x{b}.txt")
    if not os.path.exists(path):
        out = subprocess.run([os.path.join(HERE, "grid_count"), str(a), str(b), "0", str(threads)],
                             check=True, capture_output=True, text=True).stdout
        with open(path, "w") as f:
            f.write(out)
    text = open(path).read()
    hist = {}
    for line in text.splitlines():
        if line and not line.startswith("#"):
            d, c = line.split()
            hist[int(d)] = int(c)
    return hist, hashlib.sha256(text.encode()).hexdigest()


def summarise(a, b, label):
    hist, digest = run(a, b)
    n = a * b
    best = max(hist, key=lambda d: (hist[d], -d))
    rows = {}
    for N in FIXED_N + [best]:
        t = hist.get(N, 0)
        pred = C6 * int(divisor_sigma(N)) / N
        rows[str(N)] = {"T": t, "T_over_n2": t / n**2, "ratio_to_6sigma_over_pi2N": t / (n**2 * pred)}
    return {
        "label": label, "a": a, "b": b, "n": n, "sha256": digest,
        "total_triangles": sum(hist.values()),
        "best_N": best, "best_T": hist[best], "best_T_over_n2": hist[best] / n**2,
        "best_T_over_n2_loglogn": hist[best] / (n**2 * math.log(math.log(n))),
        "top5_N": sorted(hist, key=lambda d: (-hist[d], d))[:5],
        "fixed_N": rows,
    }


def main(argv):
    os.makedirs(RAW, exist_ok=True)
    squares = [25, 50, 100, 200, 400] + ([] if "--quick" in argv else [800])
    runs = [summarise(m, m, "square") for m in squares]
    for a in [200, 100, 40, 20, 10, 5, 3]:
        runs.append(summarise(a, 40000 // a, "rectangle, n about 40000"))
    for n in [10**4, 4 * 10**4, 16 * 10**4]:
        a = int(math.floor(math.sqrt(math.log(n))))  # a = [sqrt(log n)] rows and N = a!, as in [ErPu71, Thm 2]
        s = summarise(a, n // a, "ErPu71 strip")
        N = math.factorial(a)
        hist, _ = run(a, n // a)
        s["ErPu71_N"] = N
        s["ErPu71_T"] = hist.get(N, 0)
        s["ErPu71_T_over_n2"] = hist.get(N, 0) / (a * (n // a))**2
        runs.append(s)
    out = {
        "constant_6_over_pi2": C6,
        "constant_6_egamma_over_pi2": C6 * math.exp(EULER_GAMMA),
        "runs": runs,
    }
    with open(os.path.join(ROOT, "results", "lattice.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"6/pi^2 = {C6:.6f}   6 e^gamma/pi^2 = {C6 * math.exp(EULER_GAMMA):.6f}")
    for r in runs:
        extra = f"  ErPu71 N={r['ErPu71_N']}: T/n^2={r['ErPu71_T_over_n2']:.4f}" if "ErPu71_N" in r else ""
        print(f"{r['label']:26s} {r['a']:4d} x {r['b']:6d}  n={r['n']:7d}  T_1/n^2={r['fixed_N']['1']['T_over_n2']:.5f}"
              f"  best N={r['best_N']:5d} T/n^2={r['best_T_over_n2']:.4f}  T/(n^2 loglog n)={r['best_T_over_n2_loglogn']:.4f}{extra}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
