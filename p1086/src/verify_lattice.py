"""Recheck results/lattice.json from the raw histograms, without importing lattice_table.py.

For every run: the SHA-256 of its raw file; the histogram total against C(ab, 3) minus the collinear triples (a gcd
sum), and against OEIS A045996 for the squares; best N, its count and every fixed-N count against a fresh scan; the
ratio to (6/pi^2) sigma(N)/N; and, for the ErPu71 strips, that the row count is floor(sqrt(log n)) and N = rows!.
Exits 1 on any mismatch.

    python src/verify_lattice.py
"""
import hashlib
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
# https://oeis.org/A045996/b045996.txt, read September 23, 2026
A045996 = {25: 40212160, 50: 2595202884, 100: 166502148020, 200: 10663697056080,
           400: 682613759714616, 800: 43689733853321032}


def collinear(a, b):
    """Collinear triples of the a x b grid: each is counted once, by its outer pair."""
    s = 0
    for dx in range(a):
        for dy in range(-(b - 1), b):
            if dx > 0 or dy > 0:
                s += (a - dx) * (b - abs(dy)) * (math.gcd(dx, abs(dy)) - 1)
    return s


def sigma(N):
    return sum(d + (N // d if d * d != N else 0) for d in range(1, math.isqrt(N) + 1) if N % d == 0)


def main():
    data = json.loads((ROOT / "results" / "lattice.json").read_text())
    bad = []
    print(f"{'run':<26}{'grid':>12}{'n':>8}{'T_1/n^2':>10}{'best N':>8}{'T/n^2':>9}{'T/(n^2 lnln n)':>16}")
    for r in data["runs"]:
        a, b, n = r["a"], r["b"], r["n"]
        raw = ROOT / "results" / "raw" / f"grid_{a}x{b}.txt"
        blob = raw.read_bytes()
        tag = f"{r['label']} {a}x{b}"
        if hashlib.sha256(blob).hexdigest() != r["sha256"]:
            bad.append(f"{tag}: sha256")
        hist, header = {}, None
        for line in blob.decode().splitlines():
            if line.startswith("#"):
                header = dict(kv.split("=") for kv in line[1:].split())
            elif line.strip():
                d, c = line.split()
                hist[int(d)] = int(c)
        if header is None or int(header["a"]) != a or int(header["b"]) != b or int(header["dmax"]) < (a - 1) * (b - 1):
            bad.append(f"{tag}: header {header} does not cover every doubled area")
        total = sum(hist.values())
        if n != a * b or total != r["total_triangles"] or total != math.comb(a * b, 3) - collinear(a, b):
            bad.append(f"{tag}: total {total}")
        if a == b and a in A045996 and total != A045996[a]:
            bad.append(f"{tag}: total differs from A045996({a})")
        top = max(hist.values())
        argmax = [d for d, c in hist.items() if c == top]
        if r["best_T"] != top or argmax != [r["best_N"]]:
            bad.append(f"{tag}: best {r['best_N']} {r['best_T']} vs {argmax} {top}")
        for N, v in r["fixed_N"].items():
            N = int(N)
            if v["T"] != hist.get(N, 0):
                bad.append(f"{tag}: T_{N}")
            ratio = hist.get(N, 0) / n**2 / (6 / math.pi**2 * sigma(N) / N)
            if not math.isclose(ratio, v["ratio_to_6sigma_over_pi2N"], rel_tol=1e-12, abs_tol=1e-15):
                bad.append(f"{tag}: ratio for N={N}")
        if "ErPu71_N" in r:
            rows = min(a, b)
            if rows != int(math.sqrt(math.log(n))):
                bad.append(f"{tag}: rows != floor(sqrt(log n))")
            if r["ErPu71_N"] != math.factorial(rows) or r["ErPu71_T"] != hist.get(r["ErPu71_N"], 0):
                bad.append(f"{tag}: ErPu71 N or T")
        print(f"{r['label']:<26}{f'{a}x{b}':>12}{n:>8}{hist.get(1, 0) / n**2:>10.5f}{argmax[0]:>8}"
              f"{top / n**2:>9.4f}{top / (n**2 * math.log(math.log(n))):>16.4f}")
    print(f"6/pi^2 = {6 / math.pi**2:.6f}, 6 e^gamma/pi^2 = {6 * math.exp(0.5772156649015329) / math.pi**2:.6f}")
    for m in bad:
        print("MISMATCH", m)
    print(f"{len(data['runs'])} runs checked, {len(bad)} mismatches")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
