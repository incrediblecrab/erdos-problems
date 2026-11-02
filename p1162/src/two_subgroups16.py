"""Every 2-subgroup of S_16 and S_17: the census of two_subgroups.py extended by the orbit types of 16 points.

two_subgroups.py counts the 2-subgroups of S_n for n <= 15 from the fixed-point-free counts
t(k), k <= 14, by T(n) = sum_k C(n, k) t(k).  From n = 16 on, t(16) is needed, and its orbit
types include (16), whose Sylow subgroup R_16 has 2^15 elements, and (8, 8), whose R has 2^14.
This driver runs src/census16.c on each orbit type of 16 points, which lists the subgroups of
R = prod_i R_{b_i} transitive on every block with their Sylow numbers s_{b_i}, and forms

    c(type) = prod_i n_2(b_i) * sum_H 1 / prod_i s_{b_i}(H|O_i),     n_2(b) = b!/2^(b-1),

exactly as two_subgroups.py does (NOTES.md §9.8).  t(k) for k <= 14 is read from
data/two_subgroups.json, so T(16) and T(17) follow; odd n adds no new orbit type.

A transitive 2-subgroup is Frattini-closed (FC) by definition (its one orbit gives N = Phi(H)),
so the FC 2-subgroups of S_16 with an orbit longer than 8 are exactly the transitive ones of
degree 16, and those of S_17 are 17 times as many.  This is the share premise (c) of
Conjecture 9.3 asks to be negligible.

Checks (exit 1 if any fails):
  * census16 reproduces, for every orbit type of k <= 14, the stored list length and the
    weighted counts c(type) and c_FC(type) of data/two_subgroups.json (the Python enumeration);
  * every weighted sum is an integer;
  * the binomial transform of the stored t(k) gives the stored rows for n <= 15;
  * n = 16, 17: the number of 2-subgroups of each order 2^j, j <= v_2(n!), is odd (Frobenius),
    and the number of order 2^{v_2(n!)} is n!/2^{v_2(n!)}; the elementary abelian counts by
    rank equal elemab.rank_counts(n); the FC counts by order with every orbit of size <= 8 equal
    block_families.order_distribution; 2-subgroups + |S_3(n)| <= f(n) = A005432(n);
  * exactly the orbit types with two or more orbits of length >= 4 have subgroups that are not FC;
  * d(P) <= 6 for every transitive P <= R_16, the bound b/4 + 2 of [RoTr25] Prop. 5.3 as quoted
    in NOTES.md §9.7.

elemab.py and block_families.py are imported only to compare final numbers.  The independent
recount is src/verify_two_subgroups16.py, which shares no code with this file or census16.c.
"""

from __future__ import annotations

import argparse
import json
import os
import pathlib
import subprocess
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction
from math import comb, factorial

ROOT = pathlib.Path(__file__).resolve().parent.parent
CHECKS = [0]
SIZES = (16, 8, 4, 2)
FIELDS_INT = ("total", "fc", "non_fc", "elementary_abelian", "fc_with_orbit_longer_than_8")
FIELDS_HIST = ("by_order", "fc_by_order", "elementary_abelian_by_order", "by_largest_orbit", "fc_by_largest_orbit")


def orbit_types(k: int, top: int = 16) -> list[tuple[int, ...]]:
    """Multisets of 2-power block sizes >= 2 summing to k, largest first."""
    if k == 0:
        return [()]
    return [(b,) + rest for b in SIZES if b <= min(k, top) for rest in orbit_types(k - b, b)]


def set_partitions(blocks: tuple[int, ...]) -> int:
    den = 1
    for b, m in Counter(blocks).items():
        den *= factorial(b) ** m * factorial(m)
    return factorial(sum(blocks)) // den


def build(binary: pathlib.Path) -> None:
    src = ROOT / "src" / "census16.c"
    if binary.exists() and binary.stat().st_mtime >= src.stat().st_mtime:
        return
    cc = os.environ.get("CC", "cc")
    subprocess.run([cc, "-O3", "-std=c99", "-Wall", "-o", str(binary), str(src)], check=True)


def run_types(binary: pathlib.Path, types: list[tuple[int, ...]], jobs: int) -> dict:
    """Run census16 on each type, at most `jobs` at a time; return the parsed output per type."""
    out, pending, running = {}, list(types), []
    while pending or running:
        while pending and len(running) < jobs:
            t = pending.pop(0)
            p = subprocess.Popen([str(binary), *map(str, t)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            running.append((t, p, time.time()))
        t, p, t0 = running.pop(0)
        so, se = p.communicate()
        if p.returncode != 0:
            raise SystemExit(f"census16 {t} exited with {p.returncode}: {se.strip()}")
        keys, layers, listed, cpu = [], {}, None, None
        for line in so.splitlines():
            f = line.split()
            if f[0] == "key":
                v = list(map(int, f[1:]))
                keys.append((v[0], bool(v[1]), bool(v[2]), v[3], tuple(v[4:4 + len(t)]), v[4 + len(t)]))
            elif f[0] == "layer":
                layers[int(f[1])] = int(f[2])
            elif f[0] == "listed":
                listed = int(f[1])
            elif f[0] == "cpu_seconds":
                cpu = float(f[1])
        out[t] = {"keys": keys, "layers": layers, "listed": listed, "cpu_seconds": cpu, "wall_seconds": round(time.time() - t0, 1)}
    return out


def weigh(t: tuple[int, ...], keys) -> tuple[dict, bool]:
    """c(type) by (log2|H|, FC, elementary abelian); and whether every sum is an integer."""
    scale = 1
    for b in t:
        scale *= factorial(b) // 2 ** (b - 1)
    c = defaultdict(Fraction)
    if any(x < 1 for j, fc, ea, d, s, cnt in keys for x in s):  # every P lies in some Sylow subgroup
        return {}, False
    for j, fc, ea, d, s, cnt in keys:
        w = Fraction(scale * cnt)
        for x in s:
            w /= x
        c[(j, fc, ea)] += w
    ok = all(v.denominator == 1 for v in c.values())
    return {k: int(v) for k, v in c.items()}, ok


def summarise(C: Counter) -> dict:
    """The fields of data/two_subgroups.json from counts keyed (log2|H|, FC, elementary abelian, largest orbit)."""

    def hist(sel):
        h = Counter()
        for (j, fc, ea, mo), v in C.items():
            if sel(fc, ea):
                h[j] += v
        return {str(2**j): h[j] for j in sorted(h)}

    mo_all, mo_fc = Counter(), Counter()
    for (j, fc, ea, mo), v in C.items():
        mo_all[mo] += v
        if fc:
            mo_fc[mo] += v
    return {
        "total": sum(C.values()),
        "fc": sum(v for (j, fc, ea, mo), v in C.items() if fc),
        "non_fc": sum(v for (j, fc, ea, mo), v in C.items() if not fc),
        "elementary_abelian": sum(v for (j, fc, ea, mo), v in C.items() if ea),
        "by_order": hist(lambda fc, ea: True),
        "fc_by_order": hist(lambda fc, ea: fc),
        "elementary_abelian_by_order": hist(lambda fc, ea: ea),
        "by_largest_orbit": {str(m): mo_all[m] for m in sorted(mo_all)},
        "fc_by_largest_orbit": {str(m): mo_fc[m] for m in sorted(mo_fc)},
        "fc_with_orbit_longer_than_8": sum(v for m, v in mo_fc.items() if m > 8),
    }


def combine(parts: list[tuple[int, dict]]) -> dict:
    """sum of coefficient * summary, field by field."""
    out = {f: 0 for f in FIELDS_INT}
    hs = {f: Counter() for f in FIELDS_HIST}
    for c, s in parts:
        for f in FIELDS_INT:
            out[f] += c * s[f]
        for f in FIELDS_HIST:
            for key, v in s[f].items():
                hs[f][key] += c * v
    for f in FIELDS_HIST:
        out[f] = {key: hs[f][key] for key in sorted(hs[f], key=int) if hs[f][key]}
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--bin", type=pathlib.Path, default=ROOT / "src" / "census16")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--census", type=pathlib.Path, default=ROOT / "data" / "census8.json")
    ap.add_argument("--small", type=pathlib.Path, default=ROOT / "data" / "two_subgroups.json")
    ap.add_argument("--out", type=pathlib.Path, default=ROOT / "data" / "two_subgroups16.json")
    args = ap.parse_args()
    t0 = time.time()
    failures: list[str] = []

    def check(ok: bool, label: str) -> None:
        CHECKS[0] += 1
        print(f"  [{'pass' if ok else 'FAIL'}] {label}")
        if not ok:
            failures.append(label)

    build(args.bin)
    small = json.loads(args.small.read_text())
    if small["n_max"] < 15:
        raise SystemExit(f"{args.small} stops at n = {small['n_max']}; this needs the census to n = 15")

    print("== census16 against data/two_subgroups.json, every orbit type of k <= 14 ==")
    old_types = [tuple(map(int, name.split(","))) for name in small["types"]]
    res = run_types(args.bin, old_types, args.jobs)
    for t in old_types:
        s = small["types"][",".join(map(str, t))]
        c, ok = weigh(t, res[t]["keys"])
        tot, fc = sum(c.values()), sum(v for (j, f, ea), v in c.items() if f)  # c is empty when not ok
        check(ok and res[t]["listed"] == s["listed_in_R"] and tot == s["per_partition"] and fc == s["per_partition_fc"], f"type {t}: {res[t]['listed']} listed, c = {tot}, c_FC = {fc}, as stored")

    print("\n== the orbit types of 16 points ==")
    new_types = orbit_types(16)
    res16 = run_types(args.bin, new_types, args.jobs)
    t16 = Counter()
    per_type = {}
    for t in new_types:
        r = res16[t]
        c, ok = weigh(t, r["keys"])
        check(ok, f"type {t}: every Sylow number is positive and every weighted sum is an integer")
        mult = set_partitions(t)
        for (j, fc, ea), v in c.items():
            t16[(j, fc, ea, max(t))] += mult * v
        per_type[",".join(map(str, t))] = {
            "listed_in_R": r["listed"],
            "per_partition": sum(c.values()),
            "per_partition_fc": sum(v for (j, f, ea), v in c.items() if f),
            "partitions": mult,
            "sylow_order": 2 ** sum(b - 1 for b in t),
            "listed_by_log2_order": {str(j): r["layers"][j] for j in sorted(r["layers"])},
            "cpu_seconds": r["cpu_seconds"],
        }
        print(f"  {str(t):26s} {r['listed']:8d} listed; c = {sum(c.values())} per partition, x {mult} partitions  ({r['cpu_seconds']:.1f}s cpu)")
    two_big = {name: sum(1 for b in name.split(",") if int(b) >= 4) >= 2 for name in per_type}
    check(all((v["per_partition_fc"] < v["per_partition"]) == two_big[name] for name, v in per_type.items()), "k = 16: exactly the types with two or more orbits of length >= 4 have non-FC subgroups")

    # the transitive subgroups of R_16 in detail
    keys16 = res16[(16,)]["keys"]
    n2 = factorial(16) // 2**15
    by_order, by_d, s_listed, d_listed = defaultdict(Fraction), defaultdict(Fraction), Counter(), Counter()
    for j, fc, ea, d, s, cnt in keys16:
        if s[0] < 1:
            continue  # already failed in weigh()
        by_order[j] += Fraction(n2 * cnt, s[0])
        by_d[d] += Fraction(n2 * cnt, s[0])
        s_listed[s[0]] += cnt
        d_listed[d] += cnt
    max_d = max(d_listed, default=0)
    check(max_d <= 6, f"every transitive P <= R_16 has d(P) <= 6 = 16/4 + 2 (largest d: {max_d})")
    transitive16 = {
        "n_2": n2,
        "listed": res16[(16,)]["listed"],
        "total": int(sum(by_order.values())),
        "by_order": {str(2**j): int(by_order[j]) for j in sorted(by_order)},
        "by_d": {str(d): int(by_d[d]) for d in sorted(by_d)},
        "listed_by_s": {str(s): s_listed[s] for s in sorted(s_listed)},
        "listed_by_d": {str(d): d_listed[d] for d in sorted(d_listed)},
        "max_d": max_d,
    }
    print(f"  transitive 2-subgroups of S_16: {transitive16['total']}  (largest d = {max_d})")

    print("\n== binomial transform: the stored t(k), k <= 14, and the new t(16) ==")
    fpf = {int(k): v for k, v in small["fixed_point_free"].items()}
    fpf[0] = summarise(Counter({(0, True, True, 1): 1}))
    fpf[16] = summarise(t16)
    for n in range(1, 16):
        row = combine([(comb(n, k), fpf[k]) for k in fpf if k <= n])
        check(row == small["two_subgroups"][str(n)], f"n={n}: the binomial transform of the stored t(k) gives the stored row")
    rows = {n: combine([(comb(n, k), fpf[k]) for k in fpf if k <= n]) for n in (16, 17)}

    print("\n== n = 16, 17: Frobenius, the Sylow count, elemab.py, block_families.py, A005432 ==")
    sys.path.insert(0, str(ROOT / "src"))
    from elemab import rank_counts  # noqa: E402
    from block_families import A005432, L, egf, load_types, order_distribution, pmul, weights  # noqa: E402

    census = json.loads(args.census.read_text())
    types_all = load_types(census)
    dist = order_distribution(types_all, 17)
    F = egf(weights(types_all), 17)
    for n, r in rows.items():
        v2 = sum(n // 2**i for i in range(1, 6))
        by = {int(k).bit_length() - 1: v for k, v in r["by_order"].items()}
        check(sorted(by) == list(range(v2 + 1)) and all(v % 2 == 1 for v in by.values()) and by.get(v2) == factorial(n) // 2**v2, f"n={n}: every order 2^0..2^{v2} occurs an odd number of times; {by.get(v2)} of order 2^{v2} = n!/2^{v2}")
        ea = {int(k).bit_length() - 1: v for k, v in r["elementary_abelian_by_order"].items()}
        rc = rank_counts(n)
        check([ea.get(j, 0) for j in range(len(rc))] == rc and sum(ea.values()) == sum(rc), f"n={n}: elementary abelian by rank = elemab.rank_counts")
        # FC with every orbit of size <= 8: subtract the C(n, 16) * (transitive of degree 16) by order
        fc8 = Counter({int(k).bit_length() - 1: v for k, v in r["fc_by_order"].items()})
        for k, v in transitive16["by_order"].items():
            fc8[int(k).bit_length() - 1] -= comb(n, 16) * v
        fc8 = {j: v for j, v in fc8.items() if v}
        check(fc8 == dist[n], f"n={n}: FC by order, orbits of size <= 8 only, = block_families.order_distribution (B(n) = {sum(dist[n].values())})")
        check(r["fc_with_orbit_longer_than_8"] == comb(n, 16) * transitive16["total"], f"n={n}: FC with an orbit longer than 8 = C(n,16) x the transitive count of degree 16")
        s3 = comb(n, 3) * L(pmul([0, 1], F[n - 3][0]))
        r["s3"] = s3
        r["fc_orbits_at_most_8"] = sum(dist[n].values())
        check(r["total"] + s3 <= A005432[n], f"n={n}: 2-subgroups + |S3(n)| = {r['total'] + s3} <= f(n) = {A005432[n]}")

    print("\n== the census ==")
    for n, r in rows.items():
        f = A005432[n]
        print(f"  n={n}: {r['total']} 2-subgroups, {r['fc']} FC, {r['non_fc']} not FC (share {r['non_fc'] / r['total']:.6f}); 2-subgroups/f = {r['total'] / f:.6f}")
        print(f"        FC with an orbit longer than 8: {r['fc_with_orbit_longer_than_8']}, share of FC {r['fc_with_orbit_longer_than_8'] / r['fc']:.6f}, of f {r['fc_with_orbit_longer_than_8'] / f:.6f}")

    out = {
        "description": "Every 2-subgroup of S_16 and S_17, from census16.c on the orbit types of 16 points and data/two_subgroups.json for k <= 14 (src/two_subgroups16.py)",
        "types": per_type,
        "transitive16": transitive16,
        "fixed_point_free_16": fpf[16],
        "two_subgroups": {str(n): r for n, r in rows.items()},
        "checks_failed": failures,
        "checks_run": CHECKS[0],
        "cpu_seconds_census16": round(sum(v["cpu_seconds"] for v in res16.values()), 1),
        "seconds": round(time.time() - t0, 1),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=1) + "\n")
    shown = args.out.resolve().relative_to(ROOT) if args.out.resolve().is_relative_to(ROOT) else args.out
    print(f"\nwrote {shown}  ({time.time() - t0:.1f}s)")
    print(f"{CHECKS[0]} checks, {len(failures)} failed" if failures else f"{CHECKS[0]} checks, all passed")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
