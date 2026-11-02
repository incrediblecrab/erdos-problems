"""Exact census of the 2-subgroups of S_n for n <= 15, split by Frattini-closure and orbit sizes.

Every 2-subgroup H <= S_n has orbits of 2-power size.  Its fixed points can be
split off: if t(k) is the number of fixed-point-free 2-subgroups of S_k, the
number of 2-subgroups of S_n is sum_k C(n, k) t(k).  A fixed-point-free H
whose orbits form a set partition of [k] of type (b_1, ..., b_t) is a 2-subgroup
of S = prod_i Sym(O_i) that is transitive on every O_i, so

    t(k) = sum over types  k! / prod_b (b!^{m_b} m_b!)  *  c(b_1, ..., b_t),

c(type) being the number of such subgroups for one fixed partition.  c(type) is
computed inside one Sylow 2-subgroup R = prod_i R_{b_i} of S, where R_b is the
iterated wreath product C2 wr ... wr C2 on b points.  The Sylow 2-subgroups of
S are the products of Sylow 2-subgroups of the factors and are all conjugate
under S, and S preserves "transitive on every block", so double counting gives

    c(type) = prod_i n_2(b_i) * sum_{H <= R transitive on every block} 1 / s(H),
    s(H) = prod_i s_{b_i}(H|O_i),

where n_2(b) = b!/|R_b| is the number of Sylow 2-subgroups of S_b and s_b(P)
the number of them that contain P.  For transitive P <= S_b, b = 2^j, s_b(P) is
the number of maximal chains in the interval [P_x, P] of the subgroup lattice
of P (P_x a point stabiliser).  The reason: a Sylow 2-subgroup of S_b is the
full stabiliser of a nest of partitions into blocks of sizes 2, 4, ..., b/2
(R_b is self-normalising, so the correspondence is a bijection); P lies in it
iff P preserves every level; a P-invariant partition is a block system of P,
and block systems of P correspond to subgroups between P_x and P.

The subgroups of R transitive on every block are listed top-down: in a 2-group
every subgroup is reached from the whole group by passing to maximal subgroups,
the maximal subgroups of M are the preimages of the hyperplanes of
M/Phi(M), Phi(M) = <squares>, and a subgroup of an intransitive group is
intransitive, so the descent is pruned at the first block on which a
candidate is not transitive.  Each H is classified by order, by whether it is
elementary abelian, and by whether it is Frattini-closed (FC): H >= Phi(H|O)
for every orbit O, tested as "every restricted square (h^2)|O lies in H".

For n <= 15 every orbit has size at most 8 (the next power of 2 is 16), so
the FC subgroups with an orbit longer than 8 -- premise (c) of NOTES.md §9.8 --
do not exist in this range; the census cannot test (c).

Checks (exit 1 if any fails):
  * sum_P n_2(b)/s_b(P) over the transitive P <= R_b equals the number of
    transitive 2-subgroups of S_b found by brute force in data/census8.json;
  * every weighted sum is an integer;
  * n <= 8: counts by order, FC counts by order and elementary abelian counts
    equal data/census8.json (brute force on permutations) and the 2-power
    orders of data/lattice_n{n}.json (src/subgroups.c);
  * every n: the number of 2-subgroups of each order 2^j, j <= v_2(n!), is odd
    (Frobenius), the number of order 2^{v_2(n!)} is n!/2^{v_2(n!)}, the
    elementary abelian counts by rank equal elemab.rank_counts(n), and the FC
    counts by order equal block_families.order_distribution (a formula on the
    transitive types that shares no code with this enumeration);
  * every n: 2-subgroups + |S_3(n)| <= f(n) = A005432(n).

elemab.py and block_families.py are imported only to compare final numbers.
The independent check is src/verify_two_subgroups.py, which shares no code
with this file.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction
from math import comb, factorial

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
BLOCK_SIZES = (2, 4, 8)
N_MAX_SUPPORTED = 15


# --------------------------------------------------------- the groups R_b


def wreath_generators(b: int) -> list[tuple[int, ...]]:
    """(0 1), (0 2)(1 3), (0 4)(1 5)(2 6)(3 7), ... restricted to b points."""
    gens, s = [], 1
    while s < b:
        p = list(range(b))
        for i in range(s):
            p[i], p[i + s] = i + s, i
        gens.append(tuple(p))
        s *= 2
    return gens


def closure(gens: list[tuple[int, ...]], b: int) -> list[tuple[int, ...]]:
    ident = tuple(range(b))
    elems, seen, frontier = [ident], {ident}, [ident]
    while frontier:
        nxt = []
        for p in frontier:
            for g in gens:
                q = tuple(p[i] for i in g)  # p after g
                if q not in seen:
                    seen.add(q)
                    elems.append(q)
                    nxt.append(q)
        frontier = nxt
    return elems


class Ambient:
    """R = prod_i R_{b_i}, elements indexed by mixed radix over the factors; index 0 is the identity."""

    def __init__(self, blocks: tuple[int, ...]):
        self.blocks = blocks
        self.factors = [closure(wreath_generators(b), b) for b in blocks]
        sizes = [len(f) for f in self.factors]
        self.N = int(np.prod(sizes))
        ar = np.arange(self.N)
        strides, s = [], 1
        for z in sizes:
            strides.append(s)
            s *= z
        self.comp = [(ar // st) % z for st, z in zip(strides, sizes)]
        self.mul = np.zeros((self.N, self.N), dtype=np.int32)
        self.sq = np.zeros(self.N, dtype=np.int64)
        self.restsq, self.base_image = [], []
        for i, els in enumerate(self.factors):
            idx = {p: k for k, p in enumerate(els)}
            table = np.array([[idx[tuple(a[j] for j in c)] for c in els] for a in els], dtype=np.int32)
            ci = self.comp[i]
            self.mul += table[ci[:, None], ci[None, :]] * strides[i]
            fsq = np.diagonal(table).astype(np.int64)
            self.sq += fsq[ci] * strides[i]
            self.restsq.append(fsq[ci] * strides[i])  # (h^2)|O_i, identity elsewhere
            self.base_image.append(np.array([p[0] for p in els])[ci])  # image of the first point of O_i


_PARITY: dict[int, np.ndarray] = {}


def parity(r: int) -> np.ndarray:
    """parity[l-1, v] = <l, v> over F_2, for l = 1..2^r-1, v = 0..2^r-1."""
    if r not in _PARITY:
        x = np.arange(1, 1 << r)[:, None] & np.arange(1 << r)[None, :]
        p = np.zeros_like(x)
        while x.any():
            p ^= x & 1
            x >>= 1
        _PARITY[r] = p
    return _PARITY[r]


def frattini_coordinates(A: Ambient, els: np.ndarray) -> tuple[np.ndarray, int, int]:
    """A map M -> F_2^r with kernel Phi(M) = <squares>, as an array over all of R (-1 off M)."""
    squares = np.unique(A.sq[els])
    phi = np.zeros(A.N, dtype=bool)
    phi[squares] = True
    gens = squares[squares != 0]
    frontier = np.nonzero(phi)[0]
    while len(frontier) and len(gens):
        new = np.unique(A.mul[np.ix_(frontier, gens)].ravel())
        new = new[~phi[new]]
        phi[new] = True
        frontier = new
    coord = np.full(A.N, -1, dtype=np.int64)
    coord[phi] = 0
    assigned, r = np.nonzero(phi)[0], 0
    for x in els:
        if coord[x] < 0:
            ys = A.mul[assigned, x]
            coord[ys] = coord[assigned] | (1 << r)
            assigned = np.concatenate([assigned, ys])
            r += 1
    assert len(assigned) == len(els) and (1 << r) * int(phi.sum()) == len(els), "M/Phi(M) is not elementary abelian"
    return coord, r, int(phi.sum())


def descend(A: Ambient, transitive_only: bool):
    """Yield the element arrays of the subgroups of R, top-down, each once.

    With transitive_only, only the subgroups transitive on every block."""
    full = np.ones(A.N, dtype=bool)
    layer = [np.packbits(full).tobytes()]
    while layer:
        nxt: dict[bytes, None] = {}
        for key in layer:
            els = np.nonzero(np.unpackbits(np.frombuffer(key, dtype=np.uint8), count=A.N))[0]
            yield els
            coord, r, _ = frattini_coordinates(A, els)
            if r == 0:
                continue
            par = parity(r)
            cv = coord[els]
            ok = np.ones(par.shape[0], dtype=bool)
            if transitive_only:
                for i, b in enumerate(A.blocks):
                    im = A.base_image[i][els]
                    for p in range(1, b):  # the identity covers p = 0 in every hyperplane
                        hit = np.zeros(1 << r, dtype=bool)
                        hit[cv[im == p]] = True
                        ok &= (par[:, hit] == 0).any(axis=1)
            for li in np.nonzero(ok)[0]:
                mem = np.zeros(A.N, dtype=bool)
                mem[els[par[li, cv] == 0]] = True
                nxt[np.packbits(mem).tobytes()] = None
        layer = list(nxt)


# ----------------------------------------------- s_b(P) by maximal chains


def sylow_weights(b: int) -> tuple[int, dict[bytes, tuple[int, int]], dict]:
    """n_2(b) and, for each transitive P <= R_b, (id, s_b(P)); plus a summary."""
    A = Ambient((b,))
    subs = []
    for els in descend(A, transitive_only=False):
        mask = 0
        for x in els.tolist():
            mask |= 1 << x
        subs.append((len(els), mask, els))
    stab0 = 0
    for x in np.nonzero(A.base_image[0] == 0)[0].tolist():
        stab0 |= 1 << x
    by_order = defaultdict(list)
    for order, mask, _ in subs:
        by_order[order].append(mask)
    weights: dict[bytes, tuple[int, int]] = {}
    chains_total = Counter()
    for order, mask, els in subs:
        if len(np.unique(A.base_image[0][els])) != b:
            continue
        bottom = mask & stab0
        ways = {bottom: 1}
        o = bottom.bit_count() * 2
        while o <= order:
            for L in by_order[o]:
                if L & ~mask == 0 and L & bottom == bottom:
                    ways[L] = sum(w for K, w in ways.items() if K.bit_count() * 2 == o and K & ~L == 0)
            o *= 2
        mem = np.zeros(A.N, dtype=bool)
        mem[els] = True
        weights[np.packbits(mem).tobytes()] = (len(weights), ways[mask])
        chains_total[ways[mask]] += 1
    n2 = factorial(b) // A.N
    return n2, weights, {"sylow_order": A.N, "subgroups_of_R_b": len(subs), "transitive_in_R_b": len(weights), "s_values": dict(sorted(chains_total.items()))}


# ------------------------------------------------------------- the census


def types_of(k: int) -> list[tuple[int, ...]]:
    """Multisets of block sizes in {2, 4, 8} summing to k, largest first."""
    out = []

    def rec(rest: int, top: int, acc: list[int]) -> None:
        if rest == 0:
            out.append(tuple(acc))
            return
        for b in BLOCK_SIZES:
            if b <= top and b <= rest:
                rec(rest - b, b, acc + [b])

    rec(k, 8, [])
    return out


def partitions_of_type(blocks: tuple[int, ...]) -> int:
    k = sum(blocks)
    den = 1
    for b, m in Counter(blocks).items():
        den *= factorial(b) ** m * factorial(m)
    return factorial(k) // den


def count_type(blocks, sw) -> tuple[dict, int]:
    """c(type) split by (log2|H|, FC, elementary abelian): exact integers."""
    A = Ambient(blocks)
    raw = Counter()
    listed = 0
    for els in descend(A, transitive_only=True):
        listed += 1
        mem = np.zeros(A.N, dtype=bool)
        mem[els] = True
        fc = all(bool(mem[A.restsq[i][els]].all()) for i, b in enumerate(blocks) if b >= 4)
        ea = bool((A.sq[els] == 0).all())
        proj = []
        for i, b in enumerate(blocks):
            pm = np.zeros(len(A.factors[i]), dtype=bool)
            pm[A.comp[i][els]] = True
            proj.append(sw[b][1][np.packbits(pm).tobytes()][1])
        raw[(len(els).bit_length() - 1, fc, ea, tuple(proj))] += 1
    scale = 1
    for b in blocks:
        scale *= sw[b][0]
    out = defaultdict(Fraction)
    for (j, fc, ea, proj), cnt in raw.items():
        w = Fraction(scale)
        for s in proj:
            w /= s
        out[(j, fc, ea)] += cnt * w
    return out, listed


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--n-max", type=int, default=N_MAX_SUPPORTED)
    ap.add_argument("--census", type=pathlib.Path, default=ROOT / "data" / "census8.json")
    ap.add_argument("--lattice-dir", type=pathlib.Path, default=ROOT / "data")
    ap.add_argument("--out", type=pathlib.Path, default=ROOT / "data" / "two_subgroups.json")
    args = ap.parse_args()
    n_max = args.n_max
    if not 1 <= n_max <= N_MAX_SUPPORTED:
        ap.error(f"--n-max must be between 1 and {N_MAX_SUPPORTED}: from n = 16 on an orbit can have 16 points")
    t0 = time.time()
    failures: list[str] = []

    def check(ok: bool, label: str) -> None:
        print(f"  [{'pass' if ok else 'FAIL'}] {label}")
        if not ok:
            failures.append(label)

    census = json.loads(args.census.read_text())

    print("== s_b(P): Sylow 2-subgroups of S_b containing a transitive P, by maximal chains in [P_x, P] ==")
    sw = {}
    for b in BLOCK_SIZES:
        n2, wts, summary = sylow_weights(b)
        sw[b] = (n2, wts, summary)
        total = sum(Fraction(n2, s) for _, s in wts.values())
        want = census["transitive"][str(b)]["transitive_total"]
        print(f"  b={b}: n_2 = {n2}, R_b has {summary['subgroups_of_R_b']} subgroups, {summary['transitive_in_R_b']} transitive; s values {summary['s_values']}")
        check(total == want, f"b={b}: sum of n_2/s_b(P) over transitive P <= R_b = {total} = census8 transitive count {want}")
        check(n2 == factorial(b) // 2 ** (b - 1), f"b={b}: n_2 = b!/2^(b-1)")

    print(f"\n== fixed-point-free 2-subgroups t(k), k <= {n_max}, by orbit type ==")
    t: dict[int, Counter] = {0: Counter({(0, True, True, 1): 1})}  # key: (log2|H|, FC, elementary abelian, largest orbit)
    per_type = {}
    for k in range(2, n_max + 1, 2):
        t[k] = Counter()
        for blocks in types_of(k):
            t1 = time.time()
            c, listed = count_type(blocks, sw)
            ok_int = all(v.denominator == 1 for v in c.values())
            check(ok_int, f"type {blocks}: every weighted sum is an integer")
            mult = partitions_of_type(blocks)
            ctot = sum(c.values())
            per_type[",".join(map(str, blocks))] = {
                "listed_in_R": listed,
                "per_partition": int(ctot),
                "per_partition_fc": int(sum(v for (j, fc, ea), v in c.items() if fc)),
                "partitions": mult,
            }
            for (j, fc, ea), v in c.items():
                t[k][(j, fc, ea, max(blocks))] += mult * int(v)
            print(f"  {str(blocks):16s} {listed:7d} subgroups of R listed; c = {int(ctot)} per partition, x {mult} partitions  ({time.time() - t1:.1f}s)")
    for k in range(1, n_max + 1, 2):
        t[k] = Counter()

    T: dict[int, Counter] = {}
    for n in range(0, n_max + 1):
        T[n] = Counter()
        for k in range(0, n + 1):
            for key, v in t[k].items():
                T[n][key] += comb(n, k) * v

    def summarise(C: Counter) -> dict:
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

    rows = {n: summarise(T[n]) for n in range(1, n_max + 1)}
    fpf = {k: summarise(t[k]) for k in range(2, n_max + 1, 2)}

    print("\n== n <= 8 against the brute force (data/census8.json) and the lattices (src/subgroups.c) ==")
    for n in range(1, min(8, n_max) + 1):
        r = rows[n]
        want = census["two_subgroups_by_order"][str(n)]
        fcw = census["frattini_closed"][str(n)]
        check(
            r["by_order"] == want and r["fc"] == fcw["count"] and r["fc_by_order"] == fcw["order_histogram"] and r["elementary_abelian"] == fcw["elementary_abelian"] and r["total"] == fcw["two_subgroups"],
            f"n={n}: {r['total']} 2-subgroups, {r['fc']} FC, {r['elementary_abelian']} elementary abelian; all orders equal census8.json",
        )
        lat = json.loads((args.lattice_dir / f"lattice_n{n}.json").read_text())
        lw = {k: v for k, v in lat["order_histogram"].items() if int(k) & (int(k) - 1) == 0}
        check(r["by_order"] == dict(sorted(lw.items(), key=lambda kv: int(kv[0]))), f"n={n}: 2-subgroups by order equal lattice_n{n}.json")

    print("\n== every n: Frobenius parity, the Sylow count, elemab.py, block_families.py, A005432 ==")
    sys.path.insert(0, str(ROOT / "src"))
    from elemab import rank_counts  # noqa: E402
    from block_families import A005432, L, egf, load_types, order_distribution, pmul, weights  # noqa: E402

    types_all = load_types(census)
    dist = order_distribution(types_all, n_max)
    F = egf(weights(types_all), n_max)
    s3 = {n: (comb(n, 3) * L(pmul([0, 1], F[n - 3][0])) if n >= 3 else 0) for n in range(n_max + 1)}
    for n in range(1, n_max + 1):
        r = rows[n]
        v2 = sum(n // 2**i for i in range(1, 5))
        by = {int(k).bit_length() - 1: v for k, v in r["by_order"].items()}
        check(
            sorted(by) == list(range(v2 + 1)) and all(v % 2 == 1 for v in by.values()) and by[v2] == factorial(n) // 2**v2,
            f"n={n}: every order 2^0..2^{v2} occurs an odd number of times; {by[v2]} subgroups of order 2^{v2} = n!/2^{v2}",
        )
        ea = {int(k).bit_length() - 1: v for k, v in r["elementary_abelian_by_order"].items()}
        rc = rank_counts(n)
        check([ea.get(j, 0) for j in range(len(rc))] == rc and sum(ea.values()) == sum(rc), f"n={n}: elementary abelian by rank = elemab.rank_counts = {rc}")
        fch = {int(k).bit_length() - 1: v for k, v in r["fc_by_order"].items()}
        check(fch == dist[n], f"n={n}: FC by order = block_families.order_distribution (B(n) = {sum(dist[n].values())})")
        check(r["fc_with_orbit_longer_than_8"] == 0, f"n={n}: no FC subgroup has an orbit longer than 8 (none can, n < 16)")
        check(r["total"] + s3[n] <= A005432[n], f"n={n}: 2-subgroups + |S3(n)| = {r['total'] + s3[n]} <= f(n) = {A005432[n]}")

    print("\n== the census ==")
    print(f"  {'n':>2}  {'2-subgroups':>14}  {'FC':>14}  {'non-FC':>12}  {'non-FC share':>12}  {'2-subgroups/f(n)':>16}")
    for n in range(1, n_max + 1):
        r = rows[n]
        print(f"  {n:2d}  {r['total']:14d}  {r['fc']:14d}  {r['non_fc']:12d}  {r['non_fc'] / r['total']:12.6f}  {r['total'] / A005432[n]:16.6f}")

    out = {
        "description": "Every 2-subgroup of S_n, n <= n_max, by orbit type and Sylow weights (src/two_subgroups.py)",
        "n_max": n_max,
        "sylow_weights": {str(b): {"n_2": sw[b][0], **sw[b][2]} for b in BLOCK_SIZES},
        "types": per_type,
        "fixed_point_free": {str(k): v for k, v in fpf.items()},
        "two_subgroups": {str(n): v for n, v in rows.items()},
        "checks_failed": failures,
        "seconds": round(time.time() - t0, 1),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=1) + "\n")
    shown = args.out.resolve().relative_to(ROOT) if args.out.resolve().is_relative_to(ROOT) else args.out
    print(f"\nwrote {shown}  ({time.time() - t0:.1f}s)")
    print(f"{len(failures)} failed" if failures else "all checks passed")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
