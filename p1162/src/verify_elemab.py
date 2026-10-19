"""Independent brute-force check of `elemab.py`.

Enumerates every elementary abelian 2-subgroup of S_n by explicit permutation
arithmetic and compares the rank-by-rank counts with the generating-function
result.  This file deliberately shares no code with `elemab.py`: it imports it
only to compare final numbers.

Method.  A group of exponent 2 is abelian, so a subgroup H of rank r+1 is
H' u H'g for a rank-r subgroup H' and an involution g centralising H' with
g not in H'.  We therefore grow subgroups one rank at a time, keeping a bitmask
of the involutions that still centralise everything chosen so far, and dedupe
on the frozen set of group elements.
"""

from __future__ import annotations

import argparse
import itertools
import json
import pathlib
import sys
import time


def symmetric_group(n: int) -> list[tuple[int, ...]]:
    return list(itertools.permutations(range(n)))


def compose(p: tuple[int, ...], q: tuple[int, ...]) -> tuple[int, ...]:
    """(p*q)(i) = p(q(i))."""
    return tuple(p[q[i]] for i in range(len(p)))


def enumerate_elementary_abelian(n: int, verbose: bool = False) -> list[int]:
    """[ count of rank-k elementary abelian 2-subgroups of S_n, k = 0.. ]."""
    identity = tuple(range(n))
    perms = symmetric_group(n)
    index = {p: i for i, p in enumerate(perms)}

    involutions = [p for p in perms if p != identity and compose(p, p) == identity]
    m = len(involutions)
    inv_index = {p: i for i, p in enumerate(involutions)}

    # commute[i] = bitmask of involutions commuting with involutions[i]
    commute = [0] * m
    for i in range(m):
        a = involutions[i]
        mask = 0
        for j in range(m):
            b = involutions[j]
            if compose(a, b) == compose(b, a):
                mask |= 1 << j
        commute[i] = mask

    full_mask = (1 << m) - 1
    # level entries: (frozenset of element indices, candidate involution mask)
    level: dict[frozenset[int], int] = {frozenset({index[identity]}): full_mask}
    counts = [1]
    rank = 0
    while level:
        rank += 1
        nxt: dict[frozenset[int], int] = {}
        for elems, mask in level.items():
            elem_perms = [perms[i] for i in elems]
            mm = mask
            while mm:
                low = mm & -mm
                j = low.bit_length() - 1
                mm ^= low
                g = involutions[j]
                gi = index[g]
                if gi in elems:
                    continue
                new_elems = frozenset(elems | {index[compose(h, g)] for h in elem_perms})
                if len(new_elems) != 2 * len(elems):
                    continue
                if new_elems in nxt:
                    continue
                nxt[new_elems] = mask & commute[j]
        if not nxt:
            break
        counts.append(len(nxt))
        if verbose:
            print(f"  rank {rank}: {len(nxt)}", file=sys.stderr, flush=True)
        level = nxt
    return counts


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-max", type=int, default=9)
    ap.add_argument("--out", type=pathlib.Path, default=None)
    args = ap.parse_args()

    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
    from elemab import rank_counts  # noqa: PLC0415

    report = {}
    ok = True
    for n in range(0, args.n_max + 1):
        t0 = time.time()
        brute = enumerate_elementary_abelian(n, verbose=args.n_max >= 10)
        formula = rank_counts(n)
        # pad to equal length: rank_counts always runs to floor(n/2)
        while len(brute) < len(formula):
            brute.append(0)
        while len(formula) < len(brute):
            formula.append(0)
        match = brute == formula
        ok &= match
        report[n] = {"brute": brute, "formula": formula, "match": match}
        print(
            f"n={n:2d}  {'OK ' if match else 'MISMATCH'}  "
            f"{brute}  ({time.time() - t0:.1f}s)",
            flush=True,
        )
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(report, indent=1))
        print(f"wrote {args.out}")
    print("ALL MATCH" if ok else "FAILURE: brute force disagrees with formula")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
