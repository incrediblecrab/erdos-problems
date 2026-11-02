"""The 2-subgroup counts of S_n, n <= 8, against OEIS A218939 (nilpotent subgroups of S_n).

A finite nilpotent group is the direct product of its Sylow subgroups, so a nilpotent
H <= S_n is H = P x Q with P its Sylow 2-subgroup and Q its odd-order part, and P is a
2-subgroup of the centraliser C(Q).  Hence

    nil(n) = t2(n) + sum over odd-order nilpotent Q != 1 of t2(C(Q)),

where t2(G) is the number of 2-subgroups of G and t2(n) = t2(S_n).  For n <= 8 the odd part
Q is a 3-, 5- or 7-group, or the product of a 3-group and a 5-group (3 + 5 = 8 points).  This
script finds every such Q by brute force on permutations, its centraliser by backtracking, and
the 2-subgroups of the centraliser by growing them one normalising element at a time.  t2(n)
is read from data/census8.json (the brute-force census of §9.1); everything else is computed
here and shares no code with any other program.  The result must equal A218939(n), which is
read from refs/oeis_A218939.txt when present and otherwise taken from the copy below.

Exit status 0 only if every n <= --n-max agrees.
"""

from __future__ import annotations

import argparse
import itertools
import json
import pathlib
import re
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parent.parent
A218939 = [1, 1, 2, 5, 24, 102, 837, 5119, 78670, 664658, 13514453, 137227213, 4919721831, 60598902665]  # n = 0..13, fetched November 2, 2026


def compose(p, q):
    return tuple(p[i] for i in q)


def inverse(p):
    out = [0] * len(p)
    for i, x in enumerate(p):
        out[x] = i
    return tuple(out)


def order(p):
    e, k, x = tuple(range(len(p))), 1, p
    while x != e:
        x, k = compose(x, p), k + 1
    return k


def close(gens, n):
    e = tuple(range(n))
    seen, todo = {e}, [e]
    while todo:
        x = todo.pop()
        for g in gens:
            y = compose(x, g)
            if y not in seen:
                seen.add(y)
                todo.append(y)
    return frozenset(seen)


def prime_power_subgroups(elements, n, is_ok):
    """Every subgroup of <elements> all of whose elements pass is_ok (a p-group test), grown from 1:
    a p-group P > K has an element of N_P(K) outside K, so P is reached one normalising element at a time."""
    e = tuple(range(n))
    start = frozenset([e])
    found, todo = {start}, [start]
    while todo:
        K = todo.pop()
        for x in elements:
            if x in K:
                continue
            xi = inverse(x)
            if all(compose(compose(x, k), xi) in K for k in K):
                L = close(list(K) + [x], n)
                if L not in found and is_ok(len(L)):
                    found.add(L)
                    todo.append(L)
    return found


def centraliser(gens, n):
    """All g in S_n with g q = q g for every q in gens, by backtracking on the image of each point."""
    gi = [inverse(q) for q in gens]
    out = []

    def extend(g, used):
        free = [x for x in range(n) if g[x] < 0]
        if not free:
            out.append(tuple(g))
            return
        x = free[0]
        for y in range(n):
            if used[y]:
                continue
            h, u = list(g), list(used)
            h[x], u[y] = y, True
            stack, ok = [x], True
            while stack and ok:
                a = stack.pop()
                for q, qi in zip(gens, gi):
                    for r in (q, qi):
                        b, c = r[a], r[h[a]]  # g(r a) = r(g a)
                        if h[b] < 0:
                            if u[c]:
                                ok = False
                                break
                            h[b], u[c] = c, True
                            stack.append(b)
                        elif h[b] != c:
                            ok = False
                            break
                    if not ok:
                        break
            if ok:
                extend(h, u)

    extend([-1] * n, [False] * n)
    return out


def is_power(m, p):
    while m % p == 0:
        m //= p
    return m == 1


def nilpotent_count(n, t2):
    e = tuple(range(n))
    perms = list(itertools.permutations(range(n)))
    by_prime = {}
    for p in (3, 5, 7):
        els = [x for x in perms if x != e and is_power(order(x), p)]
        by_prime[p] = [Q for Q in prime_power_subgroups(els, n, lambda m, p=p: is_power(m, p)) if len(Q) > 1]
    odd = [Q for p in (3, 5, 7) for Q in by_prime[p]]
    for Q3 in by_prime[3]:  # 3-group x 5-group, commuting elementwise
        for Q5 in by_prime[5]:
            if all(compose(a, b) == compose(b, a) for a in Q3 for b in Q5):
                odd.append(frozenset(compose(a, b) for a in Q3 for b in Q5))
    total = t2
    for Q in odd:
        C = centraliser(sorted(Q), n)
        twos = [x for x in C if x != e and is_power(order(x), 2)]
        total += len(prime_power_subgroups(twos, n, lambda m: is_power(m, 2)))
    return total, len(odd)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--n-max", type=int, default=8)
    ap.add_argument("--census", type=pathlib.Path, default=ROOT / "data" / "census8.json")
    ap.add_argument("--oeis", type=pathlib.Path, default=ROOT / "refs" / "oeis_A218939.txt")
    ap.add_argument("--out", type=pathlib.Path, default=ROOT / "data" / "nilpotent8.json")
    args = ap.parse_args()
    t0 = time.time()
    ref, source = A218939, "the copy in this script"
    if args.oeis.exists():
        terms = "".join(re.findall(r"^%[STU] A218939 (.*)$", args.oeis.read_text(), re.M))
        ref, source = [int(x) for x in terms.replace(" ", "").split(",") if x], args.oeis.name
        if ref[: len(A218939)] != A218939:
            print(f"  [FAIL] {args.oeis.name} differs from the copy in this script")
            return 1
    census = json.loads(args.census.read_text())
    failures, rows = [], {}
    print(f"== nilpotent subgroups of S_n = 2-subgroups + those with an odd part, against A218939 ({source}) ==")
    for n in range(1, args.n_max + 1):
        t2 = census["frattini_closed"][str(n)]["two_subgroups"]
        nil, nq = nilpotent_count(n, t2)
        ok = nil == ref[n]
        rows[str(n)] = {"two_subgroups": t2, "odd_parts": nq, "nilpotent": nil, "A218939": ref[n]}
        print(f"  [{'pass' if ok else 'FAIL'}] n={n}: {t2} 2-subgroups + {nil - t2} with one of {nq} odd parts = {nil} = A218939({n}) = {ref[n]}", flush=True)
        if not ok:
            failures.append(n)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"description": "2-subgroup counts of S_n against A218939 (src/nilpotent8.py)", "rows": rows, "match": not failures}, indent=1) + "\n")
    print(f"({time.time() - t0:.1f}s)")
    print("all checks passed" if not failures else f"{len(failures)} failed")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
