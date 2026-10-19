"""Rule out 21 equal-area triangles on 7 points.

This proves the improved upper bound g(7) <= 20, using the already checked g(5) = 7 and g(6) = 12 inputs.

If 7 points had 21 triangles of one nonzero area, then every 6-subset would have exactly 12 such triangles, because
each triangle is contained in four 6-subsets and g(6) = 12. For every 5-subset, g(5) = 7 gives at most 7; if the
trace has 7, the g(5) classification says it is K_5^(3) with exactly the three triples through one pair missing.

The script first enumerates, up to isomorphism, all 3-uniform hypergraphs on 7 vertices satisfying these necessary
incidence conditions and having 21 edges. Then, for each hypergraph with each possible marked edge affinely normalized
to (0,0), (1,0), (0,1), it checks the polynomial system det(t)^2 = 1 for all marked equal-area triples. A Groebner
basis [1] over QQ rules out even complex realizations, hence real ones.
"""
import itertools
import sys
import time

import sympy as sp

V = tuple(range(7))
TRIPLES = list(itertools.combinations(V, 3))
TRI_INDEX = {t: i for i, t in enumerate(TRIPLES)}
ALL_TRI_MASK = (1 << len(TRIPLES)) - 1


def five_allowed_masks(S):
    local = list(itertools.combinations(S, 3))
    masks = {m for m in range(1 << 10) if m.bit_count() <= 6}
    all_mask = (1 << 10) - 1
    for pair in itertools.combinations(S, 2):
        missing = 0
        for j, t in enumerate(local):
            if pair[0] in t and pair[1] in t:
                missing |= 1 << j
        masks.add(all_mask ^ missing)
    return masks


FIVE_SETS = list(itertools.combinations(V, 5))
FIVE_DATA = []
for S in FIVE_SETS:
    local = list(itertools.combinations(S, 3))
    FIVE_DATA.append(([TRI_INDEX[t] for t in local], five_allowed_masks(S)))


def canonical(mask):
    best = None
    for p in itertools.permutations(V):
        m = 0
        for i, t in enumerate(TRIPLES):
            if mask >> i & 1:
                tt = tuple(sorted(p[x] for x in t))
                m |= 1 << TRI_INDEX[tt]
        if best is None or m < best:
            best = m
    return best


def valid_five_traces(H_mask):
    for inds, allowed in FIVE_DATA:
        local_mask = 0
        for j, ti in enumerate(inds):
            if H_mask >> ti & 1:
                local_mask |= 1 << j
        if local_mask not in allowed:
            return False
    return True


def enumerate_hypergraphs():
    """Return one labelled representative of each isomorphism class that survives the combinatorial reductions."""
    order = list(range(len(TRIPLES)))
    remaining_incidence = []
    for pos in range(len(order) + 1):
        rem = [0] * 7
        for ti in order[pos:]:
            for v in TRIPLES[ti]:
                rem[v] += 1
        remaining_incidence.append(tuple(rem))

    seen = set()
    representatives = []

    # Work in the complement C. A 21-edge H has 14 missing triples. Since every 6-subset of H has 12 edges, every
    # vertex has degree 6 in C.
    def rec(pos, C_mask, deg, chosen):
        if chosen > 14 or chosen + len(order) - pos < 14:
            return
        for v in V:
            if deg[v] > 6 or deg[v] + remaining_incidence[pos][v] < 6:
                return
        if pos == len(order):
            if chosen == 14 and all(d == 6 for d in deg):
                H_mask = ALL_TRI_MASK ^ C_mask
                if valid_five_traces(H_mask):
                    c = canonical(H_mask)
                    if c not in seen:
                        seen.add(c)
                        representatives.append(H_mask)
            return
        ti = order[pos]
        new_deg = list(deg)
        for v in TRIPLES[ti]:
            new_deg[v] += 1
        rec(pos + 1, C_mask | (1 << ti), tuple(new_deg), chosen + 1)
        rec(pos + 1, C_mask, deg, chosen)

    rec(0, 0, (0,) * 7, 0)
    return representatives


XS = sp.symbols("x3 y3 x4 y4 x5 y5 x6 y6")


def point(i):
    if i == 0:
        return sp.Integer(0), sp.Integer(0)
    if i == 1:
        return sp.Integer(1), sp.Integer(0)
    if i == 2:
        return sp.Integer(0), sp.Integer(1)
    return XS[2 * (i - 3)], XS[2 * (i - 3) + 1]


def det_poly(t):
    a, b, c = [point(i) for i in t]
    return sp.expand((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))


DET = {t: det_poly(t) for t in TRIPLES}


def relabel_to_base(H_edges, base):
    rest = [v for v in V if v not in base]
    mp = {base[0]: 0, base[1]: 1, base[2]: 2}
    for i, v in enumerate(rest, 3):
        mp[v] = i
    return frozenset(tuple(sorted(mp[v] for v in t)) for t in H_edges)


def algebraically_impossible(H_edges):
    marked = set()
    for base in H_edges:
        H = relabel_to_base(H_edges, base)
        if H in marked:
            continue
        marked.add(H)
        equations = [DET[t] ** 2 - 1 for t in H if t != (0, 1, 2)]
        G = sp.groebner(equations, *XS, order="grevlex", domain=sp.QQ)
        if list(G.exprs) != [1]:
            return False, base, len(G.exprs)
    return True, None, len(marked)


def main():
    start = time.time()
    reps = enumerate_hypergraphs()
    print(f"combinatorial 21-edge hypergraphs: {len(reps)} isomorphism classes")
    all_impossible = True
    for i, mask in enumerate(reps, 1):
        H_edges = [TRIPLES[j] for j in range(len(TRIPLES)) if mask >> j & 1]
        impossible, witness, detail = algebraically_impossible(H_edges)
        print(f"class {i}: {'no complex realization' if impossible else 'unresolved'} ({detail} marked cases)")
        if not impossible:
            print(f"first unresolved marked edge: {witness}")
            all_impossible = False
    print(f"elapsed {time.time() - start:.2f} s")
    if all_impossible and len(reps) == 2:
        print("g(7) <= 20")
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
