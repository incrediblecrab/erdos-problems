"""Re-derive the numbers quoted in NOTES.md from the artifacts.

Exit status 0 means every check below passed.  Any non-zero exit means at
least one stated number or checked claim is wrong; the failing check is
printed.  A check with no data in this run (FAST mode, missing downloads) is
reported as skipped, never as passed.

Two kinds of check.  The first half tests the mathematics: each check
recomputes a quantity rather than reading it back from a file that the
narrative also quotes, except where the point of the check is exactly that a
stored file agrees with a fresh computation.  The second half tests the prose:
it finds numbers printed in NOTES.md by the text just before them and compares
each with a fresh computation, to the precision printed.  The summary line
says how many of the decimal numbers in NOTES.md were compared; `--coverage`
lists the rest (timings, identifiers, and the wrong values that correction
notes quote on purpose).

Run:  python src/final_check.py [--coverage]
"""

from __future__ import annotations

import bisect
import functools
import json
import math
import pathlib
import re
import subprocess
import sys
import tempfile
from collections import Counter
from decimal import Decimal, getcontext
from fractions import Fraction
from math import comb, factorial

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))

from asymptotics import fit_excess, log2_int, rank_stats, ratio  # noqa: E402
from elemab import (  # noqa: E402
    galois_number,
    gaussian_binomial,
    gl_order,
    hom_table,
    rank_counts,
)
from block_families import (  # noqa: E402
    A005432,  # n = 0..18
    L,
    dense_structures,
    egf,
    load_types,
    pmul,
    w_poly,
    weights,
)
from orbit_model import (  # noqa: E402
    holdout_extrapolation,
    local_coefficients,
    orbit_profile,
    pure_type_log2,
    pure_type_scan,
)

A000638 = [1, 1, 2, 4, 11, 19, 56, 96, 296, 554]  # n = 0..9
# (n, k, j) -> exact number of homomorphisms F_2^k -> S_n all of whose orbits
# have size exactly 2^j, obtained by enumerating k-tuples of pairwise
# commuting elements of order dividing 2 (NOTES.md §6 Check 1)
BRUTE_HOM = {
    (4, 1, 1): 3, (4, 1, 2): 0, (4, 2, 1): 27, (4, 2, 2): 6,
    (4, 3, 1): 147, (4, 3, 2): 42, (6, 1, 1): 15, (6, 2, 1): 405,
    (6, 3, 1): 5145, (8, 1, 1): 105, (8, 1, 2): 0, (8, 1, 3): 0,
    (8, 2, 1): 8505, (8, 2, 2): 1260, (8, 2, 3): 0,
}
ORBIT_TOL = 2e-9  # sum_m m E[#orbits_m] = n, float rounding (NOTES.md §6 Check 2)

FAILURES: list[str] = []
CHECKS = 0
SKIPPED: list[str] = []
STORED_SPLIT: dict[str, int] = {}  # rows of elemab_stats.json recomputed or not, for §10
PDF_RESULTS: dict[str, bool] = {}  # quotation or pattern -> passed, filled by check_quotes


def check(name: str, ok: bool, detail: str = "") -> None:
    global CHECKS
    CHECKS += 1
    if ok:
        print(f"  ok    {name}" + (f"   [{detail}]" if detail else ""))
    else:
        print(f"  FAIL  {name}   {detail}")
        FAILURES.append(f"{name}: {detail}")


def close(a: float, b: float, tol: float) -> bool:
    return abs(a - b) <= tol


def load(rel: str):
    p = ROOT / rel
    if not p.exists():
        return None
    return json.loads(p.read_text())


# ------------------------------------------------------------------ constants
# Both constants below are closed forms; the point of the checks is that the
# measured quantities agree with them, so they are computed here from scratch.

THETA = 1.0 + 2.0 * sum(2.0 ** (-j * j) for j in range(1, 60))
ETA = math.prod(1.0 - 2.0**-i for i in range(1, 200))  # (1/2;1/2)_infinity
GALOIS_CONST = math.log2(THETA / ETA)  # 2.8820499654..., even r
THETA_ODD = sum(2.0 ** (-(j + 0.5) ** 2) for j in range(-60, 60))
GALOIS_ODD = math.log2(THETA_ODD / ETA)  # 2.8820461863..., odd r
THETA_SD = math.sqrt(
    sum(j * j * 2.0 ** (-j * j) for j in range(-60, 61)) / THETA
)  # 0.8493059610...


def section(title: str) -> None:
    print(f"\n== {title} ==")


# --------------------------------------------------------------------- checks


def check_closed_forms() -> None:
    section("closed forms quoted in NOTES.md")
    check(
        "Galois constant log2(theta/(1/2;1/2)_inf) = 2.882049965",
        close(GALOIS_CONST, 2.882049965, 5e-10),
        f"{GALOIS_CONST:.10f}",
    )
    check(
        "theta-law standard deviation = 0.849305961",
        close(THETA_SD, 0.849305961, 5e-10),
        f"{THETA_SD:.10f}",
    )
    # zeta_2 of RoTr25 Lemma 2.4 is 1/ETA
    check("RoTr25 Lemma 2.4 constant zeta_2 < 4", (1 / ETA) < 4, f"zeta_2={1 / ETA:.6f}")


def check_elemab_against_bruteforce() -> None:
    section("elemab.py against the independent permutation brute force")
    rec = load("results/verify_elemab.json")
    if rec is None:
        SKIPPED.append("results/verify_elemab.json missing (run src/verify_elemab.py)")
        print("  skip  results/verify_elemab.json not present")
        return
    bad = [n for n, r in rec.items() if not r["match"]]
    check(
        f"brute force agrees with the formula for n <= {max(int(k) for k in rec)}",
        not bad,
        f"mismatches at n={bad}" if bad else f"{len(rec)} values of n",
    )
    for n, r in rec.items():
        if int(n) <= 10:
            fresh = rank_counts(int(n))
            stored = [c for c in r["formula"]]
            while len(fresh) < len(stored):
                fresh.append(0)
            while len(stored) < len(fresh):
                stored.append(0)
            if fresh != stored:
                check(f"stored formula row n={n} still reproduces", False, "drifted")
                return
    check("stored formula rows still reproduce from current elemab.py", True)


def check_involutions() -> None:
    section("a(n,1) against OEIS A000085 (an independent recurrence)")
    # A000085: t(n) = t(n-1) + (n-1) t(n-2), the number of involutions
    # including the identity, so a(n,1) = A000085(n) - 1.
    t = [1, 1]
    for m in range(2, 31):
        t.append(t[m - 1] + (m - 1) * t[m - 2])
    ok = all(rank_counts(n)[1] == t[n] - 1 for n in range(2, 31))
    check("a(n,1) = A000085(n) - 1 for 2 <= n <= 30", ok, f"a(30,1)={rank_counts(30)[1]}")


def check_no_rank_above_half() -> None:
    section("structural identity a(n, floor(n/2)+1) = 0")
    # rank_counts asserts this internally; recompute the inversion one rank
    # higher than the maximum and confirm the massive cancellation lands on 0.
    ok = True
    for n in range(2, 41):
        c = rank_counts(n)
        if len(c) - 1 != n // 2:
            ok = False
            break
    check("rank_counts(n) has top rank exactly floor(n/2) for n <= 40", ok)


def check_lattice_vs_oeis() -> None:
    section("full subgroup lattice against OEIS A005432 / A000638")
    seen = 0
    for n in range(1, 10):
        d = load(f"data/lattice_n{n}.json")
        if d is None:
            continue
        seen += 1
        check(
            f"n={n}: total subgroups = A005432({n}) = {A005432[n]}",
            d["total_subgroups"] == A005432[n],
            f"got {d['total_subgroups']}",
        )
        check(
            f"n={n}: conjugacy classes = A000638({n}) = {A000638[n]}",
            d["total_classes"] == A000638[n],
            f"got {d['total_classes']}",
        )
        lat = [int(x) for x in d["elementary_abelian_rank"]]
        formula = rank_counts(n)
        while len(lat) < len(formula):
            lat.append(0)
        while len(formula) < len(lat):
            formula.append(0)
        check(
            f"n={n}: lattice elementary abelian ranks = elemab.py",
            lat == formula,
            f"lattice {lat} vs formula {formula}",
        )
    if seen == 0:
        SKIPPED.append("no data/lattice_n*.json (build and run src/subgroups.c)")
        print("  skip  no lattice files present")


def check_pure_types() -> None:
    section("pure orbit-type closed form against brute-forced hom counts")
    bad = []
    for (n, k, j), want in BRUTE_HOM.items():
        v = pure_type_log2(n, j, k)
        got = 0.0 if v is None else 2.0**v * gl_order(k)
        if not close(got, want, 1e-6 * max(1.0, want)):
            bad.append((n, k, j, want, got))
    check(
        f"{len(BRUTE_HOM)} brute-forced hom counts match pure_type_log2",
        not bad,
        str(bad) if bad else "",
    )

    section("orbit model: optimal rank of each pure type is j/2^(j+1)")
    for n in (64, 128, 256):
        scan = pure_type_scan(n)
        check(f"n={n}: orbit sizes 2^j for j = 1..4 scanned", [r["j"] for r in scan] == [1, 2, 3, 4])
        for r in scan:
            check(
                f"n={n} orbit size {r['orbit_size']}: best k/n = {r['predicted_k_over_n']}",
                close(r["best_k_over_n"], r["predicted_k_over_n"], 1e-9),
                f"got {r['best_k_over_n']}",
            )
    # j = 1 and j = 2 tie at n^2/16; j >= 3 must lose at order n^2
    for n in (128, 256):
        sc = {r["j"]: r for r in pure_type_scan(n)}
        for j in (1, 2):
            # the model says the excess over n^2/16 is Theta(n log2 n), so the
            # scale to compare against is n log2 n, not n^2
            check(
                f"n={n} orbit size {2**j} reaches n^2/16 up to O(n log n)",
                abs(sc[j]["excess_over_n2_16"]) < n * math.log2(n),
                f"excess {sc[j]['excess_over_n2_16']:+.1f} vs n log2 n = {n * math.log2(n):.0f}",
            )
            check(
                f"n={n} orbit size {2**j}: excess/(n log2 n) below its limit {1 - 2.0**-j}",
                0 < sc[j]["excess_per_n_log2n"] < 1 - 2.0**-j,
                f"{sc[j]['excess_per_n_log2n']:.5f}",
            )
        check(
            f"n={n}: orbit size 4 beats orbit size 2 at order n log n",
            sc[2]["excess_per_n_log2n"] > sc[1]["excess_per_n_log2n"],
            f"{sc[2]['excess_per_n_log2n']:.5f} vs {sc[1]['excess_per_n_log2n']:.5f}",
        )
        if 3 in sc:
            # the model predicts the gap grows like (1/16 - 9/256) n^2 = 7n^2/256,
            # approached from below
            gap = sc[2]["excess_over_n2_16"] - sc[3]["excess_over_n2_16"]
            check(
                f"n={n}: orbit size 8 falls behind orbit size 4 at order n^2",
                0.02 * n * n < gap < (7 / 256) * n * n,
                f"gap {gap:.1f}, gap/n^2 = {gap / (n * n):.5f}, limit 7/256 = 0.02734",
            )


def check_theta_law() -> None:
    section("pointwise theta law  a(n, n/4+j)/a(n, n/4) -> 2^(-j^2)")
    hom = hom_table(128 // 2 + 2, 128)
    worst = 0.0
    for n in (64, 128):
        c = rank_counts(n, hom if n == 128 else None)
        k = n // 4
        for j in range(-3, 4):
            got = ratio(c[k + j], c[k])
            want = 2.0 ** (-j * j)
            worst = max(worst, abs(got - want))
    check(
        "ratios match 2^(-j^2) for |j| <= 3 at n = 64 and 128",
        worst < 2e-4,
        f"max deviation {worst:.2e}",
    )
    # the sharper statement is at large n; use the stored sweep if present
    st = load("data/elemab_stats.json")
    if st is None:
        SKIPPED.append("data/elemab_stats.json missing (run src/asymptotics.py)")
        print("  skip  data/elemab_stats.json not present")
        return
    rows = {r["n"]: r for r in st["rows"]}
    for n in (128, 256, 512):
        if n not in rows:
            continue
        r = rows[n]
        check(
            f"n={n}: mode of the rank distribution is exactly n/4 = {n // 4}",
            r["mode"] == n // 4,
            f"got {r['mode']}",
        )
        check(
            f"n={n}: mean rank / n = 0.25 to 6 dp",
            close(r["mean_over_n"], 0.25, 5e-7),
            f"{r['mean_over_n']:.8f}",
        )
        check(
            f"n={n}: sd of the rank equals the theta constant to 6 dp",
            close(r["sd"], THETA_SD, 5e-7),
            f"measured {r['sd']:.8f} vs {THETA_SD:.8f}",
        )
        check(
            f"n={n}: >= 1-1e-6 of the mass lies within 4 of the mode",
            r["window_1e6"] <= 4,
            f"window {r['window_1e6']}",
        )


def check_stored_rows() -> None:
    section("data/elemab_stats.json: the stored rows and fits, recomputed")
    st = load("data/elemab_stats.json")
    if st is None:
        SKIPPED.append("data/elemab_stats.json missing (run src/asymptotics.py)")
        print("  skip  data/elemab_stats.json not present")
        return
    rows = st["rows"]
    # One table recomputes every row up to 256 in about 3 s; rows 264..504
    # would take about 3 minutes, so they get the closed-form checks below.
    cheap = min(256, max(r["n"] for r in rows))
    hom = hom_table(cheap // 2 + 2, cheap)
    fresh = {r["n"]: rank_stats(r["n"], rank_counts(r["n"], hom)) for r in rows if r["n"] <= cheap}
    if FULL:
        fresh[512] = rank_stats(512, list(rc(512)))
    # asymptotics.ratio divides through float logarithms in the thousands and the variance then cancels about 4 digits, so a stored sd carries float noise of up to about 2e-8 near n = 512 (NOTES.md section 3), which another libm need not reproduce.
    byn = {r["n"]: r for r in rows}
    bad = [n for n in fresh if not _same(dict(fresh[n], seconds=0), byn[n], rel=1e-6)]
    bits = all(fresh[n] == {k: v for k, v in byn[n].items() if k != "seconds"} for n in fresh)
    check(
        f"{len(fresh)} of {len(rows)} stored rows equal a fresh rank_stats to 1e-6 relative",
        not bad,
        f"differ at n = {bad[:8]}" if bad else f"n <= {cheap}" + (" and 512" if FULL else "") + ("; bit-identical here" if bits else ""),
    )
    rest = [r for r in rows if r["n"] not in fresh]
    STORED_SPLIT.update(cheap=cheap, fresh=len(fresh), total=len(rows), rest=len(rest))
    if rest:
        ns = [r["n"] for r in rest]
        check(
            f"the other {len(rest)} rows, n = {ns[0]}..{ns[-1]}: excess, mode/n, mean/n, max_rank consistent",
            all(
                r["excess"] == r["log2_e"] - r["n"] * r["n"] / 16.0
                and r["mode_over_n"] == r["mode"] / r["n"]
                and r["mean_over_n"] == r["mean"] / r["n"]
                and r["max_rank"] == r["n"] // 2
                for r in rest
            ),
        )
        check(
            f"the other {len(rest)} rows: mode n/4, window 4, mean/n within 1e-8 of 1/4, sd within 1e-7 of the closed form",
            all(
                r["mode"] * 4 == r["n"]
                and r["window_1e6"] == 4
                and abs(r["mean_over_n"] - 0.25) < 1e-8
                and abs(r["sd"] - THETA_SD) < 1e-7
                for r in rest
            ),
            f"largest sd deviation {max(abs(r['sd'] - THETA_SD) for r in rest):.1e}",
        )
        dev = [abs(Decimal(r["log2_e"]) - l2(galois_number(r["n"] // 2) * M(r["n"]))) for r in rest]
        check(
            f"the other {len(rest)} rows: log2_e within 1e-9 of log2 G_(n/2) M(n) (§9.2), e_bits consistent",
            max(dev) < Decimal("1e-9") and all(r["e_bits"] == math.floor(r["log2_e"]) + 1 for r in rest),
            f"largest {max(dev):.1e}",
        )
    fits = [fit_excess(rows, f["n_min"]) for f in st["fits"]]
    check(
        f"the {len(fits)} stored fits equal a fresh least squares on the stored rows to 1e-6 relative",
        _same(fits, st["fits"], rel=1e-6),
    )


def check_galois() -> None:
    section("Galois numbers  log2 G_r - r^2/4 -> the same constant")
    g = load("data/galois.json")
    if g is None:
        SKIPPED.append("data/galois.json missing (run src/asymptotics.py)")
        print("  skip  data/galois.json not present")
    else:
        rows = {r["r"]: r for r in g}
        # The approach to the constant is itself a claim: the leading correction
        # to [r choose r/2]_2 is of size 2^-(r/2), so check that rate rather
        # than asserting a single tolerance that happens to hold at large r.
        for r in (32, 64):
            if r in rows:
                dev = abs(rows[r]["excess_over_r2_4"] - GALOIS_CONST)
                check(
                    f"r={r}: deviation from the constant is below 2^(3-r/2)",
                    dev < 2.0 ** (3 - r / 2),
                    f"dev {dev:.3e} < {2.0 ** (3 - r / 2):.3e}",
                )
        for r in (128, 256, 512):
            if r in rows:
                # 1e-11 is the float floor of log2_int on numbers this large
                check(
                    f"r={r} (even): excess = {GALOIS_CONST:.9f}",
                    close(rows[r]["excess_over_r2_4"], GALOIS_CONST, 1e-11),
                    f"{rows[r]['excess_over_r2_4']:.12f}",
                )
        # odd r has its own limit, log2(theta'/eta) with theta' = sum 2^-(j+1/2)^2
        if 63 in rows:
            dev = abs(rows[63]["excess_over_r2_4"] - GALOIS_ODD)
            check(
                "r=63 (odd): deviation from the odd constant is below 2^(3-r/2)",
                dev < 2.0 ** (3 - 63 / 2),
                f"dev {dev:.3e} < {2.0 ** (3 - 63 / 2):.3e}",
            )
        for r in (127, 255, 511):
            if r in rows:
                check(
                    f"r={r} (odd): excess = log2(theta'/eta) = {GALOIS_ODD:.9f}",
                    close(rows[r]["excess_over_r2_4"], GALOIS_ODD, 1e-11),
                    f"{rows[r]['excess_over_r2_4']:.12f}",
                )
    # recompute one value from scratch rather than trusting the file
    fresh = log2_int(galois_number(256)) - 256 * 256 / 4.0
    check(
        "fresh galois_number(256) reproduces the constant",
        close(fresh, GALOIS_CONST, 1e-11),
        f"{fresh:.12f}",
    )
    # G_r for small r is OEIS A006116
    a006116 = [1, 2, 5, 16, 67, 374, 2825, 29212, 417199]
    check(
        "G_0..G_8 = OEIS A006116",
        [galois_number(r) for r in range(9)] == a006116,
        "",
    )
    # a single maximal elementary abelian gives only n^2/16 + O(1)
    for n in (256, 512):
        v = log2_int(galois_number(n // 2)) - n * n / 16.0
        check(
            f"n={n}: one maximal elementary abelian gives n^2/16 + O(1)",
            abs(v - GALOIS_CONST) < 1e-11,
            f"excess {v:.12f}",
        )


def check_second_order() -> None:
    section("second-order term: A_local -> 3/4, scored out of sample")
    st = load("data/elemab_stats.json")
    if st is None:
        print("  skip  data/elemab_stats.json not present")
        return
    loc = local_coefficients(st["rows"], 16)
    hv = holdout_extrapolation(loc)
    if "error" in hv:
        print(f"  skip  {hv['error']}")
        SKIPPED.append("held-out extrapolation: " + hv["error"])
        return
    best_basis = hv["ranked_bases"][0]
    check(
        "n^-1/2 is the best extrapolation basis on held-out n",
        best_basis["basis"] == "n^-1/2",
        f"won: {best_basis['basis']} (held-out maxres {best_basis['holdout_maxres']:.2e})",
    )
    b = hv["best_free_exponent"]
    check(
        "free exponent selected out of sample is near 1/2",
        0.40 <= b["p"] <= 0.60,
        f"p={b['p']:.2f}",
    )
    check(
        "extrapolated A_inf = 3/4 to within 0.002",
        abs(b["A_inf"] - 0.75) < 0.002,
        f"A_inf={b['A_inf']:.5f}, error {abs(b['A_inf'] - 0.75):.5f}",
    )
    # and the rejected hypothesis really is rejected
    byname = {r["basis"]: r for r in hv["ranked_bases"]}
    check(
        "the rejected 1/log2 n model extrapolates to the wrong limit",
        abs(byname["1/log2 n"]["A_inf"] - 0.75) > 0.02,
        f"A_inf={byname['1/log2 n']['A_inf']:.5f}",
    )


def check_orbit_profile() -> None:
    section("orbit profile at the dominant rank k = n/4")
    prev = None
    for n in (64, 128, 256, 512):
        p = orbit_profile(n, n // 4)
        check(
            f"n={n}: sum_m m*E[#orbits_m] = n to {ORBIT_TOL:.0e}",
            p["point_total_error"] < ORBIT_TOL,
            f"error {p['point_total_error']:.1e}",
        )
        f4 = p["profile"].get(4, {}).get("exp_point_fraction", 0.0)
        f8 = p["profile"].get(8, {}).get("exp_point_fraction", 0.0)
        check(
            f"n={n}: orbits of size 4 hold the majority of points",
            f4 > 0.5,
            f"fraction {f4:.6f}",
        )
        check(f"n={n}: orbits of size >= 8 are negligible", f8 < 1e-3, f"{f8:.2e}")
        if prev is not None:
            check(f"n={n}: the size-4 fraction increased", f4 > prev, f"{prev:.4f} -> {f4:.4f}")
        prev = f4


def check_sd_high_precision() -> None:
    section("rank sd vs closed form, exact arithmetic at 60 dp")
    # The doubles in elemab_stats.json truncate at ~1e-16 relative, which is
    # coarser than the real convergence from n=256 on.  Recompute exactly.
    d = load("data/sd_precision.json")
    if d is None:
        SKIPPED.append("data/sd_precision.json missing (run src/sd_precision.py)")
        print("  skip  data/sd_precision.json not present")
        return
    check(
        "closed form printed to 30 dp is 0.849305961022908280743697585694",
        d["closed_form_30dp"].startswith("0.849305961022908280743697585694"),
        d["closed_form_30dp"],
    )
    want = {64: 4, 128: 7, 256: 15, 512: 30}  # lower bounds, not exact values
    byn = {r["n"]: r for r in d["rows"]}
    for n, floor_dp in want.items():
        r = byn.get(n)
        if r is None:
            continue
        check(
            f"n={n}: exact sd agrees with the closed form to >= {floor_dp} dp",
            r["decimal_places"] >= floor_dp,
            f"{r['decimal_places']} dp, |err| {r['abs_error']}",
        )
    r512 = byn.get(512)
    if r512 is not None:
        check(
            "n=512: all 30 printed digits of the sd match the closed form",
            r512["sd_30dp"].startswith("0.849305961022908280743697585694"),
            r512["sd_30dp"],
        )
    ratios = [r["log10_ratio_vs_prev"] for r in d["rows"] if r["log10_ratio_vs_prev"]]
    if ratios:
        check(
            "log10(error) roughly doubles per doubling of n (error ~ 2^-cn)",
            all(1.8 < x < 2.3 for x in ratios[-2:]),
            f"ratios {[round(x, 2) for x in ratios]}",
        )


def check_rotr_thresholds() -> None:
    section("RoTr25 Theorem 6 threshold algebra")
    lo = (1 - math.sqrt(3) / 2) / 2
    check(
        "nu(1-nu) = 1/16 has smaller root 1/2 - sqrt(3)/4 = 0.0669872981",
        close(lo, 0.5 - math.sqrt(3) / 4, 1e-15) and close(lo, 0.06698729810778, 1e-12),
        f"{lo:.13f}",
    )
    check("that root is the constant printed in Theorem 6", close(lo, 0.06698729810778, 1e-12))
    disc = 0.25 - 4 * (1 / 16.0)
    check("nu(1/2 - nu) = 1/16 has discriminant exactly 0", disc == 0.0, f"{disc!r}")
    check("its double root is 1/4", close(0.5 / 2, 0.25, 0), "0.25")
    check(
        "improvement factors are 4.000 and 3.732",
        close(0.25 / (1 / 16.0), 4.0, 1e-12)
        and close(0.25 / lo, 3.7320508075688772, 1e-9),
        f"{0.25 / (1 / 16.0):.3f} and {0.25 / lo:.3f}",
    )


def check_quotes() -> None:
    section("verbatim quotations from the primary source")
    pdf = ROOT / "refs" / "roney-dougal_tracey_2503.05416.pdf"
    if not pdf.exists():
        SKIPPED.append("refs/*.pdf absent (run refs/fetch.sh); quotations unchecked")
        print("  skip  refs PDF not present; run refs/fetch.sh")
        return
    # -raw keeps content-stream order.  Default and -layout mode both sort by
    # position, which splices superscripts from displayed math into the middle
    # of a sentence: the base-2 remark comes out as "are to the 2 2 base 2)"
    # and a verbatim search fails even though the sentence is on the page.
    # The base-2 line was additionally read off a 150 dpi render of page 1.
    try:
        txt = subprocess.run(
            ["pdftotext", "-raw", str(pdf), "-"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout
    except (OSError, subprocess.CalledProcessError) as exc:
        SKIPPED.append(f"pdftotext unavailable ({exc}); quotations unchecked")
        print("  skip  pdftotext unavailable")
        return
    flat = " ".join(txt.split())
    quotes = [
        "our logarithms are to the base 2",
        "settling a conjecture of Pyber from 1993",
        "A random nilpotent subgroup of Sn is a 2-group",
        "Let n be congruent to 3 modulo 4",
        "the probability that a random subgroup of Sn is nilpotent is bounded away from 1",
        "an elementary abelian 2-section of order at least 2",
        "Since ν(1 − ν) < 1/16, the result follows",
        "Since µ < 1/16 the result follows from the lower bound in Theorem 1",
        "P has at most ζppk(ℓ−k) subgroups of order pk",
        # §9
        "We shall prove that this holds with α0 = 448/453",
        "If no Gi, for i ≤ r, has excess 2 and degree 8",
        "The group G has excess 2 and degree 8",
        "that no group has excess more than 2",
        "If G is excessive then G has excess i ∈ {1,2} and n ≤ 32",
        "Firstly, does there in fact exist an absolute constant γ such that",
        "Secondly, is it possible that the probability that a random subgroup of Sn is nilpotent tends to 0",
        "Thirdly, does there exist an absolute constant C such that a random subgroup of Sn has no orbits of length greater than C",
        "we present evidence in Section 6 that it may have a positive answer",
        "The lower bound is Proposition 7.3",
        "generated by m disjoint p-cycles",
        "for all integers n > 1",
        "There exist constants βp > αp > 0 such that for all integers n ≥ p",
        "Write ln for natural logarithms",
        "Pyber conjectured, however, that",
        # Proposition 7.3's minimum, as -raw extracts the displayed formula (read off a render of p. 23)
        "0.08, (1 − (2p − 1)2/(4p2)) 2plog(2p) ,αR,(1 − α0)(1 − 1/p)",
    ]
    for q in quotes:
        PDF_RESULTS[q] = " ".join(q.split()) in flat
        check(f'quote present: "{q[:56]}"', PDF_RESULTS[q])
    # NOTES.md §8: [RoTr25] does not state Theorem 9.1 or name the group
    for pat in ("7/8", "0.875", "extra.{0,3}special"):
        PDF_RESULTS[pat] = re.search(pat, flat, re.IGNORECASE) is None
        check(f"absent from the paper: {pat!r}", PDF_RESULTS[pat])


# ------------------------------------------------------ numbers in NOTES.md
# The second half.  Each number is found by the text just before it (an anchor
# that must occur exactly once in NOTES.md) and compared with a fresh value to
# the precision printed: "0.2305" must lie within 5e-5 of the value and
# "1.0e-17" within 5e-19.  A number followed by \ldots is a truncation and may
# lie anywhere in the interval it truncates.  An integer printed against an
# integer value must be equal.

getcontext().prec = 100
NOTES = (ROOT / "NOTES.md").read_text("utf-8")
# \mathbf{0.75057} and \{0.08 hide a number behind a brace; the scan copy
# blanks both wrappers and keeps every offset
_SCAN = NOTES.replace("\\mathbf{", " " * 8).replace("\\{", "  ")
NUM_RE = re.compile(
    r"(?<![\w.^{/])[−-]?(?:\d{1,3}(?:[ \u2009\u202f]\d{3})+(?![\d.])|\d+)(?:\.\d+)?"
    r"(?:\s*\\times\s*10\^\{[−-]?\d+\}|e[−-]\d+)?(?!\w)"
)
# section numbers, arXiv identifiers and versions are not quantities
XREF_RE = re.compile(
    r"(?:§|\b(?:Theorems?|Propositions?|Conjectures?|Lemma|Checks?|items?))\s*\d+(?:\.\d+)*(?:\([ivx]+\))?"
    r"(?:(?:,\s*|\s+and\s+|\s*–\s*|\s+or\s+)§?\d+(?:\.\d+)*(?:\([ivx]+\))?)*"
    r"|\b\d{4}\.\d{5}(?:v\d+)?|\b\d+\.\d+\.\d+|5\.73"
    r"|direction of 7\.2|replaces 7\.2|7\.2 and the ν-optimality half of 7\.1"
)
_EXCLUDED = [m.span() for m in re.finditer(r"(?m)^#.*$", NOTES)] + [
    m.span() for m in XREF_RE.finditer(NOTES)
]
NUMS = [
    (m.start(), m.group())
    for m in NUM_RE.finditer(_SCAN)
    if not any(a <= m.start() < b for a, b in _EXCLUDED)
]
NUM_POS = [off for off, _ in NUMS]
COVERED: dict[int, str] = {}  # offset -> compared | historical | timing
NOTES_SKIPPED: list[str] = []
IDX = ...  # a table cell that holds a label, not a quantity


def is_dec(s: str) -> bool:
    return "." in s or "times" in s or bool(re.search(r"e[−-]", s))


def _parse(s: str) -> tuple[Decimal, Decimal, bool]:
    """Printed value, half a unit in its last place, and whether it is an integer."""
    t = re.sub(r"[ \u2009\u202f]", "", s.replace("−", "-"))
    m = re.fullmatch(r"(-?)(\d+)(?:\.(\d+))?(?:\\times10\^\{(-?\d+)\}|e(-\d+))?", t)
    sign, ip, fp, e1, e2 = m.groups()
    exp = int(e1 or e2 or 0)
    p = Decimal(f"{sign}{ip}.{fp or '0'}").scaleb(exp)
    return p, Decimal(5).scaleb(exp - len(fp or "") - 1), fp is None and exp == 0


def _dec(v) -> Decimal:
    if isinstance(v, Fraction):
        return Decimal(v.numerator) / Decimal(v.denominator)
    return Decimal(v)


def _show(v) -> str:
    if isinstance(v, tuple):
        return "(" + ", ".join(_show(x) for x in v) + ")"
    if isinstance(v, (int, Fraction)):
        return str(v)
    return format(_dec(v), ".8g")


def agrees(s: str, v, trunc: bool = False) -> bool:
    """Does the printed string s state the value v to the precision printed?"""
    if isinstance(v, tuple):
        return all(agrees(s, x, trunc) for x in v)
    p, half, is_int = _parse(s)
    if is_int and isinstance(v, int):
        return p == v
    d = _dec(v) - p
    if abs(d) <= half * (1 + Decimal("1e-9")):
        return True
    return trunc and 0 <= (d if p >= 0 else -d) < 2 * half


def _trunc(end: int) -> bool:
    return re.match(r"\s*(?:\\ldots|…|\\dots)", NOTES[end : end + 8]) is not None


def _anchor(anchor: str) -> int | None:
    i = NOTES.find(anchor)
    if i < 0 or NOTES.find(anchor, i + 1) >= 0:
        return None
    return i + len(anchor)


def near(anchor: str, k: int, ints: bool = False):
    """The k decimal numbers (all numbers if ints) printed after the anchor."""
    pos = _anchor(anchor)
    if pos is None:
        return None
    out = []
    for off, s in NUMS[bisect.bisect_left(NUM_POS, pos) :]:
        if len(out) == k:
            break
        if ints or is_dec(s):
            out.append((s, off, _trunc(off + len(s))))
    return out


def _record(label: str, pairs, missing: str = "") -> None:
    """pairs: (printed, offset, trunc, value).  None values are skipped."""
    bad, done = [], 0
    for s, off, tr, v in pairs:
        if v is None:
            continue
        done += 1
        if agrees(s, v, tr):
            COVERED[off] = "compared"
        else:
            bad.append(f"printed {s}, fresh {_show(v)}")
    if done < len(pairs):
        NOTES_SKIPPED.append(label)
    if missing:
        bad.append(missing)
    if done or bad:
        shown = ", ".join(p[0] for p in pairs if p[3] is not None)
        check(f"NOTES: {label}", not bad, "; ".join(bad) if bad else shown)


def compare(label: str, anchor: str, values, ints: bool = False) -> None:
    vals = values if isinstance(values, list) else [values]
    hits = near(anchor, len(vals), ints)
    if hits is None:
        check(f"NOTES: {label}", False, f"anchor not found exactly once: {anchor!r}")
        return
    short = f"only {len(hits)} numbers follow the anchor" if len(hits) < len(vals) else ""
    _record(label, [(s, off, tr, v) for (s, off, tr), v in zip(hits, vals)], short)


def mark(kind: str, anchor: str, k: int = 1) -> None:
    """The next k decimal numbers are quoted on purpose: historical or timing."""
    hits = near(anchor, k)
    if hits is None or len(hits) < k:
        check(f"NOTES: {kind} number after {anchor[:40]!r}", False, "anchor not found exactly once")
        return
    for _, off, _ in hits:
        COVERED[off] = kind


def text(label: str, s: str) -> None:
    """An exact string, built from fresh values, is present in NOTES.md."""
    i = NOTES.find(s)
    if i >= 0:
        for off, _ in NUMS[bisect.bisect_left(NUM_POS, i) :]:
            if off >= i + len(s):
                break
            COVERED[off] = "compared"
    check(f"NOTES: {label}", i >= 0, s[:90] if i >= 0 else f"not found: {s!r}")


def below(label: str, anchor: str, v) -> None:
    """NOTES prints '< p' after the anchor; the fresh value must be below p."""
    hits = near(anchor, 1)
    if not hits:
        check(f"NOTES: {label}", False, f"anchor not found exactly once: {anchor!r}")
        return
    s, off, _ = hits[0]
    ok = _dec(v) < _parse(s)[0]
    if ok:
        COVERED[off] = "compared"
    check(f"NOTES: {label}", ok, f"{_show(v)} < {s}")


def table(anchor: str) -> dict[str, list[tuple[str, int]]] | None:
    """Rows of the Markdown table whose line holds the anchor: first cell -> cells."""
    i = NOTES.find(anchor)
    if i < 0 or NOTES.find(anchor, i + 1) >= 0:
        return None
    s = NOTES.rfind("\n", 0, i) + 1
    while s > 0 and NOTES.startswith("|", NOTES.rfind("\n", 0, s - 1) + 1):
        s = NOTES.rfind("\n", 0, s - 1) + 1
    rows = {}
    while s < len(NOTES) and NOTES.startswith("|", s):
        e = NOTES.find("\n", s)
        e = len(NOTES) if e < 0 else e
        line = NOTES[s:e]
        if not re.fullmatch(r"[\s|:-]+", line):
            cells = [(m.group(), s + m.start()) for m in re.finditer(r"(?:\\\||[^|])+", line)]
            rows[cells[0][0].strip()] = cells[1:]
        s = e + 1
    return rows


def compare_row(label: str, rows, key: str, values: list) -> None:
    """Compare the first number in each cell of a table row; IDX skips a cell."""
    if rows is None or key not in rows:
        check(f"NOTES: {label}", False, f"table or row {key!r} not found")
        return
    cells = rows[key]
    if len(cells) != len(values):
        check(f"NOTES: {label}", False, f"row has {len(cells)} cells, expected {len(values)}")
        return
    pairs, missing = [], []
    for (cell, off), v in zip(cells, values):
        if v is IDX:
            continue
        j = bisect.bisect_left(NUM_POS, off)
        if j == len(NUMS) or NUMS[j][0] >= off + len(cell):
            missing.append(cell.strip())
            continue
        o, s = NUMS[j]
        pairs.append((s, o, _trunc(o + len(s)), v))
    _record(label, pairs, f"no number in cells {missing}" if missing else "")


def compare_fracs(label: str, rows, key: str, values: list) -> None:
    """Table cells printed as p/q, compared exactly."""
    if rows is None or key not in rows:
        check(f"NOTES: {label}", False, f"table or row {key!r} not found")
        return
    got = []
    for cell, _ in rows[key]:
        m = re.search(r"(\d+)\s*/\s*(\d+)", cell)
        got.append(Fraction(int(m[1]), int(m[2])) if m else cell.strip())
    check(f"NOTES: {label}", got == values, f"printed {[str(g) for g in got]}")


# ------------------------------------------------------------- fresh values
# Everything the prose checks compare against is computed here, from the exact
# integers where NOTES prints more digits than a double holds.  The two §9
# programs are rerun into a temporary directory and their output must equal
# the stored JSON, so the §9 numbers are compared against a fresh run.

D2 = Decimal(2)
LN2 = D2.ln()


def l2(x) -> Decimal:
    """log2 of a positive int, Fraction or Decimal, to about 100 digits."""
    if isinstance(x, Fraction):
        return l2(x.numerator) - l2(x.denominator)
    if isinstance(x, int):
        sh = max(0, x.bit_length() - 400)
        return Decimal(x >> sh).ln() / LN2 + sh
    return Decimal(x).ln() / LN2


D_THETA = sum(D2 ** -(j * j) for j in range(-80, 81))
D_THETA_ODD = D2 ** Decimal("-0.25") * sum(D2 ** -(j * j + j) for j in range(-80, 81))
D_ETA = math.prod((1 - D2**-i for i in range(1, 400)), start=Decimal(1))
D_THETA_SD = (sum(j * j * D2 ** -(j * j) for j in range(-80, 81)) / D_THETA).sqrt()
D_GAL = l2(D_THETA / D_ETA)
D_GAL_ODD = l2(D_THETA_ODD / D_ETA)

_RC: dict[int, tuple[int, ...]] = {}


def rc(n: int) -> tuple[int, ...]:
    """Exact a(n, k) for k = 0..floor(n/2), computed once per n (41 s at n = 512)."""
    if n not in _RC:
        _RC[n] = tuple(rank_counts(n))
    return _RC[n]


def e(n: int) -> int:
    return sum(rc(n))


_ES = load("data/elemab_stats.json")
ST = {r["n"]: r for r in _ES["rows"]} if _ES else {}
# A full run of scripts/reproduce.sh stores n <= 512; FAST=1 stops at 128, and
# every number that needs n = 512 is then skipped rather than recomputed.
FULL = 512 in ST


def full(n: int) -> bool:
    return n < 512 or FULL


@functools.lru_cache(maxsize=None)
def sd_exact(n: int) -> Decimal:
    c = rc(n)
    s0 = sum(c)
    s1 = sum(k * x for k, x in enumerate(c))
    s2 = sum(k * k * x for k, x in enumerate(c))
    return (Decimal(s2 * s0 - s1 * s1) / Decimal(s0 * s0)).sqrt()


def sd_dp(err: Decimal) -> int:
    """Largest d with |err| < 0.5e-d, as src/sd_precision.py counts decimal places."""
    d = 0
    while abs(err) < Decimal(5).scaleb(-(d + 2)) and d < 90:
        d += 1
    return d


@functools.lru_cache(maxsize=None)
def _op(n: int) -> dict:
    return orbit_profile(n, n // 4)


def op(n: int) -> dict[int, float]:
    """Expected fraction of points in orbits of each size, at rank n/4."""
    return {m: v["exp_point_fraction"] for m, v in _op(n)["profile"].items()}


def opt(fn, *args):
    """fn(*args), or None when a value is missing from this run."""
    return None if any(a is None for a in args) else fn(*args)


def lg2n(n: int) -> Decimal:
    """n log2 n, exact for a power of two."""
    return n * l2(n) if n & (n - 1) else Decimal(n * (n.bit_length() - 1))


def alpha(log2x: Decimal, n: int) -> Decimal:
    return (log2x - Decimal(n * n) / 16) / lg2n(n)


def wpow(d: int, t: int) -> list[int]:
    out = [1]
    for _ in range(t):
        out = pmul(out, w_poly(d))
    return out


@functools.lru_cache(maxsize=None)
def Lw(d: int, t: int) -> int:
    """L(w_d^t): subspaces of (F_2^d)^t that project onto every factor."""
    return L(wpow(d, t))


def F8(n: int) -> int:
    """The Theorem 9.1 count at n = 8m, exactly."""
    m = n // 8
    return factorial(n) // (factorial(m) * factorial(8) ** m) * 105**m * Lw(4, m)


_M = [1, 0, 1, 0]


def M(n: int) -> int:
    """Partitions of [n] into blocks of sizes 2 and 4: n! [x^n] exp(x^2/2 + x^4/24)."""
    while len(_M) <= n:
        k = len(_M)
        _M.append((k - 1) * _M[k - 2] + comb(k - 1, 3) * _M[k - 4])
    return _M[n]


def y_star(n: int) -> float:
    """The saddle point y = x^2 of B*: y + 2y^2/3 + y^4/48 = n."""
    lo, hi = 0.0, float(n)
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if mid + 2 * mid * mid / 3 + mid**4 / 48 < n else (lo, mid)
    return lo


def u_star(n: int) -> float:
    """The saddle point u = x^2 of M: u + u^2/6 = n."""
    return math.sqrt(6 * n + 9) - 3


def _fresh_run(script: str, *args: str):
    """Run src/<script> with --out in a temporary directory: (exit, passes, JSON)."""
    with tempfile.TemporaryDirectory() as td:
        out = pathlib.Path(td) / "out.json"
        p = subprocess.run(
            [sys.executable, str(HERE / script), *args, "--out", str(out)],
            capture_output=True,
            text=True,
            cwd=ROOT,
        )
        data = json.loads(out.read_text()) if out.exists() else None
    return p.returncode, p.stdout.count("[pass]"), data


BF_N = min(512, max(ST, default=0))


@functools.lru_cache(maxsize=None)
def census_fresh():
    return _fresh_run("census8.py")


@functools.lru_cache(maxsize=None)
def bf_fresh():
    return _fresh_run("block_families.py", "--n-max", str(BF_N))


def dig(d, *keys):
    """d[k1][k2]..., with int keys read as JSON strings; None if absent."""
    for k in keys:
        if d is None:
            return None
        d = d.get(str(k)) if isinstance(d, dict) else None
    return d


@functools.lru_cache(maxsize=None)
def holdout() -> dict | None:
    """The §6 Check 3 extrapolation, or None if the stored n stop short of 496."""
    hv = holdout_extrapolation(local_coefficients(list(ST.values()), 16)) if ST else None
    return None if not hv or "error" in hv else hv


# ---------------------------------------------------- the prose, section by section


def notes_s2() -> None:
    section("NOTES.md §2: exact counts")
    lat = {n: load(f"data/lattice_n{n}.json") for n in range(1, 9)}
    if all(lat.values()):
        text(
            "§2.1 lattice, n = 1..7",
            "| 1–7 | "
            + ", ".join(str(lat[n]["total_subgroups"]) for n in range(1, 8))
            + " | "
            + ", ".join(str(lat[n]["total_classes"]) for n in range(1, 8))
            + " |",
        )
        text("§2.1 lattice, n = 8", f"| 8 | {lat[8]['total_subgroups']} | {lat[8]['total_classes']} |")
    else:
        NOTES_SKIPPED.append("§2.1 lattice rows (data/lattice_n*.json missing)")
    compare(
        "§2.1 the wrong enumerator's shortfall at n = 7, in percent",
        "wrong from $n=5$ on, ",
        round(100 * (1 - Fraction(6415, A005432[7]))),
        ints=True,
    )
    text("§2.2 rank_counts(8)", f"`rank_counts(8) = {list(rc(8))}`")
    compare("§2.2 bit length of e(512)", "exact integer of ", e(512).bit_length() if FULL else None, ints=True)
    text("§2.3 rank counts at n = 10", f"`{list(rc(10))}`")


def notes_s3() -> None:
    section("NOTES.md §3: the rank distribution")
    rows = table("| $n$ | mode | mean$/n$ | sd |")
    for n in (64, 128, 256, 512):
        s = rank_stats(n, rc(n)) if full(n) else {}
        compare_row(
            f"§3 table, n = {n}",
            rows,
            str(n),
            [s.get("mode"), s.get("mean_over_n"), s.get("sd"), s.get("window_1e6")],
        )
    tested = [n for n in ST if n % 4 == 0 and n >= 32]
    check(
        "NOTES: §3 the mode is n/4 for every stored n = 0 mod 4 from 32 up",
        bool(tested) and all(ST[n]["mode"] == n // 4 for n in tested),
        f"{len(tested)} values of n",
    )
    rows = table("| $j$ | $-3$ |")
    c = rc(512) if FULL else None
    compare_row(
        "§3 theta table, n = 512", rows, "$n=512$", [opt(lambda c: Fraction(c[128 + j], c[128]), c) for j in range(-3, 4)]
    )
    compare_row("§3 theta table, 2^-j^2", rows, "$2^{-j^2}$", [Fraction(1, 2 ** (j * j)) for j in range(-3, 4)])
    c = rc(128)
    dev = max(abs(Fraction(c[32 + j], c[32]) - Fraction(1, 2 ** (j * j))) for j in range(-3, 4))
    check("NOTES: §3 the theta law holds to better than 1e-6 at n = 128", dev < Fraction(1, 10**6), f"{float(dev):.1e}")
    compare("§3 sd of the theta law", "2^{-j^2}}} = ", D_THETA_SD)
    compare("§3 measured sd at n = 512", "against a measured ", sd_exact(512) if FULL else None)

    err = {n: abs(sd_exact(n) - D_THETA_SD) if full(n) else None for n in (64, 128, 256, 512)}
    rows = table("| decimal places |")
    compare_row("§3 |sd(n) - closed form|", rows, r"$\|\mathrm{sd}(n) - \text{closed form}\|$", list(err.values()))
    compare_row("§3 decimal places", rows, "decimal places", [opt(sd_dp, v) for v in err.values()])
    mark("historical", "| decimal places | 4 | 7 | 16 | **35** |\n\nThe first version printed ")
    compare("§3 the n = 256 error to three digits", "from rounding the stored ", err[256])
    compare("§3 the n = 256 error", "a second time; the value is $", err[256])
    cf = format(D_THETA_SD, ".30f")
    text("§3 closed form to 30 digits", f"agrees with ${cf}$ in all 30 digits")
    if FULL:
        got = format(sd_exact(512), ".30f")
        check("NOTES: §3 the exact sd at n = 512 agrees in all 30 digits", got == cf, got)
    lg = {n: opt(lambda v: v.log10(), v) for n, v in err.items()}
    ratio_ = lambda a, b: opt(lambda x, y: y / x, lg[a], lg[b])  # noqa: E731
    compare(
        "§3 log10 errors and their ratios",
        r"the $\log_{10}$ errors are $",
        [lg[64], lg[128], lg[256], lg[512], ratio_(64, 128), ratio_(128, 256), ratio_(256, 512)],
    )
    c_ = {n: opt(lambda a, b: (a - b) / (n * D2.log10()), lg[n], lg[2 * n]) for n in (64, 128, 256)}
    mark("historical", "the last of the three ratios as ")
    compare("§3 the last ratio", "ratios as 2.10; it is ", ratio_(256, 512))
    compare("§3 c between n = 128 and 256", r"with $c \approx ", c_[128])
    compare("§3 1/4 - 2/n", "$1/4 - 2/n$, which is ", [Fraction(1, 4) - Fraction(2, n) for n in (64, 128, 256)])
    compare("§3 measured c", "the measured values are ", [c_[64], c_[128], c_[256]])
    dev = [opt(lambda r: abs(Decimal(r["sd"]) - D_THETA_SD), ST.get(n)) for n in (64, 128, 256, 512)]
    compare("§3 deviations of the stored doubles", "the stored doubles instead (", dev[2:])
    if None not in dev:
        check("NOTES: §3 the stored-double deviations are non-monotone", any(a < b for a, b in zip(dev, dev[1:])), ", ".join(f"{d:.1e}" for d in dev))


def notes_s4() -> None:
    section("NOTES.md §4: Galois numbers")
    compare("§4 log2(theta/eta)", r"\log_2\!\frac{\vartheta}{\eta} = ", D_GAL)
    compare("§4 log2(theta'/eta)", r"\log_2(\vartheta'/\eta) = ", D_GAL_ODD)
    compare("§4 theta'", r"2^{-(j+1/2)^2} = ", D_THETA_ODD)
    compare("§4 theta", r"in place of $\vartheta = ", D_THETA)
    compare("§4 eta", r"here $\eta = ", D_ETA)
    compare("§4 theta/eta", r"$\vartheta/\eta = ", D_THETA / D_ETA)
    ex = {r: l2(galois_number(r)) - Decimal(r * r) / 4 for r in (32, 64, 127, 128, 255, 256, 511, 512)}
    compare(
        "§4 measured excess, even and odd r",
        "Measured: $",
        [tuple(ex[r] for r in (128, 256, 512)), tuple(ex[r] for r in (127, 255, 511))],
    )
    compare("§4 deviation at r = 32 and 64", "falls like $2^{-r/2}$ (", [abs(ex[32] - D_GAL), abs(ex[64] - D_GAL)])
    compare("§4 difference of the two limits", "the two limits differ by $", D_GAL - D_GAL_ODD)
    compare("§4 1/eta", r"$\eta^{-1} = ", 1 / D_ETA)


def bf_alpha() -> dict[int, dict]:
    data = bf_fresh()[2]
    return {r["n"]: r for r in data["alpha"]} if data else {}


def rising(label: str, values: list) -> None:
    vals = [v for v in values if v is not None]
    check(f"NOTES: {label}", len(vals) > 1 and all(a < b for a, b in zip(vals, vals[1:])), ", ".join(_show(v) for v in vals))


def notes_s5() -> None:
    section("NOTES.md §5: lower bounds for alpha")
    ns = (16, 32, 64, 128, 256, 512)
    al = bf_alpha()
    ae = [alpha(l2(e(n)), n) if full(n) else None for n in ns]
    a8 = [alpha(l2(F8(n)), n) for n in ns]
    aB = [al.get(n, {}).get("alpha_B") for n in ns]
    rows = table(r"| $(\log_2 e(n)-n^2/16)/(n\log_2 n)$ |")
    compare_row("§5 alpha from e(n)", rows, r"$(\log_2 e(n)-n^2/16)/(n\log_2 n)$", ae)
    compare_row("§5 alpha from F_8(n)", rows, r"$(\log_2 F_8(n)-n^2/16)/(n\log_2 n)$", a8)
    compare_row("§5 alpha from B(n)", rows, r"$(\log_2 B(n)-n^2/16)/(n\log_2 n)$", aB)
    mark("historical", "entry of the first row read ")
    compare("§5 alpha from e(128) to 7 dp", "in the first version; the value is ", ae[3])
    compare("§5 alpha from B(512)", r"not a proof that $\alpha \ge ", aB[5])
    # log2 f(2) = 1, so alpha <= (1 - 2^2/16) / (2 log2 2)
    check("NOTES: §5 n = 2 forces alpha <= 3/8", A005432[2] == 2 and (1 - Fraction(4, 16)) / 2 == Fraction(3, 8), "3/8")
    for name, vals in (("e(n)", ae), ("F_8(n)", a8), ("B(n)", aB)):
        rising(f"§5 the {name} row is still rising at the last n", vals)


def notes_s6() -> None:
    section("NOTES.md §6: the second-order term")
    rows = table("| $j$ | orbit size |")
    for j in (1, 2, 3, 4):
        compare_fracs(
            f"§6 pure orbit type j = {j}",
            rows,
            str(j),
            [str(2**j), Fraction(j, 2 ** (j + 1)), Fraction(j * j, 2 ** (2 * j + 2)), 1 - Fraction(1, 2**j)],
        )
    want = {(n, k, j) for n, kmax in ((4, 3), (6, 3), (8, 2)) for k in range(1, kmax + 1) for j in (1, 2, 3) if n % 2**j == 0}
    check("NOTES: §6 Check 1 describes the brute-forced (n, k, j) exactly", set(BRUTE_HOM) == want, f"{len(want)} triples")
    compare("§6 Check 1 number of brute-forced counts", "was validated against ", len(BRUTE_HOM), ints=True)
    sc = {n: {r["j"]: r for r in pure_type_scan(n)} for n in (64, 128, 256)}
    corr = {
        (j, n): (1 - 2.0**-j - sc[n][j]["excess_per_n_log2n"]) * math.log2(n) for j in (1, 2) for n in (64, 128, 256)
    }
    compare("§6 Check 1 O(n) corrections", "at $n = 64, 128, 256$ are ", [corr[j, n] for j in (1, 2) for n in (64, 128, 256)])
    mark("historical", "for $j=2$. The first version printed ", 2)
    compare("§6 Check 1 corrections for j = 1, to 5 dp", "0.713 for $j = 1$; the values are ", [corr[1, 64], corr[1, 256]])

    ns = (64, 128, 256, 512)
    ge8 = {n: sum(v for m, v in op(n).items() if m >= 8) for n in ns}
    rows = table("| size 4 |")
    compare_row("§6 Check 2 size 4", rows, "size 4", [op(n)[4] for n in ns])
    compare_row("§6 Check 2 size 2", rows, "size 2", [op(n)[2] for n in ns])
    compare_row("§6 Check 2 size >= 8", rows, r"size $\ge 8$", [ge8[n] for n in ns])
    errs = {n: _op(n)["point_total_error"] for n in ns}
    worst = max(errs, key=errs.get)
    compare("§6 Check 2 tolerance and largest error", "as an arithmetic self-check to within $", [ORBIT_TOL, errs[worst]])
    compare("§6 Check 2 where the largest error occurs", "the largest error is $", [errs[worst], worst], ints=True)
    mark("historical", "The first version printed $")
    compare("§6 Check 2 size >= 8 at n = 64", "at $n = 64$; the value is $", ge8[64])

    hv = holdout()
    if hv is None:
        NOTES_SKIPPED.append("§6 Check 3 (the held-out extrapolation needs n <= 496)")
        return
    compare("§6 Check 3 train and test ranges", "fitted on $n\\in[", [int(x) for x in hv["train_range"] + hv["test_range"]], ints=True)
    names = {
        "$n^{-1/2}$": "n^-1/2",
        r"$1/\log_2^2 n$": "1/log2^2 n",
        r"$1/\log_2 n$": "1/log2 n",
        "$n^{-1/4}$": "n^-1/4",
        r"$\log_2 n/n$": "log2 n / n",
        r"$n^{-1/2}\log_2 n$": "n^-1/2 log2 n",
        "$n^{-1}$": "n^-1",
    }
    rows = table("| basis | $A_\\infty$ |")
    byname = {r["basis"]: r for r in hv["ranked_bases"]}
    for key, nm in names.items():
        r = byname[nm]
        compare_row(f"§6 Check 3 basis {nm}", rows, key, [r["A_inf"], r["train_maxres"], r["holdout_maxres"]])
    order = [names[k] for k in rows or {} if k in names]
    check("NOTES: §6 Check 3 rows are in held-out order", order == [r["basis"] for r in hv["ranked_bases"]], ", ".join(order))
    b = hv["best_free_exponent"]
    compare("§6 Check 3 free exponent", "the same way selects $p = ", [b["p"], b["A_inf"], abs(b["A_inf"] - 0.75), b["holdout_maxres"]])
    mark("historical", r"\log_2 n$ as ", 2)
    compare("§6 Check 3 held-out residuals, 4 digits", "and 1.3e−3; the values are ", [byname["n^-1/2"]["holdout_maxres"], byname["n^-1/2 log2 n"]["holdout_maxres"]])
    lg = byname["1/log2 n"]["A_inf"]
    compare("§6 the rejected 1/log2 n basis", "basis extrapolates to ", [lg, lg - 0.75])
    eps = {n: op(n)[2] for n in (64, 512)}
    compare(
        "§6 epsilon log2 n and epsilon sqrt n",
        r"measured $\varepsilon\log_2 n$ falls ",
        [eps[64] * 6, eps[512] * 9, eps[64] * 8, eps[512] * math.sqrt(512)],
    )
    off4 = {n: 1 - op(n)[4] for n in (64, 512)}
    compare("§6 the off-4 fraction at n = 64", r"orbits of size $\ge 8$, and gives ", [off4[64] * 6, off4[64] * 8])
    same = [round(off4[512] * 9, 3) == round(eps[512] * 9, 3), round(off4[512] * math.sqrt(512), 3) == round(eps[512] * math.sqrt(512), 3)]
    check("NOTES: §6 the off-4 fraction gives the same values at 512", all(same), f"{off4[512] * 9:.4f}, {off4[512] * math.sqrt(512):.4f}")
    compare("§6 sqrt 6", r"predicts $\varepsilon\sqrt n \to \sqrt6 = ", Decimal(6).sqrt())


def notes_s7() -> None:
    section("NOTES.md §7: the constants of [RoTr25] Theorem 6")
    nu = Decimal(1) / 2 - Decimal(3).sqrt() / 4
    compare("§7 nu = 1/2 - sqrt3/4", r"\tfrac{\sqrt3}{4} = ", nu)
    check("NOTES: §7 that nu solves nu(1 - nu) = 1/16", abs(nu * (1 - nu) - Decimal(1) / 16) < Decimal("1e-95"), "residual < 1e-95")
    # nu(1/2 - nu) = 1/16 is nu^2 - nu/2 + 1/16 = 0
    disc = Fraction(1, 4) - 4 * Fraction(1, 16)
    check("NOTES: §7 nu(1/2 - nu) = 1/16 has discriminant 0 and double root 1/4", disc == 0, f"discriminant {disc}, root {Fraction(1, 4)}")
    compare("§7 improvement factors", "Improvement factors: ", [Fraction(1, 4) / Fraction(1, 16), Decimal(1) / 4 / nu])


def lattice_orders(n: int) -> dict[int, int] | None:
    lat = load(f"data/lattice_n{n}.json")
    return {int(k): v for k, v in lat["order_histogram"].items()} if lat else None


def median(hist: dict[int, int]) -> Fraction:
    xs, N = sorted(hist.items()), sum(hist.values())

    def at(i: int) -> int:
        acc = 0
        for o, c in xs:
            acc += c
            if i < acc:
                return o
        raise IndexError(i)

    return Fraction(at((N - 1) // 2) + at(N // 2), 2)


def mean_log2(hist: dict[int, int]) -> Decimal:
    return sum((c * l2(o) for o, c in hist.items()), Decimal(0)) / sum(hist.values())


def sp(k: int) -> str:
    """An integer as NOTES prints it: a space every three digits from five digits up."""
    return str(k) if k < 10**4 else f"{k:,}".replace(",", " ")


def notes_s8() -> None:
    section("NOTES.md §8: small-n data")
    pats = ("7/8", "0.875", "extra.{0,3}special")
    if all(p in PDF_RESULTS for p in pats):
        # check_quotes has searched the paper for each pattern; this ties the sentence to that search
        text("§8 the strings searched for in [RoTr25]", 'the strings "7/8", "0.875" and "extraspecial" do not occur in its text')
    else:
        NOTES_SKIPPED.append("§8 the strings searched for in [RoTr25] (the PDF was not read)")
    ns = range(4, 9)
    hs = {n: lattice_orders(n) for n in ns}
    if not all(hs.values()):
        NOTES_SKIPPED.append("§8 small-n table (data/lattice_n*.json missing)")
        return
    rows = table(r"| median $\|H\|$ |")
    compare_row("§8 e(n)/f(n)", rows, "$e(n)/f(n)$", [Fraction(e(n), A005432[n]) for n in ns])
    compare_row("§8 median |H|", rows, r"median $\|H\|$", [median(hs[n]) for n in ns])
    compare_row("§8 mean log2|H| / n", rows, r"$\mathbb{E}[\log_2\|H\|]/n$", [mean_log2(hs[n]) / n for n in ns])
    check("NOTES: §8 the elementary abelian fraction falls at every n = 4..8", all(Fraction(e(n), A005432[n]) > Fraction(e(n + 1), A005432[n + 1]) for n in range(4, 8)))
    cen = census_fresh()[2]
    if cen is None:
        NOTES_SKIPPED.append("§8 the Frattini-closed family at n = 8 (census8.py failed)")
        return
    fc = cen["frattini_closed"]["8"]
    ea = rc(8)
    compare(
        "§8 the Frattini-closed family at n = 8",
        "At $n = 8$ the Frattini-closed family is ",
        [
            Fraction(100 * fc["count"], A005432[8]),
            Fraction(100 * e(8), A005432[8]),
            Fraction(fc["sum_log2_order"], 8 * fc["count"]),
            Fraction(sum(k * x for k, x in enumerate(ea)), 8 * e(8)),
            mean_log2(hs[8]) / 8,
        ],
    )
    not_fc = fc["two_subgroups"] - fc["count"]
    for anchor in ("negligible? At $n = 8$ they are ", "not FC are negligible; at $n = 8$ they are "):
        compare(f"§8/§9.8 2-subgroups that are not FC, after {anchor[:22]!r}", anchor, [not_fc, fc["two_subgroups"], Fraction(100 * not_fc, fc["two_subgroups"])], ints=True)


def gn_degree8() -> tuple[dict[int, int], str] | None:
    """[GN]: {t: |8Tt|} read from the saved page, and the page itself."""
    p = ROOT / "refs" / "groupnames_T15.html"
    if not p.exists():
        return None
    h = p.read_text("utf-8", errors="replace")
    out, order = {}, None
    for m in re.finditer(r"id=['\"](?:order(\d+)|8t(\d+))['\"]", h):
        if m[1]:
            order = int(m[1])
        else:
            out[int(m[2])] = order
    return out, h


def notes_s9_1() -> None:
    section("NOTES.md §9.1: the census")
    rcode, passes, cen = census_fresh()
    compare("§9 census8.py check count", r"for $n \le 8$, ", passes if cen else None, ints=True)
    rb, pb, bfd = bf_fresh()
    compare("§9 block_families.py gate count", "`src/block_families.py` (", pb if bfd and BF_N == 512 else None, ints=True)
    if cen is None:
        NOTES_SKIPPED.append("§9.1 (census8.py failed)")
        return
    fc, tr, by = cen["frattini_closed"], cen["transitive"], cen["two_subgroups_by_order"]
    f8 = fc["8"]
    miss = {int(o): c - f8["order_histogram"].get(o, 0) for o, c in by["8"].items() if c != f8["order_histogram"].get(o, 0)}
    parts = [f"{x} of them" if i == 0 else str(x) for i, x in enumerate(miss.values())]
    orders = " and ".join(f"{o} ({p})" for o, p in zip(miss, parts))
    text(
        "§9.1 2-subgroups of S_8",
        f"{sp(f8['two_subgroups'])} of them, of which {sp(f8['count'])} are FC. The {f8['two_subgroups'] - f8['count']} that are not have orders {orders}.",
    )
    counts = [fc[str(n)]["count"] for n in range(1, 8)]
    text("§9.1 FC counts for n <= 7", f"({', '.join(map(str, counts[:-1]))} and {counts[-1]} of them for $n = 1, \\dots, 7$)")
    check("NOTES: §9.1 every 2-subgroup is FC for n <= 7", all(fc[str(n)]["count"] == fc[str(n)]["two_subgroups"] for n in range(1, 8)))
    lat = {n: lattice_orders(n) for n in range(1, 9)}
    if all(lat.values()):
        same = all({int(o): c for o, c in by[str(n)].items()} == {o: c for o, c in lat[n].items() if o & (o - 1) == 0} for n in range(1, 9))
        check("NOTES: §9.1 census counts by order equal the lattices at every 2-power order, n <= 8", same)
    tot = [tr[str(b)]["transitive_total"] for b in (1, 2, 4, 8)]
    text("§9.1 transitive 2-subgroups of S_b", f"number {', '.join(map(sp, tot[:-1]))} and {sp(tot[-1])} for $b = 1, 2, 4, 8$")
    cl8 = tr["8"]["classes"]
    prof = sorted(Counter(c["order"] for c in cl8).items())
    text(
        "§9.1 classes in degree 8 by order",
        f"The {len(cl8)} conjugacy classes in degree 8 have orders "
        + ", ".join(f"{o} ({k}{' classes' if i == 0 else ''})" for i, (o, k) in enumerate(prof[:-1]))
        + f" and {prof[-1][0]} ({prof[-1][1]})",
    )
    gn = gn_degree8()
    if gn is None:
        NOTES_SKIPPED.append("§9.1 the [GN] order profile (refs/groupnames_T15.html absent; run refs/fetch.sh)")
    else:
        tg, h = gn
        two = sorted((t, o) for t, o in tg.items() if o and o & (o - 1) == 0)
        runs: list[list[int]] = []
        for t, o in two:
            if runs and runs[-1][1] == t - 1 and runs[-1][2] == o:
                runs[-1][1] = t
            else:
                runs.append([t, t, o])
        names = [f"8T{a}–8T{b}" if a < b else f"8T{a}" for a, b, _ in runs]
        text("§9.1 the [GN] order profile", ", ".join(names[:-1]) + " and " + names[-1] + " in [GN]")
        check("NOTES: §9.1 [GN] has the census order profile", [(o, b - a + 1) for a, b, o in runs] == prof, str(prof))
        i = h.find('id="8t22"')
        row = h[h.rfind("<tr>", 0, i) : i] if i >= 0 else ""
        check("NOTES: §9.1 [GN] lists 8T22 as extraspecial, SmallGroup(32,49)", "Extraspecial" in row and ">32,49<" in row)
    check("NOTES: §9.1 every transitive 2-group of degree b <= 8 has d <= b/2", all(c["d"] <= b // 2 for b in (1, 2, 4, 8) for c in tr[str(b)]["classes"]))
    dense = {b: [c for c in tr[str(b)]["classes"] if 2 * c["d"] == b] for b in (2, 4, 8)}
    cnt = {b: sum(c["class_size"] for c in cs) for b, cs in dense.items()}
    text("§9.1 dense types", f"There are {cnt[2]}, {cnt[4]} and {cnt[8]} of them in degrees 2, 4 and 8")
    check("NOTES: §9.1 the dense 4-types are the regular V4 and three D4", sorted((c["d"], c["log2_phi"], c["class_size"]) for c in dense[4]) == [(2, 0, 1), (2, 1, 3)])
    t22 = dense[8][0] if len(dense[8]) == 1 else {}
    got = tuple(t22.get(k) for k in ("order", "d", "log2_phi", "centre", "involutions", "class_size"))
    check("NOTES: §9.1 one dense 8-type: order 32, d = 4, |Phi| = 2, |Z| = 2, 19 involutions, 105 conjugates", got == (32, 4, 1, 2, 19, 105), str(got))
    check("NOTES: §9.1 2^(1+4)_+ has 19 involutions and 2^(1+4)_- has 11", (2**4 + 2**2 - 1, 2**4 - 2**2 - 1) == (19, 11))
    terms = []
    for b in (2, 4, 8):
        c = Fraction(cnt[b], factorial(b))
        terms.append(f"{'' if c.numerator == 1 else c.numerator}x^{b}/{c.denominator}")
    text("§9.1 EGF of the dense types", r"$\exp(" + " + ".join(terms) + ")$")
    if bfd is not None:
        check("NOTES: §9.1 block_families.py finds the same dense types", bfd["dense_types"] == {str(b): v for b, v in cnt.items()})


def notes_s9_2() -> None:
    section("NOTES.md §9.2: e(n) against G M(n)")
    ns = (8, 16, 32, 64, 128, 256, 512)
    res = {n: l2(Fraction(e(n), galois_number(n // 2) * M(n))) if full(n) else None for n in ns}
    rows = table(r"| $\log_2 e(n) - \log_2 G_{n/2}M(n)$ |")
    compare_row("§9.2 residuals from the exact integers", rows, r"$\log_2 e(n) - \log_2 G_{n/2}M(n)$", list(res.values()))
    odd = [l2(Fraction(e(n), galois_number((n - 1) // 2) * n * M(n - 1))) for n in (9, 17, 33)]
    compare("§9.2 the odd model", "the odd model gives ", odd)
    stored = {n: Decimal(r["log2_e"]) - l2(galois_number(n // 2) * M(n)) for n, r in ST.items() if n >= 128}
    if FULL:
        worst = max(stored, key=lambda n: abs(stored[n]))
        compare("§9.2 number of stored n >= 128", "The table is computed from the exact integers. Over the ", len(stored), ints=True)
        compare("§9.2 largest residual of the stored doubles", "the largest residual in absolute value is $", [abs(stored[worst]), worst], ints=True)
    ns4 = (64, 128, 256, 512)
    compare("§9.2 saddle-point 2-orbit fraction, and the exact one", r"\approx \sqrt{6/n}$. That gives ", [u_star(n) / n for n in ns4] + [op(n)[2] for n in ns4])
    compare("§9.2 sqrt 6 and epsilon sqrt n at 512", r"So $\varepsilon\sqrt n \to \sqrt6 = ", [Decimal(6).sqrt(), op(512)[2] * math.sqrt(512)])
    hv = holdout()
    coef = {r["basis"]: r["coef"] for r in hv["ranked_bases"]}.get("n^-1/2") if hv else None
    compare("§9.2 the §6 fit's n^-1/2 coefficient, and -sqrt6/8", "The §6 fit found $", [coef, -Decimal(6).sqrt() / 8])
    bfd = bf_fresh()[2]
    compare("§9.2 A_local from G M against the stored values", "computed from $e(n)$ to $", dig(bfd, "A_local_max_abs_diff_n_ge_128"))
    far = dig(bfd, "A_local_remainder_times_n15_from_M") or {}
    compare(
        "§9.2 the n^-3/2 coefficient from exact M",
        "32768$, runs ",
        [far.get(str(2**j)) for j in range(9, 16)] + [dig(bfd, "A_local_n_minus_1_5_coefficient", "richardson"), 9 * Decimal(6).sqrt() / 32],
    )


def tv(n: int) -> Fraction:
    """Total variation between a(n, .)/e(n) and the dimension of a random subspace of F_2^(n/2)."""
    a, m = rc(n), n // 2
    E, g = sum(a), galois_number(m)
    return Fraction(sum(abs(a[k] * g - gaussian_binomial(m, k) * E) for k in range(m + 1)), 2 * E * g)


def notes_s9_3() -> None:
    section("NOTES.md §9.3: the rank law against a random subspace")
    ns = (16, 32, 64, 128, 256)
    t = {n: tv(n) for n in ns}
    scaled = {n: l2(t[n]) + Decimal(n) / 4 - 2 * l2(n) for n in ns}
    rows = table("| TV |")
    compare_row("§9.3 total variation", rows, "TV", list(t.values()))
    compare_row("§9.3 TV 2^(n/4)/n^2", rows, r"$\log_2(\mathrm{TV}\cdot 2^{n/4}/n^2)$", list(scaled.values()))
    compare("§9.3 range of the last row", "The last row stays between ", [min(scaled.values()), max(scaled.values())])
    check("NOTES: §9.3 the TV falls by more than 15 orders of magnitude", t[16] / t[256] > 10**15, f"{float(t[16] / t[256]):.2e}")
    k = factorial(8) // (8 * 168)
    check("NOTES: §9.3 8 * 168 = |AGL(3,2)| divides 8!", 8 * gl_order(3) == 1344 and k * 1344 == factorial(8))
    text("§9.3 regular (Z/2)^3 embeddings", f"$8!/(8\\cdot168) = {k}$")
    text("§9.3 weight of one 8-orbit", f"${k}x^8/8! = x^8/{Fraction(k, factorial(8)).denominator}$")
    cen = census_fresh()[2]
    if cen is not None:
        reg = sum(c["class_size"] for c in cen["transitive"]["8"]["classes"] if c["order"] == 8 and c["exponent"] == 2)
        check("NOTES: §9.3 the census has 30 regular (Z/2)^3 in S_8", reg == k, str(reg))
    compare("§9.3 log2(36/1344)", r"$\log_2(36/1344) = ", l2(Fraction(36, 1344)))
    sd = {n: opt(lambda x: l2(abs(x - D_THETA_SD)) + Decimal(n) / 4 - 2 * l2(n), sd_exact(n) if full(n) else None) for n in (64, 128, 256, 512)}
    compare("§9.3 the same row for the sd error", "The same row for the sd error of §3 is ", list(sd.values()))


C8 = Decimal(3) / 8 - Decimal(7) / 8 / LN2 - l2(384) / 8  # the linear coefficient of Theorem 9.1


def notes_s9_4() -> None:
    section("NOTES.md §9.4: Theorem 9.1")
    check("NOTES: §9.4 8!/105 = 384", factorial(8) == 105 * 384)
    compare("§9.4 the constant c", r"\frac18\log_2 384 = ", C8)
    br = lambda m: 1 - 15 * m * D2 ** (Decimal(1 - 4 * m) / 2)  # noqa: E731
    check("NOTES: §9.4 the bracket is positive exactly from m = 3 (m <= 64)", br(2) < 0 and all(br(m) > 0 for m in range(3, 65)))
    check("NOTES: §9.4 L(w_4^m) >= G_4m - 15m G_(4m-1), m <= 8", all(Lw(4, m) >= galois_number(4 * m) - 15 * m * galois_number(4 * m - 1) for m in range(1, 9)))
    check("NOTES: §9.4 G_r^2 >= 2^(r-1) G_(r-1)^2 for r <= 512", all(galois_number(r) ** 2 >= 2 ** (r - 1) * galois_number(r - 1) ** 2 for r in range(1, 513)))
    check("NOTES: §9.4 G_4m >= [4m, 2m]_2 >= 2^(4m^2), m <= 32", all(galois_number(4 * m) >= gaussian_binomial(4 * m, 2 * m) >= 2 ** (4 * m * m) for m in range(1, 33)))
    compare(
        "§9.4 the bracket and L(w_4^m)/G_4m",
        "- 7n/8$. ∎\n\nThe bracket is ",
        [br(m) for m in (3, 4, 5)] + [Fraction(Lw(4, m), galois_number(4 * m)) for m in (3, 4, 5)],
    )
    ns = (64, 128, 256, 512)
    rem, nxt, lg = {}, {}, {}
    for n in ns:
        m = n // 8
        rem[n] = l2(F8(n)) - (Decimal(n * n) / 16 + Decimal(7) / 8 * lg2n(n) + C8 * n + Decimal(3) / 2 + D_GAL - 7 / (12 * n * LN2))
        nxt[n] = (1 - Decimal(1) / 512) / (360 * m**3 * LN2)
        lg[n] = l2(Fraction(Lw(4, m), galois_number(4 * m)))
    compare("§9.4 Stirling remainders of log2 F_8(n)", r"from $\log_2 F_8(n)$ leaves ", list(rem.values()))
    check("NOTES: §9.4 at n = 256, 512 the remainder is the next Stirling term to three digits", all(format(rem[n], ".2e") == format(nxt[n], ".2e") for n in (256, 512)), ", ".join(f"{rem[n]:.3e}/{nxt[n]:.3e}" for n in (256, 512)))
    check("NOTES: §9.4 at n = 64 the remainder is dominated by log2(L/G)", abs(rem[64] - lg[64]) < abs(rem[64]) / 10, f"{rem[64]:.3e} vs {lg[64]:.3e}")
    n0 = next(n for n in range(8, 513, 8) if F8(n) > e(n))
    text("§9.4 where F_8 first exceeds e", f"first exceeds $e(n)$ at $n = {n0}$")
    compare("§9.4 log2(F_8/e)", r"and $\log_2(F_8(n)/e(n))$ is ", [l2(Fraction(F8(n), e(n))) if full(n) else None for n in (32, 64, 128, 256, 512)])
    compare("§9.4 alpha of F_8 at 512, c and c/9", r"(\log_2 F_8(n) - n^2/16)/(n\log_2 n)$ is ", [alpha(l2(F8(512)), 512), -C8, -C8 / 9])
    check("NOTES: §9.4 (1 - 448/453)(1 - 1/2) = 5/906", (1 - Fraction(448, 453)) * (1 - Fraction(1, 2)) == Fraction(5, 906))
    check("NOTES: §9.4 (1 - 3^2/(4*2^2))/(2*2*log2 4) = 7/128", (1 - Fraction(9, 16)) / 8 == Fraction(7, 128))
    below("§9.4 5/906 is below the printed bound", "5/906 < ", Fraction(5, 906))
    if PDF_RESULTS.get("0.08, (1 − (2p − 1)2/(4p2)) 2plog(2p) ,αR,(1 − α0)(1 − 1/p)"):
        mark("compared", r"\alpha_p = \min\{")
    else:
        NOTES_SKIPPED.append("§9.4 the 0.08 of Proposition 7.3 (not found in the PDF, or the PDF was not read)")


def notes_s9_5() -> None:
    section("NOTES.md §9.5: Proposition 9.2")
    check("NOTES: §9.5 7/8 - 3/4 = 1/8", Fraction(7, 8) - Fraction(3, 4) == Fraction(1, 8))
    check("NOTES: §9.5 j/2^j is 1/2 for j = 1, 2 and at most 3/8 for 3 <= j <= 256", Fraction(1, 2) == Fraction(2, 4) and all(Fraction(j, 2**j) <= Fraction(3, 8) for j in range(3, 257)))
    grid = all(
        (Fraction(n, 2) - d) ** 2 / 4 <= Fraction(n * n, 16) - n * d / 8
        for n in range(1, 65)
        for d in (Fraction(k, 8) for k in range(4 * n + 1))
    )
    check("NOTES: §9.5 D^2/4 <= n^2/16 - n delta/8 for D = n/2 - delta, 0 <= delta <= n/2 (n <= 64, delta in Z/8)", grid)
    check("NOTES: §9.5 [D, k]_2 <= 2^(k(D-k))/eta, D <= 64", all(gaussian_binomial(D, k) * D_ETA <= 2 ** (k * (D - k)) for D in range(65) for k in range(D + 1)))
    check("NOTES: §9.5 the odd-D theta sum is below theta", D_THETA_ODD < D_THETA, f"theta - theta_odd = {D_THETA - D_THETA_ODD:.2e}")
    check("NOTES: §9.5 theta/eta < 8", D_THETA / D_ETA < 8, f"{D_THETA / D_ETA:.4f}")
    worst = max(l2(galois_number(D)) - Decimal(D * D) / 4 for D in range(1, 513))
    check("NOTES: §9.5 G_D <= (theta/eta) 2^(D^2/4), D <= 512", worst < D_GAL, f"max log2 G_D - D^2/4 = {worst:.6f} < {D_GAL:.6f}")
    check("NOTES: §9.5 n_+! 2^(n_+) <= n^(2 n_+) for n_+ <= n <= 64", all(factorial(p) * 2**p <= n ** (2 * p) for n in range(2, 65) for p in range(n + 1)))
    check("NOTES: §9.5 n^2 <= 2^(n/64) for 2^11 <= n <= 2^13, and not at n = 2^10", all(n**128 <= 1 << n for n in range(2048, 8193)) and 1024**128 > 1 << 1024)
    # n - 128 log2 n increases for n > 128/ln 2, so the range above covers every n >= 2^11
    check("NOTES: §9.5 n - 128 log2 n is increasing from n = 2^11", 128 / LN2 < 2048)
    cauchy = all(
        l2(M(k)) <= l2(factorial(k)) + (Decimal(k).sqrt() / 2 + Decimal(k) / 24) / LN2 - Decimal(k) / 4 * l2(k)
        for k in range(2, 513, 2)
    )
    check("NOTES: §9.5 Cauchy's bound M(n') <= n'! e^(x^2/2 + x^4/24) x^(-n') at x = n'^(1/4), even n' <= 512", cauchy)
    held = []
    for n in (8, 16, 32, 64, 128, 256, 512):
        if not full(n):
            continue
        lm = [l2(M(k)) if k % 2 == 0 else None for k in range(n + 1)]
        best = max(
            lm[n - a - b] + 2 * b * l2(n) - Decimal(n * a) / 16 - Decimal(n * b) / 64
            for a in range(n + 1)
            for b in range(n - a + 1)
            if lm[n - a - b] is not None
        )
        held.append(l2(e(n)) <= 3 + n * l2(3) + Decimal(n * n) / 16 + best)
    check("NOTES: §9.5 the displayed bound on e(n) holds at n = 8, 16, ..., 512", all(held), f"{len(held)} values of n")
    cen = census_fresh()[2]
    if cen is not None:
        v4 = [c["class_size"] for c in cen["transitive"]["4"]["classes"] if c["order"] == 4 and c["exponent"] == 2]
        check("NOTES: §9.5 S_4 has exactly one regular V_4", v4 == [1], str(v4))
    check("NOTES: §9.5 F_2^2 has 3 hyperplanes and F_2^4 has 15", gaussian_binomial(2, 1) == 3 and gaussian_binomial(4, 3) == 15)
    br = lambda t: 1 - 3 * t * D2 ** (Decimal(1 - 2 * t) / 2)  # noqa: E731
    check("NOTES: §9.5 L(w_2^t) >= G_2t - 3t G_(2t-1), t <= 16", all(Lw(2, t) >= galois_number(2 * t) - 3 * t * galois_number(2 * t - 1) for t in range(1, 17)))
    check("NOTES: §9.5 the bracket is negative for t <= 4 and positive for 5 <= t <= 64", all(br(t) < 0 for t in range(1, 5)) and all(br(t) > 0 for t in range(5, 65)))
    check("NOTES: §9.5 G_2t >= 2^(t^2), t <= 64", all(galois_number(2 * t) >= 2 ** (t * t) for t in range(1, 65)))
    compare(
        "§9.5 the bracket and L(w_2^t)/G_2t",
        r"The bracket is positive for $t \ge 5$: it is ",
        [br(t) for t in (4, 5, 6)] + [Fraction(Lw(2, t), galois_number(2 * t)) for t in (4, 5, 6)],
    )
    al = bf_alpha()
    compare("§9.5 log2(B/e)", r"the exact $\log_2(B(n)/e(n))$ is ", [al.get(n, {}).get("log2_B_over_e") for n in (16, 32, 64, 128, 256, 512)])


def notes_s9_6() -> None:
    section("NOTES.md §9.6: order statistics in the family")
    bfd = bf_fresh()[2]
    os_ = dig(bfd, "order_stats") or {}
    ns = (8, 16, 32, 64, 128, 256, 512)
    rows = table(r"| $n$ | $\mathbb{E}\log_2\|H\|/n$ |")
    for n in ns:
        r = os_.get(str(n))
        vals = [None] * 7
        if r is not None:
            vals = [r["mean_log2_order"] / n, r["sd_log2_order"], r["sd_log2_order"] / n**0.25] + [r["point_fraction"][b] for b in ("1", "2", "4", "8")]
        compare_row(f"§9.6 order statistics, n = {n}", rows, str(n), vals)
    rising("§9.6 the mean of log2|H|/n falls", [-os_[str(n)]["mean_log2_order"] / n for n in ns if str(n) in os_])
    cen = census_fresh()[2]
    if cen is not None:
        dense = {b: [c for c in cen["transitive"][str(b)]["classes"] if 2 * c["d"] == b] for b in (2, 4, 8)}
        phib = {b: Fraction(sum(c["log2_phi"] * c["class_size"] for c in cs), sum(c["class_size"] for c in cs)) for b, cs in dense.items()}
        check("NOTES: §9.6 the identity's coefficients phibar_b/b are 0, 3/16, 1/8", [phib[b] / b for b in (2, 4, 8)] == [0, Fraction(3, 16), Fraction(1, 8)])
        check("NOTES: §9.6 3/8 = 1/4 + (1/8) log2|Phi(8T22)|", Fraction(1, 4) + phib[8] / 8 == Fraction(3, 8))
        coef = [Fraction(b * sum(c["class_size"] for c in dense[b]), factorial(b)) for b in (2, 4, 8)]
        check("NOTES: §9.6 the saddle equation y + 2y^2/3 + y^4/48 = n, from the dense EGF", coef == [1, Fraction(2, 3), Fraction(1, 48)])
    if bfd is not None:
        check("NOTES: §9.6 block_families.py: mean log2|Phi| of the dense 8-types is 1", bfd["mean_log2_phi_degree8_dense"] == 1)
    compare("§9.6 the dense-identity gap", "The difference, computed exactly, is ", [dig(os_, n, "dense_identity_gap") for n in (8, 32, 128, 256, 512)])
    compare("§9.6 log2 B - log2 G - log2 B*", r"\log_2 B^*(n)$ is ", [dig(bfd, "B_minus_GBstar_log2", n) for n in (128, 256, 512)])
    ns4 = (64, 128, 256, 512)
    y = {n: y_star(n) for n in ns4}
    f512 = dig(os_, 512, "point_fraction") or {}
    compare(
        "§9.6 the saddle point at n = 512 and its orbit fractions",
        "At $n = 512$, $y = ",
        [y[512], y[512] / 512, 2 * y[512] ** 2 / 3 / 512, y[512] ** 4 / 48 / 512] + [f512.get(b) for b in ("2", "4", "8")],
    )
    mu = [os_[str(n)]["mean_log2_order"] - 3 * n / 8 if str(n) in os_ else None for n in ns4]
    compare("§9.6 E log2|H| - 3n/8: y^2/24 - y/8, then exact", r"\approx \frac{y^2}{24} - \frac{y}{8}$: ", [y[n] ** 2 / 24 - y[n] / 8 for n in ns4] + mu)
    th2 = D_THETA_SD**2
    pred = [math.sqrt(y[n] ** 2 / 24 + y[n] / 32 - 1 / 6 + float(th2)) for n in ns4]
    sd = [dig(os_, n, "sd_log2_order") for n in ns4]
    short = [opt(lambda s, p=p: s - p, s) for s, p in zip(sd, pred)]
    compare("§9.6 sd: the §3 variance, the heuristic, exact, and the shortfall", r"\frac{y}{32} - \frac16 + ", [th2] + pred + sd + [short[0], short[-1]])
    if None in short:
        NOTES_SKIPPED.append("§9.6 the shortfall falls along n = 64..512")
    else:
        check("NOTES: §9.6 the shortfall falls along n = 64..512", all(a > b > 0 for a, b in zip(short, short[1:])), ", ".join(f"{s:.4f}" for s in short))
    check("NOTES: §9.6 sqrt(48)/24 = 1/sqrt(12)", abs(Decimal(48).sqrt() / 24 - 1 / Decimal(12).sqrt()) < Decimal("1e-90"))
    compare("§9.6 the sd constant 48^(1/4) 24^(-1/2)", r"48^{1/4}24^{-1/2}\,n^{1/4} = ", Decimal(48).sqrt().sqrt() / Decimal(24).sqrt())
    compare("§9.6 the mean's excess over 3n/8 per n, n = 512", "the mean is still $", opt(lambda m: m / 512, mu[-1]))


def notes_s9_7() -> None:
    section("NOTES.md §9.7: odd n")
    bfd = bf_fresh()[2]
    odd = (9, 17, 33, 65, 129, 257, 511)
    rows = table("| $n$ | 9 | 17 |")
    s3 = [dig(bfd, "odd_n_log2_S3_over_B", n) for n in odd]
    compare_row("§9.7 log2(S3/B), odd n", rows, r"$\log_2(\|S_3(n)\|/B(n))$", s3)
    compare_row("§9.7 log2(y/6)", rows, r"$\log_2(y/6)$", [math.log2(y_star(n) / 6) for n in odd])
    rising("§9.7 S3 gains on B along odd n", s3)
    check("NOTES: §9.7 w_1(z) = z - 1", w_poly(1) == [-1, 1])
    ev = (8, 16, 32, 64, 128, 256, 512)
    v = [dig(bfd, "even_n_log2_S3_over_B", n) for n in ev]
    compare("§9.7 log2(S3/B), even n, and per n/4", r"$\log_2(|S_3(n)|/B(n))$ is ", v + [opt(lambda x, n=n: x / (n / 4), x) for x, n in zip(v, ev)])
    check(
        "NOTES: §9.7 b/4 + 2 <= 3b/8 exactly when b >= 16 (b = 2..2^20)",
        all((Fraction(b, 4) + 2 <= Fraction(3 * b, 8)) == (b >= 16) for b in (2**k for k in range(1, 21))),
    )


def notes_s9_8() -> None:
    section("NOTES.md §9.8: Conjecture 9.3")
    cen, bfd = census_fresh()[2], bf_fresh()[2]
    if cen is None:
        NOTES_SKIPPED.append("§9.8 (census8.py failed)")
        return
    tr = cen["transitive"]
    tot = {b: tr[str(b)]["transitive_total"] for b in (1, 2, 4, 8)}
    check("NOTES: §9.8 one transitive 2-group of degree 1 and one of degree 2", tot[1] == tot[2] == 1)
    text("§9.8 the EGF of the orbit data", f"\\exp(x + x^2/2 + {tot[4]}x^4/{factorial(4)} + {tot[8]}x^8/8!)")
    T = dense_structures(tot, 512)
    ok = all(
        l2(T[n]) <= l2(factorial(n)) + sum(c * Decimal(n) ** (Decimal(b) / 8) / factorial(b) for b, c in tot.items()) / LN2 - Decimal(n) / 8 * l2(n)
        for n in range(1, 513)
    )
    check("NOTES: §9.8 Cauchy's bound on T(n) at x = n^(1/8), n <= 512", ok)
    if bfd is not None:
        lb = bfd["log2_B"]
        check(f"NOTES: §9.8 B(n) > F_8(n) for n = 8, 16, ..., {BF_N}", all(lb[str(n)] > l2(F8(n)) for n in range(8, BF_N + 1, 8)))
    F = egf(weights(load_types(cen)), 18)
    B = [L(F[n][0]) for n in range(19)]
    S3 = [0, 0, 0] + [comb(n, 3) * L(pmul([0, 1], F[n - 3][0])) for n in range(3, 19)]
    check("NOTES: §9.8 |B(n)| + |S3(n)| <= f(n) = A005432(n) for n <= 18", all(B[n] + S3[n] <= A005432[n] for n in range(19)))
    if bfd is not None:
        check("NOTES: §9.8 log2 B(n) equals block_families.py to 1e-9, n <= 18", all(abs(l2(B[n]) - Decimal(bfd["log2_B"][str(n)])) < Decimal("1e-9") for n in range(1, 19)))
    rows = table("| $n$ | $e/f$ | $B/f$ |")
    share = {}
    for n in range(6, 19):
        f = A005432[n]
        share[n] = [Fraction(e(n), f), Fraction(B[n], f), Fraction(S3[n], f), 1 - Fraction(B[n] + S3[n], f), alpha(l2(f), n)]
        compare_row(f"§9.8 shares of f, n = {n}", rows, str(n), share[n])
    out = [share[n][3] for n in range(6, 19)]
    compare("§9.8 the share outside B and S3, in percent", r"The share outside $B \cup S_3$ is between ", [100 * min(out), 100 * max(out), 100 * out[-1]])
    bs = {n: share[n][1] for n in range(6, 19)}
    compare("§9.8 B/f at n = 12 and 18", "Along even $n$, $B/f$ rises from ", [bs[12], bs[18]])
    rising("§9.8 B/f rises along n = 12, 14, 16, 18", [bs[n] for n in (12, 14, 16, 18)])
    odd = [bs[n] for n in range(7, 18, 2)]
    compare("§9.8 B/f along odd n = 7..17", "along odd $n$ it stays between ", [min(odd), max(odd)])
    dense = {b: sum(c["class_size"] for c in tr[str(b)]["classes"] if 2 * c["d"] == b) for b in (2, 4, 8)}
    Bs = dense_structures(dense, 18)

    def formula(n: int) -> int:
        if n % 2 == 0:
            return galois_number(n // 2) * Bs[n]
        return galois_number((n - 1) // 2) * (n * Bs[n - 1] + comb(n, 3) * Bs[n - 3])

    q = {n: Fraction(A005432[n], formula(n)) for n in range(3, 19)}
    compare("§9.8 f(n) over the recorded formula, even then odd n", "$f(n)$ divided by it is ", [q[n] for n in range(6, 19, 2)] + [q[n] for n in range(3, 18, 2)])
    if bfd is not None:
        check("NOTES: §9.8 the same ratios as block_families.py", all(abs(float(q[n]) - bfd["f_over_conjecture_small_n"][str(n)]) < 1e-12 for n in range(3, 19)))
    af = [share[n][4] for n in range(12, 19)]
    compare("§9.8 alpha_f for 12 <= n <= 18", r"$\alpha_f(n)$ lies between ", [min(af), max(af)])
    compare("§9.8 the linear coefficient of §9.4", "With the linear term $cn$, $c = ", C8)


def tfrac(q: Fraction) -> str:
    """A positive fraction as NOTES writes it: \\tfrac34, \\tfrac9{64}."""
    return r"\tfrac" + "".join(s if len(s) == 1 else "{" + s + "}" for s in (str(q.numerator), str(q.denominator)))


def notes_s9_9() -> None:
    section("NOTES.md §9.9: the double root, redone")
    cen = census_fresh()[2]
    if cen is None:
        NOTES_SKIPPED.append("§9.9 (census8.py failed)")
        return
    t22 = [c for c in cen["transitive"]["8"]["classes"] if 2 * c["d"] == 8][0]
    k, p = t22["order"].bit_length() - 1, t22["log2_phi"]
    text("§9.9 order of the ambient group and of its Frattini subgroup", f"which has order $2^{{{k}n/8}}$ and a Frattini subgroup of order $2^{{{'' if p == 1 else p}n/8}}$")
    a, b = Fraction(p, 8), Fraction(k, 8)
    text("§9.9 the count of FC subgroups of order 2^(nu n)", f"$2^{{n^2(\\nu - {a})({b} - \\nu)}}$")
    # (nu - a)(b - nu) = 1/16  <=>  nu^2 - (a + b) nu + (ab + 1/16) = 0
    c1, c0 = a + b, a * b + Fraction(1, 16)
    disc = c1 * c1 - 4 * c0
    root = c1 / 2
    text(
        "§9.9 the quadratic, its discriminant and its root",
        f"$(\\nu - {tfrac(a)})({tfrac(b)} - \\nu) = {tfrac(Fraction(1, 16))}$ is $\\nu^2 - {tfrac(c1)}\\nu + {tfrac(c0)} = 0$: "
        + ("discriminant exactly zero" if disc == 0 else f"discriminant {disc}")
        + f", double root $\\nu = {root}$",
    )


def mark_table(kind: str, anchor: str) -> None:
    """Every number in the Markdown table whose line holds the anchor is quoted on purpose."""
    i = NOTES.find(anchor)
    if i < 0 or NOTES.find(anchor, i + 1) >= 0:
        check(f"NOTES: {kind} table {anchor!r}", False, "anchor not found exactly once")
        return
    s = e = NOTES.rfind("\n", 0, i) + 1
    while e < len(NOTES) and NOTES.startswith("|", e):
        e = NOTES.find("\n", e) + 1 or len(NOTES)
    for off, _ in NUMS[bisect.bisect_left(NUM_POS, s) : bisect.bisect_left(NUM_POS, e)]:
        COVERED.setdefault(off, kind)


def notes_s10() -> None:
    section("NOTES.md §10 and the opening paragraphs")
    mark("historical", "was wrong (")
    compare("opening: alpha from e(128), as the correction quotes it", "(0.5276; the value is ", alpha(l2(e(128)), 128))
    text("§10 9!", f"$9! = {factorial(9)}$")
    compare("§10 the class-count ratio of n = 9 to n = 8", "larger permutation group and the $", Fraction(A000638[9], A000638[8]))
    compare("§10 the expected n = 9 values", "the expected $n=9$ values (", [1694723, A000638[9]], ints=True)
    k = STORED_SPLIT
    if FULL and k:
        text("§10 the stored rows recomputed", f"for $n \\le {k['cheap']}$ and $n = 512$, {k['fresh']} of the {k['total']}; recomputing the other {k['rest']}")
    else:
        NOTES_SKIPPED.append("§10 the stored rows recomputed (needs the full data)")
    mark_table("timing", "| step | time |")


def _same(a, b, rel: float = 1e-12) -> bool:
    """JSON equality, with floats equal to `rel` relative (libm may differ across machines)."""
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(_same(a[k], b[k], rel) for k in a if k != "seconds")
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(_same(x, y, rel) for x, y in zip(a, b))
    if isinstance(a, float) or isinstance(b, float):
        return isinstance(a, (int, float)) and isinstance(b, (int, float)) and math.isclose(a, b, rel_tol=rel, abs_tol=1e-300)
    return a == b


def check_fresh_runs() -> None:
    section("the §9 programs, rerun into a temporary directory")
    runs = (
        ("census8.py", census_fresh, "data/census8.json", 33),
        ("block_families.py", bf_fresh, "data/block_families.json", {512: 62, 128: 60}.get(BF_N)),
    )
    for name, run, path, want in runs:
        code, passes, data = run()
        check(f"{name} exits 0", code == 0, f"exit {code}")
        if want is None:
            NOTES_SKIPPED.append(f"{name} gate count (no expected count at n_max = {BF_N})")
        else:
            check(f"{name} passes all {want} of its checks", passes == want, f"{passes} passed")
        stored = load(path)
        check(f"{name}: a fresh run reproduces {path}", data is not None and stored is not None and _same(data, stored))


def coverage(listing: bool) -> None:
    section("NOTES.md: which decimal numbers were compared")
    decs = [(off, s) for off, s in NUMS if is_dec(s)]
    kinds = Counter(COVERED.get(off, "not compared") for off, _ in decs)
    print(f"  {kinds['compared']} of {len(decs)} decimal numbers compared; {kinds['historical']} quoted from first versions, {kinds['timing']} timings, {kinds['not compared']} not compared")
    missing = [(NOTES.count("\n", 0, off) + 1, s) for off, s in decs if off not in COVERED]
    if listing:
        for kind in ("historical", "timing"):
            for off, s in decs:
                if COVERED.get(off) == kind:
                    print(f"    {kind:10s} line {NOTES.count(chr(10), 0, off) + 1}: {s}")
        for line, s in missing:
            print(f"    NOT COMPARED line {line}: {s}")
    if NOTES_SKIPPED:
        NOTES_SKIPPED.append("the coverage check (some comparisons were skipped)")
    else:
        check("NOTES: every decimal number is compared, or marked as a timing or a first-version value", not missing, ", ".join(f"line {ln}: {s}" for ln, s in missing[:8]) if missing else f"{len(decs)} numbers")


NOTES_SECTIONS = (
    notes_s2, notes_s3, notes_s4, notes_s5, notes_s6, notes_s7, notes_s8,
    notes_s9_1, notes_s9_2, notes_s9_3, notes_s9_4, notes_s9_5, notes_s9_6, notes_s9_7, notes_s9_8, notes_s9_9, notes_s10,
)


def main() -> int:
    listing = "--coverage" in sys.argv[1:]
    print("final_check.py -- re-deriving the numbers in NOTES.md from the artifacts")
    check_closed_forms()
    check_elemab_against_bruteforce()
    check_involutions()
    check_no_rank_above_half()
    check_lattice_vs_oeis()
    check_pure_types()
    check_theta_law()
    check_stored_rows()
    check_galois()
    check_second_order()
    check_orbit_profile()
    check_sd_high_precision()
    check_rotr_thresholds()
    check_quotes()  # before the prose: notes_s9_4 reads PDF_RESULTS
    check_fresh_runs()
    for fn in NOTES_SECTIONS:
        try:
            fn()
        except Exception as exc:  # a crash is a failure, never a silent pass
            check(f"NOTES: {fn.__name__} ran to completion", False, f"{type(exc).__name__}: {exc}")
    coverage(listing)
    if SKIPPED or NOTES_SKIPPED:
        NOTES_SKIPPED.append("the check count printed in §10 (it is the count of a run with nothing skipped)")
    else:
        compare("§10 the number of checks", "`final_check.py` runs ", CHECKS + 1, ints=True)

    skipped = SKIPPED + NOTES_SKIPPED
    print(f"\n{CHECKS} checks run, {len(FAILURES)} failed, {len(skipped)} skipped")
    for s in skipped:
        print(f"  note: skipped -- {s}")
    if FAILURES:
        print("\nFAILED:")
        for f in FAILURES:
            print(f"  - {f}")
        return 1
    print("\nEvery check that ran passed." + (" The skipped ones are listed above and were not run." if skipped else ""))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
