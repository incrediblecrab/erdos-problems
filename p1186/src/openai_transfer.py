"""Check cyclic-coloring transfers, including the direct periodic route that needs no wildcard coset."""

import argparse
from fractions import Fraction
import json
from pathlib import Path

from cosetprobe import divisors, safe_coset
from exact import psi_zm


SOURCE = "https://github.com/openai/math/blob/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/Quantitative-Superexponential-Bounds-for-van-der-Waerden-Numbers-September-23-2026/build/sections/07-transfers.tex"
CYCLIC_SOURCE = SOURCE.rsplit("/", 1)[0] + "/06-perturbation.tex"
CYCLIC_EXPONENT = Fraction(1, 100000)


def monochromatic(word, k, cyclic=False):
    n = len(word)
    if k < 2 or not word:
        raise ValueError("A nonempty word and k >= 2 are required")
    count = tested = 0
    for a in range(n):
        differences = range(1, n) if cyclic else range(1, (n-1-a)//(k-1) + 1)
        for d in differences:
            tested += 1
            count += len({word[(a+j*d) % n] for j in range(k)}) == 1
    return {"progressions_tested": tested, "monochromatic": count}


def digit_product(base, digits):
    if not base or digits < 1:
        raise ValueError("A nonempty base and positive digit count are required")
    n = len(base)
    return [tuple(base[(value // n**place) % n] for place in range(digits))
            for value in range(n**digits)]


def periodic_count(n, k, modulus, include_degenerate=True):
    """Count interval APs when the base is cyclic AP-free; this hypothesis is checked by the caller."""
    if type(n) is not int or n < 0 or type(k) is not int or k < 2 or type(modulus) is not int or modulus < 1 or type(include_degenerate) is not bool:
        raise ValueError("Require integer n >= 0, k >= 2, modulus >= 1, and Boolean include_degenerate")
    spacing = (k-1)*modulus
    quotient = n//spacing
    return n*(quotient + int(include_degenerate)) - spacing*quotient*(quotient+1)//2


def periodic_probe(base, k, plant=None):
    if not base or any(type(color) not in (int, bool) or color not in (0, 1) for color in base):
        raise ValueError("Require a nonempty binary coloring")
    if monochromatic(base, k, cyclic=True)["monochromatic"]:
        raise ValueError("The periodic counting formula requires a cyclic AP-free base")
    modulus = len(base)
    spacing = (k-1)*modulus
    density = Fraction(1, (1 if plant == "normalization" else 2)*spacing)
    failures = []
    limits = sorted(set(range(2*spacing+2)) | {10*spacing-1, 10*spacing, 10*spacing+1})
    for n in limits:
        word = [base[i % modulus] for i in range(n)]
        increasing = monochromatic(word, k)["monochromatic"] if n else 0
        measured = increasing + n
        predicted = periodic_count(n, k, modulus, include_degenerate=plant != "degenerate")
        remainder = n % spacing
        error = Fraction(remainder*(spacing-remainder), 2*spacing)
        if measured != predicted or predicted != density*n*n + Fraction(n, 2) + error:
            failures.append({"n": n, "measured": measured, "formula": predicted})
        if not 0 <= error <= Fraction(spacing, 8):
            failures.append({"n": n, "error": str(error)})
        if periodic_count(n, k, modulus, include_degenerate=False) != increasing:
            failures.append({"n": n, "increasing_count": increasing})
    return {
        "word": base, "k": k, "modulus": modulus, "spacing": spacing,
        "limiting_density": str(density), "existing_evaluator_density": str(psi_zm(base, k)),
        "tested_interval_lengths": limits,
        "formula": "n*(floor(n/D)+1) - D*floor(n/D)*(floor(n/D)+1)/2; D=(k-1)*m",
        "error_identity": "M(n)/n^2 = 1/(2D) + 1/(2n) + r*(D-r)/(2D*n^2), r=n mod D, n>0",
        "counts_include_zero_difference": True,
        "failures": failures,
        "passed": not failures and density == psi_zm(base, k),
    }


def check(plant=None):
    base, k = [0, 0, 1, 1], 3
    if plant == "coloring":
        base[2] = 0
    cyclic = monochromatic(base, k, cyclic=True)
    cosets = []
    for size in divisors(len(base)):
        if size == len(base):
            continue
        step = len(base)//size
        for residue in range(step):
            free = {residue + step*j for j in range(size)}
            safe, dangerous = safe_coset(base, k, free)
            if plant == "interface":
                safe = True
            cosets.append({"free_coset": sorted(free), "safe": safe, "dangerous_progressions": dangerous})

    products = []
    for digits in (1, 2, 3):
        word = digit_product(base, digits)
        products.append({"digits": digits, "length": len(word), "colors": len(set(word)),
                         **monochromatic(word, k)})

    qr = [int(pow(x, 5, 11) == 1) for x in range(11)]
    qr_cyclic = monochromatic(qr, 4, cyclic=True)
    qr_safe, qr_dangerous = safe_coset(qr, 4, {0})
    periodic = []
    if cyclic["monochromatic"] == 0:
        periodic = [periodic_probe(base, k, plant), periodic_probe(qr, 4, plant)]
    checks = [
        {"name": "small base is cyclic AP-free", "passed": cyclic["monochromatic"] == 0},
        {"name": "digit products are interval AP-free", "passed": all(item["monochromatic"] == 0 for item in products)},
        {"name": "cyclic AP-freeness does not imply a safe proper wildcard coset",
         "passed": cyclic["monochromatic"] == 0 and not any(item["safe"] for item in cosets)},
        {"name": "known quadratic-residue wildcard base still passes",
         "passed": qr_cyclic["monochromatic"] == 0 and qr_safe},
        {"name": "direct periodic counts and normalization need no safe coset",
         "passed": len(periodic) == 2 and all(row["passed"] for row in periodic)},
    ]
    return {
        "source": SOURCE,
        "scope": "Exact finite checks of the periodic and digit-product interfaces. The separate periodic-lean.json records the formal bridge. The superexponential consequence is conditional on the cited OpenAI cyclic-coloring theorem, which was not independently verified here. No small-k record improved.",
        "plant": plant,
        "base": base,
        "k": k,
        "cyclic": cyclic,
        "proper_cosets": cosets,
        "digit_products": products,
        "positive_control": {"base": qr, "k": 4, "free_coset": [0], "cyclic": qr_cyclic,
                             "safe": qr_safe, "dangerous_progressions": qr_dangerous},
        "direct_periodic_transfer": periodic,
        "conditional_consequence": {
            "source": CYCLIC_SOURCE,
            "source_statement_label": "perturb:cyclic",
            "source_exponent": str(CYCLIC_EXPONENT),
            "premise": "For every sufficiently large k, a cyclic AP-free two-coloring exists with modulus m >= k^(c*k).",
            "conclusion": "delta_k <= 1/(2*(k-1)*k^(c*k)); hence delta_k^(1/k) tends to 0.",
            "status": "Conditional on the imported cyclic-coloring theorem; not inferred from the interval-only Comparator challenge.",
            "effective_small_k_threshold_obtained": False,
        },
        "checks": checks,
        "passed": all(item["passed"] for item in checks),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--plant", choices=("coloring", "interface", "degenerate", "normalization"))
    args = parser.parse_args()
    report = check(args.plant)
    for item in report["checks"]:
        print(f"{'PASS' if item['passed'] else 'FAIL'} {item['name']}")
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
