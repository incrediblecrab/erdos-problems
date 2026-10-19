#!/usr/bin/env python3
"""Re-derive every claim this project makes, from the artifacts on disk.

Exit code 0 means every check passed.  Nothing here trusts a stored number: the
words in results/*.json are re-evaluated from scratch in exact rational
arithmetic, and the headline records are additionally re-measured by brute-force
C programs (src/verify_zp, src/verify_n) that share no code with src/exact.py.

Checks, in order:

  A  evaluator self-consistency  -- carry-region areas, all-ones normalisation,
                                    and the Fourier identities for k = 3, 4
                                    checked exhaustively at small m
  B  literature reproduction     -- Lu-Peng's B_20 must give Phi_4 = 17/150 and
                                    PRS08's 12-block word must give Psi_3 =
                                    117/2192, both exactly
  C  certification               -- every word in results/*.json re-certifies to
                                    its stored value
  D  independent cross-check     -- the record words re-measured by brute force
                                    over Z_p and over {1,...,n}
  E  record table                -- the best certified bound per (problem, k),
                                    compared against the published values
  F  second evaluator            -- src/blocks.py computes Psi_k by a completely
                                    different decomposition (continuous block
                                    boundaries, its own polygon clipper) and
                                    must agree with src/exact.py
  G  AP-free certificates        -- every SAT model in results/apfree_*.json is
                                    re-checked by a plain loop
  H  wildcard recursion bases    -- every base in results/safebase_*.json has
                                    its coset re-tested by cosetprobe.safe_coset
                                    and its word confirmed cyclic AP-free, and
                                    the recursion identity is re-measured on
                                    explicit lifts by src/verify_zm
  J  OpenAI transfer             -- exact periodic counts, stored transfer data, and the core-only Lean bridge; not OpenAI's upstream theorem

Usage:  python src/final_check.py [--fast] [--no-brute]
"""

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
from fractions import Fraction
from itertools import product

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import exact  # noqa: E402
import blocks  # noqa: E402  -- second, independent evaluator (check F)

try:
    import numpy as np
except ImportError:                                    # brute force needs it
    np = None

# ---------------------------------------------------------------- literature

# [LuPe12] section 3: "B_20 = (1,1,1,0,1,1,0,1,1,1,0,0,0,1,0,0,1,0,0,0)".
LUPE_B20 = "11101101110001001000"
LUPE_PHI4 = Fraction(17, 150)

# [PRS08] Theorem 5: 0^28 1^6 0^28 1^37 0^59 1^116 0^116 1^59 0^37 1^28 0^6 1^28,
# block lengths in units of n/548.
PRS08_LENS = [28, 6, 28, 37, 59, 116, 116, 59, 37, 28, 6, 28]
PRS08_WORD = "".join(str(j % 2) * L for j, L in enumerate(PRS08_LENS))
PRS08_PSI3 = Fraction(117, 2192)

# Published values, in this project's normalisation (delta_k / delta~_k).
# delta~_k = Phi_k / 2 because Wolf's and Lu-Peng's m_k is twice the unoriented
# density; delta_k = Psi_k with no rescaling, matching PRS08's V(n)/n^2.
PUBLISHED = {
    ("zp", 3): (Fraction(1, 8), "= 1/8 exactly [Wo10]"),
    ("zp", 4): (Fraction(17, 300), "<= 17/300 [LuPe12 Thm 1, prime p]"),
    ("zp", 5): (Fraction(3629, 131424), "<= 3629/131424 [LuPe12 Thm 4, odd n]"),
    ("n", 3): (Fraction(117, 2192), "<= 117/2192 [PRS08 Thm 5]"),
    ("n", 4): (Fraction(1, 72), "<= 1/72 [LuPe12]"),
    ("n", 5): (Fraction(1, 304), "<= 1/304 [LuPe12]"),
}
RANDOM = {
    "zp": lambda k: Fraction(1, 2 ** k),
    "n": lambda k: Fraction(1, (k - 1) * 2 ** k),
}

# Which problem each search family bounds.  zp_* bounds delta~_k over F_p;
# the other three all bound delta_k over {1,...,n} and compete with each other.
FAMILY = {"zp": "zp", "n": "n", "zm": "n", "two": "n", "base": "n"}
FAMILY_NAME = {
    "zp": "block colouring of F_p",
    "n": "block colouring of {1..n}",
    "zm": "periodic colouring mod m",
    "two": "two-scale colouring",
    "base": "wildcard chain limit",
}

# The claims actually made in README.md and NOTES.md.  Every one is asserted
# below against the independently re-certified records, so a claim that drifts
# ahead of the artifacts fails the check rather than being quietly published.
# "beats" is the published value the claim is asserted to strictly improve on;
# None means no competing published value was found.
CLAIMS = {}


def load_claims():
    """Claims are stored next to the prose that makes them, so that editing
    README.md without re-running a search cannot silently change a claim."""
    path = os.path.join(ROOT, "results", "claims.json")
    if not os.path.exists(path):
        return {}
    out = {}
    for rec in json.load(open(path)):
        out[(rec["problem"], rec["k"])] = {
            "bound": Fraction(rec["bound"]),
            "beats": Fraction(rec["beats"]) if rec.get("beats") else None,
            "beats_random": rec.get("beats_random", False),
        }
    return out


class Checker:
    def __init__(self):
        self.failed = 0
        self.passed = 0
        self.skipped = 0

    def check(self, ok, label, detail=""):
        self.passed += bool(ok)
        self.failed += not ok
        mark = "PASS" if ok else "FAIL"
        print(f"  [{mark}] {label}" + (f"   {detail}" if detail else ""), flush=True)
        return ok

    def skip(self, label, detail=""):
        """A check that cannot run here. It is neither a pass nor a failure, and the summary says so."""
        self.skipped += 1
        print(f"  [SKIP] {label}" + (f"   {detail}" if detail else ""), flush=True)


def bits(word):
    return [1 if ch == "1" else 0 for ch in word]


# ------------------------------------------------------------------ check A

def check_evaluator(ck, fast):
    print("\nA. evaluator self-consistency")
    kmax = 6 if fast else 7
    for k in range(2, kmax + 1):
        regs = exact.carry_regions(k)
        tot = sum(a for _, a in regs)
        ck.check(tot == 1, f"carry-region areas sum to 1 (k={k})",
                 f"{len(regs)} cells")
    for k in range(3, kmax + 1):
        m = 5
        ones = [1] * m
        ck.check(exact.phi(ones, k) == 1, f"all-ones gives Phi_{k} = 1")
        ck.check(exact.psi(ones, k) == Fraction(1, 2 * (k - 1)),
                 f"all-ones gives Psi_{k} = 1/(2(k-1))")

    # Phi_3 = (1 + 3 mu^2)/4 and Phi_4 = (1 + 6 mu^2 + Lambda_4)/8.  These come
    # from expanding the monochromatic indicator into even-subset correlations;
    # every pair (x+ad, x+bd) with a != b is uniform on the torus squared, so
    # every 2-subset contributes exactly mu^2.  Checking them exhaustively is
    # the strongest available test that the cell decomposition is right.
    mrange = range(2, 8) if fast else range(2, 10)
    bad3 = bad4 = 0
    for m in mrange:
        for c in product((0, 1), repeat=m):
            mu = exact.mean(c)
            if exact.phi(c, 3) != (1 + 3 * mu * mu) / 4:
                bad3 += 1
            if m <= 7:
                lam = exact.lambda_form(c, 4)
                if exact.phi(c, 4) != (1 + 6 * mu * mu + lam) / 8:
                    bad4 += 1
    ck.check(bad3 == 0, "identity Phi_3 = (1+3mu^2)/4 holds exhaustively",
             f"m in {mrange.start}..{mrange.stop-1}")
    ck.check(bad4 == 0, "identity Phi_4 = (1+6mu^2+Lambda_4)/8 holds exhaustively",
             "m in 2..7")

    # Scale invariance.  Two different operations, with two different scopes:
    #  * stretch (subdivide every block r ways) denotes the *same* subset of
    #    [0,1), so it must preserve both Phi and Psi exactly;
    #  * repeat (tile the word r times) is the pullback under x -> r*x, which is
    #    a measure-preserving cover of the (x,d) torus and so preserves Phi, but
    #    does *not* preserve Psi, because x -> r*x does not map the triangle
    #    x + (k-1)t <= 1 to itself.
    bad_s = bad_r = 0
    for word in ("101", "1100", "10110", LUPE_B20):
        c = bits(word)
        for r in (2, 3):
            cs = bits("".join(ch * r for ch in word))
            cr = bits(word * r)
            for k in (3, 4, 5):
                bad_s += exact.phi(c, k) != exact.phi(cs, k)
                bad_s += exact.psi(c, k) != exact.psi(cs, k)
                bad_r += exact.phi(c, k) != exact.phi(cr, k)
    ck.check(bad_s == 0, "stretch invariance: Phi and Psi unchanged by subdivision")
    ck.check(bad_r == 0, "repeat invariance: Phi unchanged by tiling (Psi is not)")


# ------------------------------------------------------------------ check B

def check_literature(ck):
    print("\nB. reproduction of published constructions")
    got = exact.phi(bits(LUPE_B20), 4)
    ck.check(got == LUPE_PHI4,
             "Lu-Peng B_20 gives Phi_4 = 17/150 exactly",
             f"got {got}")
    ck.check(got / 2 == Fraction(17, 300),
             "  hence delta~_4 <= 17/300, the published bound")

    ck.check(len(PRS08_WORD) == 548 and sum(PRS08_LENS) == 548,
             "PRS08 block lengths sum to 548")
    got = exact.psi(bits(PRS08_WORD), 3)
    ck.check(got == PRS08_PSI3,
             "PRS08 12-block word gives Psi_3 = 117/2192 exactly",
             f"got {got}")

    # Lambda_4 = -7/75 is the content of 17/150 once the mean term is removed.
    lam = exact.lambda_form(bits(LUPE_B20), 4)
    ck.check(exact.mean(bits(LUPE_B20)) == 0 and lam == Fraction(-7, 75),
             "  B_20 is balanced with Lambda_4 = -7/75")


# ------------------------------------------------------------------ check C

def normalise(fam, k, key, entry):
    """One stored record -> a common shape carrying an upper bound on the constant."""
    if fam in ("zp", "n"):
        if "value" not in entry or "word" not in entry:
            return None
        val = Fraction(entry["value"])
        return {"fam": fam, "k": k, "key": key, "word": entry["word"],
                "value": val, "bound": val / 2 if fam == "zp" else val,
                "label": "m=%s" % key}
    if fam == "zm":
        if "psi" not in entry or "word" not in entry:
            return None
        val = Fraction(entry["psi"])
        return {"fam": fam, "k": k, "key": key, "word": entry["word"],
                "value": val, "bound": val, "m": entry["m"],
                "label": "m=%s%s" % (key, ", AP-free" if entry.get("ap_free") else "")}
    if fam == "two":
        if "psi" not in entry or "colour" not in entry:
            return None
        val = Fraction(entry["psi"])
        return {"fam": fam, "k": k, "key": key, "word": entry["colour"],
                "value": val, "bound": val, "m": entry["m"], "B": entry["B"],
                "label": "m=%d B=%d" % (entry["m"], entry["B"])}
    if fam == "base":
        if not all(x in entry for x in ("b", "f", "r", "word")):
            return None
        b, f = entry["b"], entry["f"]
        # the value is the recursion's limit 1/(b+f), proved in NOTES.md
        # section 3.1 from the two conditions recertify() checks; the stored
        # "bound" field is never read
        return {"fam": fam, "k": k, "key": key, "word": entry["word"],
                "b": b, "f": f, "r": entry["r"], "value": Fraction(1, b + f),
                "bound": Fraction(1, 2 * (k - 1) * (b + f)),
                "label": "b=%d f=%d" % (b, f)}
    return None


def base_coset(rec):
    b, f = rec["b"], rec["f"]
    return {(rec["r"] + (b // f) * y) % b for y in range(f)}


def base_violations(c, k, F):
    """Progressions of Z_b with nonzero difference whose terms outside F are
    all one colour: an H1 failure if no term is in F, an H2 failure otherwise.
    Either kind breaks the recursion, so the base is valid iff this is 0."""
    b = len(c)
    bad = 0
    for a in range(b):
        for d in range(1, b):
            out = {(a + j * d) % b for j in range(k)} - F
            if out and len({c[x] for x in out}) == 1:
                bad += 1
    return bad


def load_results():
    """Every certified record, keyed by (family, k), merged over search shards.

    Concurrent sweeps rewrite their whole json file on each save, so a single
    shard is not guaranteed to hold the best record for a given modulus; the
    merge here takes the minimum over all shards of a family.

        zp_*   block colourings over F_p       bound = Phi_k / 2  -> delta~_k
        n_*    block colourings of {1,...,n}   bound = Psi_k      -> delta_k
        zm_*   periodic colourings mod m       bound = Psi_k      -> delta_k
        two_*  two-scale colourings            bound = Psi_k      -> delta_k
    """
    out = {}
    for path in sorted(glob.glob(os.path.join(ROOT, "results", "*_k*.json"))):
        base = os.path.basename(path)[:-5]
        fam = None
        # wildcard_* is the restricted search inside the periodic family, so it
        # certifies the same quantity and is merged into zm.
        for pre, f in (("zp_", "zp"), ("n_", "n"), ("zm_", "zm"),
                       ("two_", "two"), ("wildcard_", "zm"),
                       ("safebase_", "base")):
            if base.startswith(pre):
                fam = f
                break
        if fam is None:
            continue
        try:
            k = int(base.split("_k")[1].split("_")[0])
        except (IndexError, ValueError):
            continue
        try:
            raw = json.load(open(path))
        except (json.JSONDecodeError, OSError):
            continue                                   # a sweep mid-write
        for key, entry in raw.items():
            rec = normalise(fam, k, key, entry)
            if rec is None:
                continue
            db = out.setdefault((fam, k), {})
            cur = db.get(key)
            if cur is None or rec["bound"] < cur["bound"]:
                db[key] = rec
    return out


def recertify(rec):
    """Recompute a record's value from its word, in exact rational arithmetic."""
    c = bits(rec["word"])
    fam, k = rec["fam"], rec["k"]
    if fam == "zp":
        if len(c) != int(rec["key"]):
            return None
        return exact.phi(c, k)
    if fam == "n":
        if len(c) != int(rec["key"]):
            return None
        return exact.psi(c, k)
    if fam == "zm":
        if len(c) != rec["m"]:
            return None
        return exact.psi_zm(c, k)
    if fam == "two":
        if len(c) != rec["m"] * rec["B"]:
            return None
        return exact.psi_two_scale(c, rec["m"], rec["B"], k)
    if fam == "base":
        b, f = rec["b"], rec["f"]
        if len(c) != b or not (1 <= f < b) or b % f:
            return None
        if base_violations(c, k, base_coset(rec)):
            return None
        return Fraction(1, b + f)
    return None


def check_certification(ck, results, fast):
    print("\nC. exact re-certification of every stored word")
    total = 0
    for (fam, k), db in sorted(results.items()):
        bad = []
        items = sorted(db.items(), key=lambda kv: kv[1]["bound"])
        if fast:
            items = items[:12]
        for key, rec in items:
            val = recertify(rec)
            if val is None:
                bad.append((key, "H1/H2 or shape" if fam == "base" else "shape"))
            elif val != rec["value"]:
                bad.append((key, f"{val} != {rec['value']}"))
            total += 1
        ck.check(not bad, f"{fam:<4} k={k}: {len(items)} words re-certify",
                 "" if not bad else f"{len(bad)} bad: {bad[:3]}")
    print(f"      ({total} words re-evaluated in exact rational arithmetic)")


# ------------------------------------------------------------------ check F

def check_second_evaluator(ck, results, fast):
    """src/blocks.py decomposes Psi_k over continuous block boundaries with its
    own polygon clipper, sharing no code with src/exact.py.  Agreement between
    the two is a real check on both."""
    print("\nF. second, independent evaluator (src/blocks.py)")
    lens = [28, 6, 28, 37, 59, 116, 116, 59, 37, 28, 6, 28]
    got = blocks.psi_blocks(blocks.bounds_of_lengths(lens), 3)
    ck.check(got == PRS08_PSI3,
             "blocks.py reproduces PRS08 from continuous block lengths",
             f"{got} vs {PRS08_PSI3}")
    for k in (3, 4, 5, 6):
        got = blocks.psi_blocks([Fraction(0), Fraction(1)], k)
        ck.check(got == Fraction(1, 2 * (k - 1)),
                 f"blocks.py: constant colouring gives the triangle area for k={k}",
                 f"{got}")
    bad, n = [], 0
    for (fam, k), db in sorted(results.items()):
        if fam != "n":
            continue
        items = sorted(db.items(), key=lambda kv: kv[1]["bound"])
        items = items[:6] if fast else items[:25]
        for key, rec in items:
            a = rec["value"]
            b = blocks.psi_blocks(blocks.bounds_of_word(rec["word"]), k)
            n += 1
            if a != b:
                bad.append((k, key, str(a), str(b)))
    ck.check(not bad, f"blocks.py agrees with exact.py on {n} stored block words",
             "" if not bad else str(bad[:3]))


# ------------------------------------------------------------------ check G

def check_apfree(ck):
    """Every SAT model claiming a k-AP-free colouring of Z_m is re-checked by a
    plain loop, so a solver bug cannot manufacture a bound."""
    print("\nG. AP-free certificates re-checked by direct enumeration")
    paths = sorted(glob.glob(os.path.join(ROOT, "results", "apfree_k*.json")))
    if not paths:
        print("      (none found)")
        return
    for path in paths:
        k = int(os.path.basename(path).split("_k")[1].split(".")[0])
        bad, good = [], 0
        for entry in json.load(open(path)).values():
            w = entry.get("word")
            if not w:
                continue
            m = entry["m"]
            c = [int(ch) for ch in w]
            viol = sum(1 for a in range(m) for b in range(1, m)
                       if len({c[(a + j * b) % m] for j in range(k)}) == 1)
            if viol or len(w) != m:
                bad.append((m, viol))
            else:
                good += 1
                cert = exact.psi_zm(c, k)
                if cert != Fraction(1, 2 * (k - 1) * m):
                    bad.append((m, "psi %s" % cert))
        ck.check(not bad, f"k={k}: {good} AP-free colourings verified "
                          f"(0 monochromatic k-APs, Psi = 1/(2(k-1)m))",
                 "" if not bad else str(bad[:3]))


# ------------------------------------------------------------------ check H

def zm_mono(word, k):
    """Monochromatic (a, d) pairs of Z_m, d = 0 included, counted by
    src/verify_zm -- a plain double loop sharing no code with this file."""
    exe = os.path.join(HERE, "verify_zm")
    if not os.path.exists(exe):
        return None
    out = subprocess.run([exe, str(k), "".join(map(str, word))],
                         capture_output=True, text=True, timeout=1800)
    for tok in out.stdout.split():
        if tok.startswith("mono="):
            return int(tok[5:])
    return None


def lift(c, f, r, h):
    """Colour Z_{bt} by the base off the lift of F and by h on it."""
    b = len(c)
    g, t = b // f, len(h) // f
    w = [c[x % b] for x in range(b * t)]
    for u in range(f * t):
        w[(r + g * u) % (b * t)] = h[u]
    return w


def check_bases(ck, results):
    """The recursion limit is a theorem (NOTES.md section 3.1) given H1 and
    H2, which check C re-verified.  Here H2 is re-tested by the older,
    independently written cosetprobe.safe_coset, the stored word (colours on
    F included) is confirmed cyclic k-AP-free by src/verify_zm, and the
    recursion identity

        #mono(Z_bt, lift) = (b - f) t^2 + #mono(Z_ft, h)

    is re-measured at t = 2 and t = 3 by src/verify_zm, for a constant and a
    non-constant inner colouring h."""
    from cosetprobe import safe_coset
    import random
    print("\nH. wildcard recursion bases: coset safety and the recursion identity")
    got = False
    for (fam, k), db in sorted(results.items()):
        if fam != "base":
            continue
        got = True
        for key, rec in sorted(db.items(), key=lambda kv: kv[1]["b"]):
            c, F = bits(rec["word"]), base_coset(rec)
            b, f, r = rec["b"], rec["f"], rec["r"]
            ok, nd = safe_coset(c, k, F)
            ck.check(ok, f"k={k} b={b} f={f} r={r}: cosetprobe.safe_coset accepts F",
                     f"{nd} dangerous progressions")
            mono = zm_mono(c, k)
            ck.check(mono == b, f"k={k} b={b} f={f}: stored word is cyclic {k}-AP-free",
                     f"{mono} monochromatic pairs; the {b} with d=0 are unavoidable")
            rng = random.Random(1186 * b + f)
            mixed = [0] * (3 * f)
            while len(set(mixed)) == 1:
                mixed = [rng.randint(0, 1) for _ in range(3 * f)]
            for t, name, h in ((2, "zeros", [0] * (2 * f)), (3, "random", mixed)):
                lhs = zm_mono(lift(c, f, r, h), k)
                inner = zm_mono(h, k)
                if lhs is None or inner is None:
                    ck.check(False, f"k={k} b={b}: verify_zm unavailable or failed")
                    continue
                want = (b - f) * t * t + inner
                ck.check(lhs == want,
                         f"k={k} b={b} f={f}: identity at t={t}, h={name}",
                         f"{lhs} == {b - f}*{t}^2 + {inner}")
    if not got:
        print("      (none found)")


# ------------------------------------------------------------------ check D


def brute(binary, args):
    exe = os.path.join(HERE, binary)
    if not os.path.exists(exe):
        return None
    try:
        out = subprocess.run([exe] + [str(a) for a in args],
                             capture_output=True, text=True, timeout=1800)
    except subprocess.TimeoutExpired:
        return None
    if out.returncode != 0:
        return None
    return json.loads(out.stdout.strip().splitlines()[-1])


def brute_count(word, k, n, m=None, B=None, periodic=False):
    """Count monochromatic k-APs in {1,...,n} directly, returning the density.

    This is an independent re-derivation: it colours each integer and counts
    progressions, sharing no code and no decomposition with exact.py.  For a
    colouring whose exact density is Psi, the count is Psi*n^2 + O(n), so the
    residual times n must stay bounded as n grows.
    """
    c = np.array([int(ch) for ch in word], dtype=np.int8)
    i = np.arange(1, n + 1)
    if periodic:
        col = c[i % len(word)]
    elif B is not None:
        col = c[(i % m) * B + np.minimum((B * i) // n, B - 1)]
    else:
        col = c[np.minimum((len(word) * i) // n, len(word) - 1)]
    tot = 0
    for d in range(0, (n - 1) // (k - 1) + 1):
        a = np.arange(0, n - (k - 1) * d)
        ok = np.ones(len(a), bool)
        v = col[a]
        for j in range(1, k):
            ok &= (col[a + j * d] == v)
        tot += int(ok.sum())
    return tot / float(n) ** 2


def check_brute(ck, results, p=60013, n=120000):
    print("\nD. independent brute force (no code shared with exact.py)")
    for (fam, k), db in sorted(results.items()):
        if fam == "base":
            continue                    # a limit, not one word: see check H
        key, rec = min(db.items(), key=lambda kv: kv[1]["bound"])
        want = float(rec["value"])
        if fam == "zp":
            res = brute("verify_zp", [p, rec["word"], k, 8])
            if res is None:
                ck.check(False, f"zp k={k}: verify_zp unavailable or failed")
                continue
            meas, tol = res["phi_measured"], 12.0 / p
            ck.check(abs(meas - want) < tol,
                     f"zp  k={k} {rec['label']}: Phi {meas:.9f} vs exact {want:.9f}",
                     f"|diff| = {abs(meas-want):.2e} < {tol:.2e}")
            continue
        if fam == "n":
            res = brute("verify_n", [n, rec["word"], k, 8])
            if res is not None:
                meas, tol = res["psi_measured"], 12.0 / n
                ck.check(abs(meas - want) < tol,
                         f"n   k={k} {rec['label']}: Psi {meas:.9f} vs exact {want:.9f}",
                         f"|diff| = {abs(meas-want):.2e} < {tol:.2e}")
                continue
        if np is None:
            ck.check(False, f"{fam} k={k}: numpy unavailable for brute force")
            continue
        kw = {"periodic": True} if fam == "zm" else (
            {"m": rec["m"], "B": rec["B"]} if fam == "two" else {})
        meas = brute_count(rec["word"], k, n, **kw)
        tol = 12.0 / n
        ck.check(abs(meas - want) < tol,
                 f"{fam:<4}k={k} {rec['label']}: Psi {meas:.9f} vs exact {want:.9f}",
                 f"|diff| = {abs(meas-want):.2e} < {tol:.2e}")


# ------------------------------------------------------------------ check E

def check_records(ck, results):
    print("\nE. certified records vs published values and vs random")
    print("  per family:")
    best = {}
    for (fam, k), db in sorted(results.items()):
        key, rec = min(db.items(), key=lambda kv: kv[1]["bound"])
        prob = FAMILY[fam]
        rnd = RANDOM[prob](k)
        print(f"    {FAMILY_NAME[fam]:<28} k={k}  {rec['label']:<14}"
              f" {str(rec['bound']):>20} = {float(rec['bound']):.10f}"
              f"  {float(rec['bound']/rnd):.4f}x random")
        cur = best.get((prob, k))
        if cur is None or rec["bound"] < cur["bound"]:
            best[(prob, k)] = rec

    print("\n  best certified bound per problem:")
    print(f"    {'':<7}{'k':<3}{'family':<22}{'where':<14}{'bound':>20}"
          f" {'':>13}  {'/random':>8}  status")
    for (prob, k), rec in sorted(best.items()):
        bound, rnd = rec["bound"], RANDOM[prob](k)
        name = "delta~_%d" % k if prob == "zp" else "delta_%d" % k
        if bound < rnd:
            status = "beats random"
        elif bound == rnd:
            status = "equals random"
        else:
            status = "WORSE THAN RANDOM (not a result)"
        pub = PUBLISHED.get((prob, k))
        if pub:
            status += "; beats" if bound < pub[0] else "; does NOT beat"
            status += f" published {pub[1]}"
        else:
            status += "; no published bound found"
        print(f"    {name:<7}{'':<3}{FAMILY_NAME[rec['fam']]:<22}{rec['label']:<14}"
              f"{str(bound):>20} {float(bound):>13.10f}"
              f"  {float(bound/rnd):>7.4f}x  {status}")

    claims = load_claims()
    if not claims:
        print("\n  (results/claims.json absent -- no claims to assert)")
        return
    print("\n  assertions on the claims made in README.md / NOTES.md:")
    for (prob, k), claim in sorted(claims.items()):
        rec = best.get((prob, k))
        name = "delta~_%d" % k if prob == "zp" else "delta_%d" % k
        if rec is None:
            ck.check(False, f"{name}: claimed {claim['bound']} but no record found")
            continue
        bound = rec["bound"]
        ck.check(bound <= claim["bound"],
                 f"{name}: certified {bound} <= claimed {claim['bound']}")
        if claim.get("beats") is not None:
            ck.check(bound < claim["beats"],
                     f"{name}: claim strictly improves on {claim['beats']}",
                     f"{float(bound):.10f} < {float(claim['beats']):.10f} "
                     f"({100*(1-float(bound/claim['beats'])):.2f}% better)")
        if claim.get("beats_random"):
            ck.check(bound < RANDOM[prob](k),
                     f"{name}: claim beats the random colouring",
                     f"{float(bound):.10f} < {float(RANDOM[prob](k)):.10f}")


def check_lrat(ck):
    """Cheap consistency check for the external LRAT proof-checking ledger.

    The expensive work is done by src/prove_unsat.py with CaDiCaL and
    drat-trim's lrat-check.  This gate only checks that the ledger matches the
    counts stated in the prose and that its validation tests have teeth."""
    print("\nI. LRAT-checked UNSAT ledger")
    path = os.path.join(ROOT, "results", "lrat_check.json")
    if not os.path.exists(path):
        print("      (results/lrat_check.json absent)")
        return
    data = json.load(open(path))
    validation = data.get("validation", {})
    ck.check(validation.get("correct_proof_accepted") is True,
             "lrat-check accepts the validation proof")
    ck.check(validation.get("corrupted_hint_or_clause_rejected") is True,
             "lrat-check rejects a proof with a corrupted hint/clause")
    ck.check(validation.get("dropped_clause_rejected") is True,
             "lrat-check rejects the same proof after dropping a CNF clause")
    expected = {
        "3": {"accepted": 89, "timeout": 0, "rejected": 0, "solver_error": 0},
        "4": {"accepted": 423, "timeout": 0, "rejected": 0, "solver_error": 0},
        "5": {"accepted": 538, "timeout": 0, "rejected": 0, "solver_error": 0},
        "6": {"accepted": 68, "timeout": 4, "rejected": 0, "solver_error": 0},
    }
    by_k = data.get("summary", {}).get("by_k", {})
    for k, want in expected.items():
        got = by_k.get(k, {})
        ok = all(got.get(field) == value for field, value in want.items())
        ck.check(ok, f"LRAT ledger k={k} counts match the prose",
                 f"got {got}, want {want}")


def check_openai_transfer(ck):
    from openai_transfer import check

    print("\nJ. cyclic-to-periodic transfer")
    current = check()
    with open(os.path.join(ROOT, "results", "openai-transfer.json")) as stream:
        stored = json.load(stream)
    ck.check(current == stored, "stored transfer artifact matches a fresh exact run")
    for item in current["checks"]:
        ck.check(item["passed"], item["name"])
    shared_gate = os.path.join(os.path.dirname(ROOT), "formal", "scripts", "verify_lean.py")
    if not os.path.exists(shared_gate):
        ck.skip("core-only Lean transfer replays and matches its recorded result",
                f"{shared_gate} is missing; the Lean check needs the sibling formal/ checkout, see README.md")
    elif shutil.which("elan") is None:
        ck.skip("core-only Lean transfer replays and matches its recorded result",
                "elan is not installed; the checker does not install toolchains")
    else:
        proof = subprocess.run(
            [sys.executable, "-B", os.path.join(HERE, "check_periodic_lean.py"), "--check"],
            capture_output=True, text=True, timeout=480,
        )
        ck.check(proof.returncode == 0, "core-only Lean transfer replays and matches its recorded result",
                 "" if proof.returncode == 0 else proof.stdout + proof.stderr)
    print("      The superexponential consequence remains conditional on the cited upstream cyclic theorem.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", action="store_true",
                    help="skip the deepest exhaustive identity checks")
    ap.add_argument("--no-brute", action="store_true",
                    help="skip the C brute-force cross-checks")
    args = ap.parse_args()

    ck = Checker()
    results = load_results()
    if not results:
        print("no results/*.json found -- nothing to check", file=sys.stderr)
        return 1
    check_evaluator(ck, args.fast)
    check_literature(ck)
    check_certification(ck, results, args.fast)
    check_second_evaluator(ck, results, args.fast)
    check_apfree(ck)
    check_bases(ck, results)
    if not args.no_brute:
        check_brute(ck, results)
    check_records(ck, results)
    check_lrat(ck)
    check_openai_transfer(ck)

    tail = f", {ck.skipped} skipped (incomplete)" if ck.skipped else ""
    print(f"\n{ck.passed} passed, {ck.failed} failed{tail}")
    return 1 if ck.failed else 0


if __name__ == "__main__":
    sys.exit(main())
