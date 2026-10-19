#!/usr/bin/env python
"""Driver for the block-colouring search (Erdos #1186).

For each modulus m in a ladder this

  1. exports the exact integer weight table for (m, k, setting),
  2. searches it with src/search -- exhaustive when m is small enough, otherwise
     restarted annealing warm-started from the blow-ups of the best words already
     found at divisors of m (a blow-up denotes the same subset of [0,1) and so
     has exactly the same objective value, and a warm start can never lose),
  3. re-certifies the winning word in exact rational arithmetic through
     src/exact.py, which shares no code with the C searcher,
  4. keeps the certified result if it improves on what is already stored.

Results accumulate in results/<setting>_k<k>.json and are reused across runs, so
a later run can only improve on an earlier one.

Settings:  zp -> minimise Phi_k, bound is delta~_k <= Phi_k/2
           n  -> minimise Psi_k, bound is delta_k  <= Psi_k

Usage:  run_search.py <k> <m-list-or-range> [--setting zp|n] [--exhaustive-max M]
                      [--restarts R] [--sweeps S] [--threads T] [--rounds N]
"""
import argparse
import json
import subprocess
import sys
import time
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from exact import phi, psi, block_weights_int, triangle_weights_int  # noqa: E402
from export_weights import export  # noqa: E402

SEARCH = HERE / "search"
WDIR = ROOT / "results" / "weights"


def objective(c, k, setting):
    return phi(c, k) if setting == "zp" else psi(c, k)


def nsets(m, k, setting):
    t = block_weights_int if setting == "zp" else triangle_weights_int
    return len(t(m, k)[1])


def divisors(m):
    return [d for d in range(2, m) if m % d == 0]


def blowup(word, m):
    """Blow a word of length d up to length m = r*d; the subset of [0,1) and
    therefore the objective value are unchanged."""
    return "".join(ch * (m // len(word)) for ch in word)


def store_path(k, setting):
    return ROOT / "results" / f"{setting}_k{k}.json"


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"{' '.join(map(str,cmd))} failed rc={r.returncode}: {r.stderr[:400]}")
    return json.loads(r.stdout.strip())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("k", type=int)
    ap.add_argument("ms", help="comma list, or A-B for an inclusive range")
    ap.add_argument("--setting", choices=("zp", "n"), default="zp")
    ap.add_argument("--exhaustive-max", type=int, default=26)
    ap.add_argument("--restarts", type=int, default=12)
    ap.add_argument("--sweeps", type=int, default=3000)
    ap.add_argument("--threads", type=int, default=11)
    ap.add_argument("--rounds", type=int, default=1)
    a = ap.parse_args()

    if "-" in a.ms and "," not in a.ms:
        lo, hi = a.ms.split("-")
        ms = list(range(int(lo), int(hi) + 1))
    else:
        ms = [int(x) for x in a.ms.split(",")]

    WDIR.mkdir(parents=True, exist_ok=True)
    sp = store_path(a.k, a.setting)
    store = json.loads(sp.read_text()) if sp.exists() else {}

    for m in ms:
        t0 = time.time()
        wfile = WDIR / f"{a.setting}_{m}_{a.k}.txt"
        if not wfile.exists():
            export(m, a.k, wfile, a.setting)
        cands = []

        if m <= a.exhaustive_max:
            cands.append(("exhaustive",
                          run([str(SEARCH), str(wfile), "exhaustive", str(a.threads)])["word"]))
        else:
            seeds, seen = [], set()
            for d in divisors(m):
                rec = store.get(str(d))
                if rec:
                    s = blowup(rec["word"], m)
                    if s not in seen:
                        seen.add(s)
                        seeds.append(s)
            prev = store.get(str(m))
            if prev and prev["word"] not in seen:
                seeds.append(prev["word"])
            for rd in range(a.rounds):
                cands.append(("anneal",
                              run([str(SEARCH), str(wfile), "anneal", str(a.restarts),
                                   str(a.sweeps), str(1000 + 7919 * rd), str(a.threads)])["word"]))
                for si, s in enumerate(seeds):
                    cands.append(("seed",
                                  run([str(SEARCH), str(wfile), "seed", s,
                                       str(max(2, a.restarts // 3)), str(a.sweeps),
                                       str(50021 + 131 * si + 7 * rd), str(a.threads)])["word"]))
            cands.extend(("blowup", s) for s in seeds)

        # exact certification of every candidate; the C output is never trusted
        best = None
        for how, w in cands:
            v = objective([int(ch) for ch in w], a.k, a.setting)
            if best is None or v < best[0]:
                best = (v, w, how)

        prev = store.get(str(m))
        if prev is None or best[0] < Fraction(prev["value"]):
            bound = best[0] / 2 if a.setting == "zp" else best[0]
            store[str(m)] = {
                "m": m, "k": a.k, "setting": a.setting, "word": best[1],
                "value": str(best[0]), "value_float": float(best[0]),
                "bound": str(bound), "bound_float": float(bound),
                "found_by": best[2], "exhaustive": bool(m <= a.exhaustive_max),
            }
            sp.write_text(json.dumps(store, indent=1, sort_keys=True))
        rec = store[str(m)]
        print(f"m={m:4d} k={a.k} {a.setting} sets={nsets(m, a.k, a.setting):7d} "
              f"val={rec['value']:>24s} = {rec['value_float']:.10f}  "
              f"bound={rec['bound_float']:.10f}  "
              f"[{rec['found_by']}{'/exh' if rec['exhaustive'] else ''}] "
              f"{time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
