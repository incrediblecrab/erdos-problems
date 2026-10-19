#!/usr/bin/env python3
"""Emit and check LRAT certificates for UNSAT wildcard-base scans.

This is a proof-checking path for the UNSAT records in results/safesat_k*.json.
It rebuilds each CNF with safesat.encode_cnf(), runs Homebrew CaDiCaL with LRAT
proof logging, checks the proof with drat-trim's lrat-check, and records the
result in results/lrat_check.json.  The script is resumable: accepted records
are skipped unless --force is passed.
"""

import argparse
import concurrent.futures
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))

from safesat import encode_cnf  # noqa: E402


DEFAULT_RESULT = ROOT / "results" / "lrat_check.json"
DEFAULT_WORK = ROOT / "refs" / "lrat-work"
DEFAULT_LRAT = ROOT / "refs" / "drat-trim" / "lrat-check"
DEFAULT_DRAT = ROOT / "refs" / "drat-trim" / "drat-trim"


def run(cmd, seconds=None):
    if seconds is not None:
        cmd = ["perl", "-e", "alarm shift; exec @ARGV", str(seconds)] + cmd
    start = time.monotonic()
    proc = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True)
    return proc, time.monotonic() - start


def tool_output(cmd):
    try:
        return subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, timeout=20)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return exc


def tool_versions(cadical, lrat_check, drat_trim):
    cad = tool_output([cadical, "--version"])
    lrat = tool_output([str(lrat_check)])
    drat = tool_output([str(drat_trim)])
    commit = subprocess.run(["git", "-C", str(lrat_check.parent), "rev-parse", "HEAD"],
                            text=True, capture_output=True)
    return {
        "cadical": (cad.stdout or cad.stderr).strip().splitlines()[0] if hasattr(cad, "stdout") else str(cad),
        "cadical_path": shutil.which(cadical) or cadical,
        "lrat_check": str(lrat_check),
        "lrat_check_usage": (lrat.stdout or lrat.stderr).strip().splitlines()[0] if hasattr(lrat, "stdout") else str(lrat),
        "drat_trim": str(drat_trim),
        "drat_trim_usage": (drat.stdout or drat.stderr).strip().splitlines()[0] if hasattr(drat, "stdout") else str(drat),
        "drat_trim_commit": commit.stdout.strip() if commit.returncode == 0 else None,
    }


def load_json(path, default):
    if not path.exists():
        return default
    with path.open() as fh:
        return json.load(fh)


def instance_key(k, b, f, symbreak):
    return f"{k}:{b}:{f}:{symbreak or 0}"


def write_cnf(cnf, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write(f"p cnf {cnf.nv} {len(cnf.clauses)}\n")
        for clause in cnf.clauses:
            fh.write(" ".join(map(str, clause)) + " 0\n")


def proof_size(path):
    if not path.exists():
        return {"bytes": 0, "lines": 0}
    lines = 0
    with path.open("rb") as fh:
        for _ in fh:
            lines += 1
    return {"bytes": path.stat().st_size, "lines": lines}


def check_proof(cnf_path, proof_path, lrat_check, timeout):
    proc, seconds = run([str(lrat_check), str(cnf_path), str(proof_path)], timeout)
    verdict = proc.returncode == 0 and "VERIFIED" in (proc.stdout + proc.stderr)
    return {
        "accepted": verdict,
        "exit_code": proc.returncode,
        "seconds": round(seconds, 3),
        "stdout_tail": proc.stdout[-500:],
        "stderr_tail": proc.stderr[-500:],
    }


def prove_one(case, args):
    k, b, f, symbreak = case
    key = instance_key(k, b, f, symbreak)
    work = args.work_dir / key.replace(":", "_")
    work.mkdir(parents=True, exist_ok=True)
    cnf_path = work / "case.cnf"
    proof_path = work / "case.lrat"
    cnf, _ = encode_cnf(k, b, f, symbreak=symbreak)
    if cnf is None:
        return {"key": key, "k": k, "b": b, "f": f, "symbreak": symbreak, "status": "FORCED"}
    write_cnf(cnf, cnf_path)
    cmd = [args.cadical, "-q", "--plain", "--no-binary", "--lrat", str(cnf_path), str(proof_path)]
    solver, solve_seconds = run(cmd, args.timeout)
    rec = {
        "key": key,
        "k": k,
        "b": b,
        "f": f,
        "symbreak": symbreak,
        "variables": cnf.nv,
        "clauses": len(cnf.clauses),
        "solver": "cadical",
        "solver_exit": solver.returncode,
        "solver_seconds": round(solve_seconds, 3),
        "status": "timeout" if solver.returncode in (142, -14) else "solver_error",
        "proof": proof_size(proof_path),
    }
    if solver.returncode == 20:
        verdict = check_proof(cnf_path, proof_path, args.lrat_check, args.timeout)
        rec.update({
            "status": "accepted" if verdict["accepted"] else "rejected",
            "checker": "lrat-check",
            "checker_exit": verdict["exit_code"],
            "checker_seconds": verdict["seconds"],
            "checker_verdict": "accepted" if verdict["accepted"] else "rejected",
            "checker_stdout_tail": verdict["stdout_tail"],
            "checker_stderr_tail": verdict["stderr_tail"],
            "proof": proof_size(proof_path),
        })
    elif solver.returncode == 10:
        rec["status"] = "sat"
    rec["solver_stdout_tail"] = solver.stdout[-500:]
    rec["solver_stderr_tail"] = solver.stderr[-500:]
    if not args.keep_proofs:
        shutil.rmtree(work, ignore_errors=True)
    return rec


def corrupt_first_hint_or_clause(src, dst):
    lines = src.read_text().splitlines()
    for i, line in enumerate(lines):
        toks = line.split()
        if not toks or toks[0] == "d" or (len(toks) > 1 and toks[1] == "d"):
            continue
        try:
            zero = toks.index("0")
        except ValueError:
            continue
        if zero + 1 < len(toks) and toks[zero + 1] != "0":
            toks[zero + 1] = "999999999"
        elif zero > 1:
            toks[1] = str(-int(toks[1]))
        else:
            continue
        lines[i] = " ".join(toks)
        dst.write_text("\n".join(lines) + "\n")
        return
    raise RuntimeError("no proof line to corrupt")


def drop_first_clause(src, dst):
    lines = src.read_text().splitlines()
    out = []
    dropped = False
    for line in lines:
        if line.startswith("p cnf "):
            parts = line.split()
            parts[3] = str(int(parts[3]) - 1)
            out.append(" ".join(parts))
        elif not dropped and line and not line.startswith("c"):
            dropped = True
            continue
        else:
            out.append(line)
    dst.write_text("\n".join(out) + "\n")


def validate_path(args):
    work = args.work_dir / "validation"
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True, exist_ok=True)
    cnf, _ = encode_cnf(4, 12, 2)
    cnf_path = work / "valid.cnf"
    proof_path = work / "valid.lrat"
    write_cnf(cnf, cnf_path)
    solver, solve_seconds = run([args.cadical, "-q", "--plain", "--no-binary", "--lrat",
                                 str(cnf_path), str(proof_path)], args.timeout)
    correct = check_proof(cnf_path, proof_path, args.lrat_check, args.timeout)
    corrupt_path = work / "corrupt.lrat"
    corrupt_first_hint_or_clause(proof_path, corrupt_path)
    corrupt = check_proof(cnf_path, corrupt_path, args.lrat_check, args.timeout)
    dropped_path = work / "dropped.cnf"
    drop_first_clause(cnf_path, dropped_path)
    dropped = check_proof(dropped_path, proof_path, args.lrat_check, args.timeout)
    if not args.keep_proofs:
        shutil.rmtree(work, ignore_errors=True)
    return {
        "instance": {"k": 4, "b": 12, "f": 2, "symbreak": 0},
        "solver_exit": solver.returncode,
        "solver_seconds": round(solve_seconds, 3),
        "correct_proof_accepted": correct["accepted"],
        "corrupted_hint_or_clause_rejected": not corrupt["accepted"],
        "dropped_clause_rejected": not dropped["accepted"],
        "correct": correct,
        "corrupted_hint_or_clause": corrupt,
        "dropped_clause": dropped,
    }


def collect_cases(ks, k6_limit=None):
    cases = []
    for k in ks:
        data = load_json(ROOT / "results" / f"safesat_k{k}.json", {})
        local = []
        for rec in data.values():
            if rec.get("status") != "UNSAT":
                continue
            local.append((rec["k"], rec["b"], rec["f"], rec.get("symbreak", 0)))
        local.sort(key=lambda x: (x[1], x[2], x[3]))
        if k == 6 and k6_limit is not None:
            local = local[:k6_limit]
        cases.extend(local)
    return cases


def summarize(instances):
    out = {"accepted": 0, "rejected": 0, "timeout": 0, "solver_error": 0, "sat": 0, "other": 0}
    by_k = {}
    for rec in instances:
        status = rec.get("status", "other")
        out[status if status in out else "other"] += 1
        k = str(rec["k"])
        by_k.setdefault(k, {"accepted": 0, "rejected": 0, "timeout": 0, "solver_error": 0, "sat": 0, "other": 0})
        by_k[k][status if status in by_k[k] else "other"] += 1
    out["by_k"] = by_k
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ks", default="3,4,5", help="comma-separated k values")
    ap.add_argument("--k6-limit", type=int, default=None, help="only try the first N sorted k=6 UNSAT cases")
    ap.add_argument("--timeout", type=int, default=600, help="per solver/checker alarm seconds")
    ap.add_argument("--max-wall", type=int, default=0, help="stop launching new work after this many seconds")
    ap.add_argument("--jobs", type=int, default=1, help="concurrent solver processes, capped at 2")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--keep-proofs", action="store_true")
    ap.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    ap.add_argument("--work-dir", type=Path, default=DEFAULT_WORK)
    ap.add_argument("--cadical", default="cadical")
    ap.add_argument("--lrat-check", type=Path, default=DEFAULT_LRAT)
    ap.add_argument("--drat-trim", type=Path, default=DEFAULT_DRAT)
    args = ap.parse_args()
    args.jobs = max(1, min(2, args.jobs))
    args.work_dir.mkdir(parents=True, exist_ok=True)

    data = load_json(args.result, {"tool_versions": {}, "validation": {}, "instances": []})
    existing = {rec.get("key"): rec for rec in data.get("instances", [])
                if rec.get("status") == "accepted" and not args.force}
    ks = [int(x) for x in args.ks.split(",") if x]
    cases = [case for case in collect_cases(ks, args.k6_limit)
             if instance_key(*case) not in existing]

    data["tool_versions"] = tool_versions(args.cadical, args.lrat_check, args.drat_trim)
    data["validation"] = validate_path(args)
    data.setdefault("instances", [])
    old = {rec.get("key"): rec for rec in data["instances"]}
    start = time.monotonic()
    launched = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        pending = []
        for case in cases:
            if args.max_wall and time.monotonic() - start >= args.max_wall:
                break
            pending.append(pool.submit(prove_one, case, args))
            launched += 1
            if len(pending) >= args.jobs:
                done, pending = concurrent.futures.wait(pending, return_when=concurrent.futures.FIRST_COMPLETED)
                pending = list(pending)
                for fut in done:
                    rec = fut.result()
                    old[rec["key"]] = rec
                    data["instances"] = sorted(old.values(), key=lambda r: (r["k"], r["b"], r["f"], r.get("symbreak", 0)))
                    data["summary"] = summarize(data["instances"])
                    args.result.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
                    print(f"{rec['key']} {rec['status']} solver={rec.get('solver_seconds')} checker={rec.get('checker_seconds')}", flush=True)
        for fut in concurrent.futures.as_completed(pending):
            rec = fut.result()
            old[rec["key"]] = rec
            data["instances"] = sorted(old.values(), key=lambda r: (r["k"], r["b"], r["f"], r.get("symbreak", 0)))
            data["summary"] = summarize(data["instances"])
            args.result.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
            print(f"{rec['key']} {rec['status']} solver={rec.get('solver_seconds')} checker={rec.get('checker_seconds')}", flush=True)

    data["instances"] = sorted(old.values(), key=lambda r: (r["k"], r["b"], r["f"], r.get("symbreak", 0)))
    data["summary"] = summarize(data["instances"])
    data["last_run"] = {
        "ks": ks,
        "k6_limit": args.k6_limit,
        "timeout": args.timeout,
        "jobs": args.jobs,
        "launched": launched,
        "seconds": round(time.monotonic() - start, 3),
    }
    args.result.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    bad_validation = not (data["validation"]["correct_proof_accepted"]
                          and data["validation"]["corrupted_hint_or_clause_rejected"]
                          and data["validation"]["dropped_clause_rejected"])
    sys.exit(1 if bad_validation else 0)


if __name__ == "__main__":
    main()
