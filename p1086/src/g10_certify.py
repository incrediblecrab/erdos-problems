"""Certify the UNSAT answers behind the hypergraph lemma of g10_exact.py with LRAT proofs; writes results/g10_lrat.json.

    python src/g10_certify.py [--only NAME,...] [--timeout S]   prove and check every case, merging into the JSON
    python src/g10_certify.py --drat ...                          use the kissat + drat-trim route instead (prove_drat)
    python src/g10_certify.py --out PART.json ...                 write to PART.json instead (parallel runs); then
    python src/g10_certify.py --merge PART.json ...               folds the parts into results/g10_lrat.json
    python src/g10_certify.py --control                           checker control on one small case (see control())
    python src/g10_certify.py --verify [PATH]                     no solving: rebuild every CNF and require that its
                                                                  SHA-256 matches an accepted record; exit 0 iff so

The cases are the lemma's CNF as g10_exact.py builds it (python-sat encoding, with the symmetry reduction), named
"exact", certified by the --drat route because CaDiCaL with proof logging had not finished it after 17 CPU minutes, and the 21 cases of g10_check.py's independent encoding, named "check-codegree-4" to "-6" and
"check-link-00" to "-17". Each CNF is written in DIMACS with the clauses in the order the builder emits them and
hashed; CaDiCaL (Homebrew) runs with the arguments recorded in SOLVER_ARGS and writes its LRAT proof into a FIFO that
is relayed to drat-trim's lrat-check (built by `sh refs/fetch.sh drat-trim` at the pinned commit), so no proof is
stored on disk; the relay counts the proof's bytes and lines. A case is accepted when CaDiCaL exits 20
(UNSATISFIABLE), lrat-check exits 0 and prints "c VERIFIED". The solver is killed if free disk space in refs/ falls
below MIN_FREE_GIB or the run exceeds the timeout.

--verify is what final_check.py runs: the CNFs cost seconds to rebuild, so the gate confirms that the certified
formulas are exactly the ones the two scripts now generate, without repeating the proofs.
"""
import hashlib
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path

import g10_check
import g10_exact

ROOT = Path(__file__).resolve().parent.parent
RESULT = ROOT / "results" / "g10_lrat.json"
WORK = ROOT / "refs" / "lrat-work"
LRAT_CHECK = ROOT / "refs" / "drat-trim" / "lrat-check"
DRAT_TRIM = ROOT / "refs" / "drat-trim" / "drat-trim"
SOLVER_ARGS = ["-q", "--no-binary", "--lrat"]
KISSAT_ARGS = ["-q", "--no-binary"]
MIN_FREE_GIB = 20
DRAT_TRIM_COMMIT = "2e3b2dc0ecf938addbd779d42877b6ed69d9a985"


def cases():
    """(name, builder) for every case; builders return (number of variables, clauses)."""
    out = [("exact", lambda: _pysat(g10_exact.regular_formula(g10_exact.LEMMA_RULES)))]
    rules = ("a3", "b5", "b4")
    for k in (4, 5, 6):
        out.append((f"check-codegree-{k}", lambda k=k: _own(g10_check.build_case(k, rules))))
    for i, E in enumerate(g10_check.link_classes()):
        out.append((f"check-link-{i:02d}", lambda E=E: _own(g10_check.build_case(set(E), rules))))
    return out


def _pysat(f):
    return max(abs(l) for c in f.clauses for l in c), f.clauses


def _own(f):
    return f.nv, f.cls


def dimacs(nv, clauses):
    return (f"p cnf {nv} {len(clauses)}\n" + "".join(" ".join(map(str, c)) + " 0\n" for c in clauses)).encode()


def free_gib():
    return shutil.disk_usage(WORK).free / 2**30


def prove(name, blob, timeout):
    work = WORK / name
    work.mkdir(parents=True, exist_ok=True)
    cnf, fifo = work / "case.cnf", work / "proof.fifo"
    cnf.write_bytes(blob)
    if fifo.exists():
        fifo.unlink()
    os.mkfifo(fifo)
    checker = subprocess.Popen([str(LRAT_CHECK), str(cnf)], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    counts, outs = {"bytes": 0, "lines": 0}, {}

    def drain():
        outs["out"], outs["err"] = checker.stdout.read(), checker.stderr.read()

    def relay():
        try:
            with open(fifo, "rb", buffering=0) as src:
                while chunk := src.read(1 << 20):
                    counts["bytes"] += len(chunk)
                    counts["lines"] += chunk.count(b"\n")
                    try:
                        checker.stdin.write(chunk)
                    except BrokenPipeError:
                        pass  # the checker stopped on a bad step; keep draining so the solver is not blocked
        finally:
            try:
                checker.stdin.close()
            except BrokenPipeError:
                pass

    threads = [threading.Thread(target=drain, daemon=True), threading.Thread(target=relay, daemon=True)]
    for t in threads:
        t.start()
    start = time.monotonic()
    solver = subprocess.Popen([shutil.which("cadical") or "/opt/homebrew/bin/cadical", *SOLVER_ARGS, str(cnf),
                               str(fifo)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    abort = None
    while solver.poll() is None:
        time.sleep(2)
        if free_gib() < MIN_FREE_GIB:
            abort = "disk"
        elif time.monotonic() - start > timeout:
            abort = "timeout"
        if abort:
            solver.send_signal(signal.SIGTERM)
            break
    sout, serr = solver.communicate()
    solve_seconds = time.monotonic() - start
    if threads[1].is_alive():  # the solver never opened the proof: open the writer once so the relay sees EOF
        try:
            os.close(os.open(fifo, os.O_WRONLY | os.O_NONBLOCK))
        except OSError:
            pass
    threads[1].join()
    checker.wait()
    threads[0].join()
    total = time.monotonic() - start
    fifo.unlink()
    shutil.rmtree(work, ignore_errors=True)
    text = (outs.get("out", b"") + outs.get("err", b"")).decode(errors="replace")
    accepted = abort is None and solver.returncode == 20 and checker.returncode == 0 and \
        re.search(r"^c VERIFIED$", text, re.M) is not None
    return {
        "name": name,
        "cnf_sha256": hashlib.sha256(blob).hexdigest(),
        "cnf_bytes": len(blob),
        "solver": "cadical",
        "solver_args": " ".join(SOLVER_ARGS),
        "solver_exit": solver.returncode,
        "solver_seconds": round(solve_seconds, 1),
        "abort": abort,
        "proof_mode": "stream",
        "proof": counts,
        "checker": "lrat-check",
        "checker_exit": checker.returncode,
        "checker_seconds_after_start": round(total, 1),
        "checker_output_tail": text[-300:],
        "status": "accepted" if accepted else "rejected",
    }


def control(name="check-link-09"):
    """Negative control for the checker: the proof of one small case, written to a file, must be accepted, and the
    same proof with one added clause deleted, one that a later step uses as a hint, must be rejected."""
    blob = dimacs(*dict(cases())[name]())
    work = WORK / "control"
    work.mkdir(parents=True, exist_ok=True)
    cnf, proof = work / "case.cnf", work / "case.lrat"
    cnf.write_bytes(blob)
    solver = subprocess.run([shutil.which("cadical") or "/opt/homebrew/bin/cadical", *SOLVER_ARGS, str(cnf),
                             str(proof)], capture_output=True, text=True)
    lines = proof.read_bytes().splitlines(keepends=True)
    used = set()  # clause ids that appear as hints of some addition step
    for l in lines:
        t = l.split()
        if len(t) > 1 and t[1] != b"d":
            used.update(t[t.index(b"0") + 1:-1])
    adds = [i for i, l in enumerate(lines) if l.split()[1:2] != [b"d"] and l.split()[0] in used]
    cut = adds[len(adds) // 2]  # an addition some later step relies on
    verdicts = {}
    for label, body in (("intact", b"".join(lines)), ("one_addition_deleted", b"".join(lines[:cut] + lines[cut + 1:]))):
        proof.write_bytes(body)
        r = subprocess.run([str(LRAT_CHECK), str(cnf), str(proof)], capture_output=True, text=True)
        verdicts[label] = r.returncode == 0 and re.search(r"^c VERIFIED$", r.stdout + r.stderr, re.M) is not None
    shutil.rmtree(work, ignore_errors=True)
    rec = {"case": name, "solver_exit": solver.returncode, "proof_lines": len(lines), "deleted_line": cut,
           "intact_accepted": verdicts["intact"], "corrupted_accepted": verdicts["one_addition_deleted"]}
    print(f"control on {name}: intact proof accepted {rec['intact_accepted']}, proof with line {cut} of "
          f"{len(lines):,} deleted accepted {rec['corrupted_accepted']}")
    return rec


def guarded(cmd, timeout):
    """Run cmd, killing it if free disk space falls below MIN_FREE_GIB or it runs past timeout."""
    start = time.monotonic()
    proc = subprocess.Popen([str(c) for c in cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    abort = None
    while proc.poll() is None:
        time.sleep(1)
        if free_gib() < MIN_FREE_GIB:
            abort = "disk"
        elif time.monotonic() - start > timeout:
            abort = "timeout"
        if abort:
            proc.send_signal(signal.SIGTERM)
            break
    out, err = proc.communicate()
    return proc.returncode, out + err, round(time.monotonic() - start, 1), abort


def size(path):
    if not path.exists():
        return {"bytes": 0, "lines": 0}
    with path.open("rb") as fh:
        return {"bytes": path.stat().st_size, "lines": sum(1 for _ in fh)}


def prove_drat(name, blob, timeout):
    """Second route, for a CNF CaDiCaL is slow on: kissat writes a DRAT proof, drat-trim checks it and writes the
    LRAT form of the checked proof, and lrat-check checks that LRAT proof. All files are deleted afterwards."""
    work = WORK / name
    work.mkdir(parents=True, exist_ok=True)
    cnf, drat, lrat = work / "case.cnf", work / "case.drat", work / "case.lrat"
    cnf.write_bytes(blob)
    k_exit, _, k_s, k_abort = guarded([shutil.which("kissat") or "/opt/homebrew/bin/kissat", *KISSAT_ARGS, cnf, drat],
                                      timeout)
    d_exit, d_out, d_s, d_abort = guarded([DRAT_TRIM, cnf, drat, "-L", lrat, "-t", int(timeout)], timeout)
    l_exit, l_out, l_s, l_abort = guarded([LRAT_CHECK, cnf, lrat], timeout)
    rec = {
        "name": name, "cnf_sha256": hashlib.sha256(blob).hexdigest(), "cnf_bytes": len(blob),
        "route": "kissat DRAT, checked by drat-trim, whose LRAT output is checked by lrat-check",
        "solver": "kissat", "solver_args": " ".join(KISSAT_ARGS), "solver_exit": k_exit, "solver_seconds": k_s,
        "proof": {"format": "DRAT", **size(drat)}, "drat_trim_exit": d_exit, "drat_trim_seconds": d_s,
        "drat_trim_output_tail": d_out[-300:], "lrat_from_drat_trim": size(lrat),
        "checker": "lrat-check", "checker_exit": l_exit, "checker_seconds": l_s, "checker_output_tail": l_out[-300:],
        "abort": k_abort or d_abort or l_abort,
    }
    shutil.rmtree(work, ignore_errors=True)
    accepted = rec["abort"] is None and k_exit == 20 and re.search(r"^s VERIFIED", d_out, re.M) is not None and \
        l_exit == 0 and re.search(r"^c VERIFIED$", l_out, re.M) is not None
    rec["status"] = "accepted" if accepted else "rejected"
    return rec


def tools():
    cad = subprocess.run([shutil.which("cadical") or "/opt/homebrew/bin/cadical", "--version"], capture_output=True,
                         text=True).stdout.strip()
    commit = subprocess.run(["git", "-C", str(LRAT_CHECK.parent), "rev-parse", "HEAD"], capture_output=True,
                            text=True).stdout.strip()
    kis = subprocess.run([shutil.which("kissat") or "/opt/homebrew/bin/kissat", "--version"], capture_output=True,
                         text=True).stdout.strip()
    return {"cadical": cad, "kissat": kis, "drat_trim_commit": commit}


def verify(path):
    data = json.loads(Path(path).read_text())
    records = {r["name"]: r for r in data["records"]}
    ok = data["tools"]["drat_trim_commit"] == DRAT_TRIM_COMMIT
    names = []
    for name, build in cases():
        names.append(name)
        digest = hashlib.sha256(dimacs(*build())).hexdigest()
        r = records.get(name)
        good = r is not None and r["status"] == "accepted" and r["solver_exit"] == 20 and r["checker_exit"] == 0 \
            and "c VERIFIED" in r["checker_output_tail"] and r["cnf_sha256"] == digest
        print(f"{name}: {'accepted LRAT proof for this CNF' if good else 'NO accepted record for this CNF'}"
              + (f", proof {r['proof']['bytes']:,} bytes" if good else ""))
        ok &= good
    ok &= set(records) == set(names)
    c = data.get("control", {})
    good_control = c.get("solver_exit") == 20 and c.get("intact_accepted") is True and c.get("corrupted_accepted") is False
    print(f"checker control: intact proof accepted and corrupted proof rejected: {good_control}")
    ok &= good_control
    print(f"{len(names)} cases, {sum(1 for n in names if n in records)} records; "
          f"{'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def main(argv):
    if "--merge" in argv:  # fold the records of parallel runs (--out files) into results/g10_lrat.json
        data = json.loads(RESULT.read_text()) if RESULT.exists() else {"records": []}
        records = {r["name"]: r for r in data["records"]}
        for path in argv[argv.index("--merge") + 1:]:
            part = json.loads(Path(path).read_text())
            records.update({r["name"]: r for r in part["records"]})
            data["tools"] = part["tools"]
        data["records"] = sorted(records.values(), key=lambda r: r["name"])
        RESULT.write_text(json.dumps(data, indent=1) + "\n")
        return 0
    if "--control" in argv:
        data = json.loads(RESULT.read_text()) if RESULT.exists() else {"records": []}
        data["control"] = control()
        RESULT.write_text(json.dumps(data, indent=1) + "\n")
        return 0 if data["control"]["intact_accepted"] and not data["control"]["corrupted_accepted"] else 1
    if "--verify" in argv:
        rest = [a for a in argv if a != "--verify"]
        return verify(rest[0] if rest else RESULT)
    timeout = float(argv[argv.index("--timeout") + 1]) if "--timeout" in argv else 7200
    only = set(argv[argv.index("--only") + 1].split(",")) if "--only" in argv else None
    result = Path(argv[argv.index("--out") + 1]) if "--out" in argv else RESULT
    drat_route = "--drat" in argv
    if not LRAT_CHECK.exists() or (drat_route and not DRAT_TRIM.exists()):
        print("build lrat-check and drat-trim first: sh refs/fetch.sh drat-trim")
        return 1
    WORK.mkdir(parents=True, exist_ok=True)
    data = json.loads(result.read_text()) if result.exists() else {"records": []}
    data["tools"] = tools()
    records = {r["name"]: r for r in data["records"]}
    for name, build in cases():
        if only is not None and name not in only:
            continue
        rec = (prove_drat if drat_route else prove)(name, dimacs(*build()), timeout)
        records[name] = rec
        print(f"{name}: {rec['status']}, solver exit {rec['solver_exit']} in {rec['solver_seconds']} s, proof "
              f"{rec['proof']['bytes']:,} bytes, {rec['proof']['lines']:,} lines, checker exit {rec['checker_exit']}",
              flush=True)
        data["records"] = sorted(records.values(), key=lambda r: r["name"])
        result.write_text(json.dumps(data, indent=1) + "\n")
    return 0 if all(r["status"] == "accepted" for r in records.values()) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
