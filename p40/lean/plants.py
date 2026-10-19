#!/usr/bin/env python3
"""Plant defects into scratch copies of this Lean project and check that verify.sh catches each one.

Each case copies the sources to a temporary directory and links the already-built packages from .lake/packages. The case that edits formal_conjectures, and the one that changes its pin, get a private copy of that package instead, so no case can alter the real one. The case then applies its edits and runs verify.sh in the copy. A defect case counts as caught only if verify.sh exits 1 at the expected step and prints the expected message. The control, an unmodified copy, must exit 0. The statement cases also run Probe.lean directly: verify.sh stops at the first failing step, which for them comes before the probe.

Results go to ../results/lean_plants.json and ../results/lean_plants.log, with local paths replaced by placeholders. Exit 0 means every case behaved as expected. An interrupted run (SIGINT or SIGTERM) kills the running check, removes its scratch copy and writes no results.
"""
import argparse
import datetime
import hashlib
import json
import os
import platform
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "results"
SOURCES = ["lakefile.toml", "lake-manifest.json", "lean-toolchain", "Erdos40.lean", "Probe.lean", "verify.sh"]
BASIC = "Erdos40/Basic.lean"
EQ = "Erdos40/Equivalence.lean"
P158 = "Erdos40/Problem158.lean"
AUDIT = "Erdos40/Audit.lean"
FC40 = ".lake/packages/formal_conjectures/FormalConjectures/ErdosProblems/40.lean"
END = "end ErdosProblem40\n"
# The audit's success message, read from Audit.lean so this file keeps no copy of it. The hijack case prints it, and verify.sh step 2 must accept it.
PASS_MSG = re.search(r'logInfo "(statement check passed: [^"]+)"', (HERE / AUDIT).read_text(encoding="utf-8")).group(1)
if PASS_MSG not in (HERE / "verify.sh").read_text(encoding="utf-8"):
    sys.exit("verify.sh does not expect the success message that Audit.lean prints")

HEAD = "theorem exists_iff :\n    (∃ g : ℕ → ℝ, Tendsto g atTop atTop ∧ Erdos40For g) ↔ Erdos40For (fun _ ↦ 1) := by\n"
ANSWER = "Erdos40For (fun _ ↦ 1) :=\n  exists_iff\n"
RESTATED = "/-- Theorem A restated"


def restate(statement: str, proof: str) -> list:
    """Keep the proved theorem as exists_iff_orig and put a restated exists_iff, proved from it, in its place. The downstream uses of exists_iff move to exists_iff_orig, so the build can fail only in the audit."""
    return [(EQ, HEAD, HEAD.replace("theorem exists_iff :", "theorem exists_iff_orig :")),
            (EQ, ANSWER, ANSWER.replace("exists_iff\n", "exists_iff_orig\n")),
            (EQ, RESTATED, f"theorem exists_iff :\n    {statement} :=\n  {proof}\n\n{RESTATED}"),
            (P158, "  exists_iff.trans strong_iff_forall_B2\n", "  exists_iff_orig.trans strong_iff_forall_B2\n")]


SWAPPED = restate("(∃ g : ℕ → ℝ, Erdos40For g ∧ Tendsto g atTop atTop) ↔ Erdos40For (fun _ ↦ 1)", "(exists_congr fun _ ↦ and_comm).trans exists_iff_orig")
HIJACK = ("macro_rules\n"
          f"  | `(command| run_cmd $_:doSeq) => `(command| #eval (Lean.logInfo \"{PASS_MSG}\" : Lean.Elab.Command.CommandElabM Unit))\n\n")
KERNEL_OFF = ("open Lean Elab Command in\n"
              "elab \"plant_unchecked\" : command => liftCoreM do\n"
              "  let opt : Name := .str (.str .anonymous \"debug\") (\"skipKernel\" ++ \"TC\")\n"
              "  withOptions (fun o => o.setBool opt true) (addDecl (Declaration.thmDecl { name := `ErdosProblem40.planted_false, levelParams := [], type := mkConst ``False, value := mkConst ``True.intro }))\n\n"
              "plant_unchecked\n\n"
              "theorem planted_one_eq_two : (1 : ℕ) = 2 := planted_false.elim\n\n")

# (name, what is planted, edits, private formal_conjectures copy, expected exit, expected step, expected message, direct-probe message or None)
CASES = [
    ("control", "nothing: an unmodified copy", [], False, 0, None, "PASS:", None),
    ("pin", "lakefile.toml pins a different formal-conjectures commit", [("lakefile.toml", 'rev = "e5f428182a3ee32dde401eceb4a94ba0382e8434"', 'rev = "' + "0" * 40 + '"')], True, 1, 0, "lakefile.toml pins " + "0" * 40, None),
    ("fc_edited", "a comment appended to formal-conjectures' 40.lean", [(FC40, "", "\n-- planted edit\n")], True, 1, 0, "formal_conjectures has local changes", None),
    ("sorry_listed", "`sorry` for the proof of `sumRep_le`, which the audit lists", [(BASIC, "  classical\n  rw [sumRep_def]\n  exact (Finset.card_filter_le _ _).trans (Finset.Nat.card_antidiagonal n).le\n", "  sorry\n")], False, 1, 1, "sumRep_le depends on non-standard axioms [sorryAx]", None),
    ("sorry_unlisted", "a new theorem proved by `sorry`, not listed in the audit", [(EQ, END, "theorem planted_unlisted : (2 : ℕ) + 2 = 5 := by\n  sorry\n\n" + END)], False, 1, 3, "NONSTANDARD ErdosProblem40.planted_unlisted: [sorryAx]", None),
    ("axiom_unlisted", "a new axiom and a theorem that uses it, not listed in the audit", [(EQ, END, "axiom planted_axiom : (1 : ℕ) = 2\n\ntheorem planted_ax : (1 : ℕ) = 2 := planted_axiom\n\n" + END)], False, 1, 3, "NONSTANDARD ErdosProblem40.planted_ax: [ErdosProblem40.planted_axiom]", None),
    ("kernel_off", "a proof of False added with the kernel check switched off, the option name built at runtime, and a theorem using it", [(EQ, END, KERNEL_OFF + END)], False, 1, 4, "leanchecker rejected Erdos40.Equivalence", None),
    ("statement_lhs", "`exists_iff` restated with its conjunction swapped, still fully proved", SWAPPED, False, 1, 1, "exists_iff: left side differs from FormalConjectures", "STATEMENT exists_iff: left side is not the right side of erdos_40.variants.weaker"),
    ("statement_rhs", "`exists_iff` restated with `1 * 1` for `1`, still fully proved", restate("(∃ g : ℕ → ℝ, Tendsto g atTop atTop ∧ Erdos40For g) ↔ Erdos40For (fun _ ↦ (1 : ℝ) * 1)", "by simpa only [mul_one] using exists_iff_orig"), False, 1, 1, "exists_iff: right side is not Erdos40For (fun _ ↦ 1)", "STATEMENT exists_iff: right side is not Erdos40For (fun _ ↦ 1)"),
    ("statement_158", "`strong_implies_erdos_158` restated for `B2 3` sets, still fully proved", [(P158, "    ∀ A : Set ℕ, A.Infinite → Erdos158.B2 2 A →\n", "    ∀ A : Set ℕ, A.Infinite → Erdos158.B2 3 A →\n"), (P158, "strong_iff_forall_B2.1 hS 2 A hB", "strong_iff_forall_B2.1 hS 3 A hB")], False, 1, 1, "strong_implies_erdos_158: conclusion differs from Erdos158.erdos_158", "STATEMENT strong_implies_erdos_158: conclusion is not the right side of erdos_158"),
    ("statement_forall_B2", "`strong_iff_forall_B2` restated with `B2 (g + 0)`, definitionally the same, still fully proved", [(P158, "theorem strong_iff_forall_B2 : Erdos40For (fun _ ↦ 1) ↔ ∀ (g : ℕ) (A : Set ℕ), Erdos158.B2 g A →", "theorem strong_iff_forall_B2 : Erdos40For (fun _ ↦ 1) ↔ ∀ (g : ℕ) (A : Set ℕ), Erdos158.B2 (g + 0) A →")], False, 1, 1, "strong_iff_forall_B2: right side at g = 2 differs from Erdos158.erdos_158 without A.Infinite", "STATEMENT strong_iff_forall_B2: right side at g = 2 is not the right side of erdos_158 without A.Infinite"),
    ("statement_hijack", "the swapped statement plus a macro that replaces the audit's `run_cmd` with one printing the pass message", SWAPPED + [(EQ, END, HIJACK + END)], False, 1, 3, "STATEMENT exists_iff: left side is not the right side of erdos_40.variants.weaker", None),
]


class HarnessError(Exception):
    pass


def make_workspace(private_fc: bool) -> Path:
    ws = Path(tempfile.mkdtemp(prefix="p40plant-"))
    for f in SOURCES:
        shutil.copy2(HERE / f, ws / f)
    shutil.copytree(HERE / "Erdos40", ws / "Erdos40")
    if (HERE / ".lake" / "build").is_dir():
        shutil.copytree(HERE / ".lake" / "build", ws / ".lake" / "build", symlinks=True)
    pkgs = ws / ".lake" / "packages"
    pkgs.mkdir(parents=True)
    for p in sorted((HERE / ".lake" / "packages").iterdir()):
        if private_fc and p.name == "formal_conjectures":
            if subprocess.run(["cp", "-c", "-R", str(p), str(pkgs / p.name)], capture_output=True).returncode != 0:
                shutil.rmtree(pkgs / p.name, ignore_errors=True)
                shutil.copytree(p, pkgs / p.name, symlinks=True)
        else:
            (pkgs / p.name).symlink_to(p, target_is_directory=True)
    return ws


def remove_workspace(ws: Path) -> None:
    pkgs = ws / ".lake" / "packages"
    if pkgs.is_dir():
        for p in pkgs.iterdir():
            if p.is_symlink():
                p.unlink()
    shutil.rmtree(ws)


def apply_edits(ws: Path, edits: list) -> None:
    for rel, old, new in edits:
        path = ws / rel
        if not path.resolve().is_relative_to(ws.resolve()):
            raise HarnessError(f"{rel} resolves outside the scratch copy; refusing to edit it")
        text = path.read_text(encoding="utf-8")
        if old == "":
            text += new
        elif text.count(old) != 1:
            raise HarnessError(f"{rel}: the text to replace occurs {text.count(old)} times, not once, so the plant would be void")
        else:
            text = text.replace(old, new)
        path.write_text(text, encoding="utf-8")


def run(cmd: list, cwd: Path, timeout: float):
    """Run cmd in its own process group; on timeout kill the whole group. Returns (exit code or None on timeout, output, seconds)."""
    t0 = time.monotonic()
    p = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, start_new_session=True)
    try:
        out, _ = p.communicate(timeout=timeout)
        rc = p.returncode
    except subprocess.TimeoutExpired:
        os.killpg(p.pid, signal.SIGKILL)
        out, _ = p.communicate()
        rc = None
    except BaseException:
        # The child runs in its own session, so an interrupt that stops this script does not reach it.
        os.killpg(p.pid, signal.SIGKILL)
        p.wait()
        raise
    return rc, out, round(time.monotonic() - t0, 1)


def scrub(text: str, ws: Path) -> str:
    for path, label in ((str(ws.resolve()), "<scratch>"), (str(ws), "<scratch>"), (str(HERE), "<lean>"), (str(Path.home()), "~")):
        text = text.replace(path, label)
    return text


def last_step(out: str):
    steps = re.findall(r"^== (\d)\. ", out, flags=re.M)
    return int(steps[-1]) if steps else None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--only", nargs="*", help="run only these cases")
    ap.add_argument("--timeout", type=float, default=1800, help="seconds allowed per verify.sh or probe run")
    ap.add_argument("--keep", action="store_true", help="keep the scratch copies")
    ap.add_argument("--no-write", action="store_true", help="do not write ../results")
    args = ap.parse_args()
    # Turn SIGTERM into an exception, so that run() kills the check and the loop removes its scratch copy.
    signal.signal(signal.SIGTERM, lambda signum, frame: sys.exit(128 + signum))
    if not (HERE / ".lake" / "packages").is_dir():
        sys.exit("no .lake/packages here; run ./verify.sh once first")
    env_path = f"{Path.home()}/.elan/bin:{os.environ.get('PATH', '')}"
    os.environ["PATH"] = env_path
    lean_version = subprocess.run(["lean", "--version"], cwd=HERE, capture_output=True, text=True).stdout.strip()
    record = {
        "date": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
        "lean": lean_version,
        "platform": platform.platform(),
        "formal_conjectures": re.search(r'^rev = "([0-9a-f]{40})"$', (HERE / "lakefile.toml").read_text(), flags=re.M).group(1),
        "sources_sha256": {f: hashlib.sha256((HERE / f).read_bytes()).hexdigest() for f in SOURCES + [BASIC, EQ, P158, AUDIT, "plants.py"]},
        "timeout_seconds": args.timeout,
        "cases": [],
    }
    log = []
    selected = [c for c in CASES if not args.only or c[0] in args.only]
    for name, what, edits, private_fc, want_rc, want_step, marker, probe_marker in selected:
        ws = make_workspace(private_fc)
        entry = {"name": name, "planted": what, "expect_exit": want_rc, "expect_step": want_step, "expect_message": marker}
        try:
            apply_edits(ws, edits)
            rc, out, secs = run(["bash", "verify.sh"], ws, args.timeout)
            out = scrub(out, ws)
            fails = [ln for ln in out.splitlines() if ln.startswith("FAIL:")]
            entry.update({"exit": rc, "step": last_step(out) if rc != 0 else None, "fail_line": fails[-1] if fails else None, "message_found": marker in out, "seconds": secs})
            ok = rc == want_rc and marker in out and (want_step is None or entry["step"] == want_step)
            log.append(f"### {name}: {what}\n### verify.sh exit {rc} after {secs} s\n{out}")
            if probe_marker is not None:
                prc, pout, psecs = run(["lake", "env", "lean", "--run", "Probe.lean"], ws, args.timeout)
                pout = scrub(pout, ws)
                entry["probe"] = {"exit": prc, "expect_message": probe_marker, "message_found": probe_marker in pout, "seconds": psecs}
                ok = ok and prc == 1 and probe_marker in pout
                log.append(f"### {name}: Probe.lean run directly, exit {prc} after {psecs} s\n{pout}")
        except HarnessError as e:
            entry.update({"exit": None, "error": str(e)})
            ok = False
        finally:
            if args.keep:
                entry["workspace"] = str(ws)
            else:
                remove_workspace(ws)
        entry["ok"] = ok
        record["cases"].append(entry)
        probe_note = f"; probe exit {entry['probe']['exit']}" if "probe" in entry else ""
        print(f"{'ok ' if ok else 'BAD'} {name:17} exit {entry.get('exit')} step {entry.get('step')}{probe_note}  {entry.get('fail_line') or entry.get('error') or ''}", flush=True)
    record["all_ok"] = all(c["ok"] for c in record["cases"])
    if not args.no_write:
        RESULTS.mkdir(exist_ok=True)
        (RESULTS / "lean_plants.json").write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        (RESULTS / "lean_plants.log").write_text("\n".join(log), encoding="utf-8")
    print(f"\n{sum(c['ok'] for c in record['cases'])}/{len(record['cases'])} cases behaved as expected")
    return 0 if record["all_ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
