"""Check the core-only periodic proof using the shared Erdős axiom probe and an already installed Lean."""

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
SHARED_GATE = ROOT.parent / "formal" / "scripts" / "verify_lean.py"
SOURCE = ROOT / "src" / "PeriodicTransfer.lean"
RESULTS = ROOT / "results" / "periodic-lean.json"
TIMEOUT = 120
REQUIRED = {
    "PeriodicTransfer.mono_iff_period_dvd",
    "PeriodicTransfer.fourColor_ap_free",
    "PeriodicTransfer.fourColor_mono_iff",
}


def check(plant=None):
    spec = importlib.util.spec_from_file_location("shared_erdos_gate", SHARED_GATE)
    gate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)
    toolchain_file = gate.ROOT / "lean-toolchain"
    toolchain = toolchain_file.read_text().strip()
    installed = subprocess.run(["elan", "toolchain", "list"], capture_output=True, text=True,
                               check=True, timeout=TIMEOUT)
    if toolchain not in {line.split()[0] for line in installed.stdout.splitlines() if line.strip()}:
        raise RuntimeError(f"{toolchain} is not installed; this checker does not install toolchains")
    prefix = subprocess.run(["elan", "run", toolchain, "lean", "--print-prefix"], capture_output=True,
                            text=True, check=True, timeout=TIMEOUT).stdout.strip()
    lean, checker = Path(prefix)/"bin"/"lean", Path(prefix)/"bin"/"leanchecker"
    source = SOURCE.read_text()
    if plant == "sorry":
        source += "\ntheorem plantedHole : False := by sorry\n"
    elif plant == "axiom":
        source += "\naxiom plantedAssumption : False\n"
    elif plant == "claim":
        before = "Mono c m k a d ↔ m ∣ d"
        if source.count(before) != 1:
            raise ValueError("Could not locate the statement for the deliberate mutation")
        source = source.replace(before, "Mono c m k a d ↔ m ∣ (d + 1)")
    refusal = gate.refusal_for(source)
    report = {
        "source_sha256": hashlib.sha256(source.encode()).hexdigest(),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "shared_gate": {"path": "formal/scripts/verify_lean.py",
                        "sha256": hashlib.sha256(SHARED_GATE.read_bytes()).hexdigest()},
        "toolchain": toolchain,
        "scope": "Own declarations compiled, runtime axiom-probed, and independently kernel-replayed; imported Std is trusted. Unisolated honest-source check. No Mathlib build or library installation.",
        "plant": plant, "declarations": {}, "passed": False,
    }
    if refusal:
        return {**report, "status": "refused", "detail": refusal}
    with tempfile.TemporaryDirectory(prefix="erdos1186-periodic-") as directory:
        work = Path(directory)
        (work/"Cand.lean").write_text(source)
        (work/"Probe.lean").write_text(gate.PROBE_SOURCE.format(module="Cand"))
        env = {**os.environ, "LEAN_PATH": str(work), "LEAN_SYSROOT": prefix}
        build = subprocess.run([str(lean), "Cand.lean", "-o", "Cand.olean"], cwd=work, env=env,
                               capture_output=True, text=True, timeout=TIMEOUT)
        if build.returncode:
            return {**report, "status": "compile_error", "detail": build.stdout+build.stderr}
        probe = subprocess.run([str(lean), "--run", "Probe.lean"], cwd=work, env=env,
                               capture_output=True, text=True, timeout=TIMEOUT)
        declarations = gate.parse_probe(probe.stdout) if probe.returncode == 0 else None
        if declarations is None:
            return {**report, "status": "probe_error", "detail": probe.stdout+probe.stderr}
        replay = subprocess.run([str(checker), "Cand"], cwd=work, env=env,
                                capture_output=True, text=True, timeout=TIMEOUT)
    report["declarations"] = {
        name: {"axioms": sorted(axioms), "verdict": gate.verdict_for(axioms)}
        for name, axioms in sorted(declarations.items())
    }
    report["kernel_replay_exit"] = replay.returncode
    report["required_declarations_present"] = REQUIRED <= declarations.keys()
    report["passed"] = (
        replay.returncode == 0 and report["required_declarations_present"]
        and all(row["verdict"] == "verified" for row in report["declarations"].values())
    )
    report["status"] = "verified" if report["passed"] else "rejected"
    if replay.returncode:
        report["detail"] = replay.stdout+replay.stderr
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--plant", choices=("sorry", "axiom", "claim"))
    args = parser.parse_args()
    report = check(args.plant)
    print(f"{'PASS' if report['passed'] else 'FAIL'} periodic Lean proof: {report['status']}")
    for name, row in report["declarations"].items():
        print(name, row["verdict"], row["axioms"])
    if report.get("detail"):
        print(report["detail"])
    if not report["passed"]:
        return 1
    if args.write:
        RESULTS.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if args.check and json.loads(RESULTS.read_text()) != report:
        print("FAIL: stored periodic Lean result differs")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
