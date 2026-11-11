#!/usr/bin/env python3
"""Reproduce the finite verifications without modifying any fixed input.

Usage:
  python3 verification/run_all.py --list
  python3 verification/run_all.py --jobs 3
  python3 verification/run_all.py --select affine,origin_t --jobs 2

Each job runs in its own temporary copy. New logs, receipts and a summary
are written beneath --output (default: verification/results).
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse
import hashlib
import importlib.metadata
import json
import shutil
import subprocess
import sys
import tempfile
import time

HERE = Path(__file__).resolve().parent


def verify_manifest():
    manifest = json.loads((HERE / "MANIFEST.json").read_text())
    for relative, expected in manifest["sha256"].items():
        path = HERE / relative
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError("Fixed input hash mismatch: " + relative)
    return len(manifest["sha256"])


def nested(data, dotted):
    for key in dotted.split("."):
        data = data[key]
    return data


def run_one(job, output):
    started = time.monotonic()
    job_id = job["id"]
    log_path = output / (job_id + ".log")
    result = {"id": job_id, "status": "FAIL"}
    try:
        with tempfile.TemporaryDirectory(prefix="schmeisser_" + job_id + "_") as tmp:
            work = Path(tmp) / "verification"
            fixed = json.loads((HERE / "MANIFEST.json").read_text())["sha256"]
            for relative in fixed:
                target = work / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(HERE / relative, target)
            command = [sys.executable, str(work / job["script"])] + job.get("args", [])
            with log_path.open("w") as stream:
                completed = subprocess.run(
                    command, cwd=work / job.get("cwd", "."),
                    stdout=stream, stderr=subprocess.STDOUT, check=False)
            if completed.returncode:
                raise RuntimeError("Verifier exit code " + str(completed.returncode))
            if job.get("receipt"):
                receipt_path = work / job["receipt"]
                receipt = json.loads(receipt_path.read_text())
            else:
                receipt = json.loads(log_path.read_text())
            for field, value in job["expect"].items():
                if nested(receipt, field) != value:
                    raise AssertionError("Unaccepted receipt field " + field)
            saved = output / (job_id + ".json")
            saved.write_text(json.dumps(receipt, indent=2) + "\n")
            result.update(status="PASS", receipt=saved.name, log=log_path.name)
    except Exception as error:
        result["error"] = repr(error)
        result["log"] = log_path.name
    result["seconds"] = round(time.monotonic() - started, 3)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true", help="List available components")
    parser.add_argument("--select", help="Comma-separated component identifiers")
    parser.add_argument("--jobs", type=int, default=1, help="Independent concurrent verifiers")
    parser.add_argument("--output", type=Path, default=HERE / "results")
    parser.add_argument("--extended", action="store_true",
                        help="Include additional independent implementations")
    args = parser.parse_args()
    index = json.loads((HERE / "INDEX.json").read_text())
    jobs = index["jobs"]
    if args.list:
        for job in jobs:
            print(job["id"] + "\t" + job["description"])
        return 0
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    if args.select:
        wanted = set(args.select.split(","))
        available = {job["id"] for job in jobs}
        if not wanted <= available:
            parser.error("Unknown component: " + ",".join(sorted(wanted-available)))
        jobs = [job for job in jobs if job["id"] in wanted]
    elif not args.extended:
        jobs = [job for job in jobs if job.get("default", True)]
    checked_files = verify_manifest()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    results = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        pending = {pool.submit(run_one, job, output): job["id"] for job in jobs}
        for future in as_completed(pending):
            result = future.result()
            results.append(result)
            print(result["id"] + ": " + result["status"] +
                  " (" + str(result["seconds"]) + " s)", flush=True)
    versions = {}
    for name in ("python-flint", "sympy"):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = None
    summary = {
        "status": "PASS" if all(r["status"] == "PASS" for r in results) else "FAIL",
        "degree_floor": 6251, "python": sys.version.split()[0],
        "packages": versions, "fixed_files_checked": checked_files,
        "components_run": len(results),
        "components": sorted(results, key=lambda r: r["id"]),
    }
    (output / "SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("Overall: " + summary["status"], flush=True)
    return 0 if summary["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
