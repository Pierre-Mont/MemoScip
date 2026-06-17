#!/usr/bin/env python3
"""Run SCIP default executable over instances and collect end-of-run metrics into a CSV.

Usage examples:
  python3 Instance/run_scip_default.py
  python3 Instance/run_scip_default.py --exec scip --inst-dir Instance/knapsack --timeout 600
"""
import argparse
import csv
import os
import re
import subprocess
from datetime import datetime
from pathlib import Path
import sys


def sanitize_for_filename(s: str) -> str:
    out = re.sub(r"[^A-Za-z0-9]+", "_", s)
    out = out.strip("_")
    return out or "opts"


PATTERNS = {
    "status": re.compile(r"SCIP Status\s*:\s*(.+?)(?:\[|$)", re.IGNORECASE),
    "time": re.compile(r"Solving Time \(sec\)\s*:\s*([\d\.Ee+-]+)", re.IGNORECASE),
    "nodes": re.compile(r"Solving Nodes\s*:\s*(\d+)", re.IGNORECASE),
    "primal": re.compile(r"Primal Bound\s*:\s*([\+\-]?[\d\.Ee+-]+)", re.IGNORECASE),
    "dual": re.compile(r"Dual Bound\s*:\s*([\+\-]?[\d\.Ee+-]+)", re.IGNORECASE),
    "gap": re.compile(r"Gap\s*:\s*([\d\.]+)\s*%", re.IGNORECASE),
}


def parse_metrics(text: str):
    result = {}
    for key, regex in PATTERNS.items():
        match = regex.search(text)
        result[key] = match.group(1) if match else ""
    return result


def run_instance(exe: Path, inst: Path, timeout: int, no_presolving: bool = False):
    # Use SCIP default executable and solve the instance. We read the instance, optimize, then quit.
    cmd = [str(exe), "-c", f"read {inst}"]
    if no_presolving:
        cmd.extend(["-c", "set presolving emphasis off"])
    cmd.extend(["-c", "opt", "-c", "quit"])
    try:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=timeout)
        return proc.stdout, proc.returncode
    except subprocess.TimeoutExpired as exc:
        return f"TIMEOUT after {timeout}s", -1
    except Exception as exc:
        return f"ERROR: {exc}", -2


def main():
    parser = argparse.ArgumentParser(description="Run SCIP over a set of instances and gather results.")
    parser.add_argument("--exec", default="/home/kisouke/local/scip-master/bin/scip", help="SCIP executable to run")
    parser.add_argument("--inst-dir", default="Instance/knapsack", help="Directory with instances")
    parser.add_argument("--exts", default=".lp", help="Comma-separated list of extensions to include (e.g. .lp,.wcsp)")
    parser.add_argument("-p", "--no-presolving", action="store_true", help="Disable presolving")
    parser.add_argument("--timeout", type=int, default=600, help="Per-run timeout in seconds")
    args = parser.parse_args()

    repo_root = Path("/home/kisouke/MemoSCIP")

    exe = Path(args.exec)
    if not exe.is_absolute():
        exe = (repo_root / exe).resolve()

    inst_dir = Path(args.inst_dir)
    if not inst_dir.is_absolute():
        inst_dir = (repo_root / inst_dir).resolve()

    out_dir = inst_dir

    exts = [ext.strip() for ext in args.exts.split(",") if ext.strip()]

    if not inst_dir.exists() or not inst_dir.is_dir():
        print(f"Instance directory not found: {inst_dir}")
        sys.exit(2)

    instances = []
    for ext in exts:
        instances.extend(sorted(inst_dir.glob(f"*{ext}")))

    if not instances:
        print("No instances found matching extensions", exts)
        sys.exit(2)

    out_dir.mkdir(parents=True, exist_ok=True)
    opts_part = "default"
    if args.no_presolving:
        opts_part += "-p"
    opts_part = sanitize_for_filename(opts_part)
    ts = datetime.now().strftime("%Y%m%d")
    csv_name = f"results__scip_{opts_part}__{ts}.csv"
    csv_path = out_dir / csv_name

    fieldnames = ["file", "returncode", "status", "time", "nodes", "primal", "dual", "gap"]
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        f.flush()
        try:
            os.fsync(f.fileno())
        except Exception:
            pass

    for inst in instances:
        inst = inst.resolve()
        print("Running", inst)
        output, returncode = run_instance(exe, inst, args.timeout, args.no_presolving)
        metrics = parse_metrics(output)
        row = {
            "file": inst.name,
            "returncode": returncode,
            "status": metrics.get("status", ""),
            "time": metrics.get("time", ""),
            "nodes": metrics.get("nodes", ""),
            "primal": metrics.get("primal", ""),
            "dual": metrics.get("dual", ""),
            "gap": metrics.get("gap", ""),
        }
        with csv_path.open("a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writerow(row)
            f.flush()
            try:
                os.fsync(f.fileno())
            except Exception:
                pass

    print(f"Wrote SCIP default results to {csv_path}")


if __name__ == "__main__":
    main()
