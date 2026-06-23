#!/usr/bin/env python3
"""Run memoscipl over instances and collect end-of-run metrics into a CSV.

Usage examples:
  python3 Instance/run_all.py -l1 -v=2
  python3 Instance/run_all.py --exec /path/to/memoscipl --inst-dir Instance/knapsack --exts .lp,.wcsp --timeout 300 -- -l1 -v=2

Everything unknown to the script's own argparse is forwarded to the executable.
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
    import re
    out = re.sub(r"[^A-Za-z0-9]+", "_", s)
    out = out.strip("_")
    return out or "opts"


PATTERNS = {
    "status": re.compile(r"Solving Status\s*:\s*(\S+)", re.IGNORECASE),
    "time": re.compile(r"Solving time\s*:\s*([\d\.Ee+-]+)", re.IGNORECASE),
    "nodes": re.compile(r"Nodes\s*:\s*(\d+)", re.IGNORECASE),
    "primal": re.compile(r"Primal bound\s*:\s*(\S+)", re.IGNORECASE),
    "dual": re.compile(r"Dual bound\s*:\s*(\S+)", re.IGNORECASE),
    "cuts": re.compile(r"Number of learned cuts\s*:\s*(\d+)", re.IGNORECASE),
    "gap": re.compile(r"Gap\s*:\s*(\S+)", re.IGNORECASE),
}


def parse_metrics(text: str):
    res = {}
    for k, rx in PATTERNS.items():
        m = rx.search(text)
        res[k] = m.group(1) if m else ""
    return res


def run_instance(exe: Path, inst: Path, forward_args, timeout: int):
    cmd = [str(exe)] + forward_args + [str(inst)]
    try:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=timeout)
        return proc.stdout, proc.returncode
    except subprocess.TimeoutExpired as e:
        return f"TIMEOUT after {timeout}s", -1
    except Exception as e:
        return f"ERROR: {e}", -2


def main():
    parser = argparse.ArgumentParser(description="Run memoscipl over a set of instances and gather results.")
    parser.add_argument("--exec", default="build/memoscipl", help="Path to memoscipl executable")
    parser.add_argument("--inst-dir", default="Instance/knapsack", help="Directory with instances")
    parser.add_argument("--exts", default=".lp", help="Comma-separated list of extensions to include (e.g. .lp,.wcsp)")
    parser.add_argument("--out-dir", default=".", help="Where to write the CSV")
    parser.add_argument("--timeout", type=int, default=3600, help="Per-run timeout in seconds")
    parser.add_argument("--parallel", action="store_true", help="Run in parallel (not implemented in this script) - reserved")

    # parse_known_args: forward unknown args to the executable (e.g. -l1 -v=2)
    args, unknown = parser.parse_known_args()

    # Prefer absolute paths under the repository root
    repo_root = Path(__file__).parent.parent.resolve()

    exe = Path(args.exec)
    if not exe.is_absolute():
        exe = (repo_root / exe).resolve()

    inst_dir = Path(args.inst_dir)
    if not inst_dir.is_absolute():
        inst_dir = (repo_root / inst_dir).resolve()

    exts = [e.strip() for e in args.exts.split(",") if e.strip()]

    forward_args = unknown

    # build csv filename including sanitized forward args and timestamp
    opts_part = sanitize_for_filename("_".join(forward_args)) if forward_args else "default"
    ts = datetime.now().strftime("%Y%m%d")
    csv_name = f"results__{opts_part}__{ts}.csv"
    csv_path = inst_dir / csv_name

    # collect instance files
    instances = []
    if inst_dir.exists() and inst_dir.is_dir():
        for ext in exts:
            instances.extend(sorted(inst_dir.glob(f"*{ext}")))
    else:
        print(f"Instance directory not found: {inst_dir}")
        sys.exit(2)

    if not instances:
        print("No instances found matching extensions", exts)
        sys.exit(2)

    fieldnames = ["file", "returncode", "status", "time", "nodes", "# cuts", "primal", "dual", "gap"]

    # Create CSV and write header before running instances so file exists even if script crashes
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        f.flush()
        try:
            os.fsync(f.fileno())
        except Exception:
            pass

    # Run instances and append each result immediately
    for inst in instances:
        inst = inst.resolve()
        print("Running", inst)
        out_text, rc = run_instance(exe, inst, forward_args, args.timeout)
        # parse metrics from combined stdout/stderr
        metrics = parse_metrics(out_text if isinstance(out_text, str) else "")
        row = {
            "file": inst.name,
            "returncode": rc,
            "status": metrics.get("status", ""),
            "time": metrics.get("time", ""),
            "nodes": metrics.get("nodes", ""),
            "# cuts": metrics.get("cuts", ""),
            "primal": metrics.get("primal", ""),
            "dual": metrics.get("dual", ""),
            "gap": metrics.get("gap", ""),
        }

        # append the row and fsync to ensure it's persisted
        with csv_path.open("a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writerow(row)
            f.flush()
            try:
                os.fsync(f.fileno())
            except Exception:
                pass

    print(f"Wrote results to {csv_path}")


if __name__ == "__main__":
    main()
