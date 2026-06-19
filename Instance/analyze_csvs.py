#!/usr/bin/env python3
"""Concatenate result CSVs in this directory into a wide table and compute averages.

Usage: python3 analyze_csvs.py

The script looks for files matching `results__*__*.csv` in the current directory,
reads columns `file`, `time`, `nodes`, `# cuts`, `gap`, `primal` when available, and
creates a combined Excel file `combined_results.xlsx` with columns for each approach:

    <approach>_gap, <approach>_primal, <approach>_time, <approach>_nodes

It prints all rows (one per instance) and then prints average time and average nodes
for each approach (averaging over available values). If `openpyxl` is not available
the script falls back to writing `combined_results.csv` instead.
"""
import argparse
import csv
import re
from pathlib import Path
from collections import defaultdict
import statistics
import sys
import math


def find_result_files(dirpath: Path):
    return sorted(p for p in dirpath.glob('results__*.csv') if p.is_file())


def extract_approach(name: str):
    # match results__<opts>__<date>.csv
    m = re.match(r'results__(.+?)__(\d{8}).*\.csv$', name)
    if m:
        opts = re.sub(r'[^A-Za-z0-9_\-]+', '_', m.group(1))
        date = m.group(2)
        return f"{opts}__{date}"
    # fallback: use stem
    return re.sub(r'[^A-Za-z0-9_\-]+', '_', Path(name).stem)


def try_float(s):
    try:
        return float(s)
    except Exception:
        return None


def try_int(s):
    try:
        return int(float(s))
    except Exception:
        return None


def parse_args():
    parser = argparse.ArgumentParser(
        description='Concatenate result CSVs in a directory into a wide table and compute averages.'
    )
    parser.add_argument(
        'directory',
        nargs='?', 
        default='.',
        help='Directory containing result CSV files (default: current directory)'
    )
    return parser.parse_args()


def main():
    args = parse_args()
    dirpath = Path(args.directory)
    if not dirpath.exists() or not dirpath.is_dir():
        print(f"Directory not found or not a directory: {dirpath}")
        sys.exit(1)

    files = find_result_files(dirpath)
    if not files:
        print("No result CSV files found (looking for results__*.csv)")
        sys.exit(1)

    # instance -> approach -> metrics
    data = defaultdict(lambda: {})
    approaches = []

    for f in files:
        approach = extract_approach(f.name)
        approaches.append(approach)
        with f.open(newline='') as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                inst = row.get('file') or row.get('File') or ''
                if not inst:
                    continue
                # collect metrics
                t = try_float(row.get('time', '').strip())
                n = try_int(row.get('nodes', '').strip())
                gap = row.get('gap', '').strip()
                p = try_float(row.get('primal', '').strip())
                cuts = try_int(row.get('# cuts', '').strip()) if '# cuts' in row else try_int(row.get('cuts', '').strip()) if 'cuts' in row else None
                data[inst][approach] = {'time': t, 'nodes': n, 'gap': gap, 'primal': p, 'cuts': cuts}

    # prepare combined CSV
    # columns: file, then for each approach: <approach>_gap, <approach>_primal, <approach>_time, <approach>_nodes
    approaches = list(dict.fromkeys(approaches))  # preserve order, unique
    fieldnames = ['file']
    for a in approaches:
        fieldnames += [f"{a}_gap", f"{a}_primal", f"{a}_time", f"{a}_nodes", f"{a}_cuts"]

    # Try to write an Excel file using openpyxl; fall back to CSV if unavailable
    xlsx_path = dirpath / 'combined_results.xlsx'
    try:
        from openpyxl import Workbook
    except Exception:
        outpath = dirpath / 'combined_results.csv'
        with outpath.open('w', newline='') as outfh:
            writer = csv.DictWriter(outfh, fieldnames=fieldnames)
            writer.writeheader()
            # print rows too
            for inst in sorted(data.keys()):
                prow = {'file': inst}
                for a in approaches:
                    m = data[inst].get(a, {})
                    prow[f"{a}_gap"] = m.get('gap', '')
                    prow[f"{a}_primal"] = '' if m.get('primal') is None else f"{m.get('primal')}"
                    prow[f"{a}_time"] = '' if m.get('time') is None else f"{m.get('time')}"
                    prow[f"{a}_nodes"] = '' if m.get('nodes') is None else f"{m.get('nodes')}"
                    prow[f"{a}_cuts"] = '' if m.get('cuts') is None else f"{m.get('cuts')}"
                writer.writerow(prow)

        print(f"openpyxl not available; wrote combined CSV to {outpath}\n")
    else:
        wb = Workbook()
        ws = wb.active
        # write header
        ws.append(fieldnames)
        for inst in sorted(data.keys()):
            prow = {'file': inst}
            for a in approaches:
                m = data[inst].get(a, {})
                prow[f"{a}_gap"] = m.get('gap', '')
                prow[f"{a}_primal"] = '' if m.get('primal') is None else f"{m.get('primal')}"
                prow[f"{a}_time"] = '' if m.get('time') is None else f"{m.get('time')}"
                prow[f"{a}_nodes"] = '' if m.get('nodes') is None else f"{m.get('nodes')}"
                prow[f"{a}_cuts"] = '' if m.get('cuts') is None else f"{m.get('cuts')}"
            row = [prow.get(fn, '') for fn in fieldnames]
            ws.append(row)
        wb.save(xlsx_path)
        print(f"Wrote combined Excel to {xlsx_path}\n")

    # Plotting: average metric per instance-size for each approach
    def plot_metric_by_size(data, approaches, metric_key, ylabel, title, out_fname="plot.png"):
        # aggregate metric by instance size
        size_groups = {}  # size -> approach -> list of metric values
        size_rx = re.compile(r"knapsack-(\d+)-\d+\.lp")
        for inst, amap in data.items():
            msz = size_rx.search(inst)
            if not msz:
                continue
            size = int(msz.group(1))
            if size not in size_groups:
                size_groups[size] = {a: [] for a in approaches}
            for a in approaches:
                entry = amap.get(a)
                if not entry:
                    continue
                value = entry.get(metric_key)
                if value is None:
                    continue
                size_groups[size].setdefault(a, []).append(value)

        if not size_groups:
            print(f"No instance names matched size pattern; skipping {metric_key} plot.")
            return {}

        sizes = sorted(size_groups.keys())

        try:
            import matplotlib.pyplot as plt
        except Exception:
            print(f"matplotlib not available; skipping {metric_key} plot.")
            return {}

        plt.figure()
        approach_avgs = {}
        for a in approaches:
            y = []
            x = []
            for s in sizes:
                values = size_groups[s].get(a, [])
                if values:
                    avg_value = statistics.mean(values)
                    approach_avgs.setdefault(a, {})[s] = avg_value
                else:
                    avg_value = math.nan
                x.append(s)
                y.append(avg_value)
            xs = [xx for xx, yy in zip(x, y) if not math.isnan(yy)]
            ys = [yy for yy in y if not math.isnan(yy)]
            if xs and ys:
                plt.plot(xs, ys, marker='o', label=a)

        plt.xlabel('Instance size')
        plt.ylabel(ylabel)
        plt.title(title)
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        outpath = Path(out_fname)
        plt.savefig(outpath)
        plt.close()
        print(f"Wrote {out_fname} plot to {outpath}")
        return approach_avgs

    # call plotting functions and get averages per size for nodes and time
    nodes_avgs = plot_metric_by_size(
        data, approaches, 'nodes', 'Average nodes developed',
        'Average nodes by instance size per approach',
        out_fname=str(dirpath / 'nodes_by_size.png')
    )
    time_avgs = plot_metric_by_size(
        data, approaches, 'time', 'Average time',
        'Average time by instance size per approach',
        out_fname=str(dirpath / 'time_by_size.png')
    )

    # print coefficient of difference between size 100 and 300 for each approach
    print("\nCoefficient de différence (relative) entre taille 100 et 300:")
    for a in approaches:
        avgs = nodes_avgs.get(a, {}) if nodes_avgs else {}
        avg100 = avgs.get(100)
        avg300 = avgs.get(300)
        if avg100 is None:
            print(f"Approach: {a}  -- pas de donnée pour taille 100")
            continue
        if avg300 is None:
            print(f"Approach: {a}  -- pas de donnée pour taille 300")
            continue
        if avg100 == 0:
            print(f"Approach: {a}  -- average nodes for size 100 is zero (cannot compute relative difference)")
            continue
        coeff = (avg300 - avg100) / avg100
        print(f"Approach: {a}  coeff = {coeff:.6f}  ({coeff*100:.2f}%)  avg100={avg100:.2f}  avg300={avg300:.2f}")

    # find inconsistent instances: gap=0 with different primal values across approaches
    inconsistent = []
    for inst, amap in data.items():
        zero_gap_primal = {}
        for a, m in amap.items():
            if m.get('gap') == '0' or m.get('gap') == '0.0':
                if m.get('primal') is not None:
                    zero_gap_primal[a] = m.get('primal')
        if len(zero_gap_primal) > 1:
            unique_primals = set(zero_gap_primal.values())
            if len(unique_primals) > 1:
                inconsistent.append((inst, zero_gap_primal))

    if inconsistent:
        print(f"Found {len(inconsistent)} inconsistent instance(s) with gap=0 and differing primal values.")
        print("Details:")
        print("file,approach,primal,gap")
        for inst, zero_gap_primal in sorted(inconsistent):
            for a, primal in sorted(zero_gap_primal.items()):
                print(f"{inst},{a},{primal},0")
        print("")
    else:
        print("No inconsistent instances found where two approaches report gap=0 with different primal values.\n")

    # compute averages per approach
    print("Averages per approach (over available instances):")
    for a in approaches:
        times = []
        nodes = []
        cuts_list = []
        for inst, amap in data.items():
            m = amap.get(a)
            if not m:
                continue
            if m.get('time') is not None:
                times.append(m.get('time'))
            if m.get('nodes') is not None:
                nodes.append(m.get('nodes'))
            if m.get('cuts') is not None:
                cuts_list.append(m.get('cuts'))
        tavg = statistics.mean(times) if times else None
        navg = statistics.mean(nodes) if nodes else None
        cavg = statistics.mean(cuts_list) if cuts_list else None
        print(f"Approach: {a}")
        if tavg is not None:
            print(f"  average_time = {tavg:.6g}   (n={len(times)})")
        else:
            print("  average_time = (no data)")
        if navg is not None:
            print(f"  average_nodes = {navg:.6g}   (n={len(nodes)})")
        else:
            print("  average_nodes = (no data)")
        if cavg is not None:
            print(f"  average_cuts = {cavg:.6g}   (n={len(cuts_list)})")
        else:
            print("  average_cuts = (no data)")


if __name__ == '__main__':
    main()