#!/usr/bin/env python3
"""
Convert KPCG / Bettinelli-style knapsack-with-conflicts instances to CPLEX LP.

Input format (AMPL-like, e.g. Instance/KPCG/R1/BPPC_1_0_1.txt_0.1):
  param n := <number of items>;
  param c := <knapsack capacity>;
  param : V : p w :=          (item index, profit, weight)
     ...
  ;
  set E :=                    (undirected conflict edges)
     i  j
     ...
  ;

MIP (minimize profit via complementary variables), with x_j encoded as x_j_0 + x_j_1 = 1:
  min  sum_j p_j * x_j_0        (equivalent to max sum_j p_j * x_j_1)
  s.t. sum_j w_j * x_j_0 >= sum_j w_j - c   (equivalent to sum_j w_j * x_j_1 <= c)
       x_i_1 + x_j_1 <= 1   for each edge {i,j} in E
       x_j_0 + x_j_1 = 1   for each item j
       x_j_0, x_j_1 in {0,1}

Note: In this benchmark, `c` is capacity, not the number of conflicts.
The suffix after the last '_' in the filename (e.g. 0.1) is the conflict-graph
density, not the knapsack capacity.

Usage:
  python3 kpcg2lp.py INPUT [OUTPUT]
  python3 kpcg2lp.py Instance/KPCG/R1/BPPC_1_0_1.txt_0.1
"""

from __future__ import annotations

import argparse
import io
import lzma
import re
import sys
from pathlib import Path
from typing import BinaryIO, Iterable, TextIO


class WidthFile:
    """Write LP lines with a maximum line length (SCIP/CPLEX friendly)."""

    maxcol = 80

    def __init__(self, path: Path) -> None:
        self.f = path.open("w", encoding="utf-8")
        self.col = 0

    def write(self, text: str) -> None:
        for chunk in text.splitlines(keepends=True):
            if not chunk:
                continue
            line = chunk.rstrip("\n")
            if self.col and self.col + len(line) >= self.maxcol:
                self.f.write("\n")
                self.col = 0
            self.f.write(line)
            if chunk.endswith("\n"):
                self.f.write("\n")
                self.col = 0
            else:
                self.col += len(line)

    def writeln(self, text: str = "") -> None:
        self.write(text + "\n")

    def close(self) -> None:
        self.f.close()


def open_instance(path: Path) -> TextIO:
    """Open plain text or .xz-compressed instance."""
    if path.suffix == ".xz":
        raw = lzma.open(path, "rb")
        return io.TextIOWrapper(raw, encoding="utf-8")
    return path.open(encoding="utf-8")


def parse_int_param(line: str, name: str) -> int | None:
    m = re.match(rf"param\s+{re.escape(name)}\s*:=\s*(\d+)\s*;", line.strip())
    return int(m.group(1)) if m else None


def parse_instance(stream: TextIO) -> tuple[int, int, list[tuple[int, int, int]], list[tuple[int, int]]]:
    n_items: int | None = None
    capacity: int | None = None
    items: list[tuple[int, int, int]] = []
    edges: list[tuple[int, int]] = []

    mode: str | None = None  # "items" | "edges"

    for raw in stream:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue

        if n_items is None:
            val = parse_int_param(line, "n")
            if val is not None:
                n_items = val
                continue

        if capacity is None:
            val = parse_int_param(line, "c")
            if val is not None:
                capacity = val
                continue

        if line.startswith("param :"):
            mode = "items"
            continue

        if line.startswith("set E"):
            mode = "edges"
            continue

        if line == ";":
            if mode == "items":
                mode = None
            elif mode == "edges":
                break
            continue

        parts = line.split()
        if mode == "items" and len(parts) >= 3:
            idx, profit, weight = int(parts[0]), int(parts[1]), int(parts[2])
            items.append((idx, profit, weight))
        elif mode == "edges" and len(parts) >= 2:
            i, j = int(parts[0]), int(parts[1])
            if i > j:
                i, j = j, i
            edges.append((i, j))

    if n_items is None or capacity is None:
        raise ValueError("Missing param n and/or param c in instance file")
    if len(items) != n_items:
        raise ValueError(f"Expected {n_items} items, found {len(items)}")

    indices = [idx for idx, _, _ in items]
    if sorted(indices) != list(range(n_items)):
        raise ValueError(f"Item indices must be 0..{n_items - 1}, got {sorted(indices)}")

    return n_items, capacity, items, edges


def lin_expr(coeffs: Iterable[tuple[int, str]], first: bool = True) -> str:
    parts: list[str] = []
    for coeff, var in coeffs:
        if coeff == 0:
            continue
        sign = "+" if coeff > 0 else "-"
        mag = abs(coeff)
        term = var if mag == 1 else f"{mag} {var}"
        if not parts:
            parts.append(f"{sign}{term}" if sign == "+" else f"-{term}" if mag == 1 else f"-{mag} {var}")
        else:
            parts.append(f" {sign} {term}" if mag > 1 else f" {sign} {var}")
    return "".join(parts) if parts else "0"


def write_lp(
    out: WidthFile,
    n_items: int,
    capacity: int,
    items: list[tuple[int, int, int]],
    edges: list[tuple[int, int]],
    instance_name: str,
) -> None:
    var_one = lambda j: f"x{j}_1"
    var_zero = lambda j: f"x{j}_0"

    out.writeln("\\ KPCG instance converted to LP")
    out.writeln(f"\\ source: {instance_name}")
    out.writeln(f"\\ n={n_items} capacity={capacity} conflicts={len(edges)}")
    out.writeln("Minimize")
    obj_terms = [(profit, var_zero(idx)) for idx, profit, _ in items]
    out.writeln(" obj: " + lin_expr(obj_terms))

    out.writeln("")
    out.writeln("Subject To")
    total_weight = sum(weight for _, _, weight in items)
    knapsack_rhs = total_weight - capacity
    knapsack_terms = [(weight, var_zero(idx)) for idx, _, weight in items]
    out.writeln(" knapsack: " + lin_expr(knapsack_terms) + f" >= {knapsack_rhs}")

    for k, (i, j) in enumerate(edges):
        out.writeln(f" conflict_{k}: {var_one(i)} + {var_one(j)} <= 1")

    for j in range(n_items):
        out.writeln(f" domain_constraint_{j}: +1 {var_zero(j)} +1 {var_one(j)} = 1")

    out.writeln("")
    out.writeln("Binary")
    binary_vars = [name for j in range(n_items) for name in (var_zero(j), var_one(j))]
    for start in range(0, len(binary_vars), 10):
        chunk = " ".join(binary_vars[start : start + 10])
        out.writeln(f" {chunk}")

    out.writeln("")
    out.writeln("End")


def default_output_path(input_path: Path) -> Path:
    name = input_path.name
    if name.endswith(".xz"):
        name = name[:-3]
    return input_path.with_name(name + ".lp")


def convert(input_path: Path, output_path: Path) -> None:
    with open_instance(input_path) as f:
        n_items, capacity, items, edges = parse_instance(f)
    out = WidthFile(output_path)
    try:
        write_lp(out, n_items, capacity, items, edges, str(input_path))
    finally:
        out.close()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Convert KPCG knapsack-with-conflicts instances to LP format."
    )
    parser.add_argument("input", type=Path, help="Input instance (.txt, .txt_0.1, .xz, ...)")
    parser.add_argument(
        "output",
        type=Path,
        nargs="?",
        help="Output .lp file (default: same basename with .lp extension)",
    )
    args = parser.parse_args(argv)

    input_path: Path = args.input
    if not input_path.is_file():
        print(f"Error: input file not found: {input_path}", file=sys.stderr)
        return 1

    output_path = args.output if args.output is not None else default_output_path(input_path)

    try:
        convert(input_path, output_path)
    except (ValueError, lzma.LZMAError, OSError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print(f"Wrote {output_path} ({output_path.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
