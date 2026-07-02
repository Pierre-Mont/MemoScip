#!/usr/bin/env python3
"""Convert Pisinger knapsack CSV bundles into one LP file per instance.

Input format (per README in hardinstances_pisinger):
  instance_name
  n <int>
  c <int>
  z <int>
  time <float>
  i,p_i,w_i,x_i        for i=1..n

Blocks are separated by a line containing only '-----'.

For an input file NAME.csv, the script creates:
  <parent>/NAME/
and writes one LP per instance in that folder.

Usage:
  python3 Instance/pisinger_csv_to_lp.py Instance/hardinstances_pisinger/knapPI_11_20_1000.csv
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Item:
    index: int
    profit: int
    weight: int


@dataclass
class KnapsackInstance:
    name: str
    n_items: int
    capacity: int
    optimum: int
    solve_time: float
    items: list[Item]


def _parse_prefixed_int(line: str, prefix: str) -> int:
    parts = line.split()
    if len(parts) != 2 or parts[0] != prefix:
        raise ValueError(f"Expected '{prefix} <value>', got: {line!r}")
    try:
        return int(parts[1])
    except ValueError as exc:
        raise ValueError(f"Invalid integer in line: {line!r}") from exc


def _parse_prefixed_float(line: str, prefix: str) -> float:
    parts = line.split()
    if len(parts) != 2 or parts[0] != prefix:
        raise ValueError(f"Expected '{prefix} <value>', got: {line!r}")
    try:
        return float(parts[1])
    except ValueError as exc:
        raise ValueError(f"Invalid float in line: {line!r}") from exc


def parse_instances(csv_path: Path) -> list[KnapsackInstance]:
    """Parse all instances from a Pisinger CSV bundle file."""
    raw_lines = csv_path.read_text(encoding="utf-8").splitlines()
    lines = [ln.strip() for ln in raw_lines]

    instances: list[KnapsackInstance] = []
    i = 0
    total = len(lines)

    while i < total:
        while i < total and (not lines[i] or lines[i] == "-----"):
            i += 1
        if i >= total:
            break

        name = lines[i]
        i += 1
        if i + 3 >= total:
            raise ValueError(f"Unexpected EOF while reading header for instance {name!r}")

        n_items = _parse_prefixed_int(lines[i], "n")
        i += 1
        capacity = _parse_prefixed_int(lines[i], "c")
        i += 1
        optimum = _parse_prefixed_int(lines[i], "z")
        i += 1
        solve_time = _parse_prefixed_float(lines[i], "time")
        i += 1

        items: list[Item] = []
        for _ in range(n_items):
            if i >= total:
                raise ValueError(f"Unexpected EOF in items of instance {name!r}")
            row = lines[i]
            i += 1
            if not row:
                raise ValueError(f"Unexpected blank item line in instance {name!r}")

            parts = [p.strip() for p in row.split(",")]
            if len(parts) != 4:
                raise ValueError(f"Expected 4 comma-separated fields, got: {row!r}")

            try:
                idx = int(parts[0])
                profit = int(parts[1])
                weight = int(parts[2])
                _x_sol = int(parts[3])
            except ValueError as exc:
                raise ValueError(f"Invalid numeric item row: {row!r}") from exc

            items.append(Item(index=idx, profit=profit, weight=weight))

        expected_indices = list(range(1, n_items + 1))
        found_indices = [it.index for it in items]
        if found_indices != expected_indices:
            raise ValueError(
                f"Item indices mismatch in {name!r}: expected 1..{n_items}, got {found_indices[:5]}..."
            )

        instances.append(
            KnapsackInstance(
                name=name,
                n_items=n_items,
                capacity=capacity,
                optimum=optimum,
                solve_time=solve_time,
                items=items,
            )
        )

    if not instances:
        raise ValueError(f"No instance found in file: {csv_path}")
    return instances


def write_lp(instance: KnapsackInstance, output_path: Path) -> None:
    """Write a 0/1 knapsack LP with split binaries x_i_0/x_i_1."""
    lines: list[str] = []
    lines.append("\\ Pisinger knapsack instance converted to LP")
    lines.append(f"\\ source_instance={instance.name}")
    lines.append(f"\\ n={instance.n_items} capacity={instance.capacity} z={instance.optimum} time={instance.solve_time}")
    lines.append("Minimize")

    # Start from min(-sum p_i x_i_1), then use x_i_1 = 1 - x_i_0 when needed
    # so all objective coefficients remain nonnegative.
    obj_terms: list[str] = []
    objective_constant_shift = 0
    for it in instance.items:
        if it.profit > 0:
            obj_terms.append(f"{it.profit} x{it.index}_0")
            objective_constant_shift -= it.profit
        elif it.profit < 0:
            obj_terms.append(f"{-it.profit} x{it.index}_1")

    if obj_terms:
        lines.append(" obj: " + " + ".join(obj_terms))
    else:
        lines.append(" obj: 0")
    lines.append(f"\\ objective_constant_shift={objective_constant_shift}")
    lines.append("")
    lines.append("Subject To")

    total_weight = sum(it.weight for it in instance.items)
    rhs = total_weight - instance.capacity
    weight_terms = [f"{it.weight} x{it.index}_0" for it in instance.items]
    lines.append(" capacity: " + " + ".join(weight_terms) + f" >= {rhs}")

    for it in instance.items:
        lines.append(f" split_{it.index}: x{it.index}_0 + x{it.index}_1 = 1")

    lines.append("")
    lines.append("Binary")

    chunk_size = 12
    vars_all = [var for it in instance.items for var in (f"x{it.index}_0", f"x{it.index}_1")]
    for start in range(0, len(vars_all), chunk_size):
        lines.append(" " + " ".join(vars_all[start : start + chunk_size]))

    lines.append("")
    lines.append("End")
    lines.append("")

    output_path.write_text("\n".join(lines), encoding="utf-8")


def default_output_dir(csv_path: Path) -> Path:
    return csv_path.with_suffix("")


def convert_bundle(csv_path: Path, output_dir: Path) -> int:
    instances = parse_instances(csv_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    for inst in instances:
        out_path = output_dir / f"{inst.name}.lp"
        write_lp(inst, out_path)

    return len(instances)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Convert Pisinger CSV bundle into one .lp file per instance."
    )
    parser.add_argument("csv", type=Path, help="Input Pisinger CSV file")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory (default: input filename without .csv)",
    )
    args = parser.parse_args(argv)

    csv_path: Path = args.csv
    if not csv_path.is_file():
        print(f"Error: input file not found: {csv_path}", file=sys.stderr)
        return 1

    output_dir = args.output_dir if args.output_dir is not None else default_output_dir(csv_path)

    try:
        count = convert_bundle(csv_path, output_dir)
    except (ValueError, OSError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print(f"Wrote {count} LP files to: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())