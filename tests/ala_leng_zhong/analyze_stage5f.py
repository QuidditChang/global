#!/usr/bin/env python3
"""Classify the observation-only Leng Stage 5F diagnostic CSV files."""

import argparse
import csv
import json
import math
from pathlib import Path


def rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def number(row, key, default=math.nan):
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return default


def last_by_cap(inner, outer, kind):
    selected = {}
    for row in inner:
        if int(row["outer"]) == outer and row["kind"] == kind:
            cap = int(row["cap"])
            if cap not in selected or int(row["iteration"]) > int(
                    selected[cap]["iteration"]):
                selected[cap] = row
    return selected


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "prefix",
        help="solver log prefix before .leng_stage5F_{outer,inner,gauge}.csv")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    prefix = Path(args.prefix)
    outer_path = Path(str(prefix) + ".leng_stage5F_outer.csv")
    inner_path = Path(str(prefix) + ".leng_stage5F_inner.csv")
    gauge_path = Path(str(prefix) + ".leng_stage5F_gauge.csv")
    outer = rows(outer_path)
    inner = rows(inner_path)
    gauge = rows(gauge_path)[0]
    if not outer or not inner:
        raise SystemExit("Stage 5F CSV files are empty")

    terminal = outer[-1]
    snapshot_outers = sorted({int(row["outer"]) for row in inner
                              if row["kind"] == "frozen_replay"})
    snapshots = {}
    cap_limited = False
    preconditioner_limited = False
    for outer_number in snapshot_outers:
        replay = last_by_cap(inner, outer_number, "frozen_replay")
        guard = last_by_cap(inner, outer_number, "frozen_guard")
        values = {str(cap): number(row, "R_inner")
                  for cap, row in replay.items()}
        if 60 in replay and 240 in replay:
            r60 = number(replay[60], "R_inner")
            r240 = number(replay[240], "R_inner")
            cap_limited |= math.isfinite(r60) and math.isfinite(r240) and (
                r240 < 0.5 * r60)
        if 240 in replay and 240 in guard:
            patch = number(replay[240], "R_inner")
            diagonal = number(guard[240], "R_inner")
            values["diagonal_240"] = diagonal
            preconditioner_limited |= (math.isfinite(patch) and
                                       math.isfinite(diagonal) and
                                       min(patch, diagonal) <
                                       0.5 * max(patch, diagonal))
        snapshots[str(outer_number)] = values

    compatibility = max(number(row, "compatibility_relative", 0.0)
                        for row in inner)
    transpose_defect = number(gauge, "transpose_defect")
    transpose_probe_valid = number(
        gauge, "transpose_shared_probe", 0.0) == 1.0
    c_mpi = number(terminal, "C_MPI")
    c_horizontal = number(terminal, "C_horizontal")
    c_radial = number(terminal, "C_radial")
    ratio = number(terminal, "norm_ratio_1")
    cosine = number(terminal, "cosine_1")
    two_cycle = number(terminal, "two_cycle_return")

    flags = {
        "transpose_probe_invalid": not transpose_probe_valid,
        "operator_transpose_inconsistent": (
            transpose_probe_valid and transpose_defect > 1.0e-8),
        # A relative constant projection grows automatically as CG removes
        # orthogonal residual.  Call it a blocking compatibility floor only
        # when it dominates the explicit residual, not merely when nonzero.
        "pressure_compatibility_failure": compatibility > 0.5,
        "inner_iteration_cap_limited": cap_limited,
        "preconditioner_sensitive": preconditioner_limited,
        "mpi_interface_concentration": max(c_mpi, c_horizontal, c_radial) > 10.0,
        "outer_fixed_point_stagnation": ratio >= 0.9,
        "outer_two_cycle": abs(cosine) >= 0.9 and two_cycle <= 0.25,
    }
    priority = [
        "transpose_probe_invalid",
        "operator_transpose_inconsistent",
        "pressure_compatibility_failure",
        "mpi_interface_concentration",
        "preconditioner_sensitive",
        "inner_iteration_cap_limited",
        "outer_two_cycle",
        "outer_fixed_point_stagnation",
    ]
    primary = next((name for name in priority if flags[name]), "unresolved")
    result = {
        "diagnostic": "Leng Stage 5F-A pressure spectroscopy",
        "classification": primary,
        "flags": flags,
        "thresholds": {
            "transpose_defect": 1.0e-8,
            "compatibility_dominant_fraction": 0.5,
            "mpi_concentration": 10.0,
            "cap_improvement_factor": 0.5,
            "outer_stagnation_ratio": 0.9,
            "two_cycle_cosine": 0.9,
            "two_cycle_return": 0.25,
        },
        "terminal_outer": {key: terminal.get(key) for key in terminal},
        "gauge": gauge,
        "frozen_snapshots": snapshots,
        "inputs": [str(outer_path), str(inner_path), str(gauge_path)],
    }
    output = args.output or Path(str(prefix) + ".leng_stage5F_decision.json")
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print("%s: %s" % (output, primary))


if __name__ == "__main__":
    main()
