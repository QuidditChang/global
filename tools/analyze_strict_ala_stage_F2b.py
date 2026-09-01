#!/usr/bin/env python3
"""Fail-closed Stage-F2b velocity-MG level localization."""

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

CHECKPOINTS = (1, 2, 4, 8, 16, 32, 64)
LEVEL_COLUMNS = (
    "POD_mode", "POD_energy_weight", "MG_cycle", "top_level", "level",
    "phase", "input_rms", "output_rms", "reduction_ratio", "alpha", "valid")


def load_json(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing or empty {path}")
    return json.loads(path.read_text())


def read_csv(path, columns):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing or empty {path}")
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != tuple(columns):
            raise ValueError(f"schema mismatch in {path.name}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"no rows in {path.name}")
    return rows


def finite(value, name):
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"non-finite {name}")
    return result


def weighted_rms(rows, field):
    numerator = denominator = 0.0
    for row in rows:
        weight = finite(row["POD_energy_weight"], "POD weight")
        value = finite(row[field], field)
        if weight < 0.0 or value < 0.0:
            raise ValueError(f"invalid weighted metric {field}")
        numerator += weight * value * value
        denominator += weight
    if denominator <= 0.0:
        raise ValueError("nonpositive POD weight sum")
    return math.sqrt(numerator / denominator)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def analyze(args):
    f2 = load_json(args.f2_decision)
    runtime = load_json(args.runtime)
    levels = read_csv(args.level_transfer, LEVEL_COLUMNS)
    trajectory = read_csv(args.trajectory, (
        "POD_mode", "POD_energy_weight", "MG_cycles",
        "momentum_residual_relative", "valid"))
    mode_rhs = read_csv(args.mode_rhs,
                        ("POD_mode", "POD_energy_weight", "rhs_norm", "valid"))
    f2_trajectory = read_csv(args.f2_trajectory, (
        "variant", "POD_mode", "POD_energy_weight", "MG_cycles",
        "momentum_residual_relative", "velocity_error_relative",
        "K_energy_error_relative", "Schur_action_error_relative",
        "velocity_cosine", "exact_tight_achieved", "valid"))
    f2_tight = read_csv(args.f2_tight, (
        "call_id", "POD_mode", "rhs_norm", "requested_relative_tolerance",
        "target_residual", "achieved_relative_residual", "cycles",
        "max_cycles", "converged"))

    selected = runtime.get("selected_modes")
    if not isinstance(selected, int) or selected <= 0:
        raise ValueError("invalid selected mode count")
    expected = {(mode, cycle) for mode in range(1, selected + 1)
                for cycle in CHECKPOINTS}
    trajectory_keys = {(int(row["POD_mode"]), int(row["MG_cycles"]))
                       for row in trajectory}
    if len(trajectory) != len(expected) or trajectory_keys != expected:
        raise ValueError("incomplete or duplicate F2b trajectory")
    if len(mode_rhs) != selected or {int(row["POD_mode"]) for row in mode_rhs} != \
            set(range(1, selected + 1)):
        raise ValueError("invalid F2b mode RHS rows")

    f2_configured = {(int(row["POD_mode"]), int(row["MG_cycles"])): row
                     for row in f2_trajectory
                     if row["variant"] == "CONFIGURED_VANKA"}
    if set(f2_configured) != expected:
        raise ValueError("incomplete configured F2 lineage trajectory")
    trajectory_match = True
    for row in trajectory:
        key = (int(row["POD_mode"]), int(row["MG_cycles"]))
        before = finite(f2_configured[key]["momentum_residual_relative"],
                        "F2 momentum")
        after = finite(row["momentum_residual_relative"], "F2b momentum")
        trajectory_match &= math.isclose(before, after, rel_tol=1.0e-10,
                                         abs_tol=1.0e-14)
    tight_by_mode = {int(row["POD_mode"]): row for row in f2_tight}
    rhs_match = len(tight_by_mode) == selected
    for row in mode_rhs:
        mode = int(row["POD_mode"])
        rhs_match &= mode in tight_by_mode and math.isclose(
            finite(row["rhs_norm"], "F2b RHS"),
            finite(tight_by_mode[mode]["rhs_norm"], "F2 RHS"),
            rel_tol=1.0e-12, abs_tol=1.0e-300)

    level_valid = all(
        row["valid"] == "1" and int(row["POD_mode"]) in range(1, selected + 1)
        and int(row["MG_cycle"]) in CHECKPOINTS
        and int(row["top_level"]) >= int(row["level"])
        and finite(row["input_rms"], "level input") >= 0.0
        and finite(row["output_rms"], "level output") >= 0.0
        and finite(row["reduction_ratio"], "level ratio") >= 0.0
        and math.isfinite(finite(row["alpha"], "level alpha"))
        for row in levels)
    max_top = max(int(row["top_level"]) for row in levels)
    final_rows = [row for row in levels
                  if int(row["top_level"]) == max_top and
                  int(row["MG_cycle"]) == 64]
    required_phases = {"DOWN_SMOOTH", "RESTRICTION", "BOTTOM_SOLVE",
                       "UP_CORRECTION"}
    phase_coverage = all(
        required_phases <= {row["phase"] for row in final_rows
                            if int(row["POD_mode"]) == mode}
        for mode in range(1, selected + 1))
    lineage = (f2.get("experiment_evidence_valid") is True and
               f2.get("decision") ==
                   "F2_VANKA_EFFECTIVE_BUT_AVV_CONVERGENCE_LIMITED" and
               f2.get("next_authorized_task") ==
                   "F2B_COARSE_GRID_TRANSFER_LOCALIZATION" and
               f2.get("production_default_change_authorized") is False)
    gates = {
        "F2_lineage": lineage,
        "runtime_schema": runtime.get("schema") ==
            "strict-ala-stage-F2b-runtime-v1",
        "runtime_mode_count": selected == len(mode_rhs),
        "runtime_selected_energy": finite(runtime.get("selected_energy"),
                                           "selected energy") >= 0.95,
        "frozen_operator": runtime.get("operator") == "frozen_current_K_gamma",
        "configured_vanka": runtime.get("velocity_smoother") ==
            "configured_element_vanka",
        "checkpoint_contract": tuple(runtime.get("checkpoints", ())) == CHECKPOINTS,
        "mode_rhs_identity": rhs_match,
        "F2_trajectory_identity": trajectory_match,
        "level_rows": level_valid and phase_coverage,
        "production_freeze":
            runtime.get("production_default_change_authorized") is False,
    }
    structural_pass = all(gates.values())

    grouped = defaultdict(list)
    for row in final_rows:
        if row["phase"] in ("DOWN_SMOOTH", "BOTTOM_SOLVE", "UP_CORRECTION"):
            grouped[(row["phase"], int(row["level"]))].append(row)
    reductions = {f"{phase}_L{level}": weighted_rms(rows, "reduction_ratio")
                  for (phase, level), rows in grouped.items()}
    if not reductions:
        raise ValueError("no F2b level reduction metrics")
    worst_key = max(reductions, key=reductions.get)
    worst_phase, worst_level_text = worst_key.rsplit("_L", 1)
    worst_level = int(worst_level_text)
    trajectory64 = [row for row in trajectory if row["MG_cycles"] == "64"]
    metrics = {
        "maximum_top_level": max_top,
        "momentum_residual_RMS_64": weighted_rms(
            trajectory64, "momentum_residual_relative"),
        "level_reduction_RMS_cycle_64": reductions,
        "worst_reduction_phase": worst_phase,
        "worst_reduction_level": worst_level,
        "worst_reduction_ratio": reductions[worst_key],
    }
    if not structural_pass:
        decision = "F2B_INVALID_STRUCTURAL_EVIDENCE"
        next_task = "REPEAT_F2B_VALID_EXPERIMENT"
    elif worst_phase == "BOTTOM_SOLVE":
        decision = "F2B_COARSE_SOLVE_BOTTLENECK"
        next_task = "F2C_COARSE_VELOCITY_SOLVER_REDESIGN"
    elif worst_phase == "DOWN_SMOOTH":
        decision = "F2B_LEVEL_SMOOTHER_BOTTLENECK"
        next_task = "F2C_LEVEL_SMOOTHER_REDESIGN"
    else:
        decision = "F2B_UPWARD_TRANSFER_BOTTLENECK"
        next_task = "F2C_PROLONGATION_GALERKIN_REVIEW"
    result = {
        "schema": "strict-ala-stage-F2b-decision-v1",
        "experiment_evidence_valid": structural_pass,
        "candidate": "configured_vanka_level_localization",
        "metrics": metrics,
        "gates": gates,
        "decision": decision,
        "production_default_change_authorized": False,
        "next_authorized_task": next_task,
    }
    write_json(args.output, result)
    return 0


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def audit(args):
    decision = load_json(args.decision)
    provenance = {}
    for name in ("binary", "inputs", "source"):
        before = Path(getattr(args, f"{name}_pre"))
        after = Path(getattr(args, f"{name}_post"))
        provenance[name] = {
            "pre_sha256": sha256(before), "post_sha256": sha256(after),
            "unchanged": before.read_bytes() == after.read_bytes()}
    valid = (decision.get("experiment_evidence_valid") is True and
             all(item["unchanged"] for item in provenance.values()))
    result = {
        "schema": "strict-ala-stage-F2b-final-audit-v1",
        "experiment_valid": valid,
        "decision": decision.get("decision") if valid else "INVALID_EXPERIMENT",
        "provenance": provenance,
        "F2b": decision,
        "production_default_change_authorized": False,
        "next_authorized_task": (decision.get("next_authorized_task") if valid
                                 else "REPEAT_F2B_VALID_EXPERIMENT"),
    }
    write_json(args.output, result)
    return 0


def parser():
    root = argparse.ArgumentParser()
    commands = root.add_subparsers(required=True)
    command = commands.add_parser("analyze")
    for name in ("f2-decision", "f2-trajectory", "f2-tight",
                 "level-transfer", "trajectory", "mode-rhs", "runtime",
                 "output"):
        command.add_argument(f"--{name}", required=True)
    command.set_defaults(fn=analyze)
    command = commands.add_parser("audit")
    command.add_argument("--decision", required=True)
    for name in ("binary-pre", "binary-post", "inputs-pre", "inputs-post",
                 "source-pre", "source-post", "output"):
        command.add_argument(f"--{name}", required=True)
    command.set_defaults(fn=audit)
    return root


def main():
    args = parser().parse_args()
    return args.fn(args)


if __name__ == "__main__":
    raise SystemExit(main())
