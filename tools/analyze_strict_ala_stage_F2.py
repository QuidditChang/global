#!/usr/bin/env python3
"""Fail-closed Stage-F2 frozen actual-mode velocity-MG review."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

CHECKPOINTS = (1, 2, 4, 8, 16, 32, 64)
VARIANTS = ("CONFIGURED_VANKA", "DIAGONAL_ROLLBACK")


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
    f1d = load_json(args.f1d_decision)
    runtime = load_json(args.runtime)
    trajectory = read_csv(args.trajectory, (
        "variant", "POD_mode", "POD_energy_weight", "MG_cycles",
        "momentum_residual_relative", "velocity_error_relative",
        "K_energy_error_relative", "Schur_action_error_relative",
        "velocity_cosine", "exact_tight_achieved", "valid"))
    tight = read_csv(args.tight, (
        "call_id", "POD_mode", "rhs_norm",
        "requested_relative_tolerance", "target_residual",
        "achieved_relative_residual", "cycles", "max_cycles", "converged"))

    selected = runtime.get("selected_modes")
    if not isinstance(selected, int) or selected <= 0:
        raise ValueError("invalid selected mode count")
    expected = {(variant, mode, cycle) for variant in VARIANTS
                for mode in range(1, selected + 1) for cycle in CHECKPOINTS}
    keys = {(row["variant"], int(row["POD_mode"]),
             int(row["MG_cycles"])) for row in trajectory}
    if len(trajectory) != len(expected) or keys != expected:
        raise ValueError("incomplete or duplicate F2 trajectory")
    if len(tight) != selected or [int(row["call_id"]) for row in tight] != \
            list(range(1, selected + 1)):
        raise ValueError("invalid F2 tight-solve count")

    by_key = {(row["variant"], int(row["POD_mode"]),
               int(row["MG_cycles"])): row for row in trajectory}
    for mode in range(1, selected + 1):
        weights = [finite(by_key[variant, mode, cycle]["POD_energy_weight"],
                          "POD weight") for variant in VARIANTS
                   for cycle in CHECKPOINTS]
        if max(weights) - min(weights) > 1.0e-12 * max(max(weights), 1.0):
            raise ValueError("POD weight changed across F2 trajectories")
    numeric_fields = (
        "momentum_residual_relative", "velocity_error_relative",
        "K_energy_error_relative", "Schur_action_error_relative")
    rows_valid = all(row["valid"] == "1" and
                     all(finite(row[field], field) >= 0.0
                         for field in numeric_fields) and
                     -1.0 - 1.0e-12 <= finite(row["velocity_cosine"],
                                              "velocity cosine") <= 1.0 + 1.0e-12
                     for row in trajectory)
    tight_pass = all(
        row["converged"] == "1" and
        finite(row["achieved_relative_residual"], "tight achieved") <=
        finite(row["requested_relative_tolerance"], "tight requested")
        for row in tight)
    lineage = (f1d.get("experiment_evidence_valid") is True and
               f1d.get("decision") == "F1D_RESIDUAL_CLOSED_COARSE_REJECTED" and
               f1d.get("next_authorized_task") == "F2_AVV_PRECONDITIONER_REVIEW" and
               f1d.get("production_default_change_authorized") is False)
    structural_gates = {
        "F1d_lineage": lineage,
        "runtime_schema": runtime.get("schema") ==
            "strict-ala-stage-F2-runtime-v1",
        "runtime_mode_count": runtime.get("selected_modes") == selected,
        "runtime_selected_energy": finite(runtime.get("selected_energy"),
                                           "selected energy") >= 0.95,
        "frozen_operator": runtime.get("operator") == "frozen_current_K_gamma",
        "one_changed_object": runtime.get("changed_object") ==
            "velocity_multigrid_smoother_only",
        "variant_contract": tuple(runtime.get("variants", ())) == VARIANTS,
        "checkpoint_contract": runtime.get("maximum_MG_cycles") == 64,
        "tight_solves": tight_pass and runtime.get("tight_solve_count") == selected,
        "trajectory_rows": rows_valid,
        "production_freeze":
            runtime.get("production_default_change_authorized") is False,
    }
    structural_pass = all(structural_gates.values())

    final_rows = {variant: [by_key[variant, mode, 64]
                            for mode in range(1, selected + 1)]
                  for variant in VARIANTS}
    metrics = {}
    for variant in VARIANTS:
        for field in numeric_fields:
            metrics[f"{variant}_{field}_RMS_64"] = weighted_rms(
                final_rows[variant], field)
        rows32 = [by_key[variant, mode, 32]
                  for mode in range(1, selected + 1)]
        r32 = weighted_rms(rows32, "momentum_residual_relative")
        r64 = metrics[f"{variant}_momentum_residual_relative_RMS_64"]
        metrics[f"{variant}_asymptotic_contraction_per_cycle_32_64"] = \
            (r64 / max(r32, 1.0e-300)) ** (1.0 / 32.0)
    configured = metrics[
        "CONFIGURED_VANKA_Schur_action_error_relative_RMS_64"]
    rollback = metrics[
        "DIAGONAL_ROLLBACK_Schur_action_error_relative_RMS_64"]
    metrics["configured_to_rollback_Schur_error_ratio_64"] = \
        configured / max(rollback, 1.0e-300)

    if not structural_pass:
        decision = "F2_INVALID_STRUCTURAL_EVIDENCE"
        next_task = "REPEAT_F2_VALID_EXPERIMENT"
    elif configured <= 0.8 * rollback:
        decision = "F2_VANKA_EFFECTIVE_BUT_AVV_CONVERGENCE_LIMITED"
        next_task = "F2B_COARSE_GRID_TRANSFER_LOCALIZATION"
    elif rollback <= 0.8 * configured:
        decision = "F2_VANKA_SMOOTHER_REGRESSION"
        next_task = "F2B_VANKA_SMOOTHER_REDESIGN"
    else:
        decision = "F2_AVV_SMOOTHER_INSENSITIVE"
        next_task = "F2B_COARSE_GRID_TRANSFER_LOCALIZATION"
    result = {
        "schema": "strict-ala-stage-F2-decision-v1",
        "experiment_evidence_valid": structural_pass,
        "candidate": "configured_vanka_vs_diagonal_rollback",
        "metrics": metrics,
        "gates": structural_gates,
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
        "schema": "strict-ala-stage-F2-final-audit-v1",
        "experiment_valid": valid,
        "decision": decision.get("decision") if valid else "INVALID_EXPERIMENT",
        "provenance": provenance,
        "F2": decision,
        "production_default_change_authorized": False,
        "next_authorized_task": (decision.get("next_authorized_task") if valid
                                 else "REPEAT_F2_VALID_EXPERIMENT"),
    }
    write_json(args.output, result)
    return 0


def parser():
    root = argparse.ArgumentParser()
    commands = root.add_subparsers(required=True)
    command = commands.add_parser("analyze")
    for name in ("f1d-decision", "trajectory", "tight", "runtime", "output"):
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
