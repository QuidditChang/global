#!/usr/bin/env python3
"""Fail-closed Stage-F2d finest-level smoother candidate analysis."""

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

CHECKPOINTS = (1, 2, 4, 8, 16, 32, 64)
CANDIDATES = ("CONFIGURED", "L4_COLORED_FORWARD",
              "L4_COLORED_ALTERNATING")
COMPOSITIONS = {"CONFIGURED": "configured",
                "L4_COLORED_FORWARD": "colored_forward",
                "L4_COLORED_ALTERNATING": "colored_alternating"}


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
    lineage = load_json(args.f2c_decision)
    runtime = load_json(args.runtime)
    trajectory = read_csv(args.trajectory, (
        "candidate", "smoother_composition", "finest_sweeps",
        "finest_damping", "POD_mode", "POD_energy_weight", "MG_cycles",
        "momentum_residual_relative", "valid"))
    levels = read_csv(args.level_transfer, (
        "candidate", "POD_mode", "POD_energy_weight", "MG_cycle",
        "top_level", "level", "phase", "input_rms", "output_rms",
        "reduction_ratio", "alpha", "valid"))
    mode_rhs = read_csv(args.mode_rhs,
                        ("POD_mode", "POD_energy_weight", "rhs_norm", "valid"))
    performance = read_csv(args.performance, (
        "candidate", "smoother_composition", "finest_sweeps",
        "finest_damping", "colors_per_sweep", "seconds_max", "valid"))
    f2c_trajectory = read_csv(args.f2c_trajectory, (
        "candidate", "finest_sweeps", "finest_damping", "POD_mode",
        "POD_energy_weight", "MG_cycles", "momentum_residual_relative",
        "valid"))
    f2c_levels = read_csv(args.f2c_levels, (
        "candidate", "POD_mode", "POD_energy_weight", "MG_cycle",
        "top_level", "level", "phase", "input_rms", "output_rms",
        "reduction_ratio", "alpha", "valid"))
    f2c_rhs = read_csv(args.f2c_rhs,
                       ("POD_mode", "POD_energy_weight", "rhs_norm", "valid"))

    selected = runtime.get("selected_modes")
    if not isinstance(selected, int) or selected <= 0:
        raise ValueError("invalid selected mode count")
    expected = {(candidate, mode, cycle) for candidate in CANDIDATES
                for mode in range(1, selected + 1) for cycle in CHECKPOINTS}
    keys = {(row["candidate"], int(row["POD_mode"]), int(row["MG_cycles"]))
            for row in trajectory}
    trajectory_complete = len(trajectory) == len(expected) and keys == expected
    rows_valid = all(row["valid"] == "1" and row["candidate"] in CANDIDATES
                     and row["smoother_composition"] ==
                         COMPOSITIONS[row["candidate"]]
                     and finite(row["momentum_residual_relative"], "residual") >= 0
                     and int(row["finest_sweeps"]) >= 1
                     and 0 < finite(row["finest_damping"], "damping") <= 1
                     for row in trajectory)

    baseline = {(int(r["POD_mode"]), int(r["MG_cycles"])): r for r in trajectory
                if r["candidate"] == "CONFIGURED"}
    prior = {(int(r["POD_mode"]), int(r["MG_cycles"])): r
             for r in f2c_trajectory if r["candidate"] == "CONFIGURED"}
    baseline_identity = set(baseline) == set(prior)
    if baseline_identity:
        baseline_identity = all(math.isclose(
            finite(baseline[key]["momentum_residual_relative"], "F2d baseline"),
            finite(prior[key]["momentum_residual_relative"], "F2c baseline"),
            rel_tol=1e-10, abs_tol=1e-14) for key in baseline)

    baseline_levels = [r for r in levels if r["candidate"] == "CONFIGURED"]
    prior_levels = [r for r in f2c_levels if r["candidate"] == "CONFIGURED"]
    level_identity = len(baseline_levels) == len(prior_levels)
    level_keys = ("POD_mode", "MG_cycle", "top_level", "level", "phase", "valid")
    level_numbers = ("POD_energy_weight", "input_rms", "output_rms",
                     "reduction_ratio", "alpha")
    if level_identity:
        for current, old in zip(baseline_levels, prior_levels):
            level_identity &= all(current[key] == old[key] for key in level_keys)
            level_identity &= all(math.isclose(
                finite(current[key], f"F2d {key}"),
                finite(old[key], f"F2c {key}"), rel_tol=1e-10, abs_tol=1e-14)
                for key in level_numbers)

    rhs_identity = len(mode_rhs) == len(f2c_rhs) == selected
    if rhs_identity:
        current = {int(r["POD_mode"]): r for r in mode_rhs}
        old = {int(r["POD_mode"]): r for r in f2c_rhs}
        rhs_identity = set(current) == set(old) and all(math.isclose(
            finite(current[k]["rhs_norm"], "F2d RHS"),
            finite(old[k]["rhs_norm"], "F2c RHS"),
            rel_tol=1e-12, abs_tol=1e-300) for k in current)

    level_valid = all(row["valid"] == "1" and row["candidate"] in CANDIDATES
                      and int(row["POD_mode"]) in range(1, selected + 1)
                      and int(row["MG_cycle"]) in CHECKPOINTS
                      and finite(row["input_rms"], "input") >= 0
                      and finite(row["output_rms"], "output") >= 0
                      and finite(row["reduction_ratio"], "ratio") >= 0
                      for row in levels)
    max_top = max(int(row["top_level"]) for row in levels)
    final_l4_down = [row for row in levels if row["phase"] == "DOWN_SMOOTH"
                     and int(row["level"]) == max_top
                     and int(row["top_level"]) == max_top
                     and int(row["MG_cycle"]) == 64]
    l4_coverage = {(r["candidate"], int(r["POD_mode"])) for r in final_l4_down} == {
        (candidate, mode) for candidate in CANDIDATES
        for mode in range(1, selected + 1)}

    lineage_ok = (lineage.get("experiment_evidence_valid") is True and
                  lineage.get("decision") ==
                  "F2C_NO_ROBUST_L4_SMOOTHER_CANDIDATE" and
                  lineage.get("next_authorized_task") ==
                  "F2C_BLOCK_STRUCTURE_REDESIGN" and
                  lineage.get("production_default_change_authorized") is False)
    configured_sweeps = runtime.get("configured_finest_sweeps")
    configured_damping = finite(runtime.get("configured_finest_damping"),
                                "configured damping")
    sweeps_frozen = (isinstance(configured_sweeps, int) and
                     configured_sweeps > 0 and all(
                         int(row["finest_sweeps"]) == configured_sweeps
                         for row in trajectory))
    damping_frozen = all(math.isclose(
        finite(row["finest_damping"], "candidate damping"),
        configured_damping, rel_tol=0.0, abs_tol=1e-15)
        for row in trajectory)
    performance_by_candidate = {row["candidate"]: row for row in performance}
    performance_valid = (len(performance) == len(CANDIDATES) and
        set(performance_by_candidate) == set(CANDIDATES) and all(
            row["valid"] == "1" and
            row["smoother_composition"] == COMPOSITIONS[row["candidate"]] and
            int(row["colors_per_sweep"]) ==
                (1 if row["candidate"] == "CONFIGURED" else 8) and
            finite(row["seconds_max"], "candidate seconds") > 0.0
            for row in performance))
    gates = {
        "F2c_lineage": lineage_ok,
        "runtime_schema": runtime.get("schema") ==
            "strict-ala-stage-F2d-runtime-v1",
        "frozen_current_viscosity": runtime.get("viscosity") == "frozen_current",
        "frozen_operator": runtime.get("operator") ==
            "frozen_current_K_gamma_D_plus_C",
        "finest_level_only": runtime.get("scope") == "finest_level_only",
        "candidate_count": runtime.get("candidate_count") == len(CANDIDATES),
        "colored_structure": runtime.get("candidate_structure") ==
            "eight_color_multiplicative_residual_update",
        "cached_blocks_frozen": runtime.get("cached_blocks") ==
            "unchanged_24x24_element_K_gamma",
        "color_count": runtime.get("colors") == 8,
        "finest_sweeps_frozen": sweeps_frozen,
        "finest_damping_frozen": damping_frozen,
        "candidate_performance_rows": performance_valid,
        "selected_energy": finite(runtime.get("selected_energy"),
                                   "selected energy") >= 0.95,
        "trajectory_complete": trajectory_complete and rows_valid,
        "level_rows": level_valid and l4_coverage,
        "baseline_identity": baseline_identity,
        "baseline_level_identity": level_identity,
        "mode_rhs_identity": rhs_identity,
        "production_freeze":
            runtime.get("production_default_change_authorized") is False,
    }
    structural_pass = all(gates.values())

    metrics = {}
    for candidate in CANDIDATES:
        final = [r for r in trajectory if r["candidate"] == candidate and
                 int(r["MG_cycles"]) == 64]
        at32 = [r for r in trajectory if r["candidate"] == candidate and
                int(r["MG_cycles"]) == 32]
        l4 = [r for r in final_l4_down if r["candidate"] == candidate]
        residual64 = weighted_rms(final, "momentum_residual_relative")
        residual32 = weighted_rms(at32, "momentum_residual_relative")
        metrics[candidate] = {
            "smoother_composition": final[0]["smoother_composition"],
            "finest_sweeps": int(final[0]["finest_sweeps"]),
            "finest_damping": finite(final[0]["finest_damping"], "damping"),
            "momentum_residual_RMS_32": residual32,
            "momentum_residual_RMS_64": residual64,
            "effective_cycle_factor_32_64":
                (residual64 / max(residual32, 1e-300)) ** (1.0 / 32.0),
            "colors_per_sweep": int(
                performance_by_candidate[candidate]["colors_per_sweep"]),
            "seconds_max": finite(
                performance_by_candidate[candidate]["seconds_max"], "seconds"),
            "L4_down_reduction_RMS_64": weighted_rms(l4, "reduction_ratio"),
            "per_mode_residual_64": {
                str(int(r["POD_mode"])):
                    finite(r["momentum_residual_relative"], "residual")
                for r in final},
        }
    base = metrics["CONFIGURED"]
    for item in metrics.values():
        item["wall_time_ratio_to_configured"] = item["seconds_max"] / max(
            base["seconds_max"], 1e-300)
    qualified = []
    for candidate in CANDIDATES[1:]:
        item = metrics[candidate]
        residual_gain = item["momentum_residual_RMS_64"] / max(
            base["momentum_residual_RMS_64"], 1e-300)
        contraction_gain = item["effective_cycle_factor_32_64"] / max(
            base["effective_cycle_factor_32_64"], 1e-300)
        no_mode_regression = all(
            value <= 1.02 * base["per_mode_residual_64"][mode]
            for mode, value in item["per_mode_residual_64"].items())
        item.update({"residual_ratio_to_configured": residual_gain,
                     "asymptotic_factor_ratio_to_configured": contraction_gain,
                     "no_mode_regression": no_mode_regression,
                     "candidate_gate_pass": residual_gain <= 0.50 and
                         contraction_gain <= 0.98 and no_mode_regression})
        if item["candidate_gate_pass"]:
            qualified.append(candidate)

    if not structural_pass:
        decision = "F2D_INVALID_STRUCTURAL_EVIDENCE"
        next_task = "REPEAT_F2D_VALID_EXPERIMENT"
        selected_candidate = None
    elif qualified:
        selected_candidate = min(
            qualified, key=lambda name: metrics[name]["momentum_residual_RMS_64"])
        decision = "F2D_L4_SMOOTHER_CANDIDATE_IDENTIFIED"
        next_task = "F2E_PRODUCTION_VELOCITY_AB"
    else:
        selected_candidate = None
        decision = "F2D_NO_ROBUST_L4_SMOOTHER_CANDIDATE"
        next_task = "F2E_FACE_PATCH_BLOCK_REDESIGN"
    write_json(args.output, {
        "schema": "strict-ala-stage-F2d-decision-v1",
        "experiment_evidence_valid": structural_pass,
        "gates": gates, "metrics": metrics,
        "qualified_candidates": qualified,
        "selected_candidate": selected_candidate,
        "decision": decision,
        "production_default_change_authorized": False,
        "next_authorized_task": next_task,
    })
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
    write_json(args.output, {
        "schema": "strict-ala-stage-F2d-final-audit-v1",
        "experiment_valid": valid,
        "decision": decision.get("decision") if valid else "INVALID_EXPERIMENT",
        "provenance": provenance, "F2d": decision,
        "production_default_change_authorized": False,
        "next_authorized_task": decision.get("next_authorized_task") if valid
            else "REPEAT_F2D_VALID_EXPERIMENT",
    })
    return 0


def parser():
    root = argparse.ArgumentParser()
    commands = root.add_subparsers(required=True)
    command = commands.add_parser("analyze")
    for name in ("f2c-decision", "f2c-trajectory", "f2c-levels", "f2c-rhs",
                 "level-transfer", "trajectory", "mode-rhs", "performance",
                 "runtime", "output"):
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
