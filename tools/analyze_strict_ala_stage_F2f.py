#!/usr/bin/env python3
"""Fail-closed Stage-F2f overlapping L4 face-patch analysis."""

import argparse
import csv
import hashlib
import importlib.util
import math
from pathlib import Path


_SPEC = importlib.util.spec_from_file_location(
    "strict_ala_f2e_common", Path(__file__).with_name(
        "analyze_strict_ala_stage_F2e.py"))
COMMON = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(COMMON)

CHECKPOINTS = COMMON.CHECKPOINTS
CANDIDATES = ("CONFIGURED", "L4_OVERLAP_FACE_Z", "L4_OVERLAP_FACE_X",
              "L4_OVERLAP_FACE_Y")
BLOCKS = {
    "CONFIGURED": "configured",
    "L4_OVERLAP_FACE_Z": "overlap_face_z",
    "L4_OVERLAP_FACE_X": "overlap_face_x",
    "L4_OVERLAP_FACE_Y": "overlap_face_y",
}


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


def analyze(args):
    lineage = COMMON.load_json(args.f2e_decision)
    runtime = COMMON.load_json(args.runtime)
    trajectory_columns = (
        "candidate", "smoother_composition", "finest_sweeps",
        "finest_damping", "POD_mode", "POD_energy_weight", "MG_cycles",
        "momentum_residual_relative", "valid")
    level_columns = (
        "candidate", "POD_mode", "POD_energy_weight", "MG_cycle",
        "top_level", "level", "phase", "input_rms", "output_rms",
        "reduction_ratio", "alpha", "valid")
    rhs_columns = ("POD_mode", "POD_energy_weight", "rhs_norm", "valid")
    trajectory = read_csv(args.trajectory, trajectory_columns)
    levels = read_csv(args.level_transfer, level_columns)
    mode_rhs = read_csv(args.mode_rhs, rhs_columns)
    performance = read_csv(args.performance, (
        "candidate", "smoother_block", "finest_sweeps",
        "finest_damping", "block_dof", "seconds_max", "valid"))
    cache = read_csv(args.cache, (
        "orientation", "patches_per_cap", "block_dof", "cache_mb_max",
        "build_seconds_max", "min_pivot_ratio", "valid"))
    prior_trajectory = read_csv(args.f2e_trajectory, trajectory_columns)
    prior_levels = read_csv(args.f2e_levels, level_columns)
    prior_rhs = read_csv(args.f2e_rhs, rhs_columns)

    selected = runtime.get("selected_modes")
    if not isinstance(selected, int) or selected <= 0:
        raise ValueError("invalid selected mode count")
    expected = {(candidate, mode, cycle) for candidate in CANDIDATES
                for mode in range(1, selected + 1) for cycle in CHECKPOINTS}
    actual = {(row["candidate"], int(row["POD_mode"]),
               int(row["MG_cycles"])) for row in trajectory}
    trajectory_valid = (
        len(trajectory) == len(expected) and actual == expected and
        all(row["valid"] == "1" and row["candidate"] in CANDIDATES and
            row["smoother_composition"] == BLOCKS[row["candidate"]] and
            int(row["finest_sweeps"]) >= 1 and
            0.0 < COMMON.finite(row["finest_damping"], "damping") <= 1.0 and
            COMMON.finite(row["momentum_residual_relative"], "residual") >= 0
            for row in trajectory))

    current_baseline = [row for row in trajectory
                        if row["candidate"] == "CONFIGURED"]
    prior_baseline = [row for row in prior_trajectory
                      if row["candidate"] == "CONFIGURED"]
    baseline_identity = COMMON.rows_identical(
        current_baseline, prior_baseline,
        ("candidate", "POD_mode", "MG_cycles", "valid"),
        ("finest_sweeps", "finest_damping", "POD_energy_weight",
         "momentum_residual_relative"), "F2f", "F2e")
    current_levels = [row for row in levels
                      if row["candidate"] == "CONFIGURED"]
    prior_baseline_levels = [row for row in prior_levels
                             if row["candidate"] == "CONFIGURED"]
    baseline_level_identity = COMMON.rows_identical(
        current_levels, prior_baseline_levels,
        ("candidate", "POD_mode", "MG_cycle", "top_level", "level",
         "phase", "valid"),
        ("POD_energy_weight", "input_rms", "output_rms",
         "reduction_ratio", "alpha"), "F2f", "F2e")
    rhs_identity = COMMON.rows_identical(
        mode_rhs, prior_rhs, ("POD_mode", "valid"),
        ("POD_energy_weight", "rhs_norm"), "F2f", "F2e")

    level_valid = all(
        row["valid"] == "1" and row["candidate"] in CANDIDATES and
        int(row["POD_mode"]) in range(1, selected + 1) and
        int(row["MG_cycle"]) in CHECKPOINTS and
        COMMON.finite(row["input_rms"], "level input") >= 0 and
        COMMON.finite(row["output_rms"], "level output") >= 0 and
        COMMON.finite(row["reduction_ratio"], "level ratio") >= 0
        for row in levels)
    max_top = max(int(row["top_level"]) for row in levels)
    final_l4 = [row for row in levels if row["phase"] == "DOWN_SMOOTH" and
                int(row["level"]) == max_top and
                int(row["top_level"]) == max_top and
                int(row["MG_cycle"]) == 64]
    l4_coverage = {
        (row["candidate"], int(row["POD_mode"])) for row in final_l4
    } == {(candidate, mode) for candidate in CANDIDATES
          for mode in range(1, selected + 1)}

    performance_by = {row["candidate"]: row for row in performance}
    performance_valid = (
        len(performance) == len(CANDIDATES) and
        set(performance_by) == set(CANDIDATES) and
        all(row["valid"] == "1" and
            row["smoother_block"] == BLOCKS[row["candidate"]] and
            int(row["block_dof"]) == (
                24 if row["candidate"] == "CONFIGURED" else 36) and
            COMMON.finite(row["seconds_max"], "seconds") > 0
            for row in performance))
    overlap_names = {"overlap_face_z", "overlap_face_x", "overlap_face_y"}
    cache_by = {row["orientation"]: row for row in cache}
    cache_valid = (
        len(cache) == 3 and set(cache_by) == overlap_names and
        all(row["valid"] == "1" and int(row["patches_per_cap"]) > 0 and
            int(row["block_dof"]) == 36 and
            0.0 < COMMON.finite(row["cache_mb_max"], "cache MB") <= 128.0 and
            COMMON.finite(row["build_seconds_max"], "build seconds") > 0 and
            COMMON.finite(row["min_pivot_ratio"], "pivot ratio") > 0
            for row in cache))

    configured_sweeps = runtime.get("configured_finest_sweeps")
    configured_damping = COMMON.finite(
        runtime.get("configured_finest_damping"), "configured damping")
    frozen_controls = (
        isinstance(configured_sweeps, int) and configured_sweeps > 0 and
        all(int(row["finest_sweeps"]) == configured_sweeps
            for row in trajectory) and
        all(math.isclose(COMMON.finite(row["finest_damping"], "damping"),
                         configured_damping, rel_tol=0.0, abs_tol=1e-15)
            for row in trajectory))
    lineage_ok = (
        lineage.get("experiment_evidence_valid") is True and
        lineage.get("decision") == "F2E_NO_ROBUST_FACE_PATCH_CANDIDATE" and
        lineage.get("next_authorized_task") ==
            "F2F_OVERLAPPING_OR_CROSS_RANK_BLOCK_REDESIGN" and
        lineage.get("production_default_change_authorized") is False)
    gates = {
        "F2e_lineage": lineage_ok,
        "runtime_schema":
            runtime.get("schema") == "strict-ala-stage-F2f-runtime-v1",
        "frozen_current_viscosity": runtime.get("viscosity") == "frozen_current",
        "frozen_operator":
            runtime.get("operator") == "frozen_current_K_gamma_D_plus_C",
        "finest_level_only": runtime.get("scope") == "finest_level_only",
        "candidate_count": runtime.get("candidate_count") == len(CANDIDATES),
        "overlapping_face_structure": runtime.get("candidate_structure") ==
            "overlapping_two_element_face_pairs_all_local_internal_faces",
        "rank_local_isolation": runtime.get("cross_rank_faces") is False,
        "single_cache_contract": runtime.get("cached_blocks") ==
            "one_orientation_at_a_time_36x36_K_gamma",
        "face_patch_dof": runtime.get("face_patch_dof") == 36,
        "frozen_sweeps_and_damping": frozen_controls,
        "cache_rows_and_memory_bound": cache_valid,
        "candidate_performance_rows": performance_valid,
        "selected_energy":
            COMMON.finite(runtime.get("selected_energy"), "energy") >= 0.95,
        "trajectory_complete": trajectory_valid,
        "level_rows": level_valid and l4_coverage,
        "baseline_identity": baseline_identity,
        "baseline_level_identity": baseline_level_identity,
        "mode_rhs_identity": rhs_identity,
        "production_freeze":
            runtime.get("production_default_change_authorized") is False,
    }
    structural_pass = all(gates.values())

    metrics = {}
    for candidate in CANDIDATES:
        final = [row for row in trajectory if row["candidate"] == candidate and
                 int(row["MG_cycles"]) == 64]
        at32 = [row for row in trajectory if row["candidate"] == candidate and
                int(row["MG_cycles"]) == 32]
        l4 = [row for row in final_l4 if row["candidate"] == candidate]
        residual64 = COMMON.weighted_rms(final, "momentum_residual_relative")
        residual32 = COMMON.weighted_rms(at32, "momentum_residual_relative")
        item = {
            "smoother_block": BLOCKS[candidate],
            "block_dof": int(performance_by[candidate]["block_dof"]),
            "finest_sweeps": int(performance_by[candidate]["finest_sweeps"]),
            "finest_damping": COMMON.finite(
                performance_by[candidate]["finest_damping"], "damping"),
            "seconds_max": COMMON.finite(
                performance_by[candidate]["seconds_max"], "seconds"),
            "momentum_residual_RMS_32": residual32,
            "momentum_residual_RMS_64": residual64,
            "effective_cycle_factor_32_64":
                (residual64 / max(residual32, 1e-300)) ** (1.0 / 32.0),
            "L4_down_reduction_RMS_64":
                COMMON.weighted_rms(l4, "reduction_ratio"),
            "per_mode_residual_64": {
                str(int(row["POD_mode"])): COMMON.finite(
                    row["momentum_residual_relative"], "residual")
                for row in final},
        }
        if candidate != "CONFIGURED":
            cache_row = cache_by[BLOCKS[candidate]]
            item.update({
                "patches_per_cap": int(cache_row["patches_per_cap"]),
                "cache_mb_max": COMMON.finite(
                    cache_row["cache_mb_max"], "cache MB"),
                "build_seconds_max": COMMON.finite(
                    cache_row["build_seconds_max"], "build seconds"),
                "min_pivot_ratio": COMMON.finite(
                    cache_row["min_pivot_ratio"], "pivot ratio"),
            })
        metrics[candidate] = item

    base = metrics["CONFIGURED"]
    qualified = []
    for candidate in CANDIDATES[1:]:
        item = metrics[candidate]
        residual_ratio = item["momentum_residual_RMS_64"] / max(
            base["momentum_residual_RMS_64"], 1e-300)
        factor_ratio = item["effective_cycle_factor_32_64"] / max(
            base["effective_cycle_factor_32_64"], 1e-300)
        no_mode_regression = all(
            value <= 1.02 * base["per_mode_residual_64"][mode]
            for mode, value in item["per_mode_residual_64"].items())
        item.update({
            "wall_time_ratio_to_configured": item["seconds_max"] / max(
                base["seconds_max"], 1e-300),
            "residual_ratio_to_configured": residual_ratio,
            "asymptotic_factor_ratio_to_configured": factor_ratio,
            "no_mode_regression": no_mode_regression,
            "candidate_gate_pass":
                residual_ratio <= 0.50 and factor_ratio <= 0.98 and
                no_mode_regression,
        })
        if item["candidate_gate_pass"]:
            qualified.append(candidate)
    base["wall_time_ratio_to_configured"] = 1.0

    if not structural_pass:
        decision = "F2F_INVALID_STRUCTURAL_EVIDENCE"
        next_task = "REPEAT_F2F_VALID_EXPERIMENT"
        selected_candidate = None
    elif qualified:
        selected_candidate = min(
            qualified, key=lambda name: metrics[name]["momentum_residual_RMS_64"])
        decision = "F2F_OVERLAPPING_FACE_PATCH_CANDIDATE_IDENTIFIED"
        next_task = "F2G_PRODUCTION_VELOCITY_AB"
    else:
        selected_candidate = None
        decision = "F2F_NO_ROBUST_OVERLAPPING_FACE_PATCH_CANDIDATE"
        next_task = "F2G_CROSS_RANK_FACE_PATCH_REDESIGN"
    COMMON.write_json(args.output, {
        "schema": "strict-ala-stage-F2f-decision-v1",
        "experiment_evidence_valid": structural_pass,
        "gates": gates,
        "metrics": metrics,
        "qualified_candidates": qualified,
        "selected_candidate": selected_candidate,
        "decision": decision,
        "production_default_change_authorized": False,
        "next_authorized_task": next_task,
    })
    return 0


def audit(args):
    decision = COMMON.load_json(args.decision)
    provenance = {}
    for name in ("binary", "inputs", "source"):
        before = Path(getattr(args, f"{name}_pre"))
        after = Path(getattr(args, f"{name}_post"))
        provenance[name] = {
            "pre_sha256": hashlib.sha256(before.read_bytes()).hexdigest(),
            "post_sha256": hashlib.sha256(after.read_bytes()).hexdigest(),
            "unchanged": before.read_bytes() == after.read_bytes(),
        }
    valid = (decision.get("experiment_evidence_valid") is True and
             all(item["unchanged"] for item in provenance.values()))
    COMMON.write_json(args.output, {
        "schema": "strict-ala-stage-F2f-final-audit-v1",
        "experiment_valid": valid,
        "decision": decision.get("decision") if valid else "INVALID_EXPERIMENT",
        "provenance": provenance,
        "F2f": decision,
        "production_default_change_authorized": False,
        "next_authorized_task": decision.get("next_authorized_task") if valid
            else "REPEAT_F2F_VALID_EXPERIMENT",
    })
    return 0


def parser():
    root = argparse.ArgumentParser()
    commands = root.add_subparsers(required=True)
    command = commands.add_parser("analyze")
    for name in (
            "f2e-decision", "f2e-trajectory", "f2e-levels", "f2e-rhs",
            "level-transfer", "trajectory", "mode-rhs", "performance",
            "cache", "runtime", "output"):
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
