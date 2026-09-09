#!/usr/bin/env python3
"""Fail-closed analysis for the bounded F3 current-plateau POD refresh."""

import argparse
import csv
import hashlib
import json
import math
import pathlib
import statistics


def load_json(path):
    return json.loads(pathlib.Path(path).read_text())


def rows(path):
    with pathlib.Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def digest(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def close(a, b, rel, absolute):
    return abs(a - b) <= absolute + rel * max(abs(a), abs(b))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--thresholds", required=True)
    parser.add_argument("--raw", required=True)
    parser.add_argument("--spectrum", required=True)
    parser.add_argument("--modes", required=True)
    parser.add_argument("--trajectory", required=True)
    parser.add_argument("--overlap", required=True)
    parser.add_argument("--current-iterations", required=True)
    parser.add_argument("--source-iterations", required=True)
    parser.add_argument("--decision", required=True)
    parser.add_argument("--audit", required=True)
    args = parser.parse_args()

    thresholds = load_json(args.thresholds)
    raw = load_json(args.raw)
    spectrum = rows(args.spectrum)
    modes = rows(args.modes)
    trajectory = rows(args.trajectory)
    overlap = rows(args.overlap)
    current = rows(args.current_iterations)
    source = rows(args.source_iterations)
    errors = []

    if thresholds.get("production_default_change_authorized") is not False:
        errors.append("threshold production freeze missing")
    if raw.get("reference_operator_actions_performed") is not False:
        errors.append("reference actions were unexpectedly performed")
    if raw.get("plateau_window") != [20, 40] or raw.get("snapshot_count") != 21:
        errors.append("plateau window contract failed")
    if len(current) != 60 or len(source) != 60:
        errors.append("trajectory length mismatch")

    integer_fields = [
        "iteration", "restart_cycle", "cumulative_inner_solves",
        "cumulative_inner_cycles", "cumulative_K_gamma_applications",
        "cumulative_schur_actions", "cumulative_preconditioner_applications",
        "restart_boundary", "best_iterate", "final_iterate",
    ]
    numeric_fields = [
        "krylov_recursive", "krylov_explicit", "krylov_drift",
        "continuity_numerator", "continuity_denominator", "continuity_relative",
        "momentum_numerator", "momentum_denominator", "momentum_relative",
        "momentum_rms",
    ]
    rel = float(thresholds["trajectory_numeric_relative_tolerance"])
    absolute = float(thresholds["trajectory_numeric_absolute_tolerance"])
    if len(current) == len(source):
        for index, (observed, expected) in enumerate(zip(current, source), 1):
            for field in integer_fields:
                if int(observed[field]) != int(expected[field]):
                    errors.append("trajectory integer mismatch row %d %s" % (index, field))
            for field in numeric_fields:
                a, b = float(observed[field]), float(expected[field])
                if not (math.isfinite(a) and math.isfinite(b) and close(a, b, rel, absolute)):
                    errors.append("trajectory numeric mismatch row %d %s" % (index, field))

    selected = [r for r in spectrum if int(r["selected"])]
    if len(selected) != int(raw.get("selected_mode_count", -1)):
        errors.append("selected mode count mismatch")
    cumulative = sum(float(r["energy_fraction"]) for r in selected)
    if not close(cumulative, float(raw.get("selected_cumulative_energy", -1)), 1e-12, 1e-14):
        errors.append("selected cumulative energy mismatch")
    if len(selected) > int(thresholds["selected_mode_hard_cap"]):
        errors.append("selected mode cap exceeded")
    if any(not r["checksum"] for r in selected):
        errors.append("selected mode checksum missing")
    if float(raw.get("orthogonality_max_abs_defect", math.inf)) > \
            thresholds["pod_orthogonality_max_abs_defect"]:
        errors.append("refreshed POD orthogonality failed")

    by_iteration = {}
    for row in trajectory:
        by_iteration.setdefault(int(row["iteration"]), []).append(row)
    poke = sorted(by_iteration)
    if poke != list(range(20, 41)):
        errors.append("refresh trajectory window incomplete")
    f_values = []
    for iteration in poke:
        values = {float(r["new_POD_f_Q"]) for r in by_iteration[iteration]}
        if len(values) != 1:
            errors.append("inconsistent f_Q at iteration %d" % iteration)
        else:
            f_values.append(values.pop())
    median_f = statistics.median(f_values) if f_values else float("nan")
    final_f = f_values[-1] if f_values else float("nan")
    if not close(median_f, float(raw.get("median_f_Q", -1)), 1e-12, 1e-14):
        errors.append("median f_Q mismatch")
    if not close(final_f, float(raw.get("f_Q_at_40", -1)), 1e-12, 1e-14):
        errors.append("iteration-40 f_Q mismatch")

    weight_sum = sum(float(r["energy_fraction"]) for r in modes)
    weighted_y = math.sqrt(sum(float(r["energy_fraction"]) *
                               float(r["E_Y_at_40"]) ** 2 for r in modes) /
                           max(weight_sum, 1e-300))
    if not close(weighted_y, float(raw.get("weighted_E_Y_at_40", -1)), 1e-12, 1e-14):
        errors.append("weighted Y reachability mismatch")
    dominant = [r for r in modes if int(r["dominant"])]
    all_dominant_reachable = bool(dominant) and all(
        float(r["E_Y_at_40"]) <= thresholds["mode_Y_reachability_error_maximum"]
        for r in dominant)
    y_reachable = all_dominant_reachable and weighted_y <= thresholds[
        "weighted_Y_reachability_error_maximum"]
    dominant_slow = any(r["mode_status"] in ("STAGNATING", "WEAK_DECAY")
                        for r in dominant)
    compact = cumulative >= thresholds["selected_cumulative_energy_minimum"]
    explained = (compact and
                 median_f >= thresholds["plateau_explanation_median_minimum"] and
                 final_f >= thresholds["plateau_explanation_at_iteration_40_minimum"])
    if not explained:
        direction = "RESIDUAL_SUBSPACE_NOT_LOW_DIMENSIONAL_ENOUGH"
    elif dominant_slow and not y_reachable:
        direction = "GLOBAL_REACHABILITY_OR_DEFLATION_DESIGN"
    else:
        direction = "SCHUR_PRECONDITIONER_OR_OPERATOR_REDESIGN"
    if raw.get("preliminary_next_direction") != direction:
        errors.append("C/Python decision mismatch")
    if bool(raw.get("well_explained")) != explained:
        errors.append("C/Python explanation mismatch")
    if bool(raw.get("Y_reachable_at_40")) != y_reachable:
        errors.append("C/Python Y gate mismatch")
    old_mode_count = len({r["old_E2_mode"] for r in overlap})
    if old_mode_count < 1 or len(overlap) != len(modes) * old_mode_count:
        errors.append("old/new overlap matrix is incomplete")
    old_max_overlap = {}
    new_max_overlap = {}
    for row in overlap:
        old_id = int(row["old_E2_mode"])
        new_id = int(row["new_F3_mode"])
        value = float(row["absolute_global_pdot_overlap"])
        old_max_overlap[old_id] = max(old_max_overlap.get(old_id, 0.0), value)
        new_max_overlap[new_id] = max(new_max_overlap.get(new_id, 0.0), value)
    surviving_old_modes = sorted(mode for mode, value in old_max_overlap.items()
        if value >= thresholds["old_mode_survival_overlap_minimum"])
    dominant_ids = {int(r["mode_id"]) for r in dominant}
    important_new_modes_missing_from_E2 = sorted(mode for mode in dominant_ids
        if new_max_overlap.get(mode, 0.0) <=
           thresholds["new_mode_missing_overlap_maximum"])
    overlap_interpretation = (
        "NEW_IMPORTANT_MODES_MISSING_FROM_E2"
        if important_new_modes_missing_from_E2 else
        "OLD_DOMINANT_MODES_SURVIVE"
        if surviving_old_modes else
        "DOMINANT_SUBSPACE_MATERIALLY_ROTATED")

    valid = not errors
    decision = {
        "schema": "strict-ala-stage-F3-pod-refresh-decision-v1",
        "experiment_valid": valid,
        "validation_errors": errors,
        "plateau_window": [20, 40],
        "pressure_inner_product": "global_pdot_euclidean_pressure",
        "selected_mode_count": len(selected),
        "selected_cumulative_energy": cumulative,
        "median_f_Q_20_40": median_f,
        "f_Q_40": final_f,
        "plateau_well_explained": explained,
        "dominant_slow_mode_present": dominant_slow,
        "weighted_E_Y_40": weighted_y,
        "Y_reachable_at_40": y_reachable,
        "old_mode_max_overlap": old_max_overlap,
        "new_mode_max_overlap": new_max_overlap,
        "surviving_old_E2_modes": surviving_old_modes,
        "important_new_modes_missing_from_E2": important_new_modes_missing_from_E2,
        "old_new_overlap_interpretation": overlap_interpretation,
        "NEXT_DIRECTION": direction if valid else "INVALID_EXPERIMENT",
        "production_trajectory_reproduced": not any(
            e.startswith("trajectory") for e in errors),
        "new_HPC_run_was_required": True,
        "reference_operator_forensics_reopened": False,
        "production_default_change_authorized": False,
    }
    pathlib.Path(args.decision).write_text(json.dumps(decision, indent=2, sort_keys=True) + "\n")
    artifacts = [args.thresholds, args.raw, args.spectrum, args.modes,
                 args.trajectory, args.overlap, args.current_iterations,
                 args.source_iterations, args.decision]
    audit = {
        "schema": "strict-ala-stage-F3-pod-refresh-final-audit-v1",
        "complete": True,
        "valid": valid,
        "decision": decision["NEXT_DIRECTION"],
        "artifact_sha256": {str(pathlib.Path(p).resolve()): digest(p) for p in artifacts},
        "production_default_change_authorized": False,
    }
    pathlib.Path(args.audit).write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    if not valid:
        raise SystemExit("F3 POD refresh invalid: " + "; ".join(errors[:8]))


if __name__ == "__main__":
    main()
