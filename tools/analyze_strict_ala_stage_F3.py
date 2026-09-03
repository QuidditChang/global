#!/usr/bin/env python3
"""Fail-closed analysis for Strict-ALA F3 outer-Krylov relocalization."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np


def load_json(path):
    return json.loads(Path(path).read_text())


def load_csv(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True,
                                     allow_nan=False) + "\n")


def json_safe(value):
    """Replace raw nonfinite diagnostics by JSON null after flagging invalid."""
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_safe(item) for item in value]
    return value


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def finite(value):
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def median(values):
    values = sorted(values)
    n = len(values)
    if not n:
        return None
    return values[n // 2] if n % 2 else .5 * (values[n // 2 - 1] + values[n // 2])


def relative_difference(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def json_contains_nonfinite(value):
    if isinstance(value, float):
        return not math.isfinite(value)
    if isinstance(value, dict):
        return any(json_contains_nonfinite(item) for item in value.values())
    if isinstance(value, list):
        return any(json_contains_nonfinite(item) for item in value)
    return False


def csv_nonfinite(rows, excluded=()):
    excluded = set(excluded)
    for row in rows:
        for key, value in row.items():
            if key in excluded or value in (None, ""):
                continue
            try:
                number = float(value)
            except ValueError:
                continue
            if not math.isfinite(number):
                return True
    return False


def validate_thresholds(t):
    required = {
        "S_ref_inner_relative_tolerance": 1e-10,
        "S_ref_max_MG_cycles": 1024,
        "repeat_action_relative_tolerance": 1e-10,
        "scaling_scalar_c": .5,
        "scaling_relative_tolerance": 1e-8,
        "additivity_relative_tolerance": 1e-8,
        "modal_significance_floor_relative": 1e-4,
        "qr_first_pass_orthogonality_threshold": 1e-10,
        "qr_rank_rejection_relative_threshold": 1e-12,
        "projected_subspace_mode_leakage_limit": .25,
        "projected_subspace_weighted_leakage_limit": .20,
        "Y_mode_reachability_error_limit": .20,
        "Y_weighted_reachability_error_limit": .20,
        "POD_late_median_explanation_min": .80,
        "POD_iteration40_explanation_min": .75,
        "mode_stagnation_ratio_min": .90,
        "mode_decay_ratio_max": .70,
        "projected_condition_number_adverse": 1e3,
        "projected_mixing_adverse": .50,
        "projected_non_normality_adverse": .25,
        "sigma_zero_relative_threshold": 1e-14,
        "H_true_symmetry_relative_tolerance": 1e-8,
        "H_true_negative_curvature_relative_tolerance": 1e-12,
        "pod_reconstruction_coefficient_relative_tolerance": 1e-10,
        "projection_monotonicity_absolute_tolerance": 1e-10,
        "pod_weight_sum_absolute_tolerance": 1e-10,
        "inner_target_relative_roundoff": 1e-12,
    }
    errors = []
    for key, expected in required.items():
        if key not in t or float(t[key]) != expected:
            errors.append(f"{key} must equal {expected!r}")
    if t.get("additivity_mode_pairs") != [[1, 2], [1, 3], [2, 3]]:
        errors.append("additivity_mode_pairs must be q1/q2, q1/q3, q2/q3")
    if t.get("numerical_floor_definition") != \
            "64*DBL_EPSILON*max(reference_action_norm_scale,1)":
        errors.append("numerical_floor_definition mismatch")
    if t.get("qr_policy") != "deterministic_always_two_pass_batched_MGS_global_pdot":
        errors.append("qr_policy mismatch")
    if t.get("production_default_change_authorized") is not False:
        errors.append("production authorization must remain false")
    if errors:
        raise ValueError("; ".join(errors))


def matrix_from_rows(rows, name):
    selected = [r for r in rows if r["matrix"] == name]
    if not selected:
        raise ValueError(f"missing projected matrix {name}")
    n = max(max(int(r["row_mode"]), int(r["column_mode"])) for r in selected)
    matrix = np.full((n, n), np.nan)
    for row in selected:
        matrix[int(row["row_mode"]) - 1, int(row["column_mode"]) - 1] = float(row["value"])
    return matrix


def projected_metrics(matrix, zero_relative):
    singular = np.linalg.svd(matrix, compute_uv=False)
    sigma_max = float(singular[0])
    sigma_min = float(singular[-1])
    infinite = sigma_max == 0.0 or sigma_min <= zero_relative * sigma_max
    kappa = None if infinite else sigma_max / sigma_min
    norm = float(np.linalg.norm(matrix, "fro"))
    off = matrix - np.diag(np.diag(matrix))
    mixing = float(np.linalg.norm(off, "fro") / max(norm, 1e-300))
    departure = matrix.T @ matrix - matrix @ matrix.T
    non_normal = float(np.linalg.norm(departure, "fro") / max(norm * norm, 1e-300))
    symmetric_defect = float(np.linalg.norm(matrix - matrix.T, "fro") / max(norm, 1e-300))
    eigenvalues = np.linalg.eigvals(matrix)
    return {
        "singular_values": [float(x) for x in singular],
        "sigma_max": sigma_max,
        "sigma_min": sigma_min,
        "kappa_2": kappa,
        "kappa_2_is_infinite": infinite,
        "kappa_2_reason": ("SIGMA_MAX_ZERO" if sigma_max == 0.0 else
                           "SIGMA_MIN_BELOW_RELATIVE_ZERO_THRESHOLD") if infinite else "FINITE",
        "eta_mix": mixing,
        "eta_non_normal": non_normal,
        "symmetry_relative_defect": symmetric_defect,
        "eigenvalues": [{"real": float(x.real), "imag": float(x.imag)}
                        for x in eigenvalues],
        "right_singular_vectors": np.linalg.svd(matrix)[2].T.tolist(),
    }


def cfg_contract(args):
    thresholds = load_json(args.thresholds)
    validate_thresholds(thresholds)
    text = Path(args.cfg).read_text()
    required = {
        "ala_inner_accuracy_max": "1e-2",
        "ala_outer_solver": "fgmres",
        "ala_pcg_restart_interval": "50",
        "piterations": "60",
    }
    observed = {}
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if "=" in line:
            key, value = [x.strip() for x in line.split("=", 1)]
            if key in required:
                observed[key] = value
    errors = [f"{key}={observed.get(key)!r}, expected {value!r}"
              for key, value in required.items() if observed.get(key) != value]
    result = {"schema": "strict-ala-stage-F3-local-contract-v1",
              "valid": not errors, "errors": errors,
              "frozen_values": thresholds,
              "observed_cfg": observed,
              "production_default_change_authorized": False}
    write_json(args.output, result)
    return 0 if not errors else 1


def analyze(args):
    t = load_json(args.thresholds)
    validate_thresholds(t)
    runtime = load_json(args.runtime)
    pod = load_json(args.pod_reconstruction)
    reference = load_json(args.reference_validation_raw)
    completion = load_json(args.completion)
    explain = load_csv(args.pod_explanation)
    modes = load_csv(args.mode_coefficients)
    subspace = load_csv(args.cumulative_subspace)
    rank_rows = load_csv(args.basis_rank)
    inner_rows = load_csv(args.inner_solves)
    matrix_rows = load_csv(args.projected_matrices)
    reconstructed = load_csv(args.reconstructed_coefficients)
    authoritative = load_csv(args.authoritative_coefficients)
    baseline = load_csv(args.baseline_iterations)
    errors = []

    numeric_tables = (explain, modes, subspace, rank_rows, inner_rows,
                      matrix_rows, reconstructed, authoritative, baseline)
    if any(csv_nonfinite(rows) for rows in numeric_tables):
        errors.append("raw CSV diagnostic contains NaN/Inf")
    if json_contains_nonfinite(runtime) or json_contains_nonfinite(pod) or \
            json_contains_nonfinite(reference) or json_contains_nonfinite(completion):
        errors.append("raw JSON diagnostic contains NaN/Inf")

    if not completion.get("complete") or not completion.get("valid") or completion.get("exit_status") != 0:
        errors.append("numerical completion invalid")
    if runtime.get("outer_preconditioning_orientation") != "RIGHT_FGMRES":
        errors.append("outer orientation mismatch")
    iterations = runtime.get("iterations")
    early_joint = (isinstance(iterations, int) and iterations < 60 and
                   runtime.get("joint_target_reached") is True)
    if (not isinstance(iterations, int) or iterations < 1 or iterations > 60 or
            (iterations < 60 and not early_joint)):
        errors.append("neither frozen 60-iteration trajectory nor genuine early joint convergence")
    if not pod.get("valid") or pod.get("selected_modes", 0) < 3:
        errors.append("POD reconstruction invalid")
    inner_targets_pass = bool(inner_rows) and all(
        row.get("status") == "CONVERGED" and
        finite(row.get("achieved_relative")) and
        finite(row.get("requested_relative_tolerance")) and
        float(row["achieved_relative"]) <=
        float(row["requested_relative_tolerance"]) *
        (1.0 + t["inner_target_relative_roundoff"])
        for row in inner_rows)
    if not inner_targets_pass:
        errors.append("required inner solve target not achieved")
    if abs(float(pod.get("pod_weight_sum", 0.0)) -
           sum(float(r["mode_energy_weight"]) for r in modes
               if int(r["iteration"]) == 0)) > t["pod_weight_sum_absolute_tolerance"]:
        errors.append("POD weight sum inconsistency")

    auth = {(r["case"], int(r["iteration"]), int(r["pod_mode"])):
            float(r["coefficient"]) for r in authoritative}
    reconstruction_max = 0.0
    for row in reconstructed:
        key = (row["case"], int(row["iteration"]), int(row["pod_mode"]))
        if key not in auth:
            errors.append(f"missing authoritative POD coefficient {key}")
            continue
        reconstruction_max = max(reconstruction_max,
            relative_difference(float(row["coefficient"]), auth[key]))
    if len(reconstructed) != len(auth):
        errors.append("POD coefficient trajectory length mismatch")
    if reconstruction_max > t["pod_reconstruction_coefficient_relative_tolerance"]:
        errors.append("POD coefficient reconstruction mismatch")

    # The established reproduction envelope is exactly zero for these fields.
    base = {int(r["iteration"]): r for r in baseline}
    reproduction_max = 0.0
    for row in explain:
        iteration = int(row["iteration"])
        if iteration == 0:
            continue
        if iteration not in base:
            errors.append(f"baseline missing iteration {iteration}")
            continue
        for field in ("continuity_relative", "momentum_relative"):
            reproduction_max = max(reproduction_max,
                relative_difference(float(row[field]), float(base[iteration][field])))
    reproduction_pass = reproduction_max <= t["base_reproduction_allowed_relative_difference"]
    if not reproduction_pass:
        errors.append("BASE reproduction gate failed")

    for row in rank_rows:
        numeric = ("pre_norm", "post_first_pass_norm", "post_second_pass_norm",
                   "first_pass_relative_orthogonality")
        if any(not finite(row[x]) for x in numeric):
            errors.append("nonfinite QR diagnostic")
            break
        if int(row["second_pass_performed"]) != 1:
            errors.append("non-deterministic QR pass count")
            break

    # All-history projection errors must not increase beyond a roundoff guard.
    for mode_id in sorted({int(r["mode_id"]) for r in subspace if r["row_type"] == "MODE"}):
        rows = sorted((r for r in subspace if r["row_type"] == "MODE" and
                       int(r["mode_id"]) == mode_id), key=lambda r: int(r["iteration"]))
        for field in ("E_V_all", "E_Z_all", "E_Y_all"):
            values = [float(r[field]) for r in rows]
            if any(b > a + t["projection_monotonicity_absolute_tolerance"]
                   for a, b in zip(values, values[1:])):
                errors.append(f"nonmonotone {field} mode {mode_id}")

    mode_by_iteration = {}
    for row in modes:
        mode_by_iteration.setdefault(int(row["mode_id"]), {})[int(row["iteration"])] = row
    residual_norm = {int(r["iteration"]): float(r["residual_global_pdot_norm"])
                     for r in explain}
    significant = []
    mode_summary = []
    for mode_id, rows in sorted(mode_by_iteration.items()):
        weight = float(next(iter(rows.values()))["mode_energy_weight"])
        required = set(range(20, 30)) | set(range(31, 41))
        if not required.issubset(rows) or not set(range(20, 30)).issubset(residual_norm):
            pre, late, rpre = None, None, None
            status, ratio = "UNAVAILABLE", None
        else:
            pre = median([float(rows[k]["absolute_modal_amplitude"])
                          for k in range(20, 30)])
            late = median([float(rows[k]["absolute_modal_amplitude"])
                           for k in range(31, 41)])
            rpre = median([residual_norm[k] for k in range(20, 30)])
        if pre is not None and pre < t["modal_significance_floor_relative"] * rpre:
            status, ratio = "BELOW_SIGNIFICANCE_FLOOR", None
        elif pre is not None:
            ratio = late / max(pre, 1e-300)
            if ratio >= t["mode_stagnation_ratio_min"]:
                status = "STAGNATING"
            elif ratio <= t["mode_decay_ratio_max"]:
                status = "DECAYING"
            else:
                status = "WEAK_DECAY"
        if weight >= t["authoritative_mode_energy_min"] and status not in \
                ("BELOW_SIGNIFICANCE_FLOOR", "UNAVAILABLE"):
            significant.append((mode_id, status))
        mode_summary.append({"mode_id": mode_id, "energy_weight": weight,
                             "A_pre": pre, "R_pre": rpre,
                             "rho_mode": ratio, "mode_status": status})

    pod_endpoint_iteration = 40 if iterations >= 40 else max(iterations - 1, 0)
    pod_window = [float(r["f_Q"]) for r in explain
                  if 20 <= int(r["iteration"]) <= pod_endpoint_iteration]
    f_late = median(pod_window)
    endpoint_row = next((r for r in explain
                         if int(r["iteration"]) == pod_endpoint_iteration), None)
    f_endpoint = float(endpoint_row["f_Q"]) if endpoint_row else None
    pod_explains = (f_late is not None and f_endpoint is not None and
                    f_late >= t["POD_late_median_explanation_min"] and
                    f_endpoint >= t["POD_iteration40_explanation_min"])

    reachability_iteration = 40 if iterations >= 40 else max(iterations - 1, 0)
    reach_rows = [r for r in subspace
                  if int(r["iteration"]) == reachability_iteration]
    reach_modes = {int(r["mode_id"]): r for r in reach_rows
                   if r["row_type"] == "MODE"}
    reach_aggregate = next((r for r in reach_rows
                            if r["row_type"] == "AGGREGATE"), None)
    important_ids = [m["mode_id"] for m in mode_summary
                     if m["energy_weight"] >= t["authoritative_mode_energy_min"]]
    y_reachable = bool(reach_aggregate) and all(mid in reach_modes and
        float(reach_modes[mid]["E_Y_all"]) <= t["Y_mode_reachability_error_limit"]
        for mid in important_ids) and \
        float(reach_aggregate["E_Y_weighted_all"]) <= t["Y_weighted_reachability_error_limit"]

    try:
        Htrue = matrix_from_rows(matrix_rows, "H_true")
        Hright = matrix_from_rows(matrix_rows, "H_right")
        if Htrue.shape != Hright.shape or Htrue.shape[0] != pod.get("selected_modes"):
            raise ValueError("projected matrix dimension mismatch")
    except (KeyError, TypeError, ValueError) as error:
        errors.append(str(error))
        Htrue = np.zeros((1, 1))
        Hright = np.zeros((1, 1))
    if not np.isfinite(Htrue).all() or not np.isfinite(Hright).all():
        errors.append("nonfinite projected matrix")
    right_metrics = projected_metrics(Hright, t["sigma_zero_relative_threshold"])
    true_metrics = projected_metrics(Htrue, t["sigma_zero_relative_threshold"])
    try:
        leakages = [float(x) for x in reference["leakage_per_mode"]]
        weighted_leakage = float(reference["leakage_weighted"])
        if (len(leakages) != pod.get("selected_modes") or
                not all(math.isfinite(x) for x in leakages) or
                not math.isfinite(weighted_leakage)):
            raise ValueError("invalid projected-subspace leakage")
    except (KeyError, TypeError, ValueError) as error:
        errors.append(str(error))
        leakages = [math.inf] * int(pod.get("selected_modes", 0))
        weighted_leakage = math.inf
    projected_closed = all(leakages[mid - 1] <= t["projected_subspace_mode_leakage_limit"]
                           for mid in important_ids) and \
        weighted_leakage <= t["projected_subspace_weighted_leakage_limit"]
    if not reference.get("S_ref_validation_pass"):
        errors.append("S_ref validation failed")
    if not reference.get("H_true_validation_pass"):
        errors.append("H_true validation failed")
    if not reference.get("H_right_fixed_operator_valid"):
        errors.append("H_right fixed-operator validation failed")
    if reference.get("raw_vector_action_nonfinite") is not False:
        errors.append("raw vector/action nonfinite or unverified")
    if (true_metrics["symmetry_relative_defect"] >
            t["H_true_symmetry_relative_tolerance"]):
        errors.append("H_true symmetry threshold failed")
    raw_minimum = reference.get("H_true_minimum_symmetric_eigenvalue")
    raw_scale = reference.get("reference_action_norm_scale")
    raw_floor = reference.get("numerical_floor")
    if not all(finite(x) for x in (raw_minimum, raw_scale, raw_floor)):
        errors.append("H_true curvature validation diagnostics nonfinite")
    elif float(raw_minimum) < -t["H_true_negative_curvature_relative_tolerance"] * \
            max(abs(float(raw_scale)), float(raw_floor)):
        errors.append("H_true negative-curvature threshold failed")

    plateau = (reproduction_pass and runtime.get("iterations") == 60 and
               not runtime.get("joint_target_reached") and
               float(runtime["best_R_cont"]) > float(runtime["R_cont_target"]))
    adverse = []
    if not projected_closed:
        adverse.append("POD_SUBSPACE_LEAKAGE")
    elif reference.get("H_right_fixed_operator_valid"):
        if right_metrics["kappa_2_is_infinite"] or \
                right_metrics["kappa_2"] >= t["projected_condition_number_adverse"]:
            adverse.append("ILL_CONDITIONED")
        if right_metrics["eta_mix"] >= t["projected_mixing_adverse"]:
            adverse.append("STRONG_MODE_MIXING")
        if right_metrics["eta_non_normal"] >= t["projected_non_normality_adverse"]:
            adverse.append("STRONGLY_NON_NORMAL")
    slow_important = any(status in ("STAGNATING", "WEAK_DECAY")
                         for _, status in significant)

    if errors:
        decision = "INVALID_EXPERIMENT"
        next_task = "REPAIR_AND_REPEAT_STRICT_STAGE_F3"
    elif plateau and not pod_explains:
        decision = "POD_SUBSPACE_NO_LONGER_EXPLAINS_PLATEAU"
        next_task = "STOP_AND_REVIEW_BEFORE_POD_REFRESH"
    elif not y_reachable and slow_important:
        decision = "OUTER_SCHUR_IMAGE_REACHABILITY_LIMITED"
        next_task = "GLOBAL_OR_DEFLATION_REACHABILITY_DESIGN"
    elif y_reachable and slow_important and reference["H_right_fixed_operator_valid"] and adverse:
        decision = "OUTER_SPACE_REACHABLE_BUT_PROJECTED_OPERATOR_ADVERSE"
        next_task = "STOP_AND_SELECT_PROJECTED_OPERATOR_REMEDY"
    else:
        decision = "OUTER_BOTTLENECK_UNRESOLVED"
        next_task = "STOP_AND_REVIEW_BEFORE_NEXT_SOLVER_DESIGN"

    projected = {
        "schema": "strict-ala-stage-F3-projected-operator-v1",
        "outer_preconditioning_orientation": "RIGHT_FGMRES",
        "H_true_definition": "Q^T S_ref Q",
        "H_right_definition": "Q^T S_ref M^-1 Q",
        "H_true": Htrue.tolist(), "H_right": Hright.tolist(),
        "H_true_metrics": true_metrics, "H_right_metrics": right_metrics,
        "projected_subspace_closed": projected_closed,
        "leakage_per_mode": json_safe(leakages),
        "leakage_weighted": json_safe(weighted_leakage),
        "spectral_interpretation_authoritative": bool(projected_closed and
            reference["H_right_fixed_operator_valid"]),
    }
    write_json(args.projected_operator_output, projected)
    reference_final = json_safe(dict(reference))
    reference_final.update({
        "schema": "strict-ala-stage-F3-reference-operator-validation-v1",
        "H_true_metrics": true_metrics,
        "H_right_metrics": right_metrics,
        "derived_infinite_condition_number_is_valid": True,
        "raw_vector_action_nonfinite_invalid": True,
    })
    write_json(args.reference_validation_output, reference_final)
    result = {
        "schema": "strict-ala-stage-F3-decision-v1",
        "decision": decision, "next_authorized_task": next_task,
        "valid": not errors, "complete": True, "validity_errors": errors,
        "base_reproduction_pass": reproduction_pass,
        "base_reproduction_max_relative_difference": reproduction_max,
        "base_plateau_confirmed": plateau,
        "base_plateau_definition": "reproduction_gate_pass && iterations==60 && !joint_target_reached && best_R_cont>target",
        "pod_reconstruction_max_relative_coefficient_error": reconstruction_max,
        "all_required_inner_solves_satisfy_requested_target": inner_targets_pass,
        "pod_explanation_late_median_20_40": f_late,
        "pod_explanation_endpoint_iteration": pod_endpoint_iteration,
        "pod_explanation_endpoint": f_endpoint,
        "shortened_window_due_to_early_joint_convergence": early_joint,
        "pod_subspace_explains_plateau": pod_explains,
        "mode_summary": mode_summary,
        "reachability_evaluation_iteration": reachability_iteration,
        "reachability_iteration_40_substituted": reachability_iteration != 40,
        "Y_REACHABLE_AT_40": y_reachable if reachability_iteration == 40 else None,
        "Y_REACHABLE_AT_EVALUATION_ITERATION": y_reachable,
        "V_Z_descriptive_geometry_only": True,
        "projected_subspace_closed": projected_closed,
        "adverse_causes": adverse,
        "restart_retention_status": "INCONCLUSIVE_SINGLE_BOUNDARY",
        "restart_boundary_observed": runtime.get("restart_boundary_observed", False),
        "recovery_metric_available": False,
        "production_default_change_authorized": False,
        "thresholds_sha256": sha256(args.thresholds),
    }
    write_json(args.output, result)
    return 0


def audit(args):
    decision = load_json(args.decision)
    pairs = [(args.binary_pre, args.binary_post),
             (args.inputs_pre, args.inputs_post),
             (args.source_pre, args.source_post)]
    provenance = all(Path(a).read_bytes() == Path(b).read_bytes() for a, b in pairs)
    value = {
        "schema": "strict-ala-stage-F3-final-audit-v1",
        "valid": bool(decision.get("valid") and provenance),
        "complete": bool(decision.get("complete") and provenance),
        "decision": decision.get("decision"),
        "next_authorized_task": decision.get("next_authorized_task"),
        "provenance_unchanged": provenance,
        "production_default_change_authorized": False,
    }
    write_json(args.output, value)
    return 0 if provenance else 1


def parser():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="command", required=True)
    q = sub.add_parser("cfg-contract")
    q.add_argument("--cfg", required=True); q.add_argument("--thresholds", required=True)
    q.add_argument("--output", required=True); q.set_defaults(fn=cfg_contract)
    q = sub.add_parser("analyze")
    for name in ("thresholds", "runtime", "pod-reconstruction",
                 "reference-validation-raw", "completion", "pod-explanation",
                 "inner-solves",
                 "mode-coefficients", "cumulative-subspace", "basis-rank",
                 "projected-matrices", "reconstructed-coefficients",
                 "authoritative-coefficients", "baseline-iterations",
                 "projected-operator-output", "reference-validation-output",
                 "output"):
        q.add_argument("--" + name, required=True)
    q.set_defaults(fn=analyze)
    q = sub.add_parser("audit")
    for name in ("decision", "binary-pre", "binary-post", "inputs-pre",
                 "inputs-post", "source-pre", "source-post", "output"):
        q.add_argument("--" + name, required=True)
    q.set_defaults(fn=audit)
    return p


if __name__ == "__main__":
    args = parser().parse_args()
    raise SystemExit(args.fn(args))
