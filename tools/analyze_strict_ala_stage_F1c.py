#!/usr/bin/env python3
"""Fail-closed Stage-F1c balanced actual-mode coarse qualification."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

RMS_RATIO_LIMIT = 0.50
SYMMETRY_LIMIT = 1.0e-8
COARSE_PROJECTION_LIMIT = 1.0e-8
MAX_TIGHT_SOLVES = 12


def load_json(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing or empty {path}")
    return json.loads(path.read_text())


def read_csv(path, columns=None):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing or empty {path}")
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if columns and tuple(reader.fieldnames or ()) != tuple(columns):
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


def weighted_rms(rows):
    numerator = denominator = 0.0
    for row in rows:
        weight = finite(row["POD_energy_weight"], "POD weight")
        error = finite(row["E_P"], "E_P")
        if weight < 0.0:
            raise ValueError("negative POD weight")
        numerator += weight * error * error
        denominator += weight
    if denominator <= 0.0:
        raise ValueError("nonpositive POD weight sum")
    return math.sqrt(numerator / denominator)


def weighted_column_rms(rows, column):
    numerator = denominator = 0.0
    for row in rows:
        weight = finite(row["POD_energy_weight"], "POD weight")
        value = finite(row[column], column)
        numerator += weight * value * value
        denominator += weight
    if denominator <= 0.0:
        raise ValueError("nonpositive decomposition weight sum")
    return math.sqrt(numerator / denominator)


def relative_symmetry(rows, matrix):
    selected = [row for row in rows if row["matrix"] == matrix]
    values = {(int(row["row_mode"]), int(row["column_mode"])):
              finite(row["value"], matrix) for row in selected}
    modes = sorted({index for key in values for index in key})
    if set(values) != {(i, j) for i in modes for j in modes}:
        raise ValueError(f"incomplete projected matrix {matrix}")
    scale = max((abs(value) for value in values.values()), default=0.0)
    defect = max((abs(values[i, j] - values[j, i])
                  for i in modes for j in modes), default=0.0)
    return defect / max(scale, 1.0e-300)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def analyze(args):
    baseline_all = read_csv(args.f1b_true_mode)
    baseline = [row for row in baseline_all
                if row.get("operator") == "operator_consistent"]
    candidate = read_csv(args.candidate_true_mode, (
        "POD_mode", "POD_energy_weight", "E_P", "cosine", "qTPq",
        "qTSPq", "coarse_residual_projection_max",
        "tight_solve_achieved", "valid"))
    runtime = load_json(args.runtime)
    projected = read_csv(args.projected,
                         ("matrix", "row_mode", "column_mode", "value"))
    tight = read_csv(args.tight, (
        "call_id", "role", "POD_mode", "rhs_norm",
        "requested_relative_tolerance", "target_residual",
        "achieved_relative_residual", "cycles", "max_cycles",
        "converged"))
    decomposition = read_csv(args.decomposition, (
        "POD_mode", "POD_energy_weight", "E_coarse_only",
        "E_conditioned_fine", "E_right_additive", "E_balanced",
        "conditioned_rhs_relative", "conditioned_fine_action_relative",
        "left_balance_action_relative", "coarse_projection_max",
        "right_projection_max", "balanced_projection_max", "valid"))
    alpha_scan = read_csv(args.alpha_scan,
                          ("alpha", "E_P_RMS", "is_clipped_optimum", "valid"))
    if not baseline:
        raise ValueError("missing operator-consistent F1b baseline")
    baseline_by_mode = {row["POD_mode"]: row for row in baseline}
    candidate_by_mode = {row["POD_mode"]: row for row in candidate}
    decomposition_by_mode = {row["POD_mode"]: row for row in decomposition}
    if set(baseline_by_mode) != set(candidate_by_mode):
        raise ValueError("F1c changed the frozen POD mode set")
    if set(candidate_by_mode) != set(decomposition_by_mode):
        raise ValueError("F1c decomposition changed the frozen POD mode set")
    for mode in baseline_by_mode:
        before = finite(baseline_by_mode[mode]["POD_energy_weight"], "weight")
        after = finite(candidate_by_mode[mode]["POD_energy_weight"], "weight")
        if not math.isclose(before, after, rel_tol=1.0e-12, abs_tol=1.0e-300):
            raise ValueError("F1c changed a frozen POD weight")
        decomposed = finite(decomposition_by_mode[mode]["POD_energy_weight"],
                            "decomposition weight")
        if not math.isclose(after, decomposed, rel_tol=1.0e-12,
                            abs_tol=1.0e-300):
            raise ValueError("F1c decomposition changed a frozen POD weight")
    baseline_rms = weighted_rms(baseline)
    candidate_rms = weighted_rms(candidate)
    ratio = candidate_rms / max(baseline_rms, 1.0e-300)
    total_weight = sum(finite(row["POD_energy_weight"], "weight")
                       for row in candidate)
    positive_weight = sum(finite(row["POD_energy_weight"], "weight")
                          for row in candidate
                          if finite(row["qTSPq"], "qTSPq") > 0.0)
    heavy_guard = all(
        finite(candidate_by_mode[mode]["E_P"], "candidate E_P") <=
        1.10 * finite(baseline_by_mode[mode]["E_P"], "baseline E_P")
        for mode in baseline_by_mode
        if finite(baseline_by_mode[mode]["POD_energy_weight"], "weight")
        >= 0.05 * total_weight)
    projection = max(finite(row["coarse_residual_projection_max"],
                            "coarse projection") for row in candidate)
    tight_pass = (runtime.get("tight_solve_count") == len(tight) and
                  len(tight) == 2 * len(candidate) and
                  len(tight) <= MAX_TIGHT_SOLVES and
                  all(row.get("converged") == "1" and
                      finite(row["achieved_relative_residual"],
                             "tight achieved residual") <=
                      finite(row["requested_relative_tolerance"],
                             "tight requested tolerance")
                      for row in tight))
    t_symmetry = relative_symmetry(projected, "T")
    m_symmetry = relative_symmetry(projected, "M_balanced")
    h_symmetry = relative_symmetry(projected, "H_balanced")
    coarse_rms = weighted_column_rms(decomposition, "E_coarse_only")
    fine_rms = weighted_column_rms(decomposition, "E_conditioned_fine")
    right_rms = weighted_column_rms(decomposition, "E_right_additive")
    decomposed_balanced_rms = weighted_column_rms(decomposition, "E_balanced")
    optimum_rows = [row for row in alpha_scan
                    if row["is_clipped_optimum"] == "1"]
    fixed_alpha = [row for row in alpha_scan
                   if row["is_clipped_optimum"] == "0"]
    alpha_zero_rows = [row for row in fixed_alpha
                       if math.isclose(finite(row["alpha"], "alpha"), 0.0,
                                       rel_tol=0.0, abs_tol=1.0e-15)]
    alpha_one_rows = [row for row in alpha_scan
                      if row["is_clipped_optimum"] == "0" and
                      math.isclose(finite(row["alpha"], "alpha"), 1.0,
                                   rel_tol=0.0, abs_tol=1.0e-15)]
    expected_alpha = {0.0, 0.125, 0.25, 0.5, 0.75, 1.0}
    observed_alpha = {finite(row["alpha"], "fixed alpha")
                      for row in fixed_alpha}
    if (len(optimum_rows) != 1 or len(alpha_zero_rows) != 1 or
            len(alpha_one_rows) != 1 or observed_alpha != expected_alpha or
            len(fixed_alpha) != len(expected_alpha)):
        raise ValueError("F1c alpha scan does not close fixed/optimum rows")
    optimum_alpha = finite(optimum_rows[0]["alpha"], "optimum alpha")
    optimum_rms = finite(optimum_rows[0]["E_P_RMS"], "optimum RMS")
    decomposition_consistent = (
        all(row["valid"] == "1" for row in decomposition) and
        all(math.isclose(
                finite(row["conditioned_rhs_relative"], "conditioned RHS"),
                finite(row["E_coarse_only"], "coarse-only E_P"),
                rel_tol=1.0e-8, abs_tol=1.0e-10) and
            math.isclose(
                finite(row["E_right_additive"], "right-additive E_P"),
                finite(row["E_coarse_only"], "coarse-only E_P") *
                finite(row["E_conditioned_fine"], "conditioned fine E_P"),
                rel_tol=1.0e-8, abs_tol=1.0e-10)
            for row in decomposition) and
        all(math.isclose(finite(decomposition_by_mode[mode]["E_balanced"],
                                "decomposed balanced E_P"),
                         finite(candidate_by_mode[mode]["E_P"],
                                "candidate E_P"),
                         rel_tol=1.0e-10, abs_tol=1.0e-12)
            for mode in candidate_by_mode) and
        math.isclose(decomposed_balanced_rms, candidate_rms,
                     rel_tol=1.0e-10, abs_tol=1.0e-12))
    alpha_consistent = (
        all(row["valid"] == "1" for row in alpha_scan) and
        0.0 <= optimum_alpha <= 1.0 and
        optimum_rms <= min(finite(row["E_P_RMS"], "fixed alpha RMS")
                           for row in fixed_alpha) + 1.0e-12 and
        math.isclose(finite(alpha_zero_rows[0]["E_P_RMS"], "alpha zero RMS"),
                     coarse_rms, rel_tol=1.0e-8, abs_tol=1.0e-10) and
        math.isclose(finite(alpha_one_rows[0]["E_P_RMS"], "alpha one RMS"),
                     candidate_rms, rel_tol=1.0e-8, abs_tol=1.0e-10))
    gates = {
        "runtime_candidate": (
            runtime.get("candidate") == "balanced_actual_mode_coarse"),
        "runtime_causal_review": (
            runtime.get("causal_review") ==
            "coarse_only_conditioned_fine_right_additive_left_balance_alpha_family" and
            runtime.get("fixed_alpha_count") == 6),
        "production_freeze": (
            runtime.get("production_default_change_authorized") is False),
        "runtime_mode_count": runtime.get("selected_modes") == len(candidate),
        "runtime_selected_energy": finite(runtime.get("selected_energy"),
                                            "selected energy") >= 0.95,
        "coarse_SPD": (finite(runtime.get("coarse_matrix_minimum_pivot_ratio"),
                               "pivot ratio") > 0.0),
        "coarse_symmetry": max(
            finite(runtime.get("coarse_matrix_symmetry_defect"),
                   "runtime symmetry"), t_symmetry) <= SYMMETRY_LIMIT,
        # P_balanced is symmetric. S*P_balanced need not be Euclidean-
        # symmetric, so H_balanced is reported but is not a symmetry gate.
        "balanced_map_symmetry": m_symmetry <= SYMMETRY_LIMIT,
        "decomposition_consistency": decomposition_consistent,
        "alpha_scan_consistency": alpha_consistent,
        "tight_solves": tight_pass,
        "coarse_residual_projection": projection <= COARSE_PROJECTION_LIMIT,
        "true_mode_RMS": ratio <= RMS_RATIO_LIMIT,
        "true_mode_heavy": heavy_guard,
        "dominant_mode_positivity": positive_weight >= 0.95 * total_weight,
        "rows_valid": all(row["valid"] == "1" for row in candidate),
    }
    passed = all(gates.values())
    causal_findings = {
        "coarse_basis_or_Galerkin_insufficient": coarse_rms > baseline_rms,
        "conditioned_fine_map_leakage": right_rms > 1.10 * coarse_rms,
        "left_balancing_amplification": candidate_rms > 1.10 * right_rms,
        "fitted_alpha_suppresses_fine_map": (
            optimum_alpha < 0.25 and optimum_rms < 0.90 * candidate_rms),
    }
    if passed:
        next_task = "G0_GENERALIZABLE_COARSE_BASIS"
    elif not decomposition_consistent or not alpha_consistent:
        next_task = "REPEAT_F1C_CAUSAL_EVIDENCE"
    elif causal_findings["coarse_basis_or_Galerkin_insufficient"]:
        next_task = "F1D_COARSE_BASIS_REDESIGN"
    elif causal_findings["conditioned_fine_map_leakage"]:
        next_task = "F1D_DEFLATED_FINE_MAP_REDESIGN"
    elif causal_findings["left_balancing_amplification"]:
        next_task = "F1D_ONE_SIDED_COARSE_TEST"
    else:
        next_task = "F1D_PRESSURE_HYPOTHESIS_UNRESOLVED"
    result = {
        "schema": "strict-ala-stage-F1c-decision-v1",
        "experiment_evidence_valid": True,
        "candidate": "balanced_actual_mode_coarse",
        "metrics": {
            "E_P_RMS_F1b": baseline_rms,
            "E_P_RMS_F1c": candidate_rms,
            "E_P_RMS_ratio": ratio,
            "coarse_projection_max": projection,
            "T_symmetry_defect": t_symmetry,
            "M_balanced_symmetry_defect": m_symmetry,
            "H_balanced_symmetry_defect": h_symmetry,
            "E_P_RMS_coarse_only": coarse_rms,
            "E_P_RMS_conditioned_fine": fine_rms,
            "E_P_RMS_right_additive": right_rms,
            "E_P_RMS_decomposed_balanced": decomposed_balanced_rms,
            "fitted_alpha": optimum_alpha,
            "E_P_RMS_fitted_alpha": optimum_rms,
        },
        "thresholds": {
            "E_P_RMS_ratio_limit": RMS_RATIO_LIMIT,
            "symmetry_limit": SYMMETRY_LIMIT,
            "coarse_projection_limit": COARSE_PROJECTION_LIMIT,
            "max_tight_solves": MAX_TIGHT_SOLVES,
        },
        "gates": gates,
        "causal_findings": causal_findings,
        "decision": ("F1C_PRESSURE_PHYSICS_QUALIFIED" if passed
                     else "F1C_BALANCED_COARSE_REJECTED"),
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
        "schema": "strict-ala-stage-F1c-final-audit-v1",
        "experiment_valid": valid,
        "decision": decision.get("decision") if valid else "INVALID_EXPERIMENT",
        "candidate_accepted": (valid and decision.get("decision") ==
                               "F1C_PRESSURE_PHYSICS_QUALIFIED"),
        "provenance": provenance,
        "F1c": decision,
        "production_default_change_authorized": False,
        "next_authorized_task": (decision.get("next_authorized_task") if valid
                                 else "REPEAT_F1C_VALID_EXPERIMENT"),
    }
    write_json(args.output, result)
    return 0


def parser():
    root = argparse.ArgumentParser()
    commands = root.add_subparsers(required=True)
    p = commands.add_parser("analyze")
    for name in ("f1b-true-mode", "candidate-true-mode", "runtime",
                 "projected", "tight", "decomposition", "alpha-scan",
                 "output"):
        p.add_argument(f"--{name}", required=True)
    p.set_defaults(fn=analyze)
    p = commands.add_parser("audit")
    p.add_argument("--decision", required=True)
    for name in ("binary-pre", "binary-post", "inputs-pre", "inputs-post",
                 "source-pre", "source-post", "output"):
        p.add_argument(f"--{name}", required=True)
    p.set_defaults(fn=audit)
    return root


def main():
    args = parser().parse_args()
    return args.fn(args)


if __name__ == "__main__":
    raise SystemExit(main())
