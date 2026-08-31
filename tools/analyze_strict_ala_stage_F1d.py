#!/usr/bin/env python3
"""Fail-closed Stage-F1d maximal-SPD residual coarse qualification."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

RMS_RATIO_LIMIT = 0.50
SYMMETRY_LIMIT = 1.0e-8
PROJECTION_LIMIT = 1.0e-8
MAX_TIGHT_SOLVES = 18


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


def relative_symmetry(rows, matrix):
    selected = [row for row in rows if row["matrix"] == matrix]
    values = {(int(row["row_mode"]), int(row["column_mode"])):
              finite(row["value"], matrix) for row in selected}
    modes = sorted({index for key in values for index in key})
    if not modes or set(values) != {(i, j) for i in modes for j in modes}:
        raise ValueError(f"incomplete projected matrix {matrix}")
    scale = max(abs(value) for value in values.values())
    defect = max(abs(values[i, j] - values[j, i])
                 for i in modes for j in modes)
    return defect / max(scale, 1.0e-300)


def complete_scope(rows, scope, size):
    selected = [row for row in rows if row["scope"] == scope]
    keys = {(int(row["row_basis"]), int(row["column_basis"]))
            for row in selected}
    expected = {(i, j) for i in range(1, size + 1)
                for j in range(1, size + 1)}
    if len(selected) != size * size or keys != expected:
        raise ValueError(f"incomplete coarse matrix scope {scope}")
    for row in selected:
        finite(row["value"], f"{scope} value")
    return True


def scope_symmetry(rows, scope):
    selected = [row for row in rows if row["scope"] == scope]
    values = {(int(row["row_basis"]), int(row["column_basis"])):
              finite(row["value"], f"{scope} value") for row in selected}
    scale = max(abs(value) for value in values.values())
    return max(abs(values[i, j] - values[j, i]) for i, j in values) / \
        max(scale, 1.0e-300)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def analyze(args):
    baseline_all = read_csv(args.f1b_true_mode)
    baseline = [row for row in baseline_all
                if row.get("operator") == "operator_consistent"]
    f1c = load_json(args.f1c_decision)
    candidate = read_csv(args.candidate_true_mode, (
        "POD_mode", "POD_energy_weight", "E_P", "cosine", "qTPq",
        "qTSPq", "original_projection_max", "enriched_projection_max",
        "tight_solve_achieved", "valid"))
    basis = read_csv(args.basis, (
        "basis_index", "basis_type", "source_mode", "preorthogonal_norm",
        "postorthogonal_norm", "accepted", "maximum_orthogonality",
        "S_energy", "valid"))
    projected = read_csv(args.projected,
                         ("matrix", "row_mode", "column_mode", "value"))
    tight = read_csv(args.tight, (
        "call_id", "role", "POD_mode", "rhs_norm",
        "requested_relative_tolerance", "target_residual",
        "achieved_relative_residual", "cycles", "max_cycles", "converged"))
    coarse_matrix = read_csv(args.coarse_matrix,
        ("scope", "row_basis", "column_basis", "value"))
    pivots = read_csv(args.pivots,
        ("scope", "step", "pivot", "pivot_ratio", "evaluated", "passed"))
    spectrum = read_csv(args.spectrum,
        ("eigen_index", "eigenvalue", "relative_to_max", "positive"))
    selected_basis = read_csv(args.selected_basis,
        ("selected_index", "full_basis_index", "basis_type", "source_mode",
         "retained"))
    runtime = load_json(args.runtime)
    if not baseline:
        raise ValueError("missing operator-consistent F1b baseline")
    baseline_by_mode = {row["POD_mode"]: row for row in baseline}
    candidate_by_mode = {row["POD_mode"]: row for row in candidate}
    if set(baseline_by_mode) != set(candidate_by_mode):
        raise ValueError("F1d changed the frozen POD mode set")
    for mode, row in candidate_by_mode.items():
        before = finite(baseline_by_mode[mode]["POD_energy_weight"], "weight")
        after = finite(row["POD_energy_weight"], "weight")
        if not math.isclose(before, after, rel_tol=1.0e-12,
                            abs_tol=1.0e-300):
            raise ValueError("F1d changed a frozen POD weight")

    selected = len(candidate)
    accepted = sum(row["accepted"] == "1" for row in basis)
    retained = int(runtime.get("retained_residual_modes", -1))
    full_dimension = selected + accepted
    active_dimension = selected + retained
    expected_tight = 2 * selected + accepted
    baseline_rms = weighted_rms(baseline)
    candidate_rms = weighted_rms(candidate)
    f1c_rms = finite(f1c.get("metrics", {}).get("E_P_RMS_F1c"), "F1c RMS")
    ratio = candidate_rms / max(baseline_rms, 1.0e-300)
    projection = max(finite(row["enriched_projection_max"], "projection")
                     for row in candidate)
    original_projection = max(
        finite(row["original_projection_max"], "original projection")
        for row in candidate)
    total_weight = sum(finite(row["POD_energy_weight"], "weight")
                       for row in candidate)
    positive_weight = sum(finite(row["POD_energy_weight"], "weight")
                          for row in candidate
                          if finite(row["qTSPq"], "qTSPq") > 0.0)
    heavy_guard = all(
        finite(candidate_by_mode[mode]["E_P"], "candidate E_P") <=
        1.10 * finite(baseline_by_mode[mode]["E_P"], "baseline E_P")
        for mode in baseline_by_mode
        if finite(baseline_by_mode[mode]["POD_energy_weight"], "weight") >=
        0.05 * total_weight)
    selected_indices = [int(row["selected_index"]) for row in selected_basis]
    selected_full_indices = [int(row["full_basis_index"])
                             for row in selected_basis]
    selected_original = [row for row in selected_basis
                         if row["basis_type"] == "POD"]
    selected_residual = [row for row in selected_basis
                         if row["basis_type"] == "coarse_residual"]
    basis_valid = (
        len(basis) == selected and accepted >= 0 and
        accepted == runtime.get("accepted_residual_modes") and
        0 <= retained <= accepted and
        runtime.get("full_enriched_basis_dimension") == full_dimension and
        runtime.get("enriched_basis_dimension") == active_dimension and
        all(row["basis_type"] == "coarse_residual" and row["valid"] == "1"
            for row in basis) and
        all(finite(row["postorthogonal_norm"], "basis norm") > 0.0 and
            finite(row["maximum_orthogonality"], "basis orthogonality") <=
            1.0e-8 and finite(row["S_energy"], "basis S energy") > 0.0
            for row in basis if row["accepted"] == "1") and
        len(selected_basis) == active_dimension and
        selected_indices == list(range(1, active_dimension + 1)) and
        len(set(selected_full_indices)) == active_dimension and
        len(selected_original) == selected and
        sorted(int(row["full_basis_index"]) for row in selected_original) ==
            list(range(1, selected + 1)) and
        len(selected_residual) == retained and
        all(row["retained"] == "1" for row in selected_basis))
    tight_pass = (
        runtime.get("tight_solve_count") == len(tight) == expected_tight and
        len(tight) <= MAX_TIGHT_SOLVES and
        all(row["converged"] == "1" and
            finite(row["achieved_relative_residual"], "tight achieved") <=
            finite(row["requested_relative_tolerance"], "tight requested")
            for row in tight))
    t_symmetry = relative_symmetry(projected, "T_enriched")
    m_symmetry = relative_symmetry(projected, "M_balanced")
    h_symmetry = relative_symmetry(projected, "H_balanced")
    matrix_valid = (
        complete_scope(coarse_matrix, "FULL_RAW", full_dimension) and
        complete_scope(coarse_matrix, "FULL_SYMMETRIZED", full_dimension) and
        complete_scope(coarse_matrix, "SELECTED_RAW", active_dimension) and
        complete_scope(coarse_matrix, "SELECTED_SYMMETRIZED", active_dimension))
    full_pivots = [row for row in pivots if row["scope"] == "FULL"]
    selected_pivots = [row for row in pivots if row["scope"] == "SELECTED"]
    pivot_valid = (
        len(full_pivots) == full_dimension and
        len(selected_pivots) == active_dimension and
        [int(row["step"]) for row in full_pivots] ==
            list(range(1, full_dimension + 1)) and
        [int(row["step"]) for row in selected_pivots] ==
            list(range(1, active_dimension + 1)) and
        all(row["evaluated"] == "1" and row["passed"] == "1" and
            finite(row["pivot"], "selected pivot") > 0.0 and
            finite(row["pivot_ratio"], "selected pivot ratio") > 0.0
            for row in selected_pivots) and
        all(row["evaluated"] in ("0", "1") and row["passed"] in ("0", "1")
            for row in full_pivots))
    spectrum_valid = (
        len(spectrum) == full_dimension and
        [int(row["eigen_index"]) for row in spectrum] ==
            list(range(1, full_dimension + 1)) and
        all(math.isfinite(finite(row["eigenvalue"], "coarse eigenvalue")) and
            row["positive"] in ("0", "1") for row in spectrum))
    full_symmetry = scope_symmetry(coarse_matrix, "FULL_RAW")
    selected_symmetry = scope_symmetry(coarse_matrix, "SELECTED_RAW")
    spectral_values = [finite(row["eigenvalue"], "coarse eigenvalue")
                       for row in spectrum]
    full_factor_pass = runtime.get("full_matrix_cholesky_pass")
    full_failure_step = runtime.get("full_matrix_failure_step")
    full_factor_consistent = (
        isinstance(full_factor_pass, bool) and
        isinstance(full_failure_step, int) and
        ((full_factor_pass and full_failure_step == -1 and
          all(row["evaluated"] == "1" and row["passed"] == "1"
              for row in full_pivots)) or
         (not full_factor_pass and 0 <= full_failure_step < full_dimension and
          full_pivots[full_failure_step]["evaluated"] == "1" and
          full_pivots[full_failure_step]["passed"] == "0")))
    runtime_diagnostics_consistent = (
        math.isclose(full_symmetry,
            finite(runtime.get("full_matrix_symmetry_defect"),
                   "full matrix symmetry"), rel_tol=1.0e-12,
            abs_tol=1.0e-15) and
        math.isclose(selected_symmetry,
            finite(runtime.get("enriched_matrix_symmetry_defect"),
                   "selected matrix symmetry"), rel_tol=1.0e-12,
            abs_tol=1.0e-15) and
        math.isclose(min(spectral_values),
            finite(runtime.get("full_matrix_eigenvalue_minimum"),
                   "minimum eigenvalue"), rel_tol=1.0e-12,
            abs_tol=1.0e-300) and
        math.isclose(max(spectral_values),
            finite(runtime.get("full_matrix_eigenvalue_maximum"),
                   "maximum eigenvalue"), rel_tol=1.0e-12,
            abs_tol=1.0e-300))
    lineage_valid = (
        f1c.get("experiment_evidence_valid") is True and
        f1c.get("decision") == "F1C_BALANCED_COARSE_REJECTED" and
        f1c.get("next_authorized_task") == "F1D_COARSE_BASIS_REDESIGN" and
        f1c.get("production_default_change_authorized") is False)
    structural_gates = {
        "F1c_lineage": lineage_valid,
        "runtime_candidate": runtime.get("candidate") ==
            "maximal_spd_residual_closed_balanced_actual_mode_coarse",
        "fine_map_frozen": runtime.get("fine_map") ==
            "F1b_operator_consistent_frozen",
        "production_freeze":
            runtime.get("production_default_change_authorized") is False,
        "runtime_mode_count": runtime.get("selected_modes") == selected,
        "runtime_selected_energy":
            finite(runtime.get("selected_energy"), "selected energy") >= 0.95,
        "basis_residual_closure": basis_valid,
        "full_matrix_diagnostics": matrix_valid and spectrum_valid and
            runtime_diagnostics_consistent,
        "cholesky_diagnostics": pivot_valid and full_factor_consistent,
        "coarse_SPD": finite(
            runtime.get("enriched_matrix_minimum_pivot_ratio"),
            "enriched pivot ratio") > 0.0,
        "coarse_symmetry": max(
            finite(runtime.get("enriched_matrix_symmetry_defect"),
                   "runtime symmetry"), t_symmetry) <= SYMMETRY_LIMIT,
        "balanced_map_symmetry": m_symmetry <= SYMMETRY_LIMIT,
        "tight_solves": tight_pass,
        "rows_valid": all(row["valid"] == "1" for row in candidate),
    }
    performance_gates = {
        "retains_residual_mode": retained > 0,
        "enriched_residual_projection": projection <= PROJECTION_LIMIT,
        "improves_F1c": candidate_rms < f1c_rms,
        "true_mode_RMS": ratio <= RMS_RATIO_LIMIT,
        "true_mode_heavy": heavy_guard,
        "dominant_mode_positivity": positive_weight >= 0.95 * total_weight,
    }
    gates = {**structural_gates, **performance_gates}
    passed = all(gates.values())
    structural_pass = all(structural_gates.values())
    if passed:
        decision = "F1D_RESIDUAL_CLOSED_COARSE_QUALIFIED"
        next_task = "G0_GENERALIZABLE_COARSE_BASIS"
    elif structural_pass:
        decision = "F1D_RESIDUAL_CLOSED_COARSE_REJECTED"
        next_task = "F2_AVV_PRECONDITIONER_REVIEW"
    else:
        decision = "F1D_INVALID_STRUCTURAL_EVIDENCE"
        next_task = "REPEAT_F1D_VALID_EXPERIMENT"
    result = {
        "schema": "strict-ala-stage-F1d-decision-v2",
        "experiment_evidence_valid": structural_pass,
        "candidate": "maximal_spd_residual_closed_balanced_actual_mode_coarse",
        "metrics": {
            "E_P_RMS_F1b": baseline_rms,
            "E_P_RMS_F1c": f1c_rms,
            "E_P_RMS_F1d": candidate_rms,
            "E_P_RMS_ratio_to_F1b": ratio,
            "E_P_RMS_ratio_to_F1c": candidate_rms / max(f1c_rms, 1.0e-300),
            "original_projection_max": original_projection,
            "enriched_projection_max": projection,
            "accepted_residual_modes": accepted,
            "retained_residual_modes": retained,
            "full_enriched_basis_dimension": full_dimension,
            "enriched_basis_dimension": active_dimension,
            "full_matrix_symmetry_defect": finite(
                runtime.get("full_matrix_symmetry_defect"),
                "full matrix symmetry"),
            "full_matrix_cholesky_pass":
                runtime.get("full_matrix_cholesky_pass"),
            "full_matrix_failure_step":
                runtime.get("full_matrix_failure_step"),
            "full_matrix_eigenvalue_minimum": finite(
                runtime.get("full_matrix_eigenvalue_minimum"),
                "minimum eigenvalue"),
            "full_matrix_eigenvalue_maximum": finite(
                runtime.get("full_matrix_eigenvalue_maximum"),
                "maximum eigenvalue"),
            "T_enriched_symmetry_defect": t_symmetry,
            "M_balanced_symmetry_defect": m_symmetry,
            "H_balanced_symmetry_defect": h_symmetry,
        },
        "thresholds": {
            "E_P_RMS_ratio_limit": RMS_RATIO_LIMIT,
            "symmetry_limit": SYMMETRY_LIMIT,
            "projection_limit": PROJECTION_LIMIT,
            "max_tight_solves": MAX_TIGHT_SOLVES,
        },
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
        "schema": "strict-ala-stage-F1d-final-audit-v1",
        "experiment_valid": valid,
        "decision": decision.get("decision") if valid else "INVALID_EXPERIMENT",
        "candidate_accepted": (valid and decision.get("decision") ==
                               "F1D_RESIDUAL_CLOSED_COARSE_QUALIFIED"),
        "provenance": provenance,
        "F1d": decision,
        "production_default_change_authorized": False,
        "next_authorized_task": (decision.get("next_authorized_task") if valid
                                 else "REPEAT_F1D_VALID_EXPERIMENT"),
    }
    write_json(args.output, result)
    return 0


def parser():
    root = argparse.ArgumentParser()
    commands = root.add_subparsers(required=True)
    p = commands.add_parser("analyze")
    for name in ("f1b-true-mode", "f1c-decision", "candidate-true-mode",
                 "basis", "projected", "tight", "coarse-matrix", "pivots",
                 "spectrum", "selected-basis", "runtime", "output"):
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
