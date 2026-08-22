#!/usr/bin/env python3
"""Deterministic offline analysis for strict-ALA Stage-E spectroscopy."""

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np


PROBES = (
    "P0_initial_continuity", "P1_fixed_random", "P2_radial_smooth",
    "P3_radial_higher", "P4_radial_alternating",
    "P5_horizontal_checkerboard", "P6_degree_1", "P6_degree_2",
    "P6_degree_4", "P6_degree_8", "P7_depth_0_200",
    "P8_depth_200_410", "P9_depth_410_660", "P10_depth_660_1000",
    "P11_depth_1000_cmb", "P12_patch_scale",
    "P13_longer_than_patch", "P14_constant", "P15_density_gauge",
)

FAMILIES = {
    "GLOBAL_LOW_DEGREE": ("P6_degree_1", "P6_degree_2", "P6_degree_4"),
    "GLOBAL_INTERMEDIATE": ("P6_degree_8", "P13_longer_than_patch"),
    "RADIAL_SMOOTH": ("P2_radial_smooth", "P3_radial_higher"),
    "RADIAL_HIGH_FREQUENCY": ("P4_radial_alternating",),
    "LOCAL_HIGH_FREQUENCY": ("P5_horizontal_checkerboard", "P12_patch_scale"),
    "DEPTH_SHALLOW": ("P7_depth_0_200", "P8_depth_200_410"),
    "DEPTH_TRANSITION": ("P9_depth_410_660",),
    "DEPTH_DEEP": ("P10_depth_660_1000", "P11_depth_1000_cmb"),
    "GAUGE_NEAR_NULL": ("P14_constant", "P15_density_gauge"),
    "RANDOM_GUARD": ("P1_fixed_random",),
    "ACTUAL_INITIAL": ("P0_initial_continuity",),
}

PRIMARY_FAMILIES = (
    "GLOBAL_LOW_DEGREE", "GLOBAL_INTERMEDIATE", "RADIAL_SMOOTH",
    "RADIAL_HIGH_FREQUENCY", "LOCAL_HIGH_FREQUENCY", "GAUGE_NEAR_NULL",
)

FAMILY_CLASS = {
    "GLOBAL_LOW_DEGREE": "GLOBAL",
    "GLOBAL_INTERMEDIATE": "GLOBAL",
    "RADIAL_SMOOTH": "RADIAL",
    "RADIAL_HIGH_FREQUENCY": "RADIAL",
    "LOCAL_HIGH_FREQUENCY": "LOCAL",
    "GAUGE_NEAR_NULL": "GAUGE_NEAR_NULL",
}

PRE_WINDOW = tuple(range(45, 50))
POST_WINDOW = tuple(range(51, 56))
DOMINANCE_FRACTION = 0.60
MIN_REPRESENTED_FRACTION = 0.50
OVERLAP_AMBIGUITY = 0.90
MODE_COMPOSITION_CHANGE = 0.10
RESTART_SUPPORTED_RATIO = 1.05
RESTART_STRONG_RATIO = 1.10
RESTART_CONTINUITY_TOLERANCE = 1.0e-10
STAGE_B_MPI_ENVELOPE = 1.0e-14
PROJECTION_RECONCILIATION_TOLERANCE = 1.0e-10


def _rows(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def _write_json(path, payload):
    Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _finite(value):
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def _prefix_file(prefix, suffix):
    return Path(str(prefix) + "." + suffix)


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_case(prefix):
    suffixes = {
        "iterations": "strict_ala_stage_E_iterations.csv",
        "correlations": "strict_ala_stage_E_correlations.csv",
        "gram": "strict_ala_stage_E_probe_gram.csv",
        "hessenberg": "strict_ala_stage_E_hessenberg.csv",
        "work": "strict_ala_stage_E_work_counters.csv",
        "restart": "strict_ala_stage_E_restart.csv",
    }
    paths = {name: _prefix_file(prefix, suffix) for name, suffix in suffixes.items()}
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise ValueError("missing Stage-E files: " + ", ".join(missing))
    return {name: _rows(path) for name, path in paths.items()}


def gram_matrix(rows):
    index = {name: i for i, name in enumerate(PROBES)}
    matrix = np.full((len(PROBES), len(PROBES)), np.nan)
    fingerprints = {}
    for row in rows:
        i = index[row["probe_i"]]
        j = index[row["probe_j"]]
        if math.isfinite(matrix[i, j]):
            raise ValueError("duplicate Gram entry")
        matrix[i, j] = float(row["value"])
        fingerprints[row["probe_i"]] = float(row["probe_i_fingerprint"])
    if not np.isfinite(matrix).all():
        raise ValueError("incomplete/nonfinite Gram matrix")
    if np.max(np.abs(matrix - matrix.T)) > 1.0e-12:
        raise ValueError("Gram matrix is not symmetric")
    return matrix, fingerprints


def svd_projection(gram, correlation, residual_energy):
    gram = np.asarray(gram, dtype=float)
    correlation = np.asarray(correlation, dtype=float)
    u, singular, vt = np.linalg.svd(gram, hermitian=True)
    scale = singular[0] if singular.size else 0.0
    tolerance = max(gram.shape) * np.finfo(float).eps * max(scale, 1.0)
    keep = singular > tolerance
    inverse = np.zeros_like(singular)
    inverse[keep] = 1.0 / singular[keep]
    pseudo = (vt.T * inverse) @ u.T
    coefficients = pseudo @ correlation
    represented = float(coefficients @ gram @ coefficients)
    unexplained = float(
        residual_energy - 2.0 * coefficients @ correlation + represented)
    reconciliation = abs(represented + unexplained - residual_energy) / max(
        residual_energy, 1.0e-300)
    return {
        "rank": int(np.count_nonzero(keep)),
        "singular_values": singular.tolist(),
        "rank_tolerance": float(tolerance),
        "condition_number": (
            float(singular[keep][0] / singular[keep][-1])
            if np.count_nonzero(keep) else None
        ),
        "coefficients": coefficients.tolist(),
        "represented_energy": represented,
        "unexplained_energy": unexplained,
        "represented_fraction": represented / max(residual_energy, 1.0e-300),
        "unexplained_fraction": unexplained / max(residual_energy, 1.0e-300),
        "reconciliation_relative": reconciliation,
        "physical_energy_pass": (
            represented >= -PROJECTION_RECONCILIATION_TOLERANCE
            * max(residual_energy, 1.0)
            and unexplained >= -PROJECTION_RECONCILIATION_TOLERANCE
            * max(residual_energy, 1.0)
            and represented <= residual_energy
            + PROJECTION_RECONCILIATION_TOLERANCE
            * max(residual_energy, 1.0)
        ),
    }


def family_projection(gram, correlation, residual_energy, names):
    indices = [PROBES.index(name) for name in names]
    result = svd_projection(
        gram[np.ix_(indices, indices)], correlation[indices], residual_energy)
    return result["represented_fraction"]


def family_overlap(gram, family_a, family_b):
    ia = [PROBES.index(name) for name in FAMILIES[family_a]]
    ib = [PROBES.index(name) for name in FAMILIES[family_b]]
    gaa = gram[np.ix_(ia, ia)]
    gbb = gram[np.ix_(ib, ib)]
    gab = gram[np.ix_(ia, ib)]

    def inverse_sqrt(matrix):
        values, vectors = np.linalg.eigh(matrix)
        tolerance = max(matrix.shape) * np.finfo(float).eps * max(
            float(np.max(values)), 1.0)
        factors = np.zeros_like(values)
        factors[values > tolerance] = 1.0 / np.sqrt(values[values > tolerance])
        return (vectors * factors) @ vectors.T

    normalised = inverse_sqrt(gaa) @ gab @ inverse_sqrt(gbb)
    singular = np.linalg.svd(normalised, compute_uv=False)
    return min(float(singular[0]) if singular.size else 0.0, 1.0)


def _geometric(values):
    values = [float(value) for value in values if float(value) > 0.0]
    if not values:
        return None
    return float(math.exp(sum(math.log(value) for value in values) / len(values)))


def _window_average(records, iterations, field):
    selected = [records[k][field] for k in iterations if k in records]
    return float(np.mean(selected)) if selected else None


def restart_verdict(iterations, restart_rows, family_records):
    pre = _geometric(iterations[k]["rho"] for k in PRE_WINDOW if k in iterations)
    post = _geometric(iterations[k]["rho"] for k in POST_WINDOW if k in iterations)
    jump = max(float(row["residual_jump_relative"]) for row in restart_rows)
    cosine = min(float(row["pre_post_cosine"]) for row in restart_rows)
    shifts = []
    for family in PRIMARY_FAMILIES:
        a = _window_average(family_records, PRE_WINDOW, family)
        b = _window_average(family_records, POST_WINDOW, family)
        if a is not None and b is not None:
            shifts.append(abs(b - a))
    mode_shift = max(shifts) if shifts else 0.0
    ratio = post / pre if pre and post else None
    if jump > RESTART_CONTINUITY_TOLERANCE or cosine < 1.0 - 1.0e-10:
        verdict = "RESTART_SENSITIVITY_STRONG"
    elif ratio is None:
        verdict = "UNRESOLVED"
    elif ratio > RESTART_STRONG_RATIO and mode_shift > MODE_COMPOSITION_CHANGE:
        verdict = "RESTART_SENSITIVITY_STRONG"
    elif ratio > RESTART_SUPPORTED_RATIO:
        verdict = "RESTART_SENSITIVITY_SUPPORTED"
    elif pre >= 0.98:
        verdict = "RESIDUAL_TAIL_PREEXISTS_RESTART"
    else:
        verdict = "NO_RESTART_EFFECT_DETECTED"
    return {
        "verdict": verdict,
        "rho_pre_restart": pre,
        "rho_post_restart": post,
        "post_to_pre_ratio": ratio,
        "restart_jump_relative": jump,
        "restart_pre_post_cosine": cosine,
        "maximum_family_fraction_shift": mode_shift,
        "pre_window": list(PRE_WINDOW),
        "post_window": list(POST_WINDOW),
    }


def _classify_tail(iterations, projections, overlaps):
    represented = _window_average(projections, range(51, 61), "represented_fraction")
    if represented is None or represented < MIN_REPRESENTED_FRACTION:
        return "UNCLASSIFIED", None
    window_scores = []
    for window in (PRE_WINDOW, tuple(range(51, 61))):
        scores = {
            family: _window_average(projections, window, family)
            for family in PRIMARY_FAMILIES
        }
        scores = {name: value / max(_window_average(
            projections, window, "represented_fraction") or 0.0, 1.0e-300)
                  for name, value in scores.items() if value is not None}
        window_scores.append(scores)
    candidates = set(window_scores[0]).intersection(window_scores[1])
    stable = [name for name in candidates
              if min(window_scores[0][name], window_scores[1][name])
              >= DOMINANCE_FRACTION]
    if len(stable) != 1:
        return ("MIXED" if stable else "UNCLASSIFIED"), None
    dominant = stable[0]
    for other in PRIMARY_FAMILIES:
        if other == dominant:
            continue
        other_score = max(window_scores[0].get(other, 0.0),
                          window_scores[1].get(other, 0.0))
        key = "|".join(sorted((dominant, other)))
        if other_score >= 0.40 and overlaps.get(key, 0.0) >= OVERLAP_AMBIGUITY:
            return "MIXED", None
    return FAMILY_CLASS[dominant], dominant


def ritz_summary(rows):
    groups = {}
    for row in rows:
        key = (row["case"], int(row["restart_cycle"]))
        groups.setdefault(key, []).append(row)
    output = []
    for (case, cycle), entries in sorted(groups.items()):
        m = max(int(row["column"]) for row in entries) + 1
        hbar = np.zeros((m + 1, m))
        beta = float(entries[0]["cycle_start_residual_norm"])
        for row in entries:
            hbar[int(row["row"]), int(row["column"])] = float(row["value"])
        h = hbar[:m, :]
        ritz = np.linalg.eigvals(h)
        harmonic = []
        if m and abs(hbar[m, m - 1]) > 0.0:
            unit = np.zeros(m)
            unit[-1] = 1.0
            try:
                correction = np.linalg.solve(h.T, unit)
                harmonic_matrix = h + hbar[m, m - 1] ** 2 * np.outer(
                    correction, unit)
                harmonic = np.linalg.eigvals(harmonic_matrix)
            except np.linalg.LinAlgError:
                harmonic = []
        output.append((case, cycle, beta, ritz, harmonic))
    return output


def analyze_case(data):
    gram, fingerprints = gram_matrix(data["gram"])
    correlations = {}
    residual_norms = {}
    for row in data["correlations"]:
        iteration = int(row["iteration"])
        correlations.setdefault(iteration, {})[row["probe"]] = float(row["correlation"])
        residual_norms[iteration] = float(row["residual_norm"])
    iteration_rows = {int(row["iteration"]): row for row in data["iterations"]}
    if set(iteration_rows) != set(range(61)):
        raise ValueError("Stage-E iterations must be exactly 0..60")
    projections = {}
    maximum_reconciliation = 0.0
    physical_energy_pass = True
    density_index = PROBES.index("P15_density_gauge")
    for iteration in range(61):
        if set(correlations.get(iteration, {})) != set(PROBES):
            raise ValueError("incomplete Stage-E correlation block")
        vector = np.array([correlations[iteration][name] for name in PROBES])
        energy = residual_norms[iteration] ** 2
        full = svd_projection(gram, vector, energy)
        maximum_reconciliation = max(maximum_reconciliation,
                                     full["reconciliation_relative"])
        physical_energy_pass = physical_energy_pass and full["physical_energy_pass"]
        record = {
            "represented_fraction": full["represented_fraction"],
            "unexplained_fraction": full["unexplained_fraction"],
            "density_gauge_fraction": vector[density_index] ** 2 / max(energy, 1.0e-300),
        }
        for family, names in FAMILIES.items():
            record[family] = family_projection(gram, vector, energy, names)
        projections[iteration] = record
    overlaps = {}
    for n, first in enumerate(PRIMARY_FAMILIES):
        for second in PRIMARY_FAMILIES[n + 1:]:
            overlaps["|".join(sorted((first, second)))] = family_overlap(
                gram, first, second)
    numeric_iterations = {
        iteration: {
            "rho": float(row["rho"]),
            "residual_direction_cosine": float(row["residual_direction_cosine"]),
            "true_continuity_relative": float(row["true_continuity_relative"]),
        }
        for iteration, row in iteration_rows.items()
    }
    restart = restart_verdict(numeric_iterations, data["restart"], projections)
    classification, dominant = _classify_tail(numeric_iterations, projections, overlaps)
    tail = tuple(range(51, 61))
    return {
        "gram": gram,
        "gram_rank": svd_projection(gram, np.zeros(len(PROBES)), 1.0)["rank"],
        "gram_singular_values": np.linalg.svd(gram, compute_uv=False).tolist(),
        "fingerprints": fingerprints,
        "projections": projections,
        "family_overlap": overlaps,
        "maximum_projection_reconciliation_relative": maximum_reconciliation,
        "projection_physical_energy_pass": physical_energy_pass,
        "tail_classification": classification,
        "dominant_family": dominant,
        "represented_fraction": _window_average(projections, tail, "represented_fraction"),
        "unexplained_fraction": _window_average(projections, tail, "unexplained_fraction"),
        "density_gauge_fraction": _window_average(projections, tail, "density_gauge_fraction"),
        "restart": restart,
        "final_continuity": numeric_iterations[60]["true_continuity_relative"],
    }


def emit_families(args):
    payload = {
        "schema": "strict-ala-stage-E-families-v1",
        "declared_before_execution": True,
        "families": {name: list(values) for name, values in FAMILIES.items()},
        "primary_modal_families": list(PRIMARY_FAMILIES),
        "thresholds": {
            "dominance_fraction": DOMINANCE_FRACTION,
            "minimum_represented_fraction": MIN_REPRESENTED_FRACTION,
            "overlap_ambiguity": OVERLAP_AMBIGUITY,
            "mode_composition_change": MODE_COMPOSITION_CHANGE,
            "restart_supported_ratio": RESTART_SUPPORTED_RATIO,
            "restart_strong_ratio": RESTART_STRONG_RATIO,
            "restart_continuity_tolerance": RESTART_CONTINUITY_TOLERANCE,
        },
    }
    _write_json(args.output, payload)
    return 0


def analyze(args):
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    data0 = load_case(args.e0_prefix)
    data1 = load_case(args.e1_prefix)
    case0 = analyze_case(data0)
    case1 = analyze_case(data1)
    gram_difference = float(np.max(np.abs(case0["gram"] - case1["gram"])))
    fingerprint_difference = max(
        abs(case0["fingerprints"][name] - case1["fingerprints"][name])
        for name in PROBES)

    family_names = list(PRIMARY_FAMILIES)
    e0_vector = np.array([_window_average(case0["projections"], range(51, 61), name)
                          for name in family_names])
    e1_vector = np.array([_window_average(case1["projections"], range(51, 61), name)
                          for name in family_names])
    family_shift = float(np.max(np.abs(e1_vector - e0_vector)))
    composition_cosine = float(e0_vector @ e1_vector / max(
        np.linalg.norm(e0_vector) * np.linalg.norm(e1_vector), 1.0e-300))
    changes = family_shift > MODE_COMPOSITION_CHANGE or composition_cosine < 0.95
    amplitude_ratio = case1["final_continuity"] / max(case0["final_continuity"], 1.0e-300)
    if changes and abs(amplitude_ratio - 1.0) > 0.10:
        primary_effect = "MIXED"
    elif changes and family_shift > MODE_COMPOSITION_CHANGE:
        primary_effect = "MODE_REDISTRIBUTION"
    elif changes:
        primary_effect = "DIRECTION"
    else:
        primary_effect = "AMPLITUDE"

    def public(case):
        return {
            "tail_classification": case["tail_classification"],
            "represented_fraction": case["represented_fraction"],
            "unexplained_fraction": case["unexplained_fraction"],
            "dominant_family": case["dominant_family"],
            "density_gauge_fraction": case["density_gauge_fraction"],
            "restart_verdict": case["restart"]["verdict"],
            "restart_evidence": case["restart"],
            "gram_rank": case["gram_rank"],
            "gram_singular_values": case["gram_singular_values"],
            "maximum_projection_reconciliation_relative":
                case["maximum_projection_reconciliation_relative"],
            "family_overlap": case["family_overlap"],
        }

    classification = {
        "schema": "strict-ala-stage-E-classification-v1",
        "case_E0": public(case0),
        "case_E1": public(case1),
        "comparison": {
            "gram_maximum_absolute_difference": gram_difference,
            "probe_fingerprint_maximum_absolute_difference":
                fingerprint_difference,
            "schwarz_changes_mode_composition": changes,
            "schwarz_primary_effect": primary_effect,
            "late_family_fraction_maximum_difference": family_shift,
            "late_family_composition_cosine": composition_cosine,
            "final_continuity_amplitude_ratio_E1_over_E0": amplitude_ratio,
        },
        "thresholds": {
            "svd_rank_rule": "max(shape)*eps*max(smax,1)",
            "dominance_fraction": DOMINANCE_FRACTION,
            "minimum_represented_fraction": MIN_REPRESENTED_FRACTION,
            "projection_reconciliation": 1.0e-10,
            "stage_B_MPI_envelope": STAGE_B_MPI_ENVELOPE,
        },
    }
    _write_json(output / "strict_ala_stage_E_classification.json", classification)

    restart_verdicts = {case0["restart"]["verdict"], case1["restart"]["verdict"]}
    if "RESTART_SENSITIVITY_STRONG" in restart_verdicts:
        path = "KRYLOV_RESTART_PATH"
    elif case0["tail_classification"] == "LOCAL":
        path = "LOCAL_PRECONDITIONER_PATH"
    elif case0["tail_classification"] == "GLOBAL":
        path = "GLOBAL_COARSE_PATH"
    elif case0["tail_classification"] == "GAUGE_NEAR_NULL":
        path = "GAUGE_NEAR_NULL_PATH"
    elif case0["tail_classification"] == "MIXED":
        path = "MIXED_PATH"
    else:
        path = "UNRESOLVED"
    next_action = {
        "schema": "strict-ala-stage-E-next-action-v1",
        "next_authorized_path": path,
        "automatic_solver_change_authorized": False,
        "basis_sufficient": min(case0["represented_fraction"],
                                case1["represented_fraction"])
                            >= MIN_REPRESENTED_FRACTION,
    }
    _write_json(output / "strict_ala_stage_E_next_action.json", next_action)

    with (output / "strict_ala_stage_E_projection.csv").open("w", newline="") as stream:
        fields = ["case", "iteration", "represented_fraction",
                  "unexplained_fraction", "density_gauge_fraction"] + list(FAMILIES)
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for name, case in (("E0", case0), ("E1", case1)):
            for iteration, values in sorted(case["projections"].items()):
                writer.writerow({"case": name, "iteration": iteration, **values})

    with (output / "strict_ala_stage_E_ritz.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("case", "restart_cycle", "kind", "index", "real", "imag",
                         "cycle_start_residual_norm"))
        for data in (data0, data1):
            for case, cycle, beta, ritz, harmonic in ritz_summary(data["hessenberg"]):
                for kind, values in (("RITZ", ritz), ("HARMONIC_RITZ", harmonic)):
                    for index, value in enumerate(values):
                        writer.writerow((case, cycle, kind, index,
                                         float(np.real(value)), float(np.imag(value)), beta))
    summary = {
        "schema": "strict-ala-stage-E-analysis-summary-v1",
        "probe_gram_pass": (
            gram_difference <= STAGE_B_MPI_ENVELOPE
            and fingerprint_difference <= STAGE_B_MPI_ENVELOPE),
        "projection_reconciliation_pass": max(
            case0["maximum_projection_reconciliation_relative"],
            case1["maximum_projection_reconciliation_relative"])
            <= PROJECTION_RECONCILIATION_TOLERANCE
            and case0["projection_physical_energy_pass"]
            and case1["projection_physical_energy_pass"],
        "gram_maximum_absolute_difference": gram_difference,
        "probe_fingerprint_maximum_absolute_difference":
            fingerprint_difference,
        "stage_B_MPI_envelope": STAGE_B_MPI_ENVELOPE,
    }
    _write_json(output / "strict_ala_stage_E_analysis_summary.json", summary)
    return 0 if all((summary["probe_gram_pass"],
                     summary["projection_reconciliation_pass"])) else 1


def _case_complete(log_text, case):
    return (f"STRICT_ALA_STAGE_E_COMPLETE case={case}" in log_text and
            f"STRICT_ALA_STAGE_C_CASE_COMPLETE case={case}" in log_text)


def validate_case_structure(data, case):
    """Fail closed on missing, duplicate, mislabelled, or inconsistent rows."""
    expected_iterations = set(range(61))
    iteration_keys = [int(row["iteration"]) for row in data["iterations"]]
    if len(iteration_keys) != 61 or set(iteration_keys) != expected_iterations:
        return False
    if any(row.get("case") != case for rows in data.values() for row in rows):
        return False
    correlation_keys = [(int(row["iteration"]), row["probe"])
                        for row in data["correlations"]]
    expected_correlations = {(iteration, probe)
                             for iteration in expected_iterations
                             for probe in PROBES}
    if len(correlation_keys) != len(expected_correlations) or \
       set(correlation_keys) != expected_correlations:
        return False
    gram_keys = [(row["probe_i"], row["probe_j"]) for row in data["gram"]]
    expected_gram = {(first, second) for first in PROBES for second in PROBES}
    if len(gram_keys) != len(expected_gram) or set(gram_keys) != expected_gram:
        return False
    work_keys = [int(row["iteration"]) for row in data["work"]]
    if len(work_keys) != 61 or set(work_keys) != expected_iterations:
        return False
    restart_keys = [(int(row["iteration"]), row["probe"])
                    for row in data["restart"]]
    if len(restart_keys) != len(PROBES) or \
       set(restart_keys) != {(50, probe) for probe in PROBES}:
        return False

    hessenberg_keys = [(int(row["restart_cycle"]), int(row["column"]),
                        int(row["row"])) for row in data["hessenberg"]]
    expected_hessenberg = set()
    for cycle, columns in ((1, 50), (2, 10)):
        for column in range(columns):
            for row in range(column + 2):
                expected_hessenberg.add((cycle, column, row))
    if len(hessenberg_keys) != len(expected_hessenberg) or \
       set(hessenberg_keys) != expected_hessenberg:
        return False

    counters = ("pressure_schur_actions", "K_gamma_rhs_solves",
                "K_gamma_operator_applications", "velocity_MG_cycles",
                "velocity_smoother_sweeps", "preconditioner_applications")
    ordered_work = sorted(data["work"], key=lambda row: int(row["iteration"]))
    for counter in counters:
        values = [int(row[counter]) for row in ordered_work]
        if any(value < 0 for value in values) or any(
                later < earlier for earlier, later in zip(values, values[1:])):
            return False

    for row in data["iterations"]:
        if row.get("residual_state_guard_pass") != "1":
            return False
        depth_sum = sum(float(row[name]) for name in (
            "depth_0_200_fraction", "depth_200_410_fraction",
            "depth_410_660_fraction", "depth_660_1000_fraction",
            "depth_1000_cmb_fraction"))
        if abs(depth_sum - 1.0) > PROJECTION_RECONCILIATION_TOLERANCE:
            return False
        ratio = float(row["recursive_residual"]) / max(
            float(row["explicit_residual"]), 1.0e-300)
        relative = abs(float(row["recursive_residual"])
                       - float(row["explicit_residual"])) / max(
                           float(row["explicit_residual"]), 1.0e-300)
        if abs(float(row["krylov_residual_ratio"]) - ratio) > 1.0e-12 \
           * max(abs(ratio), 1.0):
            return False
        if abs(float(row["krylov_residual_rel_difference"]) - relative) \
           > 1.0e-12 * max(abs(relative), 1.0):
            return False
    return True


def audit(args):
    analysis_dir = Path(args.analysis_dir)
    classification = json.loads(
        (analysis_dir / "strict_ala_stage_E_classification.json").read_text())
    next_action = json.loads(
        (analysis_dir / "strict_ala_stage_E_next_action.json").read_text())
    summary = json.loads(
        (analysis_dir / "strict_ala_stage_E_analysis_summary.json").read_text())
    manifest = json.loads(Path(args.manifest).read_text())
    data0 = load_case(args.e0_prefix)
    data1 = load_case(args.e1_prefix)
    structure = validate_case_structure(data0, "E0") \
        and validate_case_structure(data1, "E1")
    for data in (data0, data1):
        structure = structure and all(
            _finite(value) for rows in data.values() for row in rows
            for key, value in row.items()
            if key not in {"case", "probe", "probe_i", "probe_j",
                           "aggregation_semantics"})
    logs = [Path(path).read_text(errors="replace") for path in args.case_log]
    e0_complete = len(logs) == 2 and _case_complete(logs[0], "E0")
    e1_complete = len(logs) == 2 and _case_complete(logs[1], "E1")
    fatal_tokens = ("STRICT_ALA_STAGE_C_HIDDEN_FALLBACK", "Fatal error",
                    "MPI_Abort", "MPI failure")
    nonfinite = re.compile(
        r"(?<![A-Za-z])(?:nan|[+-]?inf(?:inity)?)(?![A-Za-z])", re.I)
    clean_logs = not any(
        any(token in text for token in fatal_tokens) or nonfinite.search(text)
        for text in logs)
    binary_unchanged = Path(args.binary_pre).read_bytes() == Path(args.binary_post).read_bytes()
    inputs_unchanged = Path(args.input_pre).read_bytes() == Path(args.input_post).read_bytes()
    source_unchanged = Path(args.source_pre).read_bytes() == Path(args.source_post).read_bytes()
    cfg_single = Path(args.cfg_diff).read_text().strip() == \
        "ala_shallow_patch_preconditioner: off -> on"
    source = manifest.get("source", {})
    identity = bool(manifest.get("provenance_complete") and
                    source.get("branch_verified") and
                    not source.get("scientific_dirty") and
                    not manifest.get("missing"))
    manifest_files_unchanged = bool(manifest.get("files")) and all(
        Path(item.get("path", "")).is_file()
        and _sha256(item["path"]) == item.get("sha256")
        for item in manifest.get("files", []))
    experiment_valid = all((structure, e0_complete, e1_complete, clean_logs,
                            binary_unchanged, inputs_unchanged,
                            source_unchanged, manifest_files_unchanged,
                            cfg_single, identity,
                            summary.get("probe_gram_pass"),
                            summary.get("projection_reconciliation_pass")))
    payload = {
        "schema": "strict-ala-stage-E-final-audit-v1",
        "experiment_valid": experiment_valid,
        "source_binary_input_identity_pass": (
            identity and source_unchanged and binary_unchanged
            and inputs_unchanged and manifest_files_unchanged),
        "diagnostics_observation_only_pass": structure,
        "probe_gram_pass": bool(summary.get("probe_gram_pass")),
        "projection_reconciliation_pass": bool(
            summary.get("projection_reconciliation_pass")),
        "E0_complete": e0_complete,
        "E1_complete": e1_complete,
        "E0_tail_classification": classification["case_E0"]["tail_classification"],
        "E1_tail_classification": classification["case_E1"]["tail_classification"],
        "E0_restart_verdict": classification["case_E0"]["restart_verdict"],
        "E1_restart_verdict": classification["case_E1"]["restart_verdict"],
        "density_near_null_dominant_E0":
            classification["case_E0"]["tail_classification"] == "GAUGE_NEAR_NULL",
        "density_near_null_dominant_E1":
            classification["case_E1"]["tail_classification"] == "GAUGE_NEAR_NULL",
        "schwarz_changes_mode_composition":
            classification["comparison"]["schwarz_changes_mode_composition"],
        "next_authorized_path": next_action["next_authorized_path"],
        "automatic_solver_change_authorized": False,
        "raw_schema_complete": structure,
        "case_logs_clean": clean_logs,
        "cfg_single_variable_verified": cfg_single,
        "binary_unchanged": binary_unchanged,
        "inputs_unchanged": inputs_unchanged,
        "source_unchanged": source_unchanged,
        "manifest_files_unchanged": manifest_files_unchanged,
    }
    _write_json(args.output, payload)
    return 0 if experiment_valid else 1


def parser():
    root = argparse.ArgumentParser()
    sub = root.add_subparsers(dest="command", required=True)
    families = sub.add_parser("families")
    families.add_argument("--output", required=True)
    families.set_defaults(function=emit_families)
    analyse = sub.add_parser("analyze")
    analyse.add_argument("--e0-prefix", required=True)
    analyse.add_argument("--e1-prefix", required=True)
    analyse.add_argument("--output-dir", required=True)
    analyse.set_defaults(function=analyze)
    final = sub.add_parser("audit")
    final.add_argument("--analysis-dir", required=True)
    final.add_argument("--e0-prefix", required=True)
    final.add_argument("--e1-prefix", required=True)
    final.add_argument("--manifest", required=True)
    final.add_argument("--cfg-diff", required=True)
    final.add_argument("--binary-pre", required=True)
    final.add_argument("--binary-post", required=True)
    final.add_argument("--input-pre", required=True)
    final.add_argument("--input-post", required=True)
    final.add_argument("--source-pre", required=True)
    final.add_argument("--source-post", required=True)
    final.add_argument("--case-log", action="append", required=True)
    final.add_argument("--output", required=True)
    final.set_defaults(function=audit)
    return root


def main():
    args = parser().parse_args()
    return args.function(args)


if __name__ == "__main__":
    raise SystemExit(main())
