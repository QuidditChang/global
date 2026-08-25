#!/usr/bin/env python3
"""Fail-closed offline analysis for strict-ALA Stage E2 outputs."""

import argparse
import csv
import hashlib
import json
import math
import re
from collections import defaultdict
from pathlib import Path

SCHEDULE = (0, 5, 10, 20, 30, 35, 40, 42, 44, 46, 48, 49, 50,
            51, 52, 54, 56, 58, 60)
LATE_COUNT = 13
PROVISIONAL = "PROVISIONAL_DIAGNOSTIC"
GRAM_RELATIVE_TOLERANCE = 1.0e-12
POD_RECONCILIATION_TOLERANCE = 1.0e-8
POD_ORTHOGONALITY_TOLERANCE = 1.0e-8
TRUE_MODE_TIGHT_SOLVE_LIMIT = 2.0e-10
DOMINANCE_FRACTION = 0.60
MAP_MATERIAL_FRACTION = 0.05
MPI_SUSPECT_RATIO = 1.5
MPI_STRONG_RATIO = 2.0

SCHEMAS = {
    "snapshots_manifest": ("case", "iteration", "residual_norm",
        "global_checksum", "storage_identifier", "ranks", "local_npno"),
    "snapshot_gram": ("case_i", "iteration_i", "case_j", "iteration_j",
        "inner_product"),
    "pod": ("case", "basis_type", "mode", "singular_value",
        "energy_fraction", "cumulative_fraction", "N50", "N80", "N90",
        "N95", "weighting"),
    "joint_pod_coefficients": ("case", "iteration", "pod_mode",
        "coefficient", "contraction", "sign_change"),
    "subspace_angles": ("energy_level", "subspace_dimension", "angle_index",
        "angle_degrees", "cosine", "threshold_status"),
    "depth_scale_spectrum": ("case", "object_type", "object_id", "depth_band",
        "horizontal_scale_band", "energy", "normalized_fraction", "method"),
    "mpi_distance": ("case", "object_type", "object_id", "mpi_distance_band",
        "element_count", "energy", "energy_per_element"),
    "true_mode_preconditioner": ("mode", "map", "mode_energy_fraction",
        "cumulative_energy", "E_P", "cosine", "alpha_opt", "Pq_norm",
        "SPq_norm", "qTPq", "qTSPq", "positive_qTPq", "positive_qTSPq",
        "tight_solve_achieved", "map_semantics"),
}


def fnv1a64(payload):
    value = 1469598103934665603
    for byte in payload:
        value ^= byte
        value = (value * 1099511628211) & ((1 << 64) - 1)
    return value


def read_csv(path, schema):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing or empty CSV: {path}")
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != tuple(schema):
            raise ValueError(f"schema mismatch in {path}: {reader.fieldnames}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"header-only CSV: {path}")
    return rows


def finite(row, fields):
    for field in fields:
        value = float(row[field])
        if not math.isfinite(value):
            raise ValueError(f"non-finite {field}: {row}")


def pod_metrics(matrix, normalized=False):
    """Small deterministic method-of-snapshots reference for local tests."""
    import numpy as np
    r = np.asarray(matrix, dtype=float)
    if r.ndim != 2 or r.shape[1] == 0 or not np.isfinite(r).all():
        raise ValueError("invalid snapshot matrix")
    if normalized:
        norms = np.linalg.norm(r, axis=0)
        if np.any(norms <= 0):
            raise ValueError("zero snapshot in direction POD")
        r = r / norms
    gram = r.T @ r
    values, vectors = np.linalg.eigh(gram)
    order = np.argsort(values)[::-1]
    values = np.maximum(values[order], 0.0)
    vectors = vectors[:, order]
    tolerance = max(gram.shape) * np.finfo(float).eps * max(values[0], 1.0)
    rank = int(np.sum(values > tolerance))
    fractions = values / max(values.sum(), np.finfo(float).tiny)
    modes = r @ vectors[:, :rank] / np.sqrt(values[:rank])
    return {"gram": gram, "values": values, "vectors": vectors,
            "fractions": fractions, "rank": rank, "modes": modes}


def principal_angles(basis_a, basis_b):
    import numpy as np
    qa, _ = np.linalg.qr(np.asarray(basis_a, dtype=float))
    qb, _ = np.linalg.qr(np.asarray(basis_b, dtype=float))
    singular = np.linalg.svd(qa.T @ qb, compute_uv=False)
    singular = np.clip(singular, 0.0, 1.0)
    return np.degrees(np.arccos(singular))


def nested_projection_bands(field, factors=(8, 4, 2)):
    """Orthogonal energy ledger for nested block-constant 2-D spaces."""
    import numpy as np
    values = np.asarray(field, dtype=float)
    if values.ndim != 2:
        raise ValueError("horizontal field must be two-dimensional")
    projections = []
    for factor in factors:
        if values.shape[0] % factor or values.shape[1] % factor:
            raise ValueError("block factor does not divide synthetic mesh")
        projected = np.empty_like(values)
        for y in range(0, values.shape[0], factor):
            for x in range(0, values.shape[1], factor):
                block = values[y:y + factor, x:x + factor]
                projected[y:y + factor, x:x + factor] = block.mean()
        projections.append(float(np.sum(projected * projected)))
    total = float(np.sum(values * values))
    bands = (projections[0], projections[1] - projections[0],
             projections[2] - projections[1], total - projections[2])
    return np.maximum(bands, 0.0)


def mpi_concentration(energy, counts):
    if counts[0] <= 0 or counts[3] <= 0:
        raise ValueError("empty MPI distance population")
    return (energy[0] / counts[0]) / max(energy[3] / counts[3], 1.0e-300)


def preconditioner_metrics(operator, preconditioner, mode):
    """Reference true-mode metric calculation for a small SPD test system."""
    import numpy as np
    q = np.asarray(mode, dtype=float)
    q = q / np.linalg.norm(q)
    z = np.asarray(preconditioner, dtype=float) @ q
    action = np.asarray(operator, dtype=float) @ z
    action2 = float(action @ action)
    q_action = float(q @ action)
    return {"E_P": float(np.linalg.norm(q - action)),
            "cosine": q_action / max(math.sqrt(action2), 1e-300),
            "alpha_opt": q_action / max(action2, 1e-300),
            "qTPq": float(q @ z), "qTSPq": q_action}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_trajectory(current, baseline, tolerance):
    fields = ("iteration", "true_continuity_relative", "true_momentum_relative",
              "recursive_residual", "explicit_residual", "velocity_MG_cycles")
    def generic(path):
        with Path(path).open(newline="") as stream:
            return list(csv.DictReader(stream))
    a, b = generic(current), generic(baseline)
    if len(a) != len(b) or not a:
        raise ValueError("trajectory row-count mismatch")
    maximum = 0.0
    for left, right in zip(a, b):
        if left["iteration"] != right["iteration"]:
            raise ValueError("trajectory iteration mismatch")
        for field in fields[1:]:
            x, y = float(left[field]), float(right[field])
            if not math.isfinite(x) or not math.isfinite(y):
                raise ValueError("non-finite trajectory value")
            difference = abs(x - y) / max(abs(y), 1.0)
            maximum = max(maximum, difference)
    if maximum > tolerance:
        raise ValueError(f"trajectory drift {maximum} exceeds {tolerance}")
    return maximum


def load_all(root):
    root = Path(root)
    result = {}
    for key, schema in SCHEMAS.items():
        result[key] = read_csv(root / f"strict_ala_stage_E2_{key}.csv", schema)
    return result


def validate(data):
    manifest = data["snapshots_manifest"]
    if len(manifest) != 38:
        raise ValueError("snapshot manifest must contain exactly 38 rows")
    labels = []
    manifest_keys = set()
    local_npno = set()
    for case in ("E0", "E1"):
        rows = [r for r in manifest if r["case"] == case]
        if tuple(int(r["iteration"]) for r in rows) != SCHEDULE:
            raise ValueError(f"snapshot schedule mismatch for {case}")
        if any(int(r["ranks"]) != 384 for r in rows):
            raise ValueError("Stage E2 frozen layout is not 384 ranks")
        for row in rows:
            finite(row, ("residual_norm",))
            if float(row["residual_norm"]) <= 0.0:
                raise ValueError("non-positive snapshot norm")
            try:
                int(row["global_checksum"], 16)
            except ValueError as error:
                raise ValueError("invalid snapshot checksum") from error
            if len(row["global_checksum"]) != 16:
                raise ValueError("invalid snapshot checksum")
            key = (case, int(row["iteration"]))
            if key in manifest_keys:
                raise ValueError("duplicate snapshot manifest key")
            manifest_keys.add(key)
            labels.append(key)
            if int(row["local_npno"]) <= 0:
                raise ValueError("invalid local pressure size")
            local_npno.add(int(row["local_npno"]))
            if "%06d" not in row["storage_identifier"]:
                raise ValueError("snapshot storage identifier is not rank-addressable")
    if len(local_npno) != 1:
        raise ValueError("inconsistent local pressure size")

    gram = data["snapshot_gram"]
    if len(gram) != 38 * 38:
        raise ValueError("joint Gram must be exactly 38 by 38")
    values = {}
    for row in gram:
        finite(row, ("inner_product",))
        key = (row["case_i"], int(row["iteration_i"]),
               row["case_j"], int(row["iteration_j"]))
        if key in values:
            raise ValueError("duplicate Gram key")
        values[key] = float(row["inner_product"])
    for key, value in values.items():
        reverse = (key[2], key[3], key[0], key[1])
        if reverse not in values or abs(value - values[reverse]) > 1e-12 * max(1.0, abs(value)):
            raise ValueError("asymmetric snapshot Gram")
    expected_gram_keys = {
        (left[0], left[1], right[0], right[1])
        for left in labels for right in labels
    }
    if set(values) != expected_gram_keys:
        raise ValueError("snapshot Gram key set is incomplete")

    expected_pod = {("E0", "full_physical"): 19,
                    ("E0", "late_physical"): 13,
                    ("E0", "late_direction"): 13,
                    ("E1", "full_physical"): 19,
                    ("E1", "late_physical"): 13,
                    ("E1", "late_direction"): 13,
                    ("JOINT", "late_physical"): 26,
                    ("JOINT", "late_direction_case_balanced"): 26}
    groups = defaultdict(list)
    for row in data["pod"]:
        finite(row, ("singular_value", "energy_fraction", "cumulative_fraction"))
        groups[(row["case"], row["basis_type"])].append(row)
    if {key: len(rows) for key, rows in groups.items()} != expected_pod:
        raise ValueError("POD group/count contract mismatch")
    for rows in groups.values():
        if [int(r["mode"]) for r in rows] != list(range(1, len(rows) + 1)):
            raise ValueError("POD modes are not consecutive")
        for row in rows:
            finite(row, ("singular_value", "energy_fraction",
                         "cumulative_fraction"))
            if min(float(row["singular_value"]),
                   float(row["energy_fraction"]),
                   float(row["cumulative_fraction"])) < 0.0:
                raise ValueError("negative POD quantity")
        if abs(sum(float(r["energy_fraction"]) for r in rows) - 1.0) > 1e-8:
            raise ValueError("POD energy does not reconcile")
        if abs(float(rows[-1]["cumulative_fraction"]) - 1.0) > 1e-8:
            raise ValueError("POD cumulative energy does not end at one")
        cumulative = 0.0
        for row in rows:
            cumulative += float(row["energy_fraction"])
            if abs(cumulative - float(row["cumulative_fraction"])) > 1e-8:
                raise ValueError("POD cumulative energy is inconsistent")

    true_rows = data["true_mode_preconditioner"]
    modes = sorted({int(r["mode"]) for r in true_rows})
    if not modes or modes != list(range(1, max(modes) + 1)) or max(modes) > 6:
        raise ValueError("invalid selected-mode sequence")
    if len(true_rows) != 3 * len(modes):
        raise ValueError("true-mode map matrix incomplete")
    for mode in modes:
        if {r["map"] for r in true_rows if int(r["mode"]) == mode} != {
                "BPI", "CONFIGURED", "PURE_SCHWARZ"}:
            raise ValueError("true-mode map set incomplete")
    for row in true_rows:
        finite(row, ("mode_energy_fraction", "cumulative_energy", "E_P",
                     "cosine", "alpha_opt", "Pq_norm", "SPq_norm", "qTPq",
                     "qTSPq", "tight_solve_achieved"))
        if row["positive_qTPq"] != ("1" if float(row["qTPq"]) > 0 else "0"):
            raise ValueError("qTPq positivity flag is inconsistent")
        if row["positive_qTSPq"] != ("1" if float(row["qTSPq"]) > 0 else "0"):
            raise ValueError("qTSPq positivity flag is inconsistent")
        if float(row["qTPq"]) <= 0.0:
            raise ValueError("pressure map is not positive on a selected mode")
        if float(row["tight_solve_achieved"]) > TRUE_MODE_TIGHT_SOLVE_LIMIT:
            raise ValueError("true-mode Schur action did not meet tight-solve contract")
        expected_semantics = ("configured_minus_BPI_exact_component"
                              if row["map"] == "PURE_SCHWARZ"
                              else "complete_map")
        if row["map_semantics"] != expected_semantics:
            raise ValueError("true-mode map semantics mismatch")

    coefficients = data["joint_pod_coefficients"]
    if len(coefficients) != 38 * len(modes):
        raise ValueError("joint POD coefficient trajectory incomplete")
    coefficient_keys = set()
    for row in coefficients:
        finite(row, ("coefficient", "contraction"))
        key = (row["case"], int(row["iteration"]), int(row["pod_mode"]))
        if key in coefficient_keys:
            raise ValueError("duplicate joint POD coefficient")
        coefficient_keys.add(key)
    expected_coefficient_keys = {
        (case, iteration, mode) for case, iteration in labels for mode in modes
    }
    if coefficient_keys != expected_coefficient_keys:
        raise ValueError("joint POD coefficient keys are incomplete")

    for key in ("depth_scale_spectrum", "mpi_distance"):
        for row in data[key]:
            numeric = ("energy", "normalized_fraction") if key.startswith("depth") \
                else ("element_count", "energy", "energy_per_element")
            finite(row, numeric)
    depth_groups = defaultdict(float)
    for row in data["depth_scale_spectrum"]:
        depth_groups[(row["case"], row["object_type"], row["object_id"])] += \
            float(row["normalized_fraction"])
        if row["method"] != "geometry_native_nested_v1":
            raise ValueError("unrecognized scale method")
    if any(abs(value - 1.0) > 1e-8 for value in depth_groups.values()):
        raise ValueError("depth-scale energy does not reconcile")
    mpi_groups = defaultdict(set)
    for row in data["mpi_distance"]:
        mpi_groups[(row["case"], row["object_type"], row["object_id"])].add(
            row["mpi_distance_band"])
    if any(bands != {"0", "1", "2", "3+"} for bands in mpi_groups.values()):
        raise ValueError("MPI-distance bands incomplete")
    return modes


def numerical_validation(data, modes):
    """Independently validate the C-side method-of-snapshots products."""
    import numpy as np

    manifest = data["snapshots_manifest"]
    labels = [(row["case"], int(row["iteration"])) for row in manifest]
    index = {label: i for i, label in enumerate(labels)}
    gram = np.zeros((len(labels), len(labels)))
    for row in data["snapshot_gram"]:
        i = index[(row["case_i"], int(row["iteration_i"]))]
        j = index[(row["case_j"], int(row["iteration_j"]))]
        gram[i, j] = float(row["inner_product"])
    symmetry = float(np.max(np.abs(gram - gram.T)))
    scale = max(float(np.max(np.abs(gram))), 1.0)
    eigenvalues = np.linalg.eigvalsh((gram + gram.T) * 0.5)
    psd_tolerance = len(labels) * np.finfo(float).eps * max(
        float(np.max(np.abs(eigenvalues))), 1.0)
    minimum_eigenvalue = float(eigenvalues[0])
    if symmetry > GRAM_RELATIVE_TOLERANCE * scale:
        raise ValueError("joint snapshot Gram fails symmetry validation")
    if minimum_eigenvalue < -psd_tolerance:
        raise ValueError("joint snapshot Gram is not positive semidefinite")

    pod_groups = defaultdict(list)
    for row in data["pod"]:
        pod_groups[(row["case"], row["basis_type"])].append(row)

    def subset(case, late=False, joint=False):
        return [i for i, (label_case, iteration) in enumerate(labels)
                if (joint or label_case == case) and (not late or iteration >= 40)]

    definitions = {
        ("E0", "full_physical"): (subset("E0"), False),
        ("E0", "late_physical"): (subset("E0", True), False),
        ("E0", "late_direction"): (subset("E0", True), True),
        ("E1", "full_physical"): (subset("E1"), False),
        ("E1", "late_physical"): (subset("E1", True), False),
        ("E1", "late_direction"): (subset("E1", True), True),
        ("JOINT", "late_physical"): (subset("", True, True), False),
        ("JOINT", "late_direction_case_balanced"):
            (subset("", True, True), True),
    }
    maximum_energy_defect = 0.0
    for key, (indices, normalized) in definitions.items():
        block = gram[np.ix_(indices, indices)]
        if normalized:
            diagonal = np.sqrt(np.maximum(np.diag(block), 0.0))
            if np.any(diagonal <= 0.0):
                raise ValueError("zero norm in direction POD")
            block = block / np.outer(diagonal, diagonal)
        expected = float(np.trace(block))
        reported = sum(float(row["singular_value"]) ** 2
                       for row in pod_groups[key])
        defect = abs(reported - expected) / max(abs(expected), 1.0)
        maximum_energy_defect = max(maximum_energy_defect, defect)
        if defect > POD_RECONCILIATION_TOLERANCE:
            raise ValueError(f"POD energy reconstruction failed for {key}")

    normalized_gram = gram.copy()
    norms = np.sqrt(np.maximum(np.diag(normalized_gram), 0.0))
    normalized_gram /= np.outer(norms, norms)
    pseudo_inverse = np.linalg.pinv(normalized_gram,
                                    rcond=len(labels) * np.finfo(float).eps)
    coefficient = np.zeros((len(labels), len(modes)))
    for row in data["joint_pod_coefficients"]:
        coefficient[index[(row["case"], int(row["iteration"]))],
                    int(row["pod_mode"]) - 1] = float(row["coefficient"])
    direction_coefficient = coefficient / norms[:, None]
    mode_gram = direction_coefficient.T @ pseudo_inverse @ direction_coefficient
    orthogonality = float(np.max(np.abs(mode_gram - np.eye(len(modes)))))
    if orthogonality > POD_ORTHOGONALITY_TOLERANCE:
        raise ValueError("selected joint POD modes fail orthonormality validation")
    return {
        "schema": "strict-ala-stage-E2-numerical-validation-v1",
        "gram_symmetry_absolute": symmetry,
        "gram_minimum_eigenvalue": minimum_eigenvalue,
        "gram_psd_tolerance": float(psd_tolerance),
        "maximum_POD_energy_reconciliation_relative": maximum_energy_defect,
        "selected_mode_orthogonality_maximum_absolute": orthogonality,
        "status": "PASS",
    }


def classify(data, modes, stage_e):
    pod_groups = defaultdict(list)
    for row in data["pod"]:
        pod_groups[(row["case"], row["basis_type"])].append(row)

    def dimensions(rows):
        return tuple(next(i + 1 for i, row in enumerate(rows)
                          if float(row["cumulative_fraction"]) >= level)
                     for level in (.8, .9, .95))

    joint = pod_groups[("JOINT", "late_direction_case_balanced")]
    n80, n90, n95 = dimensions(joint)
    dimensionality = "LOW_DIMENSIONAL" if n90 <= 6 else "BROAD_SUBSPACE"
    angles = data["subspace_angles"]
    maximum_angle = max(float(r["angle_degrees"]) for r in angles)
    subspace = ("SAME_SUBSPACE" if maximum_angle <= 10 else
                "MOSTLY_SAME" if maximum_angle <= 30 else
                "MATERIAL_SUBSPACE_CHANGE")

    def spatial(case):
        basis = ("late_direction_case_balanced" if case == "JOINT"
                 else "late_direction")
        fractions = {int(row["mode"]): float(row["energy_fraction"])
                     for row in pod_groups[(case, basis)]}
        rows = [row for row in data["depth_scale_spectrum"]
                if row["case"] == case and row["object_type"] == "pod_mode"]
        present = sorted({int(row["object_id"]) for row in rows})
        if not present:
            raise ValueError(f"missing POD spatial rows for {case}")
        covered = sum(fractions[mode] for mode in present)
        scale = defaultdict(float)
        depth = defaultdict(float)
        for row in rows:
            weight = fractions[int(row["object_id"])] / max(covered, 1e-300)
            contribution = weight * float(row["normalized_fraction"])
            scale[row["horizontal_scale_band"]] += contribution
            depth[row["depth_band"]] += contribution
        dominant_scale, scale_fraction = max(scale.items(), key=lambda item: item[1])
        scale_verdict = (dominant_scale if scale_fraction >= DOMINANCE_FRACTION
                         else "BROAD_HORIZONTAL_SPECTRUM")
        shallow = depth["0_200"] + depth["200_410"]
        transition = depth["410_660"]
        deep = depth["660_1000"] + depth["1000_cmb"]
        depth_values = {"SHALLOW": shallow, "TRANSITION": transition,
                        "DEEP": deep}
        depth_name, depth_fraction = max(depth_values.items(),
                                         key=lambda item: item[1])
        depth_verdict = depth_name if depth_fraction >= DOMINANCE_FRACTION else "MIXED"
        return {"horizontal_scale": scale_verdict,
                "horizontal_scale_fractions": dict(sorted(scale.items())),
                "depth_location": depth_verdict,
                "depth_location_fractions": depth_values,
                "selected_mode_energy_coverage": covered}

    def mpi_metrics(case):
        basis = ("late_direction_case_balanced" if case == "JOINT"
                 else "late_direction")
        fractions = {int(row["mode"]): float(row["energy_fraction"])
                     for row in pod_groups[(case, basis)]}
        rows = [row for row in data["mpi_distance"]
                if row["case"] == case and row["object_type"] == "pod_mode"]
        by_mode = defaultdict(dict)
        for row in rows:
            by_mode[int(row["object_id"])][row["mpi_distance_band"]] = row
        if not by_mode:
            raise ValueError(f"missing POD MPI-distance rows for {case}")
        covered = sum(fractions[mode] for mode in by_mode)
        boundary = interior = 0.0
        concentrations = []
        for mode, mode_rows in by_mode.items():
            weight = fractions[mode] / max(covered, 1e-300)
            b = float(mode_rows["0"]["energy_per_element"])
            i = float(mode_rows["3+"]["energy_per_element"])
            boundary += weight * b
            interior += weight * i
            concentrations.append(b / max(i, 1e-300))
        concentration = boundary / max(interior, 1e-300)
        verdict = ("MPI_INTERFACE_STRONGLY_SUSPECT"
                   if concentration >= MPI_STRONG_RATIO else
                   "MPI_INTERFACE_SUSPECT"
                   if concentration >= MPI_SUSPECT_RATIO else
                   "MPI_INTERFACE_NOT_SUPPORTED")
        return {"MPI_interface_verdict": verdict,
                "MPI_interface_energy_density_ratio": concentration,
                "maximum_individual_mode_MPI_ratio": max(concentrations)}

    spatial_results = {case: spatial(case) for case in ("E0", "E1", "JOINT")}
    mpi_results = {case: mpi_metrics(case) for case in ("E0", "E1", "JOINT")}
    scale_verdict = spatial_results["JOINT"]["horizontal_scale"]
    mpi_verdict = mpi_results["JOINT"]["MPI_interface_verdict"]
    maximum_concentration = mpi_results["JOINT"][
        "maximum_individual_mode_MPI_ratio"]

    maps = defaultdict(dict)
    for row in data["true_mode_preconditioner"]:
        maps[int(row["mode"])][row["map"]] = float(row["E_P"])
    weighted_bpi = sum(maps[m]["BPI"] * float(joint[m - 1]["energy_fraction"])
                       for m in modes)
    weighted_configured = sum(maps[m]["CONFIGURED"] *
                              float(joint[m - 1]["energy_fraction"]) for m in modes)
    mode_ratios = [maps[m]["CONFIGURED"] / max(maps[m]["BPI"], 1e-300)
                   for m in modes]
    improves = [ratio < 1.0 - MAP_MATERIAL_FRACTION for ratio in mode_ratios]
    degrades = [ratio > 1.0 + MAP_MATERIAL_FRACTION for ratio in mode_ratios]
    if any(improves) and any(degrades):
        configured_effect = "MODE_DEPENDENT"
    elif weighted_configured < (1.0 - MAP_MATERIAL_FRACTION) * weighted_bpi:
        configured_effect = "USEFUL_ON_DOMINANT_MODES"
    elif weighted_configured > (1.0 + MAP_MATERIAL_FRACTION) * weighted_bpi:
        configured_effect = "DESTRUCTIVE_ON_DOMINANT_MODES"
    else:
        configured_effect = "WEAK_ON_DOMINANT_MODES"

    coefficient_lookup = defaultdict(dict)
    for row in data["joint_pod_coefficients"]:
        coefficient_lookup[(row["case"], int(row["iteration"]))][
            int(row["pod_mode"])] = float(row["coefficient"])
    restart_metrics = {}
    for case in ("E0", "E1"):
        vectors = {iteration: [coefficient_lookup[(case, iteration)][mode]
                               for mode in modes] for iteration in (49, 50, 51)}
        def jump(left, right):
            numerator = math.sqrt(sum((a - b) ** 2 for a, b in zip(left, right)))
            denominator = max(math.sqrt(sum(a * a for a in left)), 1e-300)
            return numerator / denominator
        restart_metrics[case] = {
            "pod_jump_49_to_50": jump(vectors[49], vectors[50]),
            "pod_jump_50_to_51": jump(vectors[50], vectors[51]),
            "prior_true_residual_restart_verdict": stage_e[f"case_{case}"]["restart_verdict"],
        }
    density_max = max(stage_e["case_E0"]["density_gauge_fraction"],
                      stage_e["case_E1"]["density_gauge_fraction"])

    local_scale = scale_verdict in ("PATCH_LOCAL", "HIGH_FREQUENCY",
                                    "LONG_INTERMEDIATE")
    if dimensionality == "BROAD_SUBSPACE":
        next_path = "MULTILEVEL_PRESSURE_PATH"
    elif local_scale and mpi_verdict == "MPI_INTERFACE_STRONGLY_SUSPECT":
        next_path = "MPI_OVERLAP_PATH"
    elif scale_verdict == "GLOBAL_VERY_LONG":
        next_path = "GLOBAL_COARSE_PATH"
    elif dimensionality == "LOW_DIMENSIONAL" and subspace in (
            "SAME_SUBSPACE", "MOSTLY_SAME"):
        next_path = "TARGETED_LOW_RANK_PATH"
    elif scale_verdict == "LONG_INTERMEDIATE":
        next_path = "PRECONDITIONER_REDESIGN_PATH"
    elif scale_verdict in ("PATCH_LOCAL", "HIGH_FREQUENCY"):
        next_path = "LOCAL_SCHWARZ_PATH"
    elif configured_effect in ("DESTRUCTIVE_ON_DOMINANT_MODES", "MODE_DEPENDENT"):
        next_path = "PRECONDITIONER_REDESIGN_PATH"
    else:
        next_path = "UNRESOLVED"
    case_results = {}
    for case in ("E0", "E1"):
        case_n80, case_n90, case_n95 = dimensions(
            pod_groups[(case, "late_direction")])
        case_results[case] = {
            "late_subspace": "LOW_DIMENSIONAL" if case_n90 <= 6 else "BROAD_SUBSPACE",
            "N80": case_n80, "N90": case_n90, "N95": case_n95,
            **spatial_results[case], **mpi_results[case],
        }
    return {
        "schema": "strict-ala-stage-E2-classification-v1",
        "threshold_status": PROVISIONAL,
        "late_joint_N80": n80, "late_joint_N90": n90, "late_joint_N95": n95,
        "late_dimensionality": dimensionality,
        "subspace_verdict": subspace,
        "maximum_principal_angle_degrees": maximum_angle,
        "horizontal_scale_verdict": scale_verdict,
        "maximum_MPI_interface_concentration": maximum_concentration,
        "MPI_interface_verdict": mpi_verdict,
        "configured_true_mode_effect": configured_effect,
        "schwarz_actual_mode_effect": configured_effect,
        "weighted_BPI_E_P": weighted_bpi,
        "weighted_configured_E_P": weighted_configured,
        "case_E0": case_results["E0"],
        "case_E1": case_results["E1"],
        "joint_actual_mode_spatial_evidence": spatial_results["JOINT"],
        "joint_actual_mode_MPI_evidence": mpi_results["JOINT"],
        "restart_in_joint_POD_coordinates": restart_metrics,
        "prior_restart_verdict_retained": all(value[
            "prior_true_residual_restart_verdict"] ==
            "RESIDUAL_TAIL_PREEXISTS_RESTART" for value in restart_metrics.values()),
        "maximum_density_gauge_fraction_from_valid_stage_E": density_max,
        "density_gauge_remains_non_dominant": density_max < 1e-3,
        "next_authorized_path": next_path,
        "forensic_path_authorized": next_path == "LOCAL_SCHWARZ_PATH",
        "production_schwarz_modification_authorized": False,
        "automatic_solver_change_authorized": False,
        "thresholds": {
            "status": PROVISIONAL,
            "low_dimensional_N90_maximum": 6,
            "spatial_dominance_fraction": DOMINANCE_FRACTION,
            "material_map_change_fraction": MAP_MATERIAL_FRACTION,
            "MPI_interface_suspect_ratio": MPI_SUSPECT_RATIO,
            "MPI_interface_strong_ratio": MPI_STRONG_RATIO,
        },
    }


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def analyze(args):
    data = load_all(args.input_dir)
    modes = validate(data)
    validation = numerical_validation(data, modes)
    drift0 = validate_trajectory(args.e0_iterations, args.baseline_e0_iterations,
                                 args.trajectory_tolerance)
    drift1 = validate_trajectory(args.e1_iterations, args.baseline_e1_iterations,
                                 args.trajectory_tolerance)
    stage_e = json.loads(Path(args.stage_e_reanalysis).read_text())
    if stage_e.get("comparison", {}).get("schwarz_mode_composition_verdict") != \
            "UNRESOLVED_LOW_REPRESENTATION":
        raise ValueError("Stage-E low-representation reanalysis is missing")
    result = classify(data, modes, stage_e)
    result["trajectory_maximum_relative_drift"] = max(drift0, drift1)
    result["observation_only_trajectory_pass"] = True
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    write_json(output / "strict_ala_stage_E2_numerical_validation.json",
               validation)
    write_json(output / "strict_ala_stage_E2_classification.json", result)
    write_json(output / "strict_ala_stage_E2_next_action.json", {
        "schema": "strict-ala-stage-E2-next-action-v1",
        "next_authorized_path": result["next_authorized_path"],
        "basis": result,
        "automatic_solver_change_authorized": False,
    })
    return 0


def audit(args):
    analysis = Path(args.analysis_dir)
    data = load_all(analysis)
    modes = validate(data)
    recomputed_validation = numerical_validation(data, modes)
    recorded_validation = json.loads((analysis /
        "strict_ala_stage_E2_numerical_validation.json").read_text())
    numerical_pass = recorded_validation == recomputed_validation and \
        recorded_validation.get("status") == "PASS"
    classification = json.loads((analysis / "strict_ala_stage_E2_classification.json").read_text())
    next_action = json.loads((analysis / "strict_ala_stage_E2_next_action.json").read_text())
    hashes_pass = all(Path(a).read_bytes() == Path(b).read_bytes()
                      for a, b in ((args.binary_pre, args.binary_post),
                                   (args.input_pre, args.input_post),
                                   (args.source_pre, args.source_post)))
    e0_log = Path(args.e0_log).read_text(errors="replace")
    e1_log = Path(args.e1_log).read_text(errors="replace")
    e0_marker = e0_log.count("STRICT_ALA_STAGE_E2_ARCHIVED case=E0 snapshots=19") == 1
    e1_marker = e1_log.count("STRICT_ALA_STAGE_E2_COMPLETE case=E1 snapshots=38") == 1
    wrong_case_marker = ("STRICT_ALA_STAGE_E2_COMPLETE case=E1" in e0_log or
                         "STRICT_ALA_STAGE_E2_ARCHIVED case=E0" in e1_log)
    combined_log = (e0_log + "\n" + e1_log).lower()
    logs_clean = not any((
        re.search(r"(?<![a-z])(?:nan|inf)(?![a-z])", combined_log),
        "fatal error" in combined_log,
        "mpi_abort" in combined_log,
        "traceback" in combined_log,
        re.search(r"(?:fallback_blocks|velocity_block_fallbacks)=[1-9][0-9]*",
                  combined_log),
    ))
    logs_pass = e0_marker and e1_marker and not wrong_case_marker and logs_clean
    manifest = json.loads(Path(args.manifest).read_text())
    source = manifest.get("source", {})
    provenance_pass = bool(manifest.get("provenance_complete") and
        source.get("branch") == "cmbhf_ALA_strict" and
        source.get("branch_verified") and not source.get("scientific_dirty"))
    stage_e_audit = json.loads(Path(args.stage_e_audit).read_text())
    stage_e_pass = all(stage_e_audit.get(key) is True for key in (
        "experiment_valid", "E0_complete", "E1_complete",
        "diagnostics_observation_only_pass", "probe_gram_pass",
        "projection_reconciliation_pass", "source_binary_input_identity_pass"))
    decision_pass = bool(
        classification.get("schema") == "strict-ala-stage-E2-classification-v1" and
        classification.get("threshold_status") == PROVISIONAL and
        classification.get("automatic_solver_change_authorized") is False and
        next_action.get("schema") == "strict-ala-stage-E2-next-action-v1" and
        next_action.get("next_authorized_path") ==
            classification.get("next_authorized_path") and
        next_action.get("automatic_solver_change_authorized") is False)
    valid = bool(hashes_pass and logs_pass and provenance_pass and numerical_pass and
                 stage_e_pass and decision_pass and
                 classification.get("observation_only_trajectory_pass") and
                 next_action.get("automatic_solver_change_authorized") is False)
    result = {"schema": "strict-ala-stage-E2-final-audit-v1",
              "experiment_valid": valid, "hashes_unchanged": hashes_pass,
              "case_completion_markers_pass": logs_pass,
              "case_logs_clean": logs_clean,
              "provenance_pass": provenance_pass,
              "valid_stage_E_baseline_pass": stage_e_pass,
              "numerical_validation_pass": numerical_pass,
              "decision_contract_pass": decision_pass,
              "observation_only_trajectory_pass":
                  classification.get("observation_only_trajectory_pass", False),
              "next_authorized_path": classification.get("next_authorized_path"),
              "automatic_solver_change_authorized": False}
    write_json(args.output, result)
    return 0 if valid else 1


def parser():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="command", required=True)
    a = sub.add_parser("analyze")
    a.add_argument("--input-dir", required=True)
    a.add_argument("--output-dir", required=True)
    a.add_argument("--e0-iterations", required=True)
    a.add_argument("--e1-iterations", required=True)
    a.add_argument("--baseline-e0-iterations", required=True)
    a.add_argument("--baseline-e1-iterations", required=True)
    a.add_argument("--stage-e-reanalysis", required=True)
    a.add_argument("--trajectory-tolerance", type=float, default=1e-12)
    a.set_defaults(function=analyze)
    f = sub.add_parser("audit")
    f.add_argument("--analysis-dir", required=True)
    f.add_argument("--manifest", required=True)
    f.add_argument("--binary-pre", required=True); f.add_argument("--binary-post", required=True)
    f.add_argument("--input-pre", required=True); f.add_argument("--input-post", required=True)
    f.add_argument("--source-pre", required=True); f.add_argument("--source-post", required=True)
    f.add_argument("--e0-log", required=True); f.add_argument("--e1-log", required=True)
    f.add_argument("--stage-e-audit", required=True)
    f.add_argument("--output", required=True)
    f.set_defaults(function=audit)
    return p


if __name__ == "__main__":
    arguments = parser().parse_args()
    raise SystemExit(arguments.function(arguments))
