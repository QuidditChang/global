#!/usr/bin/env python3
"""Fail-closed offline analysis for strict-ALA Stage E2 outputs."""

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

SCHEDULE = (0, 5, 10, 20, 30, 35, 40, 42, 44, 46, 48, 49, 50,
            51, 52, 54, 56, 58, 60)
LATE_COUNT = 13
PROVISIONAL = "PROVISIONAL_DIAGNOSTIC"

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
    for case in ("E0", "E1"):
        rows = [r for r in manifest if r["case"] == case]
        if tuple(int(r["iteration"]) for r in rows) != SCHEDULE:
            raise ValueError(f"snapshot schedule mismatch for {case}")
        if any(int(r["ranks"]) != 384 for r in rows):
            raise ValueError("Stage E2 frozen layout is not 384 ranks")
        for row in rows:
            finite(row, ("residual_norm",))
            if len(row["global_checksum"]) != 16:
                raise ValueError("invalid snapshot checksum")

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
        if abs(sum(float(r["energy_fraction"]) for r in rows) - 1.0) > 1e-8:
            raise ValueError("POD energy does not reconcile")
        if abs(float(rows[-1]["cumulative_fraction"]) - 1.0) > 1e-8:
            raise ValueError("POD cumulative energy does not end at one")

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

    coefficients = data["joint_pod_coefficients"]
    if len(coefficients) != 38 * len(modes):
        raise ValueError("joint POD coefficient trajectory incomplete")
    for row in coefficients:
        finite(row, ("coefficient", "contraction"))

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


def classify(data, modes, stage_e):
    joint = [r for r in data["pod"] if r["case"] == "JOINT" and
             r["basis_type"] == "late_direction_case_balanced"]
    n80 = next(i + 1 for i, r in enumerate(joint)
               if float(r["cumulative_fraction"]) >= .8)
    n90 = next(i + 1 for i, r in enumerate(joint)
               if float(r["cumulative_fraction"]) >= .9)
    n95 = next(i + 1 for i, r in enumerate(joint)
               if float(r["cumulative_fraction"]) >= .95)
    dimensionality = "LOW_DIMENSIONAL" if n90 <= 6 else "BROAD_SUBSPACE"
    angles = data["subspace_angles"]
    maximum_angle = max(float(r["angle_degrees"]) for r in angles)
    subspace = ("SAME_SUBSPACE" if maximum_angle <= 10 else
                "MOSTLY_SAME" if maximum_angle <= 30 else
                "MATERIAL_SUBSPACE_CHANGE")

    scale_energy = defaultdict(float)
    for row in data["depth_scale_spectrum"]:
        if row["case"] == "JOINT" and row["object_type"] == "pod_mode":
            scale_energy[row["horizontal_scale_band"]] += float(row["energy"])
    scale_total = sum(scale_energy.values())
    dominant_scale, dominant_energy = max(scale_energy.items(), key=lambda item: item[1])
    scale_verdict = dominant_scale if dominant_energy / max(scale_total, 1e-300) >= .6 \
        else "BROAD_HORIZONTAL_SPECTRUM"

    mpi_rows = [r for r in data["mpi_distance"]
                if r["case"] == "JOINT" and r["object_type"] == "pod_mode"]
    by_mode = defaultdict(dict)
    for row in mpi_rows:
        by_mode[int(row["object_id"])][row["mpi_distance_band"]] = row
    concentrations = []
    for rows in by_mode.values():
        concentrations.append((float(rows["0"]["energy_per_element"]) /
            max(float(rows["3+"]["energy_per_element"]), 1e-300)))
    maximum_concentration = max(concentrations)
    mpi_verdict = ("MPI_INTERFACE_STRONGLY_SUSPECT" if maximum_concentration >= 2 else
                   "MPI_INTERFACE_SUSPECT" if maximum_concentration >= 1.5 else
                   "MPI_INTERFACE_NOT_SUPPORTED")

    maps = defaultdict(dict)
    for row in data["true_mode_preconditioner"]:
        maps[int(row["mode"])][row["map"]] = float(row["E_P"])
    weighted_bpi = sum(maps[m]["BPI"] * float(joint[m - 1]["energy_fraction"])
                       for m in modes)
    weighted_configured = sum(maps[m]["CONFIGURED"] *
                              float(joint[m - 1]["energy_fraction"]) for m in modes)
    configured_effect = ("IMPROVES" if weighted_configured < .95 * weighted_bpi else
                         "DEGRADES" if weighted_configured > 1.05 * weighted_bpi else
                         "NO_MATERIAL_CHANGE")

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

    if dimensionality == "BROAD_SUBSPACE":
        next_path = "MULTILEVEL_PRESSURE_PATH"
    elif mpi_verdict != "MPI_INTERFACE_NOT_SUPPORTED":
        next_path = "MPI_OVERLAP_PATH"
    elif scale_verdict == "GLOBAL_VERY_LONG":
        next_path = "GLOBAL_COARSE_PATH"
    elif scale_verdict == "LONG_INTERMEDIATE":
        next_path = "TARGETED_LOW_RANK_PATH"
    elif scale_verdict in ("PATCH_LOCAL", "HIGH_FREQUENCY"):
        next_path = "LOCAL_SCHWARZ_PATH"
    else:
        next_path = "UNRESOLVED"
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
        "weighted_BPI_E_P": weighted_bpi,
        "weighted_configured_E_P": weighted_configured,
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
    }


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def analyze(args):
    data = load_all(args.input_dir)
    modes = validate(data)
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
    classification = json.loads((analysis / "strict_ala_stage_E2_classification.json").read_text())
    next_action = json.loads((analysis / "strict_ala_stage_E2_next_action.json").read_text())
    hashes_pass = all(Path(a).read_bytes() == Path(b).read_bytes()
                      for a, b in ((args.binary_pre, args.binary_post),
                                   (args.input_pre, args.input_post),
                                   (args.source_pre, args.source_post)))
    logs_pass = all("STRICT_ALA_STAGE_E2_COMPLETE case=E1" in Path(path).read_text(
        errors="replace") or "STRICT_ALA_STAGE_E2_ARCHIVED case=E0" in Path(path).read_text(
        errors="replace") for path in args.case_log)
    manifest = json.loads(Path(args.manifest).read_text())
    source = manifest.get("source", {})
    provenance_pass = bool(manifest.get("provenance_complete") and
        source.get("branch") == "cmbhf_ALA_strict" and
        source.get("branch_verified") and not source.get("scientific_dirty"))
    valid = bool(hashes_pass and logs_pass and provenance_pass and
                 classification.get("observation_only_trajectory_pass") and
                 next_action.get("automatic_solver_change_authorized") is False)
    result = {"schema": "strict-ala-stage-E2-final-audit-v1",
              "experiment_valid": valid, "hashes_unchanged": hashes_pass,
              "case_completion_markers_pass": logs_pass,
              "provenance_pass": provenance_pass,
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
    f.add_argument("--case-log", action="append", required=True)
    f.add_argument("--output", required=True)
    f.set_defaults(function=audit)
    return p


if __name__ == "__main__":
    arguments = parser().parse_args()
    raise SystemExit(arguments.function(arguments))
