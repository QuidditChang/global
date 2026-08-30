#!/usr/bin/env python3
"""Fail-closed Stage-F1b operator-consistent Schwarz audit."""

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path

SCALE_FLOOR = 1.0e-300
MAX_TIGHT_SOLVES = 24
K_REPLAY = 1.0e-10
K_REPLAY_REVIEW = 1.0e-8
SYMMETRY = 1.0e-10
LOCAL_SOLVE = 1.0e-10
OUTPUT_REPLICA = 1.0e-10
MOMENTUM_GATE = 1.0e-3
MOMENTUM_PRECISION = 1.0e-12

ASSEMBLY_COLUMNS = (
    "row_type", "patch_category", "rank", "patch_id", "pressure_dim",
    "velocity_dim", "symmetry_defect", "regularization",
    "factorization_status", "pivot_metric", "local_solve_residual",
    "fallback", "candidate_matrix_bytes", "candidate_factor_bytes",
    "Kgamma_replay_relative_defect",
    "production_output_replica_relative_difference",
    "B_restriction_relative_defect", "Bt_restriction_relative_defect",
    "B_Bt_bilinear_relative_defect")
SHORT_COLUMNS = (
    "case", "fgmres_candidate_iteration", "global_cancellation_L2",
    "unaugmented_momentum_relative", "recursive_krylov_residual",
    "explicit_krylov_residual", "cumulative_MG_cycles",
    "cumulative_Kgamma_applications", "elapsed_solver_seconds",
    "elapsed_wall_seconds", "terminal", "terminal_reason")
COST_COLUMNS = (
    "case", "terminal_iteration", "total_MG_cycles",
    "total_Kgamma_applications", "solver_wall_seconds", "case_wall_seconds",
    "rank_rss_max_kib", "rank_rss_sum_kib", "ranks",
    "retained_cache_bytes_per_rank", "temporary_peak_bytes_per_rank")


def finite(value, name):
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"non-finite {name}")
    return number


def load_json(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing or empty {path}")
    return json.loads(path.read_text())


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path, expected=None):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing or empty {path}")
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if expected is not None and tuple(reader.fieldnames or ()) != tuple(expected):
            raise ValueError(f"schema mismatch in {path.name}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"no rows in {path.name}")
    return rows


def cfg_values(path):
    values = {}
    for line in Path(path).read_text().splitlines():
        match = re.match(r"^\s*([A-Za-z0-9_]+)\s*(?:=\s*|\s+)(\S.*?)\s*$", line)
        if match and not match.group(1).startswith("#"):
            values[match.group(1)] = match.group(2)
    return values


def cfg_entries(path):
    entries = {}
    section = None
    for number, line in enumerate(Path(path).read_text().splitlines(), 1):
        header = re.match(r"^\s*\[([^]]+)\]\s*$", line)
        if header:
            section = header.group(1)
            continue
        match = re.match(r"^\s*([A-Za-z0-9_]+)\s*(?:=\s*|\s+)(\S.*?)\s*$", line)
        if match and not match.group(1).startswith("#"):
            if match.group(1) in entries:
                raise ValueError(f"duplicate cfg key {match.group(1)}")
            entries[match.group(1)] = (match.group(2), section, number)
    return entries


def set_cfg_value(args):
    path = Path(args.cfg)
    lines = path.read_text().splitlines(keepends=True)
    section = None
    section_start = section_end = None
    matches = []
    key_pattern = re.compile(r"^\s*" + re.escape(args.key) + r"\s*(?:=|\s)")
    for index, line in enumerate(lines):
        header = re.match(r"^\s*\[([^]]+)\]\s*$", line)
        if header:
            if section == args.section and section_end is None:
                section_end = index
            section = header.group(1)
            if section == args.section:
                if section_start is not None:
                    raise ValueError(f"duplicate cfg section {args.section}")
                section_start = index
            continue
        if key_pattern.match(line) and not line.lstrip().startswith("#"):
            matches.append((index, section))
    if section == args.section and section_end is None:
        section_end = len(lines)
    if section_start is None:
        raise ValueError(f"missing cfg section {args.section}")
    if len(matches) > 1:
        raise ValueError(f"duplicate cfg key {args.key}")
    replacement = f"{args.key} = {args.value}\n"
    if matches:
        index, owner = matches[0]
        if owner != args.section:
            raise ValueError(f"cfg key {args.key} belongs to [{owner}], not [{args.section}]")
        lines[index] = replacement
    else:
        lines.insert(section_end, replacement)
    path.write_text("".join(lines))
    return 0


def check_cfg(args):
    base_entries, cand_entries = cfg_entries(args.base), cfg_entries(args.candidate)
    base = {key: entry[0] for key, entry in base_entries.items()}
    cand = {key: entry[0] for key, entry in cand_entries.items()}
    differences = [(key, base.get(key), cand.get(key))
                   for key in sorted(set(base) | set(cand))
                   if base.get(key) != cand.get(key)]
    Path(args.diff).write_text("".join(
        f"{key}: {old} -> {new}\n" for key, old, new in differences))
    expected = [("ala_shallow_patch_local_operator", "legacy",
                 "operator_consistent")]
    if differences != expected:
        raise ValueError(f"F1b cfg scientific diff is not closed: {differences!r}")
    required = {
        "ala_shallow_patch_preconditioner": "on",
        "ala_shallow_patch_velocity_solver": "element_vanka",
        "ala_augmented_lagrangian_gamma": "10.0",
        "ala_outer_solver": "fgmres",
        "ala_pcg_restart_interval": "50",
        "ala_stage_abc_production_logging": "on",
        "ala_stage_e_diagnostic": "on",
        "piterations": "30", "nprocx": "4", "nprocy": "4", "nprocz": "2",
        "nodex": "129", "nodey": "129", "nodez": "65", "steps": "1"}
    required_sections = {
        **{key: "CitcomS.solver.vsolver" for key in (
            "ala_shallow_patch_preconditioner",
            "ala_shallow_patch_velocity_solver",
            "ala_augmented_lagrangian_gamma", "ala_outer_solver",
            "ala_pcg_restart_interval", "ala_stage_abc_production_logging",
            "ala_stage_e_diagnostic", "piterations")},
        **{key: "CitcomS.solver.mesher" for key in (
            "nprocx", "nprocy", "nprocz", "nodex", "nodey", "nodez")},
        "steps": "CitcomS"}
    for key, value in required.items():
        if base.get(key) != value or cand.get(key) != value:
            raise ValueError(f"F1b cfg requires {key}={value}")
        expected_section = required_sections[key]
        if (base_entries[key][1] != expected_section or
                cand_entries[key][1] != expected_section):
            raise ValueError(f"F1b cfg key {key} is outside [{expected_section}]")
    for entries in (base_entries, cand_entries):
        if entries["ala_shallow_patch_local_operator"][1] != \
                "CitcomS.solver.vsolver":
            raise ValueError("F1b selector is outside [CitcomS.solver.vsolver]")
    return 0


def thresholds(layout=None):
    layout = dict(layout or {})
    node_memory_kib = int(layout.pop("allocated_node_memory_kib", 0))
    wall_limit_seconds = int(layout.pop("allocated_wall_limit_seconds", 0))
    startup_budget_seconds = int(layout.pop("startup_budget_seconds", 0))
    return {
        "schema": "strict-ala-stage-F1b-thresholds-v1",
        "threshold_status": "FROZEN_PRELAUNCH",
        "MAX_F1B_TIGHT_SGAMMA_SOLVES": MAX_TIGHT_SOLVES,
        "Kgamma_replay_relative_tolerance": K_REPLAY,
        "Kgamma_replay_stop_review_upper_bound": K_REPLAY_REVIEW,
        "candidate_symmetry_tolerance": SYMMETRY,
        "local_solve_relative_tolerance": LOCAL_SOLVE,
        "production_output_replica_relative_tolerance": OUTPUT_REPLICA,
        "projected_peak_memory_ratio_limit": 1.25,
        "local_action_RMS_ratio_limit": 0.70,
        "true_mode_EP_RMS_ratio_limit": 0.80,
        "H_RAW_ratio_limit": 0.50,
        "short_route_A_residual_ratio": 0.85,
        "short_route_B_residual_ratio": 0.90,
        "short_cost_ratio_limit": 1.25,
        "momentum_gate": MOMENTUM_GATE,
        "momentum_ratio_limit": 1.10,
        "momentum_precision_absolute": MOMENTUM_PRECISION,
        "weighted_metric": "sqrt(sum(weight*error^2)/sum(weight))",
        "operator_action": "Q_i solve(A_raw + epsilon diag(A_raw), Q_i rhs)",
        "mpi_layout": layout,
        "allocated_node_memory_kib": node_memory_kib,
        "node_peak_fraction_limit": 0.80,
        "allocated_wall_limit_seconds": wall_limit_seconds,
        "declared_project_startup_budget_seconds": startup_budget_seconds,
        "startup_wall_fraction_limit": 0.25,
        "production_freeze_contract": {
            "F1c_is_final_production_qualification": False,
            "production_representative_smoke_soak_required": True,
            "minimum_scope": "one complete coupled Stokes/thermal timestep",
            "independent_repeat_required": True,
        },
    }


def schema_manifest():
    return {"schema": "strict-ala-stage-F1b-schema-v1",
            "collective_manifest": [
                {"audit_sequence_id": 0, "patch_category": "LOCAL",
                 "vector_id": "canonical_trigonometric_v1",
                 "collective_buffer_count": 1, "collectives_per_velocity_dof": 2},
                {"audit_sequence_id": 1, "patch_category": "INTERFACE",
                 "vector_id": "canonical_trigonometric_v1",
                 "collective_buffer_count": 1, "collectives_per_velocity_dof": 2}],
            "files": [
        {"filename": "strict_ala_stage_F1b_candidate_assembly.csv",
         "columns": list(ASSEMBLY_COLUMNS),
         "unique_key": ["row_type", "patch_category", "rank", "patch_id"]},
        {"filename": "strict_ala_stage_F1b_local_action.csv",
         "unique_key": ["patch_or_bin_ID", "POD_mode", "operator"]},
        {"filename": "strict_ala_stage_F1b_true_mode.csv",
         "unique_key": ["POD_mode", "operator"]},
        {"filename": "strict_ala_stage_F1b_telescoping.csv",
         "unique_key": ["POD_mode", "stage", "operator"]},
        {"filename": "strict_ala_stage_F1b_short_iterations.csv",
         "columns": list(SHORT_COLUMNS),
         "unique_key": ["case", "fgmres_candidate_iteration"]},
        {"filename": "strict_ala_stage_F1b_short_cost.csv",
         "columns": list(COST_COLUMNS), "unique_key": ["case"]},
        {"filename": "strict_ala_stage_F1b_mode_manifest.csv",
         "unique_key": ["POD_mode"]},
    ]}


def freeze_contract(args):
    output = Path(args.output_dir); output.mkdir(parents=True, exist_ok=True)
    layout = {
        "allocated_tasks": int(args.allocated_tasks),
        "solver_world_size": int(args.solver_world_size),
        "nprocx": int(args.nprocx), "nprocy": int(args.nprocy),
        "nprocz": int(args.nprocz), "caps": int(args.caps),
        "ranks_per_node": int(args.ranks_per_node),
        "allocated_node_memory_kib": int(args.node_memory_kib),
        "allocated_wall_limit_seconds": int(args.allocated_wall_seconds),
        "startup_budget_seconds": int(args.startup_budget_seconds),
    }
    if min(layout.values()) <= 0:
        raise ValueError("F1b MPI layout values must be positive")
    expected = layout["nprocx"] * layout["nprocy"] * layout["nprocz"] * layout["caps"]
    if expected != layout["solver_world_size"]:
        raise ValueError("F1b solver world does not match topology times caps")
    if layout["allocated_tasks"] < layout["solver_world_size"]:
        raise ValueError("F1b allocation is smaller than solver world")
    write_json(output / "strict_ala_stage_F1b_thresholds.json", thresholds(layout))
    write_json(output / "strict_ala_stage_F1b_schema.json", schema_manifest())
    return 0


def verify_layout(args):
    frozen = load_json(args.thresholds)["mpi_layout"]
    entries = cfg_entries(args.cfg)
    observed = {key: int(entries[key][0]) for key in ("nprocx", "nprocy", "nprocz")}
    observed.update({
        "caps": int(args.caps),
        "solver_world_size": int(args.solver_world_size),
        "allocated_tasks": int(args.allocated_tasks),
        "ranks_per_node": int(args.ranks_per_node),
    })
    expected_world = observed["nprocx"] * observed["nprocy"] * observed["nprocz"] * observed["caps"]
    checks = {
        "frozen_layout_matches_runtime": observed == frozen,
        "cfg_topology_matches_solver_world": expected_world == observed["solver_world_size"],
        "allocation_covers_solver_world": observed["allocated_tasks"] >= observed["solver_world_size"],
    }
    result = {"schema": "strict-ala-stage-F1b-layout-v1", "valid": all(checks.values()),
              "checks": checks, "layout": observed,
              "unused_allocation_tasks": observed["allocated_tasks"] - observed["solver_world_size"]}
    write_json(args.output, result)
    if not result["valid"]:
        raise ValueError(f"F1b MPI layout contract failed: {checks}")
    return 0


def verify_lineage(args):
    final = load_json(args.final_audit)
    checks = {
        "experiment_valid": final.get("experiment_valid") is True,
        "primary_defect": final.get("primary_defect_class") == "LOCAL_OPERATOR_MISMATCH",
        "dominant_stage": final.get("dominant_damage_stage") == "RAW_LOCAL",
        "F1b_authorized": final.get("F1b_authorized") is True,
        "next_task": final.get("next_authorized_task") == "F1B_LOCAL_CANDIDATE",
        "evidence_root": final.get("numerical_evidence_root") == args.expected_hpc_root,
        "projected_matrices": Path(args.evidence_root,
            "01_operator_replay/strict_ala_stage_F1a_projected_matrices.csv").is_file(),
    }
    if not all(checks.values()):
        raise ValueError(f"F1a dual-root lineage failed: {checks}")
    hraw = finite(final["weighted_telescoping_damage"]["RAW_LOCAL"], "H_RAW")
    if hraw <= 0:
        raise ValueError("F1a H_RAW is not positive")
    result = {"schema": "strict-ala-stage-F1b-lineage-v1", "valid": True,
              "checks": checks, "H_RAW_legacy": hraw,
              "postprocess_final_sha256": sha256(args.final_audit)}
    write_json(args.output, result)
    return 0


def weighted_rms(rows, error_field, weight_field):
    weighted = total = 0.0
    for row in rows:
        error = finite(row[error_field], error_field)
        weight = finite(row[weight_field], weight_field)
        if weight < 0:
            raise ValueError("negative weight")
        weighted += weight * error * error
        total += weight
    if total <= 0:
        raise ValueError("zero total weight")
    return math.sqrt(weighted / total)


def uniform_scalar_rescaling_bound(rows, legacy_rms, ratio_limit):
    """Return the best possible true-mode RMS from one scalar rescaling.

    The frozen mode normalization is ||q||=1.  With a=S_gamma P q,
    E_P^2=1+||a||^2-2 q^T a reconstructs the action norm without another
    tight solve.  This is a diagnostic impossibility bound only; F1b does not
    authorize changing W or any other scaling.
    """
    total_weight = sum(finite(row["POD_energy_weight"], "POD weight")
                       for row in rows)
    if total_weight <= 0.0:
        raise ValueError("nonpositive POD weight sum")
    weighted_dot = 0.0
    weighted_norm2 = 0.0
    modes = []
    for row in rows:
        weight = finite(row["POD_energy_weight"], "POD weight")
        error = finite(row["E_P"], "E_P")
        dot = finite(row["qTSPq"], "qTSPq")
        cosine = finite(row["cosine"], "cosine")
        action_norm2 = error * error - 1.0 + 2.0 * dot
        roundoff = 1.0e-12 * max(1.0, error * error, abs(2.0 * dot))
        if action_norm2 < -roundoff:
            raise ValueError("true-mode metrics imply negative action norm")
        action_norm2 = max(action_norm2, 0.0)
        action_norm = math.sqrt(action_norm2)
        if action_norm > 0.0:
            reconstructed_cosine = dot / action_norm
            if not math.isclose(reconstructed_cosine, cosine,
                                rel_tol=1.0e-8, abs_tol=1.0e-10):
                raise ValueError("true-mode E_P/cosine/qTSPq are inconsistent")
            individual_floor = math.sqrt(max(
                0.0, 1.0 - dot * dot / action_norm2))
            individual_alpha = max(dot / action_norm2, 0.0)
        else:
            individual_floor = 1.0
            individual_alpha = 0.0
        weighted_dot += weight * dot
        weighted_norm2 += weight * action_norm2
        modes.append({
            "POD_mode": row["POD_mode"],
            "POD_energy_weight": weight,
            "action_norm": action_norm,
            "qTSPq": dot,
            "individual_optimal_nonnegative_scalar": individual_alpha,
            "individual_residual_lower_bound": individual_floor,
        })
    alpha = max(weighted_dot / weighted_norm2, 0.0) \
        if weighted_norm2 > 0.0 else 0.0
    minimum_squared = sum(
        mode["POD_energy_weight"] * (
            1.0 + alpha * alpha * mode["action_norm"] ** 2
            - 2.0 * alpha * mode["qTSPq"])
        for mode in modes) / total_weight
    minimum_rms = math.sqrt(max(minimum_squared, 0.0))
    minimum_ratio = minimum_rms / max(legacy_rms, SCALE_FLOOR)
    return {
        "diagnostic_only": True,
        "changes_frozen_scaling": False,
        "optimal_common_nonnegative_scalar": alpha,
        "minimum_candidate_E_P_RMS": minimum_rms,
        "minimum_candidate_to_legacy_E_P_RMS_ratio": minimum_ratio,
        "frozen_ratio_limit": ratio_limit,
        "can_reach_frozen_gate": minimum_ratio <= ratio_limit,
        "modes": modes,
    }


def keyed(rows, fields):
    result = {}
    for row in rows:
        key = tuple(row[field] for field in fields)
        if key in result:
            raise ValueError(f"duplicate key {key}")
        result[key] = row
    return result


def write_tagged_csv(path, groups, tag="operator"):
    first = next(iter(groups.values()))
    fields = list(first[0].keys()) + [tag]
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for label, rows in groups.items():
            for row in rows:
                tagged = dict(row); tagged[tag] = label
                writer.writerow(tagged)


def identity(args):
    """Gate distributed restriction evidence before scientific candidate use."""
    feasibility = load_json(args.feasibility)
    thresholds = load_json(args.thresholds)
    if feasibility.get("collective_manifest_sha256") != sha256(args.schema):
        raise ValueError("F1b collective manifest hash mismatch")
    if feasibility.get("factorization_phase") != "diagnostic_factor_attempt":
        raise ValueError("identity preflight used a non-diagnostic factor phase")
    tolerance = finite(thresholds["Kgamma_replay_relative_tolerance"],
                       "frozen restriction tolerance")
    review = finite(thresholds["Kgamma_replay_stop_review_upper_bound"],
                    "frozen restriction review bound")
    if tolerance != K_REPLAY or review != K_REPLAY_REVIEW:
        raise ValueError("F1b identity thresholds differ from source contract")
    rows = read_csv(args.assembly, ASSEMBLY_COLUMNS)
    if {row["patch_category"] for row in rows} != {"LOCAL", "INTERFACE"}:
        raise ValueError("identity evidence lacks local/interface representatives")
    identity_fields = ("Kgamma_replay_relative_defect",
                       "B_restriction_relative_defect",
                       "Bt_restriction_relative_defect",
                       "B_Bt_bilinear_relative_defect")
    values = [finite(row[field], field) for row in rows
              for field in identity_fields]
    values.extend(finite(feasibility[field], f"feasibility {field}")
                  for field in identity_fields)
    replicas = [finite(row["production_output_replica_relative_difference"],
                       "output replica difference") for row in rows]
    replicas.append(finite(
        feasibility["production_output_replica_relative_difference"],
        "feasibility output replica difference"))
    maximum = max(values)
    maximum_replica = max(replicas)
    passed = maximum <= tolerance and maximum_replica <= OUTPUT_REPLICA
    review_required = tolerance < maximum <= review
    gates = {"lineage": True, "candidate_isolatable": passed,
             "memory": True, "startup": True, "assembly": passed,
             "tight_solves": False, "local_RMS": False,
             "local_heavy_mode": False, "true_mode_RMS": False,
             "true_mode_heavy": False, "dominant_mode_positivity": False,
             "H_RAW": False}
    result = {
        "schema": "strict-ala-stage-F1b-identity-gate-v1",
        "experiment_evidence_valid": True,
        "identity_gate_pass": passed,
        "offline_gate_pass": False,
        "gates": gates,
        "restriction_tolerance_review_required": review_required,
        "maximum_restriction_identity_relative_defect": maximum,
        "maximum_production_output_replica_relative_difference": maximum_replica,
        "factorization_phase": "diagnostic_factor_attempt",
        "authorized_candidate_factorization": False,
        "thresholds_sha256": sha256(args.thresholds),
        "schema_sha256": sha256(args.schema)}
    write_json(args.output, result)
    return 0


def offline(args):
    output = Path(args.output_dir); output.mkdir(parents=True, exist_ok=True)
    lineage = load_json(args.lineage)
    identity_result = load_json(args.identity_gate)
    feasibility = load_json(args.feasibility)
    frozen_thresholds = load_json(args.thresholds)
    if feasibility.get("collective_manifest_sha256") != sha256(args.schema):
        raise ValueError("F1b collective manifest hash mismatch")
    if feasibility.get("factorization_phase") != "authorized_candidate_factorization":
        raise ValueError("offline candidate was not authorized after identity gate")
    if identity_result.get("identity_gate_pass") is not True:
        raise ValueError("authorized candidate lacks a passing identity gate")
    replay_tolerance = finite(
        frozen_thresholds["Kgamma_replay_relative_tolerance"],
        "frozen Kgamma replay tolerance")
    replay_review = finite(
        frozen_thresholds["Kgamma_replay_stop_review_upper_bound"],
        "frozen Kgamma replay review bound")
    if replay_tolerance != K_REPLAY or replay_review != K_REPLAY_REVIEW:
        raise ValueError("F1b runtime replay thresholds differ from source contract")
    base_memory = memory_row(args.base_memory, "BASE_MEMORY_PREFLIGHT")
    cand_memory = memory_row(args.cand_memory, "CAND_MEMORY_PREFLIGHT")
    if (int(base_memory["ranks"]) != frozen_thresholds["mpi_layout"]["solver_world_size"] or
            int(cand_memory["ranks"]) != frozen_thresholds["mpi_layout"]["solver_world_size"]):
        raise ValueError("F1b memory preflight MPI world mismatch")
    for field in ("factorization_count", "RHS_solve_count",
                  "estimated_dense_factor_and_solve_work",
                  "temporary_peak_bytes_per_rank",
                  "retained_cache_bytes_per_rank",
                  "legacy_retained_cache_bytes_per_rank",
                  "retained_cache_ratio_vs_legacy",
                  "projected_peak_ratio_vs_legacy_cache",
                  "new_mpi_payload_bytes_max_per_rank",
                  "diagonal_completion_relative_max",
                  "Kgamma_replay_relative_defect",
                  "production_output_replica_relative_difference",
                  "B_restriction_relative_defect",
                  "Bt_restriction_relative_defect",
                  "B_Bt_bilinear_relative_defect",
                  "startup_seconds_max"):
        if finite(feasibility[field], field) < 0:
            raise ValueError(f"negative feasibility field {field}")
    m1 = (cand_memory["retained_cache_bytes_per_rank"] <=
          1.25 * base_memory["retained_cache_bytes_per_rank"])
    m2 = (cand_memory["rank_rss_max_kib"] <=
          1.25 * base_memory["rank_rss_max_kib"])
    projected_node_peak = (cand_memory["rank_rss_max_kib"] *
                           frozen_thresholds["mpi_layout"]["ranks_per_node"])
    m3 = projected_node_peak <= (frozen_thresholds["node_peak_fraction_limit"] *
                                 frozen_thresholds["allocated_node_memory_kib"])
    startup_limit = min(
        frozen_thresholds["startup_wall_fraction_limit"] *
        frozen_thresholds["allocated_wall_limit_seconds"],
        frozen_thresholds["declared_project_startup_budget_seconds"])
    startup_pass = finite(feasibility["startup_seconds_max"],
                          "startup seconds") <= startup_limit
    if (int(feasibility["factorization_count"]) <= 0 or
            int(feasibility["RHS_solve_count"]) <= 0):
        raise ValueError("candidate did not factor or solve any real patch")
    assembly = read_csv(args.assembly, ASSEMBLY_COLUMNS)
    if {row["patch_category"] for row in assembly} != {"LOCAL", "INTERFACE"}:
        raise ValueError("candidate assembly lacks local/interface aggregates")
    for row in assembly:
        for field in ("pressure_dim", "velocity_dim", "symmetry_defect",
                      "regularization", "pivot_metric", "local_solve_residual",
                      "candidate_matrix_bytes", "candidate_factor_bytes",
                      "Kgamma_replay_relative_defect",
                      "production_output_replica_relative_difference",
                      "B_restriction_relative_defect",
                      "Bt_restriction_relative_defect",
                      "B_Bt_bilinear_relative_defect"):
            finite(row[field], field)
        if int(row["pressure_dim"]) <= 0 or int(row["velocity_dim"]) <= 0:
            raise ValueError("candidate assembly has an empty aggregate")
        if (row["row_type"] != "AGGREGATE" or
                row["factorization_status"] != "PASS" or row["fallback"] != "0"):
            raise ValueError("candidate assembly status failed")
    identity_fields = ("Kgamma_replay_relative_defect",
                       "B_restriction_relative_defect",
                       "Bt_restriction_relative_defect",
                       "B_Bt_bilinear_relative_defect")
    maximum_replay = max(
        *(finite(r[field], field) for r in assembly for field in identity_fields),
        *(finite(feasibility[field], f"feasibility {field}")
          for field in identity_fields))
    replay_review_required = replay_tolerance < maximum_replay <= replay_review
    assembly_pass = (max(finite(r["symmetry_defect"], "symmetry") for r in assembly)
                     <= SYMMETRY and
                     max(finite(r["local_solve_residual"], "solve") for r in assembly)
                     <= LOCAL_SOLVE and
                     max(finite(r["Kgamma_replay_relative_defect"], "replay")
                         for r in assembly) <= replay_tolerance)
    assembly_pass = (assembly_pass and
                     all(finite(r[field], field) <= replay_tolerance
                         for r in assembly for field in identity_fields) and
                     all(finite(feasibility[field], f"feasibility {field}")
                         <= replay_tolerance for field in identity_fields) and
                     max(finite(r["production_output_replica_relative_difference"],
                                "output replica difference") for r in assembly)
                     <= OUTPUT_REPLICA and
                     finite(feasibility["production_output_replica_relative_difference"],
                            "feasibility output replica difference")
                     <= OUTPUT_REPLICA)

    legacy_local = read_csv(args.legacy_local)
    cand_local = read_csv(args.candidate_local)
    write_tagged_csv(output / "strict_ala_stage_F1b_local_action.csv",
                     {"legacy": legacy_local, "operator_consistent": cand_local})
    lk = keyed(legacy_local, ("patch_or_bin_ID", "POD_mode"))
    ck = keyed(cand_local, ("patch_or_bin_ID", "POD_mode"))
    if set(lk) != set(ck):
        raise ValueError("candidate local representative set changed")
    legacy_rms = weighted_rms(legacy_local, "local_action_error", "contribution_weight")
    cand_rms = weighted_rms(cand_local, "local_action_error", "contribution_weight")
    local_heavy_guard = True
    total_weight = sum(finite(r["contribution_weight"], "weight") for r in legacy_local)
    for key in lk:
        weight = finite(lk[key]["contribution_weight"], "weight")
        if not math.isclose(weight, finite(ck[key]["contribution_weight"], "weight"),
                            rel_tol=1e-12, abs_tol=1e-300):
            raise ValueError("candidate changed local contribution weights")
        if weight >= 0.05 * total_weight:
            old = finite(lk[key]["local_action_error"], "legacy local error")
            new = finite(ck[key]["local_action_error"], "candidate local error")
            local_heavy_guard &= new <= 1.10 * old

    legacy_tel = read_csv(args.legacy_telescoping)
    cand_tel = read_csv(args.candidate_telescoping)
    write_tagged_csv(output / "strict_ala_stage_F1b_telescoping.csv",
                     {"legacy": legacy_tel, "operator_consistent": cand_tel})
    legacy_cfg = [r for r in legacy_tel if r["stage"] == "CONFIGURED"]
    cand_cfg = [r for r in cand_tel if r["stage"] == "CONFIGURED"]
    lcfg = keyed(legacy_cfg, ("POD_mode",)); ccfg = keyed(cand_cfg, ("POD_mode",))
    if set(lcfg) != set(ccfg):
        raise ValueError("candidate changed selected POD modes")
    true_fields = ("POD_mode", "POD_energy_weight", "E_P", "cosine",
                   "qTPq", "qTSPq", "tight_solve_achieved")
    with (output / "strict_ala_stage_F1b_true_mode.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=true_fields + ("operator",))
        writer.writeheader()
        for label, rows in (("legacy", legacy_cfg),
                            ("operator_consistent", cand_cfg)):
            for row in rows:
                writer.writerow({**{field: row[field] for field in true_fields},
                                 "operator": label})
    ep_legacy = weighted_rms(legacy_cfg, "E_P", "POD_energy_weight")
    ep_cand = weighted_rms(cand_cfg, "E_P", "POD_energy_weight")
    scalar_bound = uniform_scalar_rescaling_bound(
        cand_cfg, ep_legacy,
        finite(frozen_thresholds["true_mode_EP_RMS_ratio_limit"],
               "true-mode RMS ratio limit"))
    ep_heavy_guard = True; positive_weight = 0.0
    total_pod = sum(finite(row["POD_energy_weight"], "POD weight")
                    for row in legacy_cfg)
    for key in lcfg:
        weight = finite(lcfg[key]["POD_energy_weight"], "POD weight")
        if not math.isclose(weight, finite(ccfg[key]["POD_energy_weight"], "POD weight"),
                            rel_tol=1e-12, abs_tol=1e-300):
            raise ValueError("candidate changed POD weights")
        if finite(ccfg[key]["qTSPq"], "qTSPq") > 0:
            positive_weight += weight
        if weight >= 0.05 * total_pod:
            ep_heavy_guard &= (finite(ccfg[key]["E_P"], "candidate E_P") <=
                               1.10 * finite(lcfg[key]["E_P"], "legacy E_P"))

    legacy_decision = load_json(args.legacy_decision)
    cand_decision = load_json(args.candidate_decision)
    hraw_legacy = finite(lineage["H_RAW_legacy"], "H_RAW_legacy")
    if not math.isclose(hraw_legacy,
        finite(legacy_decision["weighted_telescoping_damage"]["RAW_LOCAL"], "legacy H_RAW"),
        rel_tol=1e-13):
        raise ValueError("legacy H_RAW lineage changed")
    hraw_cand = finite(cand_decision["weighted_telescoping_damage"]["RAW_LOCAL"],
                       "candidate H_RAW")
    tight = read_csv(args.candidate_tight)
    tight_pass = (len(tight) <= MAX_TIGHT_SOLVES and
                  all(row.get("converged") == "1" for row in tight))
    snapshot_hash = sha256(args.snapshot_manifest)
    mode_manifest = output / "strict_ala_stage_F1b_mode_manifest.csv"
    with mode_manifest.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("POD_mode", "POD_energy_weight", "mode_norm_contract",
                         "snapshot_manifest_sha256", "lineage"))
        for key in sorted(lcfg, key=lambda value: int(value[0])):
            writer.writerow((key[0], lcfg[key]["POD_energy_weight"], "unit_global_pdot",
                             snapshot_hash, "frozen_F1a_selected_joint_POD"))

    metrics = {
        "H_RAW_legacy": hraw_legacy, "H_RAW_candidate": hraw_cand,
        "H_RAW_ratio": hraw_cand / hraw_legacy,
        "E_local_RMS_legacy": legacy_rms, "E_local_RMS_candidate": cand_rms,
        "E_local_RMS_ratio": cand_rms / max(legacy_rms, SCALE_FLOOR),
        "E_P_RMS_legacy": ep_legacy, "E_P_RMS_candidate": ep_cand,
        "E_P_RMS_ratio": ep_cand / max(ep_legacy, SCALE_FLOOR)}
    gates = {
        "lineage": lineage.get("valid") is True,
        "candidate_isolatable": feasibility.get("candidate_isolatable") is True,
        "memory": m1 and m2 and m3, "startup": startup_pass,
        "assembly": assembly_pass, "tight_solves": tight_pass,
        "local_RMS": metrics["E_local_RMS_ratio"] <= 0.70,
        "local_heavy_mode": local_heavy_guard,
        "true_mode_RMS": metrics["E_P_RMS_ratio"] <= 0.80,
        "true_mode_heavy": ep_heavy_guard,
        "dominant_mode_positivity": positive_weight >= 0.95 * total_pod,
        "H_RAW": metrics["H_RAW_ratio"] <= 0.50}
    failed_gates = [name for name, passed in gates.items() if not passed]
    result = {"schema": "strict-ala-stage-F1b-offline-gate-v1",
              "experiment_evidence_valid": True, "metrics": metrics,
              "gates": gates, "offline_gate_pass": all(gates.values()),
              "independent_failed_gates": failed_gates,
              "uniform_scalar_rescaling_bound": scalar_bound,
              "candidate_rejected_even_if_memory_repaired":
                  (not gates["memory"] and any(
                      not gates[name] for name in gates if name != "memory")),
              "restriction_tolerance_review_required": replay_review_required,
              "maximum_restriction_identity_relative_defect": maximum_replay,
              "thresholds_sha256": sha256(args.thresholds),
              "memory_gates": {
                  "M1_retained_cache": m1,
                  "M2_total_process_peak": m2,
                  "M3_node_safety": m3,
                  "BASE_rank_rss_max_kib": base_memory["rank_rss_max_kib"],
                  "CAND_rank_rss_max_kib": cand_memory["rank_rss_max_kib"],
                  "BASE_retained_cache_bytes_per_rank":
                      base_memory["retained_cache_bytes_per_rank"],
                  "CAND_retained_cache_bytes_per_rank":
                      cand_memory["retained_cache_bytes_per_rank"],
                  "projected_candidate_node_peak_kib": projected_node_peak,
                  "allocated_node_memory_kib": frozen_thresholds["allocated_node_memory_kib"],
              },
              "startup_gate": {"pass": startup_pass,
                  "startup_seconds_max": feasibility["startup_seconds_max"],
                  "frozen_limit_seconds": startup_limit},
              "tight_solve_count": len(tight),
              "mode_manifest_sha256": sha256(mode_manifest)}
    write_json(output / "strict_ala_stage_F1b_offline_gate.json", result)
    return 0


def stage_rows(path, case, case_wall):
    source = read_csv(path)
    rows = []
    final_solver = finite(source[-1]["elapsed_seconds"], "elapsed_seconds")
    startup = max(case_wall - final_solver, 0.0)
    seen = set()
    for row in source:
        iteration = int(row["iteration"])
        if iteration in seen or iteration < 0 or iteration > 30:
            raise ValueError("invalid Stage-E iteration index")
        seen.add(iteration)
        values = {name: finite(row[name], name) for name in (
            "true_continuity_relative", "true_momentum_relative",
            "recursive_residual", "explicit_residual", "velocity_MG_cycles",
            "K_gamma_operator_applications", "elapsed_seconds")}
        rows.append({
            "case": case, "fgmres_candidate_iteration": iteration,
            "global_cancellation_L2": values["true_continuity_relative"],
            "unaugmented_momentum_relative": values["true_momentum_relative"],
            "recursive_krylov_residual": values["recursive_residual"],
            "explicit_krylov_residual": values["explicit_residual"],
            "cumulative_MG_cycles": int(values["velocity_MG_cycles"]),
            "cumulative_Kgamma_applications": int(values["K_gamma_operator_applications"]),
            "elapsed_solver_seconds": values["elapsed_seconds"],
            "elapsed_wall_seconds": startup + values["elapsed_seconds"],
            "terminal": 0, "terminal_reason": ""})
    if not rows or rows[0]["fgmres_candidate_iteration"] != 0:
        raise ValueError("short case lacks real iteration zero")
    last = rows[-1]
    joint = (last["global_cancellation_L2"] < 1.0e-2 and
             last["unaugmented_momentum_relative"] <= MOMENTUM_GATE)
    if joint and last["fgmres_candidate_iteration"] < 30:
        reason = "EARLY_JOINT_CONVERGENCE"
    elif last["fgmres_candidate_iteration"] == 30:
        reason = "ITERATION_30"
    else:
        raise ValueError("short case terminated without a declared terminal state")
    last["terminal"] = 1; last["terminal_reason"] = reason
    return rows


def best_at(rows, field, budget):
    eligible = [row for row in rows if row[field] <= budget + 1e-12]
    if not eligible:
        raise ValueError("no achieved point at common budget")
    # Earliest actual row wins an exact residual tie. Momentum is taken from
    # this same row, never from another state.
    return min(eligible, key=lambda row: (
        row["global_cancellation_L2"], row["fgmres_candidate_iteration"]))


def momentum_guard(base, cand):
    if base["unaugmented_momentum_relative"] <= MOMENTUM_GATE and \
            cand["unaugmented_momentum_relative"] > MOMENTUM_GATE:
        return False
    if cand["unaugmented_momentum_relative"] <= \
            1.10 * base["unaugmented_momentum_relative"]:
        return True
    return (base["unaugmented_momentum_relative"] < 1.0e-6 and
            cand["unaugmented_momentum_relative"] < 1.0e-6 and
            abs(cand["unaugmented_momentum_relative"] -
                base["unaugmented_momentum_relative"]) <= MOMENTUM_PRECISION)


def memory_row(path, expected_case):
    rows = read_csv(path)
    if len(rows) != 1 or rows[0]["case"] != expected_case:
        raise ValueError("rank memory evidence mismatch")
    return {key: finite(rows[0][key], key) for key in
            ("rank_rss_max_kib", "rank_rss_sum_kib", "ranks",
             "retained_cache_bytes_per_rank", "temporary_peak_bytes_per_rank")}


def short(args):
    base_wall = finite(Path(args.base_wall).read_text().strip(), "BASE wall")
    cand_wall = finite(Path(args.cand_wall).read_text().strip(), "CAND wall")
    base = stage_rows(args.base_iterations, "BASE", base_wall)
    cand = stage_rows(args.cand_iterations, "CAND", cand_wall)
    base_mem = memory_row(args.base_memory, "BASE")
    cand_mem = memory_row(args.cand_memory, "CAND")
    output = Path(args.output_dir); output.mkdir(parents=True, exist_ok=True)
    with (output / "strict_ala_stage_F1b_short_iterations.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=SHORT_COLUMNS)
        writer.writeheader(); writer.writerows(base + cand)
    costs = []
    for case, rows, wall, mem in (("BASE", base, base_wall, base_mem),
                                  ("CAND", cand, cand_wall, cand_mem)):
        last = rows[-1]
        costs.append({"case": case,
            "terminal_iteration": last["fgmres_candidate_iteration"],
            "total_MG_cycles": last["cumulative_MG_cycles"],
            "total_Kgamma_applications": last["cumulative_Kgamma_applications"],
            "solver_wall_seconds": last["elapsed_solver_seconds"],
            "case_wall_seconds": wall, **mem})
    with (output / "strict_ala_stage_F1b_short_cost.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=COST_COLUMNS)
        writer.writeheader(); writer.writerows(costs)
    bt, ct = base[-1], cand[-1]
    base_early = bt["terminal_reason"] == "EARLY_JOINT_CONVERGENCE"
    cand_early = ct["terminal_reason"] == "EARLY_JOINT_CONVERGENCE"
    early = "NOT_APPLICABLE"
    early_evidence = {}
    if cand_early and not base_early:
        reached = [row for row in base if row["global_cancellation_L2"] <=
                   ct["global_cancellation_L2"]]
        comparator = reached[0] if reached else bt
        ratios = (ct["fgmres_candidate_iteration"] /
                  max(comparator["fgmres_candidate_iteration"], 1),
                  ct["cumulative_MG_cycles"] / max(comparator["cumulative_MG_cycles"], 1),
                  ct["elapsed_wall_seconds"] / max(comparator["elapsed_wall_seconds"], 1e-300))
        early = "PASS" if max(ratios) <= 1.25 and momentum_guard(comparator, ct) else "FAIL"
        early_evidence = {"comparator_iteration": comparator["fgmres_candidate_iteration"],
                          "work_ratios": ratios,
                          "base_reached_candidate_target": bool(reached)}
    elif base_early and not cand_early:
        early = "FAIL"
    elif base_early and cand_early:
        ratios = (ct["fgmres_candidate_iteration"] / max(bt["fgmres_candidate_iteration"], 1),
                  ct["cumulative_MG_cycles"] / max(bt["cumulative_MG_cycles"], 1),
                  ct["elapsed_wall_seconds"] / max(bt["elapsed_wall_seconds"], 1e-300))
        early = "PASS" if max(ratios) <= 1.10 and min(ratios) <= 0.90 and \
            momentum_guard(bt, ct) else "FAIL"
        early_evidence = {"work_ratios": ratios}

    route_a = None
    if bt["fgmres_candidate_iteration"] == ct["fgmres_candidate_iteration"] == 30:
        route_a = (ct["global_cancellation_L2"] <= 0.85 * bt["global_cancellation_L2"] and
                   cand_wall <= 1.25 * base_wall and
                   ct["cumulative_MG_cycles"] <= 1.25 * bt["cumulative_MG_cycles"] and
                   momentum_guard(bt, ct))
    common_wall = min(bt["elapsed_wall_seconds"], ct["elapsed_wall_seconds"])
    bw, cw = best_at(base, "elapsed_wall_seconds", common_wall), \
             best_at(cand, "elapsed_wall_seconds", common_wall)
    route_bw = (cw["global_cancellation_L2"] <= 0.90 * bw["global_cancellation_L2"] and
                momentum_guard(bw, cw))
    common_mg = min(bt["cumulative_MG_cycles"], ct["cumulative_MG_cycles"])
    bm, cm = best_at(base, "cumulative_MG_cycles", common_mg), \
             best_at(cand, "cumulative_MG_cycles", common_mg)
    route_bm = (cm["global_cancellation_L2"] <= 0.90 * bm["global_cancellation_L2"] and
                momentum_guard(bm, cm))
    late = None
    if route_a is not None:
        bmap = {r["fgmres_candidate_iteration"]: r for r in base}
        cmap = {r["fgmres_candidate_iteration"]: r for r in cand}
        late = sum(cmap[k]["global_cancellation_L2"] < bmap[k]["global_cancellation_L2"]
                   for k in range(20, 31)) >= 8
    memory_pass = cand_mem["rank_rss_max_kib"] <= 1.25 * base_mem["rank_rss_max_kib"]
    performance = (early == "PASS" or route_a is True or (route_bw and route_bm))
    passed = performance and (late is not False) and memory_pass
    result = {"schema": "strict-ala-stage-F1b-short-decision-v1",
              "experiment_valid": True,
              "base_terminal_reason": bt["terminal_reason"],
              "cand_terminal_reason": ct["terminal_reason"],
              "route_A_pass": route_a, "route_B_wall_pass": route_bw,
              "route_B_MG_pass": route_bm,
              "early_joint_convergence_rule": early,
              "early_evidence": early_evidence,
              "late_window_stable": late, "memory_gate_pass": memory_pass,
              "matched_wall_rows": {"BASE": bw["fgmres_candidate_iteration"],
                                    "CAND": cw["fgmres_candidate_iteration"]},
              "matched_MG_rows": {"BASE": bm["fgmres_candidate_iteration"],
                                  "CAND": cm["fgmres_candidate_iteration"]},
              "short_AB_pass": passed}
    write_json(output / "strict_ala_stage_F1b_short_decision.json", result)
    return 0


def final_audit(args):
    offline_result = load_json(args.offline)
    short_result = load_json(args.short) if args.short else None
    # Evidence integrity is independent of whether the scientific candidate
    # passes. A deterministic negative result remains a valid experiment.
    provenance = {}
    for label in ("binary", "inputs", "source"):
        before = Path(getattr(args, f"{label}_pre"))
        after = Path(getattr(args, f"{label}_post"))
        provenance[label] = {
            "pre_sha256": sha256(before), "post_sha256": sha256(after),
            "unchanged": before.read_bytes() == after.read_bytes()}
    provenance_valid = all(item["unchanged"] for item in provenance.values())
    experiment_valid = (offline_result.get("experiment_evidence_valid") is True
                        and provenance_valid)
    gates = offline_result.get("gates", {})
    independent_rejections = []
    if experiment_valid and not offline_result.get("offline_gate_pass"):
        rejection_map = (
            ("candidate_isolatable", "CANDIDATE_NOT_ISOLATABLE"),
            ("memory", "CANDIDATE_REJECTED_MEMORY"),
            ("startup", "F1B_STARTUP_COST_REVIEW"),
            ("assembly", "CANDIDATE_REJECTED_ASSEMBLY"),
            ("tight_solves", "CANDIDATE_REJECTED_ASSEMBLY"),
            ("local_RMS", "CANDIDATE_REJECTED_LOCAL_OPERATOR"),
            ("local_heavy_mode", "CANDIDATE_REJECTED_LOCAL_OPERATOR"),
            ("H_RAW", "CANDIDATE_REJECTED_RAW_LOCAL"),
            ("true_mode_RMS", "CANDIDATE_REJECTED_ACTUAL_MODES"),
            ("true_mode_heavy", "CANDIDATE_REJECTED_ACTUAL_MODES"),
            ("dominant_mode_positivity", "CANDIDATE_REJECTED_ACTUAL_MODES"),
        )
        for gate, rejection in rejection_map:
            if gates.get(gate, True) is False and rejection not in independent_rejections:
                independent_rejections.append(rejection)
    if not experiment_valid:
        decision = "INVALID_EXPERIMENT"
    elif not offline_result.get("offline_gate_pass"):
        if offline_result.get("restriction_tolerance_review_required"):
            decision = "F1B_RESTRICTION_TOLERANCE_REVIEW"
        elif not gates.get("candidate_isolatable"): decision = "CANDIDATE_NOT_ISOLATABLE"
        elif not gates.get("memory"): decision = "CANDIDATE_REJECTED_MEMORY"
        elif not gates.get("startup", True): decision = "F1B_STARTUP_COST_REVIEW"
        elif not gates.get("assembly"): decision = "CANDIDATE_REJECTED_ASSEMBLY"
        elif not gates.get("local_RMS") or not gates.get("local_heavy_mode"):
            decision = "CANDIDATE_REJECTED_LOCAL_OPERATOR"
        elif not gates.get("H_RAW"): decision = "CANDIDATE_REJECTED_RAW_LOCAL"
        else: decision = "CANDIDATE_REJECTED_ACTUAL_MODES"
    elif short_result is None or not short_result.get("experiment_valid"):
        decision = "INVALID_EXPERIMENT"
        experiment_valid = False
    elif not short_result.get("short_AB_pass"):
        decision = ("SHORT_AB_UNSTABLE" if short_result.get("late_window_stable") is False
                    else "CANDIDATE_REJECTED_SHORT_PRODUCTION")
    else:
        decision = "CANDIDATE_PASSES_SHORT_AB"
    if decision == "CANDIDATE_PASSES_SHORT_AB":
        next_task = "F1C_FULL_PRODUCTION_QUALIFICATION"
    elif decision == "INVALID_EXPERIMENT":
        next_task = "REPEAT_F1B_VALID_EXPERIMENT"
    elif decision == "F1B_RESTRICTION_TOLERANCE_REVIEW":
        next_task = "REVIEW_RESTRICTION_IMPLEMENTATION_NOT_THRESHOLD"
    elif decision == "F1B_STARTUP_COST_REVIEW":
        next_task = "REVIEW_F1B_STARTUP_COST"
    else:
        next_task = "F1B_LOCAL_OPERATOR_REDESIGN"
    result = {"schema": "strict-ala-stage-F1b-final-audit-v1",
              "experiment_valid": experiment_valid,
              "candidate_valid": offline_result.get("offline_gate_pass", False),
              "candidate_accepted": decision == "CANDIDATE_PASSES_SHORT_AB",
              "provenance": provenance,
              "offline": offline_result, "short_AB": short_result,
              "decision": decision,
              "independent_rejections": independent_rejections,
              "candidate_rejected_even_if_memory_repaired":
                  offline_result.get(
                      "candidate_rejected_even_if_memory_repaired", False),
              "uniform_scalar_rescaling_can_reach_true_mode_gate":
                  offline_result.get("uniform_scalar_rescaling_bound", {}).get(
                      "can_reach_frozen_gate"),
              "production_default_change_authorized": False,
              "production_freeze_authorized": False,
              "production_freeze_requires": [
                  "F1c pressure physics qualification",
                  "G0 qualified-solver profile",
                  "production-representative coupled smoke/soak",
                  "independent repeat"],
              "next_authorized_task": next_task}
    write_json(args.output, result)
    return 0


def parser():
    root = argparse.ArgumentParser(); commands = root.add_subparsers(required=True)
    p = commands.add_parser("freeze-contract")
    p.add_argument("--output-dir", required=True)
    for name in ("allocated-tasks", "solver-world-size", "nprocx", "nprocy",
                 "nprocz", "caps", "ranks-per-node", "node-memory-kib",
                 "allocated-wall-seconds", "startup-budget-seconds"):
        p.add_argument(f"--{name}", required=True)
    p.set_defaults(fn=freeze_contract)
    p = commands.add_parser("verify-layout")
    for name in ("cfg", "thresholds", "allocated-tasks", "solver-world-size",
                 "caps", "ranks-per-node", "output"):
        p.add_argument(f"--{name}", required=True)
    p.set_defaults(fn=verify_layout)
    p = commands.add_parser("check-cfg")
    for name in ("base", "candidate", "diff"): p.add_argument(f"--{name}", required=True)
    p.set_defaults(fn=check_cfg)
    p = commands.add_parser("set-cfg")
    for name in ("cfg", "section", "key", "value"):
        p.add_argument(f"--{name}", required=True)
    p.set_defaults(fn=set_cfg_value)
    p = commands.add_parser("verify-lineage")
    for name in ("final-audit", "evidence-root", "expected-hpc-root", "output"):
        p.add_argument(f"--{name}", required=True)
    p.set_defaults(fn=verify_lineage)
    p = commands.add_parser("identity")
    for name in ("feasibility", "assembly", "thresholds", "schema", "output"):
        p.add_argument(f"--{name}", required=True)
    p.set_defaults(fn=identity)
    p = commands.add_parser("offline")
    for name in ("lineage", "identity-gate", "feasibility", "assembly", "legacy-local",
                 "candidate-local", "legacy-telescoping", "candidate-telescoping",
                 "legacy-decision", "candidate-decision", "candidate-tight",
                 "snapshot-manifest", "thresholds", "schema", "base-memory",
                 "cand-memory", "output-dir"):
        p.add_argument(f"--{name}", required=True)
    p.set_defaults(fn=offline)
    p = commands.add_parser("short")
    for name in ("base-iterations", "cand-iterations", "base-memory", "cand-memory",
                 "base-wall", "cand-wall", "output-dir"):
        p.add_argument(f"--{name}", required=True)
    p.set_defaults(fn=short)
    p = commands.add_parser("audit")
    p.add_argument("--offline", required=True); p.add_argument("--short")
    for name in ("binary-pre", "binary-post", "inputs-pre", "inputs-post",
                 "source-pre", "source-post"):
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--output", required=True); p.set_defaults(fn=final_audit)
    return root


if __name__ == "__main__":
    arguments = parser().parse_args()
    raise SystemExit(arguments.fn(arguments))
