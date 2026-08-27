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
K_REPLAY = 5.0e-6
SYMMETRY = 1.0e-10
LOCAL_SOLVE = 1.0e-10
MOMENTUM_GATE = 1.0e-3
MOMENTUM_PRECISION = 1.0e-12

ASSEMBLY_COLUMNS = (
    "row_type", "patch_category", "rank", "patch_id", "pressure_dim",
    "velocity_dim", "symmetry_defect", "regularization",
    "factorization_status", "pivot_metric", "local_solve_residual",
    "fallback", "candidate_matrix_bytes", "candidate_factor_bytes",
    "Kgamma_replay_relative_defect")
SHORT_COLUMNS = (
    "case", "fgmres_candidate_iteration", "global_cancellation_L2",
    "unaugmented_momentum_relative", "recursive_krylov_residual",
    "explicit_krylov_residual", "cumulative_MG_cycles",
    "cumulative_Kgamma_applications", "elapsed_solver_seconds",
    "elapsed_wall_seconds", "terminal", "terminal_reason")
COST_COLUMNS = (
    "case", "terminal_iteration", "total_MG_cycles",
    "total_Kgamma_applications", "solver_wall_seconds", "case_wall_seconds",
    "rank_rss_max_kib", "rank_rss_sum_kib", "ranks")


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


def thresholds():
    return {
        "schema": "strict-ala-stage-F1b-thresholds-v1",
        "threshold_status": "PROVISIONAL_DIAGNOSTIC",
        "MAX_F1B_TIGHT_SGAMMA_SOLVES": MAX_TIGHT_SOLVES,
        "Kgamma_replay_relative_tolerance": K_REPLAY,
        "candidate_symmetry_tolerance": SYMMETRY,
        "local_solve_relative_tolerance": LOCAL_SOLVE,
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
    }


def schema_manifest():
    return {"schema": "strict-ala-stage-F1b-schema-v1", "files": [
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
    write_json(output / "strict_ala_stage_F1b_thresholds.json", thresholds())
    write_json(output / "strict_ala_stage_F1b_schema.json", schema_manifest())
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


def offline(args):
    output = Path(args.output_dir); output.mkdir(parents=True, exist_ok=True)
    lineage = load_json(args.lineage)
    feasibility = load_json(args.feasibility)
    for field in ("factorization_count", "RHS_solve_count",
                  "estimated_dense_factor_and_solve_work",
                  "temporary_peak_bytes_per_rank",
                  "retained_cache_bytes_per_rank",
                  "projected_peak_ratio_vs_legacy_cache",
                  "new_mpi_payload_bytes_max_per_rank",
                  "diagonal_completion_relative_max",
                  "Kgamma_replay_relative_defect", "startup_seconds_max"):
        if finite(feasibility[field], field) < 0:
            raise ValueError(f"negative feasibility field {field}")
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
                      "Kgamma_replay_relative_defect"):
            finite(row[field], field)
        if int(row["pressure_dim"]) <= 0 or int(row["velocity_dim"]) <= 0:
            raise ValueError("candidate assembly has an empty aggregate")
        if (row["row_type"] != "AGGREGATE" or
                row["factorization_status"] != "PASS" or row["fallback"] != "0"):
            raise ValueError("candidate assembly status failed")
    assembly_pass = (max(finite(r["symmetry_defect"], "symmetry") for r in assembly)
                     <= SYMMETRY and
                     max(finite(r["local_solve_residual"], "solve") for r in assembly)
                     <= LOCAL_SOLVE and
                     max(finite(r["Kgamma_replay_relative_defect"], "replay")
                         for r in assembly) <= K_REPLAY)
    assembly_pass = (assembly_pass and
                     finite(feasibility["Kgamma_replay_relative_defect"],
                            "feasibility replay") <= K_REPLAY)

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
        "memory": feasibility.get("memory_gate_pass") is True,
        "assembly": assembly_pass, "tight_solves": tight_pass,
        "local_RMS": metrics["E_local_RMS_ratio"] <= 0.70,
        "local_heavy_mode": local_heavy_guard,
        "true_mode_RMS": metrics["E_P_RMS_ratio"] <= 0.80,
        "true_mode_heavy": ep_heavy_guard,
        "dominant_mode_positivity": positive_weight >= 0.95 * total_pod,
        "H_RAW": metrics["H_RAW_ratio"] <= 0.50}
    result = {"schema": "strict-ala-stage-F1b-offline-gate-v1",
              "experiment_evidence_valid": True, "metrics": metrics,
              "gates": gates, "offline_gate_pass": all(gates.values()),
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
            ("rank_rss_max_kib", "rank_rss_sum_kib", "ranks")}


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
    if not experiment_valid:
        decision = "INVALID_EXPERIMENT"
    elif not offline_result.get("offline_gate_pass"):
        gates = offline_result["gates"]
        if not gates.get("candidate_isolatable"): decision = "CANDIDATE_NOT_ISOLATABLE"
        elif not gates.get("memory"): decision = "CANDIDATE_REJECTED_MEMORY"
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
    else:
        next_task = "F1B_LOCAL_OPERATOR_REDESIGN"
    result = {"schema": "strict-ala-stage-F1b-final-audit-v1",
              "experiment_valid": experiment_valid,
              "candidate_valid": offline_result.get("offline_gate_pass", False),
              "candidate_accepted": decision == "CANDIDATE_PASSES_SHORT_AB",
              "provenance": provenance,
              "offline": offline_result, "short_AB": short_result,
              "decision": decision,
              "production_default_change_authorized": False,
              "next_authorized_task": next_task}
    write_json(args.output, result)
    return 0


def parser():
    root = argparse.ArgumentParser(); commands = root.add_subparsers(required=True)
    p = commands.add_parser("freeze-contract"); p.add_argument("--output-dir", required=True); p.set_defaults(fn=freeze_contract)
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
    p = commands.add_parser("offline")
    for name in ("lineage", "feasibility", "assembly", "legacy-local",
                 "candidate-local", "legacy-telescoping", "candidate-telescoping",
                 "legacy-decision", "candidate-decision", "candidate-tight",
                 "snapshot-manifest", "output-dir"):
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
