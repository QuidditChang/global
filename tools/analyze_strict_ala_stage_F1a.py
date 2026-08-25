#!/usr/bin/env python3
"""Fail-closed analysis for the strict-ALA Stage-F1a Schwarz audit.

The numerical executable writes primitive inner products and energies.  This
tool owns all provisional scientific thresholds and makes no solver changes.
"""

import argparse
import csv
import json
import math
import re
from pathlib import Path

PROVISIONAL = "PROVISIONAL_DIAGNOSTIC"
MATERIAL = 0.10
PRIMARY = 0.60
CONSISTENCY = 0.80
MIXED_MIN = 0.25
MIXED_TOTAL = 0.75
LOCAL_SOLVE = 1.0e-10
LOCAL_ACTION = 0.30
NONLOCAL = 0.30
MAX_PATCHES = 128
MAX_TIGHT_SOLVES = 40
REPLAY_RELATIVE = 1.0e-8
SCALE_FLOOR = 1.0e-300

STAGES = ("BPI", "S0_RAW_LOCAL", "S1_RIGHT_SCALING",
          "S2_PARTITION_WEIGHTING", "S3_LEFT_SCALING", "CONFIGURED")
CHAIN = ("RAW_LOCAL", "RIGHT_SCALING", "PARTITION_WEIGHTING", "LEFT_SCALING")
PRIMARY_CLASSES = {"LOCAL_SOLVE_FAILURE", "LOCAL_OPERATOR_MISMATCH",
                   "DEPTH_SCALING_DAMAGE", "PARTITION_WEIGHTING_DAMAGE",
                   "ADDITIVE_INTERFERENCE", "PATCH_TOPOLOGY_MISMATCH",
                   "NONLOCAL_SCHUR_MISMATCH"}

TELESCOPING_COLUMNS = (
    "POD_mode", "stage", "POD_energy_weight", "x_norm", "Sx_norm",
    "residual_energy", "standalone_residual_energy", "E_P", "cosine",
    "qTPq", "qTSPq", "eB_yS", "yB_yS", "self_term",
    "useful_cross_term", "delta_energy", "delta_identity_error",
    "tight_solve_achieved", "valid")
TIGHT_COLUMNS = ("call_id", "role", "POD_mode", "rhs_norm",
                 "requested_relative_tolerance", "target_residual",
                 "achieved_relative_residual", "cycles", "max_cycles",
                 "converged")
PATCH_COLUMNS = ("patch_ID", "selection_rule", "POD_mode", "depth_stratum",
                 "contribution_quantile", "deterministic_order",
                 "outlier_guard")
LOCAL_COLUMNS = ("patch_or_bin_ID", "POD_mode", "local_RHS_norm",
                 "local_solution_norm", "local_solve_relative_residual",
                 "Ai_v_norm", "Ti_v_norm", "local_action_error",
                 "local_action_cosine", "contribution_weight",
                 "depth_category", "scale_category", "valid")
NONLOCAL_COLUMNS = ("source_vector_or_bin", "inside_patch_energy",
                    "distance_1_energy", "distance_2_energy",
                    "distance_3_energy", "farther_energy",
                    "outside_patch_plus_one_fraction", "valid")
PROJECTED_COLUMNS = ("matrix", "row_mode", "column_mode", "value")
STRUCTURE_COLUMNS = ("supported_DOFs", "multiplicity_min", "multiplicity_max",
                     "invalid_multiplicity", "top_W", "mid_W", "transition_W",
                     "partition_formula", "valid")
SUPPORT_COLUMNS = ("POD_mode", "q_energy_in_support_fraction",
                   "eB_energy_in_support_fraction", "valid")
LINEARITY_COLUMNS = ("role", "difference_norm", "combined_norm", "sum_norm",
                     "relative_defect", "solve_floor", "warn_limit",
                     "hard_limit", "status", "valid")
F1A_REQUIRED_OFF = ("ala_radial_line_preconditioner", "ala_two_level_preconditioner",
                    "ala_pressure_multigrid", "ala_pressure_multigrid_galerkin",
                    "ala_global_coarse_preconditioner", "ala_geneo_preconditioner")
F1A_REQUIRED_VALUES = {
    "ala_pressure_bpi_weight": "1.0",
    "ala_augmented_lagrangian_gamma": "10.0",
    "ala_outer_solver": "fgmres",
    "ala_pcg_restart_interval": "50",
    "ala_shallow_patch_velocity_solver": "element_vanka",
    "nprocx": "4", "nprocy": "4", "nprocz": "2",
    "nodex": "129", "nodey": "129", "nodez": "65",
    "levels": "5", "steps": "1", "piterations": "60",
}


def finite(value, name):
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"non-finite {name}")
    return number


def nonlocal_mode(identifier):
    match = re.search(r"_mode_([1-9][0-9]*)$", identifier)
    if match is None:
        raise ValueError("nonlocality row does not identify its POD mode")
    return int(match.group(1))


def read_csv(path, columns, allow_empty=False):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing or empty {path.name}")
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != tuple(columns):
            raise ValueError(f"schema mismatch in {path.name}")
        rows = list(reader)
    if not rows and not allow_empty:
        raise ValueError(f"no data rows in {path.name}")
    return rows


def cfg_values(path):
    values = {}
    for line in Path(path).read_text().splitlines():
        match = re.match(r"^\s*([A-Za-z0-9_]+)\s*(?:=\s*|\s+)(\S.*?)\s*$", line)
        if match and not match.group(1).startswith("#"):
            values[match.group(1)] = match.group(2)
    return values


def check_cfg_contract(c0, c1, diff_path):
    baseline, configured = cfg_values(c0), cfg_values(c1)
    keys = sorted(set(baseline) | set(configured))
    differences = [(key, baseline.get(key), configured.get(key)) for key in keys
                   if baseline.get(key) != configured.get(key)]
    Path(diff_path).write_text("\n".join(
        "%s: %s -> %s" % difference for difference in differences) + "\n")
    errors = []
    if differences != [("ala_shallow_patch_preconditioner", "off", "on")]:
        errors.append("unexpected scientific diff %r" % (differences,))
    for key in F1A_REQUIRED_OFF:
        if baseline.get(key) != "off" or configured.get(key) != "off":
            errors.append("%s must remain off" % key)
    for key, expected in F1A_REQUIRED_VALUES.items():
        if baseline.get(key) != expected or configured.get(key) != expected:
            errors.append("%s must equal %s" % (key, expected))
    for key in ("ala_stage_abc_adjoint_diagnostic", "ala_stage_abc_production_logging"):
        if baseline.get(key, "off") != "off" or configured.get(key, "off") != "off":
            errors.append("%s must be absent or off for F1a" % key)
    if errors:
        raise ValueError("F1a cfg contract failed: " + "; ".join(errors))
    return 0


def validate_linearity(rows):
    if len(rows) != 1 or rows[0]["role"] != "CONFIGURED":
        raise ValueError("linearity guard completeness violation")
    row = rows[0]
    for column in LINEARITY_COLUMNS[1:8]:
        value = finite(row[column], column)
        if value < 0.0:
            raise ValueError("negative linearity metric")
    if (row["valid"] != "1" or row["status"] not in ("PASS", "WARN") or
            float(row["relative_defect"]) > float(row["hard_limit"]) or
            float(row["warn_limit"]) > float(row["hard_limit"])):
        raise ValueError("normalized linearity guard failed")
    return row


def thresholds():
    return {
        "schema": "strict-ala-stage-F1a-thresholds-v1",
        "threshold_status": PROVISIONAL,
        "MATERIAL_EFFECT_FRACTION": MATERIAL,
        "PRIMARY_EXPLAINED_HARM_FRACTION": PRIMARY,
        "MODE_CONSISTENCY_ENERGY_FRACTION": CONSISTENCY,
        "MIXED_CLASS_MIN_FRACTION": MIXED_MIN,
        "MIXED_CLASSES_COMBINED_FRACTION": MIXED_TOTAL,
        "LOCAL_SOLVE_RELATIVE_RESIDUAL_LIMIT": LOCAL_SOLVE,
        "LOCAL_ACTION_ERROR_THRESHOLD": LOCAL_ACTION,
        "NONLOCAL_OUTSIDE_PATCH_PLUS_ONE_THRESHOLD": NONLOCAL,
        "MPI_INTERFACE_SUSPECT_THRESHOLD": 1.5,
        "MAX_REPRESENTATIVE_PATCHES": MAX_PATCHES,
        "MAX_TIGHT_SGAMMA_SOLVES": MAX_TIGHT_SOLVES,
        "OPERATOR_REPLAY_RELATIVE_TOLERANCE": REPLAY_RELATIVE,
        "LINEARITY_RELATIVE_WARN_FLOOR": 1.0e-8,
        "LINEARITY_RELATIVE_HARD_FLOOR": 1.0e-5,
        "scale_floor": SCALE_FLOOR,
        "normalization": "all residual-energy changes divided by ||e_B||^2",
        "additive_rule": "standalone Schwarz useful, configured map harmful, interaction explains >=60%",
        "priority": list(CHAIN) + ["ADDITION_TO_BPI"],
    }


def schema_manifest():
    entries = (
        ("strict_ala_stage_F1a_telescoping.csv", TELESCOPING_COLUMNS,
         ["POD_mode", "stage"], "6 rows per selected mode"),
        ("strict_ala_stage_F1a_tight_solves.csv", TIGHT_COLUMNS,
         ["call_id"], "1 BPI + 4 Schwarz actions per mode, plus one guard"),
        ("strict_ala_stage_F1a_patch_selection.csv", PATCH_COLUMNS,
         ["patch_ID", "POD_mode", "selection_rule"], "1..128 rows"),
        ("strict_ala_stage_F1a_local_operator.csv", LOCAL_COLUMNS,
         ["patch_or_bin_ID", "POD_mode"], "one row per audited patch/bin and mode"),
        ("strict_ala_stage_F1a_nonlocality.csv", NONLOCAL_COLUMNS,
         ["source_vector_or_bin"], "at most 8 bin-summed actions"),
        ("strict_ala_stage_F1a_projected_matrices.csv", PROJECTED_COLUMNS,
         ["matrix", "row_mode", "column_mode"], "6*Nmode*Nmode rows"),
        ("strict_ala_stage_F1a_patch_structure.csv", STRUCTURE_COLUMNS,
         ["partition_formula"], "exactly one row"),
        ("strict_ala_stage_F1a_support.csv", SUPPORT_COLUMNS,
         ["POD_mode"], "one row per selected mode"),
        ("strict_ala_stage_F1a_linearity.csv", LINEARITY_COLUMNS,
         ["role"], "exactly one normalized configured-map guard"),
    )
    return {"schema": "strict-ala-stage-F1a-schema-v1", "files": [
        {"filename": name, "required_columns": list(columns),
         "unique_key": key, "row_count_rule": rule,
         "nullable_columns": [], "finite_numeric_fields_required": True}
        for name, columns, key, rule in entries]}


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def freeze_contract(output):
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    write_json(output / "strict_ala_stage_F1a_thresholds.json", thresholds())
    write_json(output / "strict_ala_stage_F1a_schema.json", schema_manifest())


def sign_coverage(values, weights, sign):
    return sum(w for value, w in zip(values, weights)
               if (value > 0) == (sign > 0) and value != 0)


def classify_telescoping(rows):
    keyed = {}
    weights = {}
    for row in rows:
        mode = int(row["POD_mode"]); stage = row["stage"]
        if stage not in STAGES or (mode, stage) in keyed:
            raise ValueError("invalid or duplicate telescoping key")
        if row["valid"] != "1":
            raise ValueError("invalid required telescoping action")
        for column in TELESCOPING_COLUMNS[2:-1]:
            finite(row[column], column)
        keyed[mode, stage] = row
        weights[mode] = finite(row["POD_energy_weight"], "POD_energy_weight")
    modes = sorted(weights)
    if not modes or modes != list(range(1, len(modes) + 1)) or len(modes) > 6:
        raise ValueError("selected POD modes are not contiguous 95%-or-6 modes")
    if any((mode, stage) not in keyed for mode in modes for stage in STAGES):
        raise ValueError("missing telescoping row")
    total_weight = sum(weights.values())
    if not (0.95 - 1e-12 <= total_weight <= 1.0 + 1e-12):
        raise ValueError("selected POD energy does not reach 95%")
    normalized_weights = [weights[m] / total_weight for m in modes]
    per_mode = {}
    increments = {name: [] for name in CHAIN}
    addition = []
    standalone = []
    interaction = []
    for mode in modes:
        base = finite(keyed[mode, "BPI"]["residual_energy"], "BPI energy")
        denom = max(base, SCALE_FLOOR)
        h = []
        for stage in STAGES[1:5]:
            value = (finite(keyed[mode, stage]["residual_energy"], stage) - base) / denom
            h.append(value)
        inc = [h[0], h[1] - h[0], h[2] - h[1], h[3] - h[2]]
        cfg = (finite(keyed[mode, "CONFIGURED"]["residual_energy"], "configured") - base) / denom
        stand = (finite(keyed[mode, "S3_LEFT_SCALING"]["standalone_residual_energy"],
                        "standalone") - 1.0) / denom
        inter = cfg - stand
        reported_inter = 2.0 * finite(keyed[mode, "S3_LEFT_SCALING"]["yB_yS"], "yB_yS") / denom
        closure = abs(inter - reported_inter)
        if closure > 1e-8 * max(1.0, abs(inter), abs(reported_inter)):
            raise ValueError("additive interference identity failed")
        for name, value in zip(CHAIN, inc): increments[name].append(value)
        addition.append(cfg); standalone.append(stand); interaction.append(inter)
        per_mode[mode] = {"H": h, "increments": dict(zip(CHAIN, inc)),
                          "configured_harm": cfg, "standalone_schwarz": stand,
                          "additive_interaction": inter}
    weighted = {name: sum(v * w for v, w in zip(values, normalized_weights))
                for name, values in increments.items()}
    weighted_cfg = sum(v * w for v, w in zip(addition, normalized_weights))
    weighted_standalone = sum(v * w for v, w in zip(standalone, normalized_weights))
    weighted_interaction = sum(v * w for v, w in zip(interaction, normalized_weights))
    consistency = {name: sign_coverage(values, normalized_weights, weighted[name])
                   for name, values in increments.items()}
    positive_sum = sum(max(value, 0.0) for value in weighted.values())
    contributions = {name: max(value, 0.0) / max(positive_sum, SCALE_FLOOR)
                     for name, value in weighted.items()}
    material = {name: abs(weighted[name]) >= MATERIAL and consistency[name] >= CONSISTENCY
                for name in CHAIN}
    first = next((name for name in CHAIN if weighted[name] > 0 and material[name]), "NONE")
    primary_stages = [name for name in CHAIN if weighted[name] > 0 and material[name]
                      and contributions[name] >= PRIMARY]
    dominant = primary_stages[0] if len(primary_stages) == 1 else "UNRESOLVED"
    mixed_members = [name for name in CHAIN if contributions[name] >= MIXED_MIN
                     and weighted[name] > 0]
    mixed = (not primary_stages and len(mixed_members) >= 2 and
             sum(contributions[name] for name in mixed_members) >= MIXED_TOTAL)

    # Addition is an independent composition diagnostic, not a telescoping
    # transformation: D_cfg-D_standalone = 2<y_B,y_S>.
    addition_consistency = sign_coverage(addition, normalized_weights, weighted_cfg)
    standalone_consistency = sum(w for value, w in zip(standalone, normalized_weights)
                                 if value <= -MATERIAL)
    interaction_fraction = max(weighted_interaction, 0.0) / max(weighted_cfg, SCALE_FLOOR)
    additive_primary = (weighted_cfg >= MATERIAL and weighted_standalone <= -MATERIAL and
                        addition_consistency >= CONSISTENCY and
                        standalone_consistency >= CONSISTENCY and
                        interaction_fraction >= PRIMARY)
    if additive_primary:
        first = "ADDITION_TO_BPI" if first == "NONE" else first
        dominant = "ADDITION_TO_BPI"
        primary_class = "ADDITIVE_INTERFERENCE"
        explained = min(interaction_fraction, 1.0)
        mode_consistency = min(addition_consistency, standalone_consistency)
    elif mixed:
        primary_class = "MIXED"; dominant = "MIXED"
        explained = max(contributions.values(), default=0.0)
        mode_consistency = max(consistency.values(), default=0.0)
    elif dominant in ("RIGHT_SCALING", "LEFT_SCALING"):
        primary_class = "DEPTH_SCALING_DAMAGE"
        explained = contributions[dominant]; mode_consistency = consistency[dominant]
    elif dominant == "PARTITION_WEIGHTING":
        primary_class = "PARTITION_WEIGHTING_DAMAGE"
        explained = contributions[dominant]; mode_consistency = consistency[dominant]
    elif dominant == "RAW_LOCAL":
        # The local secondary tree is resolved only after bounded A_i/T_i and
        # nonlocality evidence is merged by final_audit.
        primary_class = "UNRESOLVED"
        explained = contributions[dominant]; mode_consistency = consistency[dominant]
    else:
        primary_class = "UNRESOLVED"
        explained = max(contributions.values(), default=0.0)
        mode_consistency = max(consistency.values(), default=0.0)
    return {"modes": modes, "selected_energy": total_weight,
            "per_mode": per_mode, "weighted_telescoping_damage": weighted,
            "configured_harm": weighted_cfg,
            "standalone_schwarz_effect": weighted_standalone,
            "additive_interaction": weighted_interaction,
            "additive_interaction_explained_fraction": interaction_fraction,
            "contributions": contributions, "consistency": consistency,
            "first_material_damage_stage": first,
            "dominant_damage_stage": dominant,
            "primary_defect_class": primary_class,
            "primary_explained_harm_fraction": explained,
            "mode_consistency_energy_fraction": mode_consistency}


def validate_operator_replay(rows, e2_true_mode):
    archived = read_csv(e2_true_mode, ("mode", "map", "mode_energy_fraction",
        "cumulative_energy", "E_P", "cosine", "alpha_opt", "Pq_norm",
        "SPq_norm", "qTPq", "qTSPq", "positive_qTPq", "positive_qTSPq",
        "tight_solve_achieved", "map_semantics"))
    current = {(int(row["POD_mode"]), row["stage"]): row for row in rows}
    checked = []
    for old in archived:
        mode = int(old["mode"])
        if mode not in (1, 2) or old["map"] not in ("BPI", "CONFIGURED"):
            continue
        row = current[mode, old["map"]]
        differences = {}
        for field in ("E_P", "cosine", "qTSPq"):
            expected = finite(old[field], f"archived {field}")
            actual = finite(row[field], f"replayed {field}")
            relative = abs(actual - expected) / max(abs(expected), 1e-300)
            differences[field] = relative
            if relative > REPLAY_RELATIVE:
                raise ValueError(f"OPERATOR_REPLAY_IDENTITY_FAIL {mode} {old['map']} {field}")
        checked.append({"mode": mode, "map": old["map"],
                        "relative_differences": differences})
    if len(checked) != 4:
        raise ValueError("operator replay did not cover E2 modes 1 and 2")
    return checked


def analyze(input_dir, output_dir, e2_true_mode, e2_reanalysis):
    root = Path(input_dir); output = Path(output_dir); output.mkdir(parents=True, exist_ok=True)
    freeze_contract(output)
    rows = read_csv(root / "strict_ala_stage_F1a_telescoping.csv", TELESCOPING_COLUMNS)
    replay = validate_operator_replay(rows, e2_true_mode)
    e2 = json.loads(Path(e2_reanalysis).read_text())
    if e2.get("next_authorized_path") != "LOCAL_SCHWARZ_PATH" or not \
            e2.get("forensic_path_authorized") or \
            e2.get("production_schwarz_modification_authorized") is not False:
        raise ValueError("corrected E2 authorization gate failed")
    tight = read_csv(root / "strict_ala_stage_F1a_tight_solves.csv", TIGHT_COLUMNS)
    if len(tight) > MAX_TIGHT_SOLVES or len({r["call_id"] for r in tight}) != len(tight):
        raise ValueError("tight solve budget or unique-key violation")
    for row in tight:
        for column in TIGHT_COLUMNS[3:9]:
            finite(row[column], column)
        requested = float(row["requested_relative_tolerance"])
        rhs_norm = float(row["rhs_norm"])
        target = float(row["target_residual"])
        achieved = float(row["achieved_relative_residual"])
        if (row["converged"] != "1" or requested <= 0.0 or rhs_norm < 0.0 or
                target < 0.0 or achieved < 0.0 or
                achieved > max(2.0 * requested, 1.0e-15) or
                int(row["cycles"]) < 0 or
                int(row["cycles"]) > int(row["max_cycles"])):
            raise ValueError("failed tight solve contract")
    patches = read_csv(root / "strict_ala_stage_F1a_patch_selection.csv", PATCH_COLUMNS)
    if (len(patches) > MAX_PATCHES or len({(r["patch_ID"], r["POD_mode"],
            r["selection_rule"]) for r in patches}) != len(patches)):
        raise ValueError("representative patch budget exceeded")
    local = read_csv(root / "strict_ala_stage_F1a_local_operator.csv", LOCAL_COLUMNS)
    nonlocal_rows = read_csv(root / "strict_ala_stage_F1a_nonlocality.csv", NONLOCAL_COLUMNS)
    if len(nonlocal_rows) > 8:
        raise ValueError("bin-summed true-Schur budget exceeded")
    if any(row["valid"] != "1" for row in local + nonlocal_rows):
        raise ValueError("invalid representative patch action")
    if len({(r["patch_or_bin_ID"], r["POD_mode"]) for r in local}) != len(local):
        raise ValueError("duplicate local operator key")
    if len({r["source_vector_or_bin"] for r in nonlocal_rows}) != len(nonlocal_rows):
        raise ValueError("duplicate nonlocality key")
    for row in local:
        for column in LOCAL_COLUMNS[2:-3]:
            finite(row[column], column)
    nonlocal_modes = {}
    for row in nonlocal_rows:
        for column in NONLOCAL_COLUMNS[1:-1]:
            finite(row[column], column)
        nonlocal_modes[row["source_vector_or_bin"]] = nonlocal_mode(
            row["source_vector_or_bin"])
    selected_modes = sorted({int(row["POD_mode"]) for row in rows})
    projected = read_csv(root / "strict_ala_stage_F1a_projected_matrices.csv", PROJECTED_COLUMNS)
    expected_projected = 6 * len(selected_modes) * len(selected_modes)
    if len(projected) != expected_projected or len({(r["matrix"], r["row_mode"],
            r["column_mode"]) for r in projected}) != expected_projected:
        raise ValueError("projected matrix completeness violation")
    for row in projected: finite(row["value"], "projected matrix value")
    structure = read_csv(root / "strict_ala_stage_F1a_patch_structure.csv", STRUCTURE_COLUMNS)
    if len(structure) != 1 or structure[0]["valid"] != "1" or \
            int(structure[0]["invalid_multiplicity"]) != 0:
        raise ValueError("patch/multiplicity structural audit failed")
    support = read_csv(root / "strict_ala_stage_F1a_support.csv", SUPPORT_COLUMNS)
    if (len(support) != len(selected_modes) or
            {int(r["POD_mode"]) for r in support} != set(selected_modes) or
            any(r["valid"] != "1" for r in support)):
        raise ValueError("support audit incomplete")
    for row in support:
        for column in SUPPORT_COLUMNS[1:-1]:
            value = finite(row[column], column)
            if value < 0.0 or value > 1.0 + 1.0e-12:
                raise ValueError("support fraction outside [0,1]")
    linearity_rows = read_csv(root / "strict_ala_stage_F1a_linearity.csv",
                              LINEARITY_COLUMNS)
    linearity = validate_linearity(linearity_rows)
    result = classify_telescoping(rows)
    if result["dominant_damage_stage"] == "RAW_LOCAL" and \
            result["primary_defect_class"] == "UNRESOLVED":
        mode_weights = {int(row["POD_mode"]): float(row["POD_energy_weight"])
                        for row in rows if row["stage"] == "BPI"}
        weight_total = sum(mode_weights.values())
        solve_modes = {int(row["POD_mode"]) for row in local
            if finite(row["local_solve_relative_residual"], "local solve") > LOCAL_SOLVE}
        operator_modes = {int(row["POD_mode"]) for row in local
            if finite(row["local_action_error"], "local action") >= LOCAL_ACTION or
               finite(row["local_action_cosine"], "local cosine") <= 0.0}
        nonlocal_failure_modes = {nonlocal_modes[row["source_vector_or_bin"]]
            for row in nonlocal_rows
            if finite(row["outside_patch_plus_one_fraction"], "nonlocality") >= NONLOCAL}
        if not (solve_modes | operator_modes | nonlocal_failure_modes) <= set(mode_weights):
            raise ValueError("secondary evidence references an unknown POD mode")
        solve_failure = sum(mode_weights[mode] for mode in solve_modes)
        operator_failure = sum(mode_weights[mode] for mode in operator_modes)
        nonlocal_failure = sum(mode_weights[mode] for mode in nonlocal_failure_modes)
        solve_fraction = solve_failure / max(weight_total, SCALE_FLOOR)
        operator_fraction = operator_failure / max(weight_total, SCALE_FLOOR)
        nonlocal_fraction = nonlocal_failure / max(weight_total, SCALE_FLOOR)
        if solve_fraction >= CONSISTENCY:
            result["primary_defect_class"] = "LOCAL_SOLVE_FAILURE"
            result["mode_consistency_energy_fraction"] = solve_fraction
        elif nonlocal_fraction >= CONSISTENCY:
            result["primary_defect_class"] = "NONLOCAL_SCHUR_MISMATCH"
            result["mode_consistency_energy_fraction"] = nonlocal_fraction
        elif operator_fraction >= CONSISTENCY:
            result["primary_defect_class"] = "LOCAL_OPERATOR_MISMATCH"
            result["mode_consistency_energy_fraction"] = operator_fraction
        result["raw_local_secondary_evidence"] = {
            "local_solve_failure_energy_fraction": solve_fraction,
            "local_operator_mismatch_energy_fraction": operator_fraction,
            "nonlocality_energy_fraction": nonlocal_fraction,
            "patch_topology_rule_satisfied": False}
    result.update({"schema": "strict-ala-stage-F1a-decision-v1",
                   "experiment_valid": True,
                   "thresholds_file": "strict_ala_stage_F1a_thresholds.json",
                   "E2_reanalysis": {"next_authorized_path": "LOCAL_SCHWARZ_PATH",
                       "forensic_path_authorized": True,
                       "production_schwarz_modification_authorized": False},
                   "operator_replay_identity": replay,
                   "normalized_linearity_guard": {
                       "status": linearity["status"],
                       "relative_defect": float(linearity["relative_defect"]),
                       "warn_limit": float(linearity["warn_limit"]),
                       "hard_limit": float(linearity["hard_limit"]),
                       "solve_floor": float(linearity["solve_floor"])},
                   "selected_POD_modes": result.pop("modes"),
                   "selected_POD_energy": result.pop("selected_energy"),
                   "raw_local": {"local_solve_quality": "MEASURED",
                                 "local_operator_quality": "MEASURED",
                                 "nonlocality": "MEASURED",
                                 "patch_topology": "GUARD_ONLY"},
                   "F1b_authorized": result["primary_defect_class"] in PRIMARY_CLASSES
                       and result["primary_explained_harm_fraction"] >= PRIMARY
                       and result["mode_consistency_energy_fraction"] >= CONSISTENCY,
                   "automatic_solver_change_authorized": False,
                   "production_schwarz_modification_authorized": False})
    result["next_authorized_task"] = ("F1B_LOCAL_CANDIDATE" if result["F1b_authorized"]
                                      else "MORE_F1A_DIAGNOSTICS")
    write_json(output / "strict_ala_stage_F1a_decision.json", result)
    return result


def audit(args):
    decision = json.loads(Path(args.decision).read_text())
    manifest = json.loads(Path(args.manifest).read_text())
    hashes_pass = all(Path(a).read_bytes() == Path(b).read_bytes()
                      for a, b in ((args.binary_pre, args.binary_post),
                                   (args.input_pre, args.input_post),
                                   (args.source_pre, args.source_post)))
    snapshots_pass = Path(args.snapshot_pre).read_bytes() == Path(args.snapshot_post).read_bytes()
    model_log = Path(args.model_log).read_text(errors="replace")
    marker_pass = model_log.count("STRICT_ALA_STAGE_F1A_COMPLETE") == 1
    lower_log = model_log.lower()
    logs_clean = not any((re.search(r"(?<![a-z])(?:nan|inf)(?![a-z])", lower_log),
                          "fatal error" in lower_log, "mpi_abort" in lower_log,
                          "traceback" in lower_log,
                          re.search(r"(?:fallback_blocks|velocity_block_fallbacks)=[1-9][0-9]*",
                                    lower_log)))
    provenance_pass = bool(manifest.get("provenance_complete") and
        manifest.get("source", {}).get("branch") == "cmbhf_ALA_strict" and
        manifest.get("source", {}).get("branch_verified") and
        not manifest.get("source", {}).get("scientific_dirty"))
    e2_audit = json.loads(Path(args.e2_audit).read_text())
    e2_pass = bool(e2_audit.get("experiment_valid") is True and
                   e2_audit.get("numerical_validation_pass") is True and
                   e2_audit.get("observation_only_trajectory_pass") is True and
                   e2_audit.get("next_authorized_path") == "LOCAL_SCHWARZ_PATH")
    authorization_closed = bool(
        decision.get("automatic_solver_change_authorized") is False and
        decision.get("production_schwarz_modification_authorized") is False and
        decision.get("next_authorized_task") in
            ("F1B_LOCAL_CANDIDATE", "MORE_F1A_DIAGNOSTICS"))
    valid = bool(decision.get("experiment_valid") and hashes_pass and snapshots_pass and
                 marker_pass and logs_clean and e2_pass
                 and provenance_pass and decision.get("automatic_solver_change_authorized") is False
                 and authorization_closed)
    result = dict(decision)
    result.update({"schema": "strict-ala-stage-F1a-final-audit-v1",
                   "experiment_valid": valid, "hashes_unchanged": hashes_pass,
                   "E2_snapshots_unchanged": snapshots_pass,
                   "valid_E2_authorization_pass": e2_pass,
                   "case_log_clean": logs_clean,
                   "authorization_closed": authorization_closed,
                   "completion_marker_pass": marker_pass,
                   "provenance_pass": provenance_pass})
    write_json(args.output, result)
    return 0 if valid else 1


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    freeze = sub.add_parser("freeze-contract"); freeze.add_argument("--output-dir", required=True)
    run = sub.add_parser("analyze"); run.add_argument("--input-dir", required=True); run.add_argument("--output-dir", required=True)
    run.add_argument("--e2-true-mode", required=True); run.add_argument("--e2-reanalysis", required=True)
    cfg = sub.add_parser("check-cfg"); cfg.add_argument("--c0", required=True)
    cfg.add_argument("--c1", required=True); cfg.add_argument("--diff", required=True)
    check = sub.add_parser("audit")
    check.add_argument("--decision", required=True); check.add_argument("--manifest", required=True)
    check.add_argument("--binary-pre", required=True); check.add_argument("--binary-post", required=True)
    check.add_argument("--input-pre", required=True); check.add_argument("--input-post", required=True)
    check.add_argument("--source-pre", required=True); check.add_argument("--source-post", required=True)
    check.add_argument("--snapshot-pre", required=True); check.add_argument("--snapshot-post", required=True)
    check.add_argument("--e2-audit", required=True)
    check.add_argument("--model-log", required=True); check.add_argument("--output", required=True)
    args = parser.parse_args()
    if args.command == "freeze-contract": freeze_contract(args.output_dir)
    elif args.command == "analyze":
        analyze(args.input_dir, args.output_dir, args.e2_true_mode, args.e2_reanalysis)
    elif args.command == "check-cfg":
        check_cfg_contract(args.c0, args.c1, args.diff)
    else:
        raise SystemExit(audit(args))


if __name__ == "__main__":
    main()
