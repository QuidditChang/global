#!/usr/bin/env python3
"""Fail-closed inner-accuracy/outer-continuity causality analysis."""

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path


CASES = ("INNER_1E2", "INNER_3E3", "INNER_1E3")
CAPS = {"INNER_1E2": 1e-2, "INNER_3E3": 3e-3, "INNER_1E3": 1e-3}
OPERATIONAL_CFG_KEYS = {"datadir", "datadir_old"}


def load_json(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing or empty {path}")
    return json.loads(path.read_text())


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def finite(value, name):
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"non-finite {name}")
    return result


def read_rows(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing or empty {path}")
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
    if not reader.fieldnames or not rows:
        raise ValueError(f"invalid CSV {path}")
    return tuple(reader.fieldnames), rows


def parse_case_paths(specs):
    result = {}
    for spec in specs:
        case, separator, path = spec.partition("=")
        if not separator or case not in CASES or case in result:
            raise ValueError(f"invalid case path {spec}")
        result[case] = Path(path)
    if set(result) != set(CASES):
        raise ValueError("all three case paths are required")
    return result


def parse_cfg(path):
    values = {}
    section = ""
    for raw in Path(path).read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1].strip()
        elif "=" in line:
            key, value = line.split("=", 1)
            values[(section, key.strip())] = value.strip()
    return values


def cfg_diff(args):
    paths = parse_case_paths(args.cfg)
    values = {case: parse_cfg(path) for case, path in paths.items()}
    all_keys = set().union(*(item.keys() for item in values.values()))
    differences = []
    unauthorized = []
    for section, key in sorted(all_keys):
        by_case = {case: values[case].get((section, key)) for case in CASES}
        if len(set(by_case.values())) > 1:
            row = {"section": section, "key": key, "values": by_case}
            differences.append(row)
            if key not in OPERATIONAL_CFG_KEYS and key != "ala_inner_accuracy_max":
                unauthorized.append(row)
    cap_ok = all(math.isclose(
        finite(values[case].get(("CitcomS.solver.vsolver",
                                 "ala_inner_accuracy_max")), "inner cap"),
        CAPS[case], rel_tol=0.0, abs_tol=1e-15) for case in CASES)
    factor_values = [values[case].get(
        ("CitcomS.solver.vsolver", "ala_inner_accuracy_factor"))
        for case in CASES]
    required_frozen = {
        ("CitcomS", "steps"): "1",
        ("CitcomS.solver.mesher", "nprocx"): "4",
        ("CitcomS.solver.mesher", "nprocy"): "4",
        ("CitcomS.solver.mesher", "nprocz"): "2",
        ("CitcomS.solver.mesher", "nodex"): "129",
        ("CitcomS.solver.mesher", "nodey"): "129",
        ("CitcomS.solver.mesher", "nodez"): "65",
        ("CitcomS.solver.mesher", "levels"): "5",
        ("CitcomS.solver.vsolver", "piterations"): "60",
        ("CitcomS.solver.vsolver", "tole_compressibility"): "1e-02",
        ("CitcomS.solver.vsolver", "ala_augmented_lagrangian_gamma"): "10.0",
        ("CitcomS.solver.vsolver", "ala_element_vanka_smoother"): "on",
        ("CitcomS.solver.vsolver", "ala_element_vanka_damping"): "0.2",
        ("CitcomS.solver.vsolver", "ala_inner_accuracy_factor"): "1e-2",
        ("CitcomS.solver.vsolver", "ala_pcg_restart_interval"): "50",
        ("CitcomS.solver.vsolver", "ala_outer_solver"): "fgmres",
        ("CitcomS.solver.vsolver", "ala_shallow_patch_preconditioner"): "on",
        ("CitcomS.solver.vsolver", "ala_shallow_patch_horizontal_elements"): "6",
        ("CitcomS.solver.vsolver", "ala_shallow_patch_horizontal_stride"): "3",
        ("CitcomS.solver.vsolver", "ala_shallow_patch_mpi_overlap"): "3",
        ("CitcomS.solver.vsolver", "ala_shallow_patch_velocity_solver"):
            "element_vanka",
        ("CitcomS.solver.vsolver", "ala_stage_abc_production_logging"): "on",
        ("CitcomS.solver.vsolver", "ala_stage_e_diagnostic"): "off",
    }
    frozen_mismatches = []
    for (section, key), expected in required_frozen.items():
        observed = {case: values[case].get((section, key)) for case in CASES}
        if any(value != expected for value in observed.values()):
            frozen_mismatches.append({"section": section, "key": key,
                                      "expected": expected,
                                      "observed": observed})
    result = {
        "schema": "strict-ala-stage-F2g-causality-cfg-diff-v1",
        "valid": not unauthorized and not frozen_mismatches and cap_ok and
            len(set(factor_values)) == 1 and factor_values[0] is not None,
        "only_scientific_difference": "ala_inner_accuracy_max",
        "authorized_scientific_differences": [row for row in differences
                                               if row["key"] ==
                                               "ala_inner_accuracy_max"],
        "operational_differences": [row for row in differences
                                    if row["key"] in OPERATIONAL_CFG_KEYS],
        "unauthorized_differences": unauthorized,
        "frozen_control_mismatches": frozen_mismatches,
        "frozen_controls_verified": not frozen_mismatches,
        "ala_inner_accuracy_factor_frozen": factor_values[0]
            if len(set(factor_values)) == 1 else None,
        "case_caps": CAPS,
        "production_default_change_authorized": False,
    }
    write_json(args.output, result)
    return 0 if result["valid"] else 1


def derive_envelope(args):
    _, a = read_rows(args.reference_a)
    _, b = read_rows(args.reference_b)
    fields = ("continuity_relative", "momentum_relative")
    by_a = {int(row["iteration"]): row for row in a}
    by_b = {int(row["iteration"]): row for row in b}
    common = sorted(set(by_a) & set(by_b))
    required = list(range(1, 61))
    deviations = {field: [] for field in fields}
    for iteration in common:
        for field in fields:
            av = finite(by_a[iteration][field], field)
            bv = finite(by_b[iteration][field], field)
            deviations[field].append(abs(av - bv) / max(abs(av), abs(bv), 1e-300))
    maxima = {field: max(values, default=math.inf)
              for field, values in deviations.items()}
    valid = common == required and all(value == 0.0 for value in maxima.values())
    result = {
        "schema": "strict-ala-stage-F2g-causality-reproduction-envelope-v1",
        "valid": valid,
        "basis": "two_independent_observed_CONFIGURED_trajectories",
        "iterations": [1, 60],
        "fields": list(fields),
        "maximum_observed_relative_difference": maxima,
        "frozen_allowed_relative_difference": 0.0,
        "reference_a": {"path": str(Path(args.reference_a).resolve()),
                        "sha256": sha256(args.reference_a)},
        "reference_b": {"path": str(Path(args.reference_b).resolve()),
                        "sha256": sha256(args.reference_b)},
        "production_default_change_authorized": False,
    }
    write_json(args.output, result)
    return 0 if valid else 1


def completion_valid(path):
    value = load_json(path)
    if not (value.get("complete") is True and value.get("valid") is True and
            value.get("exit_status") == 0):
        return False
    for name, digest in value.get("artifacts", {}).items():
        candidate = Path(name)
        if not candidate.is_file() or sha256(candidate) != digest:
            return False
    return bool(value.get("artifacts"))


def median_ratio(numerator, denominator, iterations, field="R_cont"):
    return statistics.median(
        numerator[k][field] / max(denominator[k][field], 1e-300)
        for k in iterations)


def analyze(args):
    thresholds = load_json(args.thresholds)
    envelope = load_json(args.reproduction_envelope)
    cfg = load_json(args.cfg_diff)
    iteration_paths = parse_case_paths(args.iterations)
    inner_paths = parse_case_paths(args.inner)
    completion_paths = parse_case_paths(args.completion)
    required_iter_fields = {
        "case", "iteration", "krylov_recursive", "krylov_explicit",
        "continuity_relative", "momentum_relative", "cumulative_inner_solves",
        "cumulative_inner_cycles", "cumulative_K_gamma_applications",
        "cumulative_schur_actions", "elapsed_seconds", "final_iterate"}
    required_inner_fields = {
        "case", "call_id", "outer_iteration", "requested_relative_tolerance",
        "target_absolute", "achieved_absolute", "achieved_relative", "cycles",
        "max_cycles", "seconds", "status"}
    rows = {}
    inner = {}
    case_structure = {}
    for case in CASES:
        fields, raw = read_rows(iteration_paths[case])
        inner_fields, raw_inner = read_rows(inner_paths[case])
        label_ok = all(row["case"] == case for row in raw + raw_inner)
        unique = len({int(row["iteration"]) for row in raw}) == len(raw)
        finite_rows = all(all(math.isfinite(finite(row[field], field)) for field in (
            "krylov_recursive", "krylov_explicit", "continuity_relative",
            "momentum_relative", "cumulative_inner_cycles",
            "cumulative_K_gamma_applications", "elapsed_seconds")) for row in raw)
        rows[case] = {int(row["iteration"]): {
            "R_cont": finite(row["continuity_relative"], "R_cont"),
            "R_mom": finite(row["momentum_relative"], "R_mom"),
            "krylov_recursive_absolute": finite(
                row["krylov_recursive"], "recursive absolute"),
            "krylov_explicit_absolute": finite(
                row["krylov_explicit"], "explicit absolute"),
            "cumulative_inner_solves": int(row["cumulative_inner_solves"]),
            "cumulative_MG_cycles": int(row["cumulative_inner_cycles"]),
            "MG_cycles_per_K_gamma_solve_cumulative":
                int(row["cumulative_inner_cycles"]) /
                max(int(row["cumulative_inner_solves"]), 1),
            "cumulative_K_gamma_applications": int(
                row["cumulative_K_gamma_applications"]),
            "cumulative_schur_actions": int(row["cumulative_schur_actions"]),
            "cumulative_wall_time": finite(row["elapsed_seconds"], "elapsed"),
            "final_iterate": row["final_iterate"] == "1",
        } for row in raw}
        inner[case] = raw_inner
        targets_ok = all(
            row["status"] == "CONVERGED" and
            finite(row["achieved_absolute"], "achieved") <=
                finite(row["target_absolute"], "target") and
            finite(row["requested_relative_tolerance"], "requested") <=
                CAPS[case] * (1.0 + 1e-12) and
            int(row["cycles"]) < int(row["max_cycles"])
            for row in raw_inner)
        case_structure[case] = {
            "completion_integrity": completion_valid(completion_paths[case]),
            "required_iteration_columns": required_iter_fields <= set(fields),
            "required_inner_columns": required_inner_fields <= set(inner_fields),
            "case_labels": label_ok,
            "unique_iterations": unique,
            "finite_iteration_diagnostics": finite_rows,
            "inner_targets_satisfied": targets_ok,
        }

    reference_fields, reference_raw = read_rows(args.reference)
    reference = {int(row["iteration"]): row for row in reference_raw}
    base_reproduced = (envelope.get("valid") is True and
                       envelope.get("frozen_allowed_relative_difference") == 0.0 and
                       all(k in rows["INNER_1E2"] and k in reference and
                           rows["INNER_1E2"][k]["R_cont"] ==
                               finite(reference[k]["continuity_relative"], "ref cont") and
                           rows["INNER_1E2"][k]["R_mom"] ==
                               finite(reference[k]["momentum_relative"], "ref mom")
                           for k in range(1, 61)))
    structural = (thresholds.get("schema") ==
                  "strict-ala-stage-F2g-causality-thresholds-v1" and
                  thresholds.get("frozen_before_hpc_launch") is True and
                  cfg.get("valid") is True and base_reproduced and
                  all(all(gates.values()) for gates in case_structure.values()))

    primary_start, primary_end = thresholds["primary_late_window"]
    common_iterations = sorted(set.intersection(*(set(rows[c]) for c in CASES)))
    late = [k for k in common_iterations if primary_start <= k <= primary_end]
    early_joint = {}
    for case in CASES:
        matches = [k for k in sorted(rows[case])
                   if rows[case][k]["R_cont"] < thresholds["joint_targets"]["R_cont"] and
                   rows[case][k]["R_mom"] <= thresholds["joint_targets"]["R_mom"]]
        k = matches[0] if matches else None
        early_joint[case] = None if k is None else {
            "outer_iterations_to_joint_target": k,
            "MG_cycles_to_joint_target": rows[case][k]["cumulative_MG_cycles"],
            "wall_time_to_joint_target": rows[case][k]["cumulative_wall_time"],
            "time_to_joint_target": rows[case][k]["cumulative_wall_time"],
        }
    if len(late) != primary_end - primary_start + 1:
        converged_before_end = any(value is not None and
                                   value["outer_iterations_to_joint_target"] < primary_end
                                   for value in early_joint.values())
        structural = structural and converged_before_end and bool(common_iterations)
        late = [k for k in common_iterations if k >= primary_start]
        if not late and converged_before_end:
            late = common_iterations

    metrics = {"late_window_used": [min(late), max(late)] if late else None,
               "common_observed_prefix_end": max(common_iterations)
                   if common_iterations else None,
               "joint_convergence": early_joint}
    if late:
        r_3_base = median_ratio(rows["INNER_3E3"], rows["INNER_1E2"], late)
        r_1_3 = median_ratio(rows["INNER_1E3"], rows["INNER_3E3"], late)
        r_1_base = median_ratio(rows["INNER_1E3"], rows["INNER_1E2"], late)
        mom_3_base = median_ratio(rows["INNER_3E3"], rows["INNER_1E2"],
                                  late, "R_mom")
        mom_1_base = median_ratio(rows["INNER_1E3"], rows["INNER_1E2"],
                                  late, "R_mom")
        momentum = {case: statistics.median(rows[case][k]["R_mom"] for k in late)
                    for case in CASES}
        best = {case: min(value["R_cont"] for value in rows[case].values())
                for case in CASES}
        work = {case: max(value["cumulative_MG_cycles"]
                          for value in rows[case].values()) for case in CASES}
        k_gamma_solves = {case: max(value["cumulative_inner_solves"]
                                    for value in rows[case].values())
                          for case in CASES}
        wall_time = {case: max(value["cumulative_wall_time"]
                               for value in rows[case].values())
                     for case in CASES}
        metrics.update({
            "late_median_Rcont_3e3_over_1e2": r_3_base,
            "late_median_Rcont_1e3_over_3e3": r_1_3,
            "late_median_Rcont_1e3_over_1e2": r_1_base,
            "late_median_Rmom_3e3_over_1e2": mom_3_base,
            "late_median_Rmom_1e3_over_1e2": mom_1_base,
            "late_median_R_mom": momentum, "best_R_cont": best,
            "total_MG_cycles": work,
            "total_K_gamma_solves": k_gamma_solves,
            "total_wall_time": wall_time,
        })
    else:
        r_3_base = r_1_3 = r_1_base = math.inf
        mom_3_base = mom_1_base = math.inf
        momentum = best = work = {case: math.inf for case in CASES}
        structural = False

    common_budget = min(work.values()) if structural else 0
    matched_work = {}
    matched_rows = []
    for case in CASES:
        eligible = [(k, value) for k, value in rows[case].items()
                    if value["cumulative_MG_cycles"] <= common_budget]
        if eligible:
            k, value = min(eligible, key=lambda item: item[1]["R_cont"])
            matched_work[case] = {
                "best_achieved_R_cont": value["R_cont"],
                "corresponding_R_mom": value["R_mom"],
                "outer_iteration": k,
                "cumulative_MG_cycles": value["cumulative_MG_cycles"],
                "cumulative_wall_time": value["cumulative_wall_time"],
            }
            matched_rows.append({"case": case, "common_MG_budget": common_budget,
                                 **matched_work[case]})
    metrics["matched_work_actual_states_only"] = matched_work

    control_thresholds = thresholds["classification"]["controls"]
    no_control_thresholds = thresholds["classification"]["does_not_control"]
    control_3 = control_thresholds["max_late_median_Rcont_3e3_over_1e2"]
    control_1 = control_thresholds["max_late_median_Rcont_1e3_over_3e3"]
    momentum_limit = 1.0 + control_thresholds[
        "max_momentum_degradation_relative_to_base"]
    no_control_low, no_control_high = no_control_thresholds[
        "late_median_ratio_interval"]
    best_improvement_floor = 1.0 - no_control_thresholds[
        "best_R_cont_max_improvement_fraction"]
    controls = (structural and r_3_base <= control_3 and r_1_3 <= control_1 and
                best["INNER_1E3"] < best["INNER_3E3"] < best["INNER_1E2"] and
                mom_3_base <= momentum_limit and mom_1_base <= momentum_limit)
    does_not = (structural and no_control_low <= r_3_base <= no_control_high and
                no_control_low <= r_1_base <= no_control_high and
                best["INNER_3E3"] >= best_improvement_floor * best["INNER_1E2"] and
                best["INNER_1E3"] >= best_improvement_floor * best["INNER_1E2"] and
                work["INNER_3E3"] > work["INNER_1E2"] and
                work["INNER_1E3"] > work["INNER_1E2"])
    if not structural:
        decision = "INVALID_EXPERIMENT"
        next_task = "REPAIR_AND_REPEAT_STRICT_STAGE_F2G_CAUSALITY"
    elif controls:
        decision = "INNER_ACCURACY_CONTROLS_OUTER_CONTINUITY"
        next_task = "F2H_CROSS_RANK_FACE_PATCH_FALSIFICATION"
    elif does_not:
        decision = "INNER_ACCURACY_DOES_NOT_CONTROL_OUTER_CONTINUITY"
        next_task = "STOP_F2_LOCAL_BLOCK_SEARCH_AND_RELOCALIZE_OUTER_BOTTLENECK"
    else:
        decision = "INNER_ACCURACY_WEAKLY_COUPLED_TO_OUTER_CONTINUITY"
        next_task = "STOP_AND_REVIEW_BEFORE_F2H"

    with Path(args.aligned_csv).open("w", newline="") as stream:
        fields = ["outer_iteration"] + [f"{case}_{field}" for case in CASES
            for field in ("R_cont", "R_mom", "krylov_recursive_absolute",
                          "krylov_explicit_absolute", "cumulative_inner_solves",
                          "cumulative_MG_cycles",
                          "MG_cycles_per_K_gamma_solve_cumulative",
                          "cumulative_K_gamma_applications",
                          "cumulative_schur_actions", "cumulative_wall_time")]
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for k in common_iterations:
            out = {"outer_iteration": k}
            for case in CASES:
                for field in ("R_cont", "R_mom", "krylov_recursive_absolute",
                              "krylov_explicit_absolute", "cumulative_inner_solves",
                              "cumulative_MG_cycles",
                              "MG_cycles_per_K_gamma_solve_cumulative",
                              "cumulative_K_gamma_applications",
                              "cumulative_schur_actions", "cumulative_wall_time"):
                    out[f"{case}_{field}"] = rows[case][k][field]
            writer.writerow(out)
    with Path(args.matched_work_csv).open("w", newline="") as stream:
        fields = ("case", "common_MG_budget", "best_achieved_R_cont",
                  "corresponding_R_mom", "outer_iteration",
                  "cumulative_MG_cycles", "cumulative_wall_time")
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader(); writer.writerows(matched_rows)
    write_json(args.output, {
        "schema": "strict-ala-stage-F2g-causality-decision-v1",
        "experiment_valid": structural,
        "decision": decision,
        "observed_systematic_improvement": decision ==
            "INNER_ACCURACY_CONTROLS_OUTER_CONTINUITY",
        "independent_repeatability_claimed": False,
        "authoritative_residual_names": {
            "R_cont": "Stage-C continuity_relative",
            "R_mom": "Stage-C momentum_relative",
            "krylov_recursive_absolute": "Stage-C krylov_recursive",
            "krylov_explicit_absolute": "Stage-C krylov_explicit"},
        "gates": {"cfg_diff": cfg.get("valid") is True,
                  "base_reproduction": base_reproduced,
                  "case_structure": case_structure},
        "metrics": metrics,
        "thresholds_sha256": sha256(args.thresholds),
        "reproduction_envelope_sha256": sha256(args.reproduction_envelope),
        "next_authorized_task": next_task,
        "production_default_change_authorized": False,
    })
    return 0


def audit(args):
    decision = load_json(args.decision)
    unchanged = all(Path(getattr(args, f"{name}_pre")).read_bytes() ==
                    Path(getattr(args, f"{name}_post")).read_bytes()
                    for name in ("binary", "inputs", "source"))
    valid = (decision.get("experiment_valid") is True and unchanged and
             decision.get("production_default_change_authorized") is False)
    write_json(args.output, {
        "schema": "strict-ala-stage-F2g-causality-final-audit-v1",
        "experiment_valid": valid,
        "hashes_unchanged": unchanged,
        "decision": decision.get("decision") if valid else "INVALID_EXPERIMENT",
        "next_authorized_task": decision.get("next_authorized_task") if valid
            else "REPAIR_AND_REPEAT_STRICT_STAGE_F2G_CAUSALITY",
        "independent_repeatability_claimed": False,
        "production_default_change_authorized": False,
    })
    return 0


def parser():
    root = argparse.ArgumentParser()
    commands = root.add_subparsers(required=True)
    command = commands.add_parser("cfg-diff")
    command.add_argument("--cfg", action="append", required=True)
    command.add_argument("--output", required=True)
    command.set_defaults(fn=cfg_diff)
    command = commands.add_parser("derive-envelope")
    command.add_argument("--reference-a", required=True)
    command.add_argument("--reference-b", required=True)
    command.add_argument("--output", required=True)
    command.set_defaults(fn=derive_envelope)
    command = commands.add_parser("analyze")
    for name in ("iterations", "inner", "completion"):
        command.add_argument(f"--{name}", action="append", required=True)
    for name in ("thresholds", "reproduction-envelope", "cfg-diff",
                 "reference", "aligned-csv", "matched-work-csv", "output"):
        command.add_argument(f"--{name}", required=True)
    command.set_defaults(fn=analyze)
    command = commands.add_parser("audit")
    command.add_argument("--decision", required=True)
    for name in ("binary-pre", "binary-post", "inputs-pre", "inputs-post",
                 "source-pre", "source-post", "output"):
        command.add_argument(f"--{name}", required=True)
    command.set_defaults(fn=audit)
    return root


def main():
    args = parser().parse_args()
    return args.fn(args)


if __name__ == "__main__":
    raise SystemExit(main())
