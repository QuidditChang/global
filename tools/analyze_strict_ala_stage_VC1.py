#!/usr/bin/env python3
"""Validate and summarize Strict-ALA VC1 viscosity continuation runs."""

import argparse
import csv
import hashlib
import json
import math
import pathlib
import re
import sys


AUTHORITATIVE = {
    "COLD_100": {"status": "NUMERICAL_NONCONVERGENCE", "iterations": 220,
                 "final_R_cont": 1.938367e-2, "final_R_mom": 1.303084e-5},
    "SEED_1": {"status": "CONVERGED", "iterations": 215,
               "final_R_cont": 9.994139e-3, "final_R_mom": 2.486057e-4},
}
ENVELOPE = {"R_cont_abs": 5.0e-7, "R_mom_abs": 5.0e-10,
            "iterations_abs": 0}
PHASE_DELTA_S_COMPAT = (
    -0.032902012641276054,
    -0.022663818742040032,
    0.017520140281939885,
)


def digest(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def parse_cfg(path):
    values, section = {}, ""
    for raw in pathlib.Path(path).read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1]
        elif "=" in line:
            key, value = (item.strip() for item in line.split("=", 1))
            values[(section, key)] = value
    return values


def parse_float_vector(value):
    return tuple(float(item.strip()) for item in value.split(","))


def prepare_runtime_cfg(args):
    """Add the current parser's required phase entropy vector to a step-0 cfg.

    The upstream cfg is retained byte-for-byte.  phase_delta_s is not consumed
    by the frozen step-0 Stokes equations, but the current single-entropy phase
    parser requires a nonzero vector before those equations can start.
    """
    upstream = pathlib.Path(args.upstream)
    output = pathlib.Path(args.output)
    audit_path = pathlib.Path(args.audit)
    failures = []
    upstream_hash = digest(upstream)
    if upstream_hash != args.expected_upstream_sha256:
        failures.append("upstream cfg sha256 mismatch")
    values = parse_cfg(upstream)
    phase_key = ("CitcomS.solver.phase", "phase_delta_s")
    if phase_key in values:
        failures.append("upstream cfg unexpectedly already defines phase_delta_s")

    lines = upstream.read_text().splitlines()
    rewritten = []
    section = ""
    inserted = 0
    vector = ", ".join(format(value, ".17g") for value in PHASE_DELTA_S_COMPAT)
    for raw in lines:
        stripped = raw.split("#", 1)[0].strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            section = stripped[1:-1]
        rewritten.append(raw)
        if (section == "CitcomS.solver.phase" and
                re.match(r"^\s*phase_delta_rho\s*=", raw)):
            rewritten.append("# VC1 step-0 parser compatibility; thermal evolution is disabled.")
            rewritten.append("phase_delta_s = " + vector)
            inserted += 1
    if inserted != 1:
        failures.append("expected one phase_delta_rho insertion point, found %d" % inserted)

    runtime_hash = None
    changed_fields = []
    if not failures:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text("\n".join(rewritten) + "\n")
        runtime_hash = digest(output)
        runtime_values = parse_cfg(output)
        changed_fields = sorted("%s.%s" % key for key in set(values) | set(runtime_values)
                                if values.get(key) != runtime_values.get(key))
        if changed_fields != ["CitcomS.solver.phase.phase_delta_s"]:
            failures.append("runtime cfg changed fields are not phase_delta_s only")

    audit = {
        "schema": "strict-ala-stage-VC1-runtime-cfg-compatibility-v1",
        "valid": not failures,
        "failures": failures,
        "upstream_cfg_sha256": upstream_hash,
        "runtime_cfg_sha256": runtime_hash,
        "changed_fields": changed_fields,
        "phase_delta_s": list(PHASE_DELTA_S_COMPAT),
        "scope": "step_0_frozen_stokes_initialization_only",
        "stokes_operator_change": False,
        "thermal_evolution_authorized": False,
        "production_default_change_authorized": False,
    }
    write_json(audit_path, audit)
    if failures:
        raise SystemExit("; ".join(failures))


def cfg_contract(args):
    base, stage = parse_cfg(args.canonical), parse_cfg(args.cfg)
    allowed = {
        ("CitcomS", "steps"),
        ("CitcomS.solver", "datadir"), ("CitcomS.solver", "datadir_old"),
        ("CitcomS.solver", "datafile"), ("CitcomS.solver", "datafile_old"),
        ("CitcomS.solver", "logfile"),
        ("CitcomS.solver.visc", "visc_max"),
    }
    changed = sorted("%s.%s" % key for key in set(base) | set(stage)
                     if base.get(key) != stage.get(key))
    unexpected = sorted("%s.%s" % key for key in set(base) | set(stage)
                        if base.get(key) != stage.get(key) and key not in allowed)
    required = {
        ("CitcomS", "steps"): "1",
        ("CitcomS.solver.vsolver", "piterations"): "220",
        ("CitcomS.solver.vsolver", "ala_outer_solver"): "coupled_fgmres",
        ("CitcomS.solver.vsolver", "ala_element_vanka_smoother"): "on",
        ("CitcomS.solver.vsolver", "ala_coupled_element_vanka"): "on",
        ("CitcomS.solver.visc", "SDEPV"): "off",
        ("CitcomS.solver.visc", "PDEPV"): "off",
    }
    failures = []
    for key, expected in required.items():
        if stage.get(key) != expected:
            failures.append("%s.%s=%r expected %r" %
                            (key[0], key[1], stage.get(key), expected))
    phase_delta_s = stage.get(("CitcomS.solver.phase", "phase_delta_s"))
    if phase_delta_s is None:
        failures.append("CitcomS.solver.phase.phase_delta_s is required")
    else:
        try:
            parsed_phase_delta_s = parse_float_vector(phase_delta_s)
        except ValueError:
            parsed_phase_delta_s = ()
        if parsed_phase_delta_s != PHASE_DELTA_S_COMPAT:
            failures.append("phase_delta_s compatibility vector mismatch")
    if float(stage.get(("CitcomS.solver.visc", "visc_max"), "nan")) != args.visc_max:
        failures.append("visc_max mismatch")
    failures.extend("unauthorized cfg change: " + item for item in unexpected)
    value = {"schema": "strict-ala-stage-VC1-cfg-contract-v1",
             "valid": not failures, "failures": failures,
             "changed_fields": changed, "authorized_changed_fields":
             sorted("%s.%s" % key for key in allowed),
             "canonical_sha256": digest(args.canonical),
             "stage_cfg_sha256": digest(args.cfg),
             "visc_max": args.visc_max}
    write_json(args.output, value)
    if failures:
        raise SystemExit("; ".join(failures))


def search_float(text, pattern, default=None):
    matches = re.findall(pattern, text, re.MULTILINE)
    return float(matches[-1]) if matches else default


def search_int(text, pattern, default=0):
    matches = re.findall(pattern, text, re.MULTILINE)
    return int(matches[-1]) if matches else default


def parse_time(path):
    if not path or not pathlib.Path(path).is_file():
        return None
    return search_int(pathlib.Path(path).read_text(),
                      r"Maximum resident set size \(kbytes\):\s*(\d+)", None)


def stage_summary(args):
    case = pathlib.Path(args.case_dir)
    log = case / "DATA" / "0" / "raw.log"
    text = log.read_text(errors="replace") if log.is_file() else ""
    iter_files = list((case / "DATA" / "0").glob("*.strict_ala_stage_C_iterations.csv"))
    rows = list(csv.DictReader(iter_files[0].open())) if len(iter_files) == 1 else []
    complete = re.findall(
        r"STRICT_ALA_STAGE_C_CASE_COMPLETE case=\S+ status=(\S+) iterations=(\d+) "
        r"continuity=([0-9eE+.-]+) momentum=([0-9eE+.-]+)", text)
    if not complete:
        coupled = re.findall(
            r"ALA_COUPLED_FEASIBILITY_SUMMARY status=(\S+) iterations=(\d+) "
            r"cancellation=([0-9eE+.-]+).*?raw_momentum_relative=([0-9eE+.-]+)",
            text)
        if coupled:
            coupled_status = ("CONVERGED" if coupled[-1][0] == "joint_target_reached"
                              else "NUMERICAL_NONCONVERGENCE")
            complete = [(coupled_status, coupled[-1][1], coupled[-1][2], coupled[-1][3])]
    status = complete[-1][0] if complete else "INFRASTRUCTURE_FAILURE"
    iterations = int(complete[-1][1]) if complete else 0
    final_cont = float(complete[-1][2]) if complete else None
    final_mom = float(complete[-1][3]) if complete else None
    visc = re.findall(
        r"STRICT_ALA_VC1_VISCOSITY configured_visc_max=([0-9eE+.-]+) "
        r"eta_min=([0-9eE+.-]+) eta_max=([0-9eE+.-]+) eta_ratio=([0-9eE+.-]+) "
        r"lower_clamp_fraction=([0-9eE+.-]+) upper_clamp_fraction=([0-9eE+.-]+) "
        r"sample_count=(\d+) global_xor_checksum=([0-9a-f]+)", text)
    warm = re.findall(r"STRICT_ALA_VC1_WARM_STATE action=read directory=(\S+) "
                      r"global_xor_checksum=([0-9a-f]+) u_checksum=([0-9a-f]+) "
                      r"p_checksum=([0-9a-f]+)", text)
    warm_written = re.findall(r"STRICT_ALA_VC1_WARM_STATE action=write directory=(\S+) "
                              r"global_xor_checksum=([0-9a-f]+) u_checksum=([0-9a-f]+) "
                              r"p_checksum=([0-9a-f]+)", text)
    frozen = re.findall(r"STRICT_ALA_VC1_FROZEN_STATE before=([0-9a-f]+) "
                        r"after=([0-9a-f]+) bitwise_equal=true", text)
    frozen_guard = re.findall(r"STRICT_ALA_VC1_FROZEN_STATE_GUARD before=([0-9a-f]+) ", text)
    timing = re.findall(r"STRICT_ALA_VC1_TIMING operator_rebuild_seconds=([0-9eE+.-]+) "
                        r"fgmres_and_preconditioner_seconds=([0-9eE+.-]+)", text)
    vanka_build = search_float(text, r"ALA full element-Vanka factors .*?build_seconds_max=([0-9eE+.-]+)", 0.0)
    last = rows[-1] if rows else {}
    rcont_traj = [{"iteration": int(row["iteration"]),
                   "R_cont": float(row["continuity_relative"])} for row in rows]
    rmom_traj = [{"iteration": int(row["iteration"]),
                  "R_mom": float(row["momentum_relative"])} for row in rows
                 if row.get("momentum_relative")]
    if not rows:
        startup = re.findall(r"ALA COUPLED FGMRES startup .*?cancellation=([0-9eE+.-]+) "
                             r"raw_momentum_relative=([0-9eE+.-]+)", text)
        coupled_rows = re.findall(r"ALA COUPLED FGMRES iteration=(\d+) .*?"
                                  r"cancellation=([0-9eE+.-]+).*?"
                                  r"raw_momentum_relative_last_audit=([0-9eE+.-]+)", text)
        if startup:
            rcont_traj.append({"iteration": 0, "R_cont": float(startup[-1][0])})
            rmom_traj.append({"iteration": 0, "R_mom": float(startup[-1][1])})
        rcont_traj.extend({"iteration": int(item[0]), "R_cont": float(item[1])}
                          for item in coupled_rows)
        rmom_traj.extend({"iteration": int(item[0]), "R_mom": float(item[2])}
                         for item in coupled_rows)
    coupled_cost = re.findall(r"STRICT_ALA_VC1_COUPLED_COST K_gamma_rhs_solves=(\d+) "
                              r"K_gamma_operator_applications=(\d+) velocity_MG_cycles=(\d+) "
                              r"preconditioner_applications=(\d+)", text)
    k_solves = (int(last.get("cumulative_inner_solves", 0) or 0) if last else
                (int(coupled_cost[-1][0]) if coupled_cost else 0))
    mg_cycles = (int(last.get("cumulative_inner_cycles", 0) or 0) if last else
                 (int(coupled_cost[-1][2]) if coupled_cost else 0))
    frozen_valid = ((len(frozen) >= 1 and frozen[-1][0] == frozen[-1][1]) or
                    (status == "NUMERICAL_NONCONVERGENCE" and len(frozen_guard) >= 1))
    numerical_output_valid = bool(rcont_traj) and (bool(rows) or len(coupled_cost) == 1)
    valid = (status in ("CONVERGED", "NUMERICAL_NONCONVERGENCE") and
             len(visc) >= 1 and frozen_valid and
             numerical_output_valid and
             (args.warm_source == "COLD" or len(warm) == 1) and
             (status != "CONVERGED" or len(warm_written) == 1))
    value = {
        "schema": "strict-ala-stage-VC1-stage-completion-v1",
        "path_id": args.path_id, "continuation_stage_index": args.stage_index,
        "visc_max": args.visc_max, "warm_start_source": args.warm_source,
        "warm_start_u_p_checksum": warm[-1][1] if warm else None,
        "warm_start_u_checksum": warm[-1][2] if warm else None,
        "warm_start_p_checksum": warm[-1][3] if warm else None,
        "output_u_p_checksum": warm_written[-1][1] if warm_written else None,
        "output_u_checksum": warm_written[-1][2] if warm_written else None,
        "output_p_checksum": warm_written[-1][3] if warm_written else None,
        "pressure_gauge": "exact_equation_point_pressure_no_regauge",
        "FGMRES_iterations": iterations, "R_cont_trajectory": rcont_traj,
        "R_mom_trajectory": rmom_traj, "final_R_cont": final_cont,
        "final_R_mom": final_mom,
        "K_gamma_solve_count": k_solves,
        "K_gamma_operator_application_count": (int(last.get("cumulative_K_gamma_applications", 0) or 0) if last else (int(coupled_cost[-1][1]) if coupled_cost else 0)),
        "total_MG_cycles": mg_cycles,
        "MG_cycles_per_K_gamma_solve": (
            float(mg_cycles) / max(float(k_solves), 1.0)),
        "operator_rebuild_time_seconds": (max(float(timing[-1][0]) - vanka_build, 0.0)
                                          if timing else None),
        "preconditioner_cache_rebuild_time_seconds": (
            float(last.get("bpi_construction_seconds", 0) or 0) +
            float(last.get("schwarz_construction_seconds", 0) or 0) +
            vanka_build) if last else vanka_build,
        "FGMRES_solve_time_seconds": (float(last.get("elapsed_seconds", 0) or 0) if last else
                                      (float(timing[-1][1]) if timing and float(timing[-1][1]) > 0 else max(args.wall_seconds - (float(timing[-1][0]) if timing else 0), 0.0))),
        "total_stage_wall_time_seconds": args.wall_seconds,
        "peak_RSS_kbytes": parse_time(case / "case_time.txt"),
        "breakdown_fallback_state": "none" if "breakdown=1" not in text else "breakdown",
        "convergence_status": status, "process_exit_status": args.exit_status,
        "valid": valid, "complete": bool(complete),
        "physical_state_fingerprint": frozen[-1][0] if frozen else (frozen_guard[-1] if frozen_guard else None),
        "viscosity": ({"configured_visc_max": float(visc[0][0]),
                       "eta_min": float(visc[0][1]), "eta_max": float(visc[0][2]),
                       "eta_ratio": float(visc[0][3]),
                       "fraction_at_lower_clamp": float(visc[0][4]),
                       "fraction_at_upper_clamp": float(visc[0][5]),
                       "sample_count": int(visc[0][6]), "checksum": visc[0][7]}
                      if visc else None),
        "artifacts": {str(log.resolve()): digest(log)} if log.is_file() else {},
    }
    write_json(args.output, value)
    if status == "INFRASTRUCTURE_FAILURE":
        raise SystemExit("stage infrastructure failure: path=%s exit_status=%d" %
                         (args.path_id, args.exit_status))


def within(actual, expected, tolerance):
    return actual is not None and abs(actual - expected) <= tolerance


def continuation_path_complete(completions, path_id, viscosity_sequence):
    rows = sorted((row for row in completions if row["path_id"] == path_id),
                  key=lambda row: row["continuation_stage_index"])
    if not rows:
        return False
    for offset, expected_viscosity in enumerate(viscosity_sequence):
        if offset >= len(rows):
            return False
        row = rows[offset]
        if (row["continuation_stage_index"] != offset + 1 or
                row["visc_max"] != float(expected_viscosity)):
            return False
        if row["convergence_status"] != "CONVERGED":
            return len(rows) == offset + 1
    return len(rows) == len(viscosity_sequence)


def aggregate(args):
    root = pathlib.Path(args.root)
    outputs = pathlib.Path(args.output_dir)
    outputs.mkdir(parents=True, exist_ok=True)
    completions = [json.loads(path.read_text()) for path in
                   sorted(root.glob("**/completion.json"))]
    by_key = {(row["path_id"], row["continuation_stage_index"]): row
              for row in completions}
    cold = by_key.get(("COLD_100", 0), {})
    seed = by_key.get(("SEED_1", 0), {})
    baseline_ok = (cold.get("convergence_status") == AUTHORITATIVE["COLD_100"]["status"] and
                   cold.get("FGMRES_iterations") == 220 and
                   within(cold.get("final_R_cont"), AUTHORITATIVE["COLD_100"]["final_R_cont"], ENVELOPE["R_cont_abs"]) and
                   within(cold.get("final_R_mom"), AUTHORITATIVE["COLD_100"]["final_R_mom"], ENVELOPE["R_mom_abs"]))
    seed_ok = (seed.get("convergence_status") == AUTHORITATIVE["SEED_1"]["status"] and
               seed.get("FGMRES_iterations") == 215 and
               within(seed.get("final_R_cont"), AUTHORITATIVE["SEED_1"]["final_R_cont"], ENVELOPE["R_cont_abs"]) and
               within(seed.get("final_R_mom"), AUTHORITATIVE["SEED_1"]["final_R_mom"], ENVELOPE["R_mom_abs"]))
    source_by_name = {"SEED_1": seed}
    for row in completions:
        if row["path_id"] == "COARSE":
            source_by_name["COARSE_V%s" % format(row["visc_max"], "g")] = row
        if row["path_id"] == "FINE":
            source_by_name["FINE_V%s" % format(row["visc_max"], "g")] = row
    warm_chain_valid = True
    for row in completions:
        if row["warm_start_source"] == "COLD":
            continue
        source = source_by_name.get(row["warm_start_source"])
        warm_chain_valid = warm_chain_valid and source is not None and all(
            row.get("warm_start_" + component + "checksum") ==
            source.get("output_" + component + "checksum")
            for component in ("u", "p"))
    physical_fingerprints = {row.get("physical_state_fingerprint") for row in completions}
    physical_lineage_valid = len(physical_fingerprints) == 1 and None not in physical_fingerprints
    cold_and_seed_present = (("COLD_100", 0) in by_key and
                             ("SEED_1", 0) in by_key)
    if seed.get("convergence_status") == "CONVERGED":
        path_execution_complete = all((
            continuation_path_complete(completions, "DIRECT_WARM", (100,)),
            continuation_path_complete(completions, "COARSE", (10, 100)),
            continuation_path_complete(completions, "FINE", (2, 3, 5, 10, 30, 100)),
        ))
    else:
        path_execution_complete = not any(
            row["path_id"] in ("DIRECT_WARM", "COARSE", "FINE")
            for row in completions)
    experiment_complete = (cold_and_seed_present and path_execution_complete and
                           all(row.get("complete") for row in completions))
    all_valid = (all(row.get("valid") for row in completions) and
                 warm_chain_valid and physical_lineage_valid and
                 experiment_complete)
    final_success = {}
    for path in ("DIRECT_WARM", "COARSE", "FINE"):
        stages = [row for row in completions if row["path_id"] == path]
        final_success[path] = bool(stages and max(stages, key=lambda x: x["continuation_stage_index"])["visc_max"] == 100.0 and
                                   max(stages, key=lambda x: x["continuation_stage_index"])["convergence_status"] == "CONVERGED")
    if not (baseline_ok and seed_ok and all_valid):
        decision = "INVALID_EXPERIMENT"
        next_task = None
    elif final_success["DIRECT_WARM"]:
        decision = "DIRECT_WARM_START_RESCUES_TNEW_TARGET"
        next_task = "VC2_TARGET_VISCOSITY_PERSISTENCE_TEST"
    elif final_success["COARSE"]:
        decision = "COARSE_CONTINUATION_RESCUES_TNEW_TARGET"
        next_task = "VC2_TARGET_VISCOSITY_PERSISTENCE_TEST"
    elif final_success["FINE"]:
        decision = "FINE_CONTINUATION_RESCUES_TNEW_TARGET"
        next_task = "VC2_TARGET_VISCOSITY_PERSISTENCE_TEST"
    else:
        decision = "CONTINUATION_DOES_NOT_RESCUE_TNEW_TARGET"
        next_task = "STRICT_STAGE_F3_OUTER_KRYLOV_RELOCALIZATION"
    stage_fields = ["path_id", "continuation_stage_index", "visc_max",
                    "warm_start_source", "FGMRES_iterations", "final_R_cont",
                    "final_R_mom", "K_gamma_solve_count", "total_MG_cycles",
                    "MG_cycles_per_K_gamma_solve", "operator_rebuild_time_seconds",
                    "preconditioner_cache_rebuild_time_seconds", "FGMRES_solve_time_seconds",
                    "total_stage_wall_time_seconds", "peak_RSS_kbytes",
                    "breakdown_fallback_state", "convergence_status", "valid"]
    with (outputs / "strict_ala_stage_VC1_stage_summary.csv").open("w", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=stage_fields); writer.writeheader()
        writer.writerows({key: row.get(key) for key in stage_fields} for row in completions)
    visc_fields = ["path_id", "continuation_stage_index", "configured_visc_max",
                   "eta_min", "eta_max", "eta_ratio", "fraction_at_lower_clamp",
                   "fraction_at_upper_clamp", "sample_count", "checksum"]
    with (outputs / "strict_ala_stage_VC1_viscosity_stats.csv").open("w", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=visc_fields); writer.writeheader()
        for row in completions:
            if row.get("viscosity"):
                writer.writerow({"path_id": row["path_id"],
                                 "continuation_stage_index": row["continuation_stage_index"],
                                 **row["viscosity"]})
    target_fields = ["path_id", "continuation_stage_index", "warm_start_source",
                     "initial_R_cont", "initial_R_mom"]
    with (outputs / "strict_ala_stage_VC1_target100_entry_residuals.csv").open("w", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=target_fields); writer.writeheader()
        for row in completions:
            if row["visc_max"] == 100.0:
                writer.writerow({"path_id": row["path_id"],
                    "continuation_stage_index": row["continuation_stage_index"],
                    "warm_start_source": row["warm_start_source"],
                    "initial_R_cont": row["R_cont_trajectory"][0]["R_cont"] if row["R_cont_trajectory"] else None,
                    "initial_R_mom": row["R_mom_trajectory"][0]["R_mom"] if row["R_mom_trajectory"] else None})
    path_rows = []
    for path in ("COLD_100", "SEED_1", "DIRECT_WARM", "COARSE", "FINE"):
        rows = [row for row in completions if row["path_id"] == path]
        path_rows.append({"path_id": path, "completed_stages": len(rows),
            "total_FGMRES_iterations": sum(row["FGMRES_iterations"] for row in rows),
            "total_K_gamma_solves": sum(row["K_gamma_solve_count"] for row in rows),
            "total_MG_cycles": sum(row["total_MG_cycles"] for row in rows),
            "total_operator_preconditioner_rebuild_seconds": sum((row["operator_rebuild_time_seconds"] or 0) + (row["preconditioner_cache_rebuild_time_seconds"] or 0) for row in rows),
            "total_solver_wall_seconds": sum(row["FGMRES_solve_time_seconds"] or 0 for row in rows),
            "total_continuation_wall_seconds": sum(row["total_stage_wall_time_seconds"] for row in rows),
            "target_100_converged": final_success.get(path, cold.get("convergence_status") == "CONVERGED" if path == "COLD_100" else False),
            "target_100_FGMRES_iterations": next((row["FGMRES_iterations"] for row in rows if row["visc_max"] == 100.0), None),
            "target_100_K_gamma_solves": next((row["K_gamma_solve_count"] for row in rows if row["visc_max"] == 100.0), None),
            "target_100_MG_cycles": next((row["total_MG_cycles"] for row in rows if row["visc_max"] == 100.0), None),
            "target_100_solver_wall_seconds": next((row["FGMRES_solve_time_seconds"] for row in rows if row["visc_max"] == 100.0), None)})
    with (outputs / "strict_ala_stage_VC1_path_summary.csv").open("w", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=list(path_rows[0])); writer.writeheader(); writer.writerows(path_rows)
    decision_value = {"schema": "strict-ala-stage-VC1-decision-v1",
                      "decision": decision, "next_authorized_task": next_task,
                      "baseline_reproduced": baseline_ok, "seed_reproduced": seed_ok,
                      "warm_chain_checksum_valid": warm_chain_valid,
                      "physical_state_lineage_valid": physical_lineage_valid,
                      "production_default_change_authorized": False}
    write_json(outputs / "strict_ala_stage_VC1_decision.json", decision_value)
    provenance_path = root / "00_manifest" / "provenance.json"
    provenance = json.loads(provenance_path.read_text()) if provenance_path.is_file() else None
    manifest = {"schema": "strict-ala-stage-VC1-manifest-v1", "root": str(root),
                "authoritative_evidence": AUTHORITATIVE,
                "reproduction_envelope": ENVELOPE,
                "provenance": provenance,
                "provenance_sha256": digest(provenance_path) if provenance else None,
                "stage_completion_hashes": {str(path): digest(path) for path in sorted(root.glob("**/completion.json"))}}
    write_json(outputs / "strict_ala_stage_VC1_manifest.json", manifest)
    comparison_path = root / "ANALYSIS" / "state_comparisons.log"
    comparisons = []
    if comparison_path.is_file():
        pattern = re.compile(r"left=(\S+) right=(\S+) relative_velocity_difference=([0-9eE+.-]+) "
                             r"pressure_difference=([0-9eE+.-]+) relative_pressure_difference=([0-9eE+.-]+)")
        for line in comparison_path.read_text().splitlines():
            match = pattern.search(line)
            if match:
                comparisons.append({"left": match.group(1), "right": match.group(2),
                    "relative_velocity_difference": float(match.group(3)),
                    "pressure_difference_global_pdot": float(match.group(4)),
                    "relative_pressure_difference_global_pdot": float(match.group(5)),
                    "pressure_gauge": "exact_no_regauge"})
    audit = {"schema": "strict-ala-stage-VC1-final-audit-v1",
             "valid": decision != "INVALID_EXPERIMENT",
             "complete": experiment_complete,
             "decision": decision, "maximum_stage_count": 11,
             "completed_stage_count": len(completions),
             "physical_state_fingerprints": sorted(physical_fingerprints, key=lambda item: "" if item is None else item),
             "physical_state_lineage_valid": physical_lineage_valid,
             "successful_target_state_comparisons": comparisons,
             "production_default_change_authorized": False,
             "exit_status": 1 if decision == "INVALID_EXPERIMENT" else 0}
    write_json(outputs / "strict_ala_stage_VC1_final_audit.json", audit)
    if decision == "INVALID_EXPERIMENT":
        raise SystemExit("VC1 experiment is invalid")


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("prepare-runtime-cfg")
    p.add_argument("--upstream", required=True); p.add_argument("--output", required=True)
    p.add_argument("--audit", required=True)
    p.add_argument("--expected-upstream-sha256", required=True)
    p.set_defaults(func=prepare_runtime_cfg)
    p = sub.add_parser("cfg-contract")
    p.add_argument("--canonical", required=True); p.add_argument("--cfg", required=True)
    p.add_argument("--visc-max", type=float, required=True); p.add_argument("--output", required=True)
    p.set_defaults(func=cfg_contract)
    p = sub.add_parser("stage")
    p.add_argument("--path-id", required=True); p.add_argument("--stage-index", type=int, required=True)
    p.add_argument("--visc-max", type=float, required=True); p.add_argument("--warm-source", required=True)
    p.add_argument("--case-dir", required=True); p.add_argument("--exit-status", type=int, required=True)
    p.add_argument("--wall-seconds", type=float, required=True); p.add_argument("--output", required=True)
    p.set_defaults(func=stage_summary)
    p = sub.add_parser("aggregate")
    p.add_argument("--root", required=True); p.add_argument("--output-dir", required=True)
    p.set_defaults(func=aggregate)
    args = parser.parse_args(); args.func(args)


if __name__ == "__main__":
    main()
