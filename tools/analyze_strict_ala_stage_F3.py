#!/usr/bin/env python3
"""Fail-closed analysis for Strict-ALA F3 outer-Krylov relocalization."""

import argparse
import csv
import hashlib
import json
import math
import struct
from pathlib import Path

import numpy as np

STAGE_ID = "STRICT_STAGE_F3_OUTER_SCHUR_KRYLOV_RELOCALIZATION"

def load_json(path):
    return json.loads(Path(path).read_text())


def load_csv(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True,
                                     allow_nan=False) + "\n")


def json_safe(value):
    """Replace raw nonfinite diagnostics by JSON null after flagging invalid."""
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_safe(item) for item in value]
    return value


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def finite(value):
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def median(values):
    values = sorted(values)
    n = len(values)
    if not n:
        return None
    return values[n // 2] if n % 2 else .5 * (values[n // 2 - 1] + values[n // 2])


def relative_difference(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def json_contains_nonfinite(value):
    if isinstance(value, float):
        return not math.isfinite(value)
    if isinstance(value, dict):
        return any(json_contains_nonfinite(item) for item in value.values())
    if isinstance(value, list):
        return any(json_contains_nonfinite(item) for item in value)
    return False


def csv_nonfinite(rows, excluded=()):
    excluded = set(excluded)
    for row in rows:
        for key, value in row.items():
            if key in excluded or value in (None, ""):
                continue
            try:
                number = float(value)
            except ValueError:
                continue
            if not math.isfinite(number):
                return True
    return False


def validate_thresholds(t):
    # Numeric values live only in the canonical JSON.  This routine validates
    # shape and semantics; it deliberately carries no duplicate gate values.
    required = ("S_ref_inner_relative_tolerance", "S_ref_max_MG_cycles",
        "repeat_action_relative_tolerance", "scaling_scalar_c",
        "scaling_relative_tolerance", "additivity_relative_tolerance",
        "modal_significance_floor_relative",
        "qr_first_pass_orthogonality_threshold",
        "qr_rank_rejection_relative_threshold",
        "projected_subspace_mode_leakage_limit",
        "projected_subspace_weighted_leakage_limit",
        "Y_mode_reachability_error_limit", "Y_weighted_reachability_error_limit",
        "POD_late_median_explanation_min", "POD_iteration40_explanation_min",
        "mode_stagnation_ratio_min", "mode_decay_ratio_max",
        "authoritative_mode_energy_min", "projected_condition_number_adverse",
        "projected_mixing_adverse", "projected_non_normality_adverse",
        "sigma_zero_relative_threshold", "H_true_symmetry_relative_tolerance",
        "H_true_negative_curvature_relative_tolerance",
        "pod_reconstruction_coefficient_relative_tolerance",
        "projection_monotonicity_absolute_tolerance",
        "pod_weight_sum_absolute_tolerance", "inner_target_relative_roundoff")
    errors = []
    for key in required:
        if key not in t or isinstance(t.get(key), bool) or not finite(t.get(key)):
            errors.append(f"{key} must be a finite numeric value")
    pairs=t.get("additivity_mode_pairs")
    if not isinstance(pairs,list) or not pairs or any(not isinstance(x,list) or len(x)!=2 for x in pairs):
        errors.append("additivity_mode_pairs must contain mode-id pairs")
    if not isinstance(t.get("qr_policy"),str) or not t["qr_policy"]:
        errors.append("qr_policy must be nonempty")
    if t.get("production_default_change_authorized") is not False:
        errors.append("production authorization must remain false")
    if errors:
        raise ValueError("; ".join(errors))


def validate_threshold_identity(canonical, canonical_hash, mirror, runtime):
    validate_thresholds(canonical)
    errors=[]
    if mirror.get("values") != canonical:
        errors.append("generated runtime threshold values differ from canonical JSON")
    if mirror.get("canonical_threshold_sha256") != canonical_hash:
        errors.append("generated runtime threshold canonical hash mismatch")
    if runtime.get("canonical_threshold_sha256") != canonical_hash:
        errors.append("runtime canonical threshold hash mismatch")
    loaded=runtime.get("runtime_threshold_values",{})
    for key,value in loaded.items():
        if key not in canonical or canonical[key]!=value:
            errors.append("runtime-loaded threshold mismatch: "+key)
    if not loaded:
        errors.append("runtime-loaded threshold values missing")
    return errors


def sref_residual_record(rhs_norm2,residual_norm2,global_neq,tolerance,cycles,max_cycles,
                         converged=True,fallback=False):
    if global_neq<=0 or rhs_norm2<0 or residual_norm2<0:
        raise ValueError("invalid S_ref residual inputs")
    rhs_rms=math.sqrt(rhs_norm2/global_neq)
    residual_rms=math.sqrt(residual_norm2/global_neq)
    zero=rhs_norm2==0.0
    if zero and residual_norm2==0.0:
        relative,status=0.0,"ZERO_RHS_EXACT"
    elif zero:
        relative,status=None,"ZERO_RHS_NONZERO_RESIDUAL"
    else:
        relative,status=residual_rms/rhs_rms,"NONZERO_RHS"
    passed=(status!="ZERO_RHS_NONZERO_RESIDUAL" and relative is not None and
            relative<=tolerance and cycles<=max_cycles and converged and not fallback and
            all(math.isfinite(x) for x in (rhs_rms,residual_rms,relative)))
    return {"global_rhs_norm2":rhs_norm2,"global_residual_norm2":residual_norm2,
            "rhs_rms":rhs_rms,"residual_rms":residual_rms,"zero_rhs":zero,
            "relative_residual":relative,"residual_status":status,
            "requested_relative_tolerance":tolerance,"achieved_MG_cycles":cycles,
            "max_MG_cycles":max_cycles,"converged":bool(converged),
            "fallback":bool(fallback),"finite":all(math.isfinite(x) for x in
            (rhs_norm2,residual_norm2,rhs_rms,residual_rms)),"residual_gate_pass":passed}


def ordered_viscosity_checksum(records,actual_world_size):
    ordered=sorted(records,key=lambda x:int(x["rank_id"]))
    if [int(x["rank_id"]) for x in ordered] != list(range(actual_world_size)):
        raise ValueError("viscosity rank records do not match actual MPI world")
    payload=b"F3-VISCOSITY-RANK-RECORDS-V1\n"
    for row in ordered:
        payload += ("%d,%d,%s\n"%(int(row["rank_id"]),int(row["local_entry_count"]),
                    row["local_viscosity_hash"])).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def resolve_relative(stage_root,relative_path):
    relative=Path(relative_path)
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError("artifact path is not root-relative")
    return Path(stage_root)/relative


def validate_artifact(path,spec):
    path=Path(path)
    if not path.is_file() or path.stat().st_size==0:
        return False,"missing_or_empty"
    try:
        if spec["artifact_type"]=="json":
            value=load_json(path)
            missing=[key for key in spec.get("required_fields",{}) if key not in value]
            if missing: return False,"missing_json_fields:"+",".join(missing)
        elif spec["artifact_type"]=="csv":
            with path.open(newline="") as stream:
                reader=csv.reader(stream); header=next(reader,[]); data=list(reader)
            if header!=spec.get("exact_ordered_header"):
                return False,"csv_header_mismatch"
            if len(data)<int(spec.get("minimum_record_count",0)):
                return False,"csv_record_count"
    except (OSError,ValueError,json.JSONDecodeError):
        return False,"parse_failure"
    return True,"complete"


def matrix_from_rows(rows, name):
    selected = [r for r in rows if r["matrix"] == name]
    if not selected:
        raise ValueError(f"missing projected matrix {name}")
    n = max(max(int(r["row_mode"]), int(r["column_mode"])) for r in selected)
    matrix = np.full((n, n), np.nan)
    for row in selected:
        matrix[int(row["row_mode"]) - 1, int(row["column_mode"]) - 1] = float(row["value"])
    return matrix


def projected_metrics(matrix, zero_relative):
    singular = np.linalg.svd(matrix, compute_uv=False)
    sigma_max = float(singular[0])
    sigma_min = float(singular[-1])
    infinite = sigma_max == 0.0 or sigma_min <= zero_relative * sigma_max
    kappa = None if infinite else sigma_max / sigma_min
    norm = float(np.linalg.norm(matrix, "fro"))
    off = matrix - np.diag(np.diag(matrix))
    mixing = float(np.linalg.norm(off, "fro") / max(norm, 1e-300))
    departure = matrix.T @ matrix - matrix @ matrix.T
    non_normal = float(np.linalg.norm(departure, "fro") / max(norm * norm, 1e-300))
    symmetric_defect = float(np.linalg.norm(matrix - matrix.T, "fro") / max(norm, 1e-300))
    eigenvalues = np.linalg.eigvals(matrix)
    return {
        "singular_values": [float(x) for x in singular],
        "sigma_max": sigma_max,
        "sigma_min": sigma_min,
        "kappa_2": kappa,
        "kappa_2_is_infinite": infinite,
        "kappa_2_reason": ("SIGMA_MAX_ZERO" if sigma_max == 0.0 else
                           "SIGMA_MIN_BELOW_RELATIVE_ZERO_THRESHOLD") if infinite else "FINITE",
        "eta_mix": mixing,
        "eta_non_normal": non_normal,
        "symmetry_relative_defect": symmetric_defect,
        "eigenvalues": [{"real": float(x.real), "imag": float(x.imag)}
                        for x in eigenvalues],
        "right_singular_vectors": np.linalg.svd(matrix)[2].T.tolist(),
    }


def cfg_contract(args):
    thresholds = load_json(args.thresholds)
    validate_thresholds(thresholds)
    text = Path(args.cfg).read_text()
    required = {
        "ala_inner_accuracy_max": "1e-2",
        "ala_outer_solver": "fgmres",
        "ala_pcg_restart_interval": "50",
        "piterations": "60",
    }
    observed = {}
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if "=" in line:
            key, value = [x.strip() for x in line.split("=", 1)]
            if key in required:
                observed[key] = value
    errors = [f"{key}={observed.get(key)!r}, expected {value!r}"
              for key, value in required.items() if observed.get(key) != value]
    result = {"schema": "strict-ala-stage-F3-local-contract-v1",
              "valid": not errors, "errors": errors,
              "frozen_values": thresholds,
              "observed_cfg": observed,
              "production_default_change_authorized": False}
    write_json(args.output, result)
    return 0 if not errors else 1


def runtime_provenance(args):
    manifest=load_json(args.manifest); expected=Path(args.manifest_hash).read_text().split()[0]
    actual=sha256(args.manifest)
    if actual!=expected: raise ValueError("manifest changed before runtime provenance")
    raw=load_json(args.runtime_raw); pod=load_json(args.pod_reconstruction)
    mirror=load_json(args.threshold_contract); canonical=load_json(args.thresholds)
    ranks=load_csv(args.viscosity_rank_records)
    world=int(raw["actual_mpi_world_size"])
    checksum=ordered_viscosity_checksum(ranks,world)
    viscosity=dict(raw["viscosity_provenance"])
    viscosity.update({"checksum_algorithm":"SHA256(F3-VISCOSITY-RANK-RECORDS-V1 ordered world-rank records)",
      "ordered_rank_record_count":len(ranks),"distributed_viscosity_checksum":checksum,
      "actual_mpi_world_size":world})
    errors=validate_threshold_identity(canonical,sha256(args.thresholds),mirror,raw)
    errors += ([] if raw.get("manifest_sha256")==expected else ["runtime manifest hash mismatch"])
    errors += ([] if pod.get("Q_E2_validation_complete") is True else ["Q_E2 validation incomplete"])
    value={"schema_id":"strict-ala-stage-F3-runtime-provenance","schema_version":"1",
      "stage_id":STAGE_ID,"manifest_sha256":expected,
      "actual_mpi_world_size":world,"actual_restart":raw.get("actual_restart"),
      "actual_outer_budget":raw.get("actual_outer_budget"),
      "canonical_threshold_sha256":sha256(args.thresholds),
      "generated_runtime_threshold_contract_sha256":sha256(args.threshold_contract),
      "runtime_threshold_values":canonical,
      "reconstructed_mode_checksums":pod.get("mode_checksums"),
      "Q_E2_validation_complete":pod.get("Q_E2_validation_complete"),
      "viscosity_provenance":viscosity,"runtime_contract_errors":errors,
      "production_default_change_authorized":False,"complete":not errors}
    write_json(args.output,value)
    return 0 if not errors else 1


def evidence_index(args):
    root=Path(args.stage_root); schema=load_json(args.artifact_schema)
    manifest_hash=Path(args.manifest_hash).read_text().split()[0]
    entries=[]; errors=[]
    for spec in schema["artifacts"]:
        rel=spec["relative_path"]
        if rel.endswith("strict_ala_stage_F3_evidence_index.json") or rel.endswith("strict_ala_stage_F3_final_audit.json"):
            continue
        path=resolve_relative(root,rel); ok,status=validate_artifact(path,spec)
        if not ok: errors.append(rel+":"+status)
        entries.append({"relative_path":rel,"byte_size":path.stat().st_size if path.is_file() else 0,
          "sha256":sha256(path) if path.is_file() else None,"schema_id":spec["schema_id"],
          "schema_version":spec["schema_version"],"complete_status":status})
    value={"schema_id":"strict-ala-stage-F3-evidence-index","schema_version":"1",
      "stage_id":STAGE_ID,"manifest_sha256":manifest_hash,
      "runtime_provenance_sha256":sha256(args.runtime_provenance),
      "threshold_canonical_sha256":sha256(args.thresholds),
      "E2_mode_contract_sha256":sha256(args.mode_contract),"artifacts":entries,
      "errors":errors,"complete":not errors}
    write_json(args.output,value)
    return 0 if not errors else 1


def analyze(args):
    t = load_json(args.thresholds)
    validate_thresholds(t)
    runtime = load_json(args.runtime)
    pod = load_json(args.pod_reconstruction)
    reference = load_json(args.reference_validation_raw)
    completion = load_json(args.completion)
    explain = load_csv(args.pod_explanation)
    modes = load_csv(args.mode_coefficients)
    subspace = load_csv(args.cumulative_subspace)
    rank_rows = load_csv(args.basis_rank)
    inner_rows = load_csv(args.inner_solves)
    matrix_rows = load_csv(args.projected_matrices)
    reconstructed = load_csv(args.reconstructed_coefficients)
    authoritative = load_csv(args.authoritative_coefficients)
    baseline = load_csv(args.baseline_iterations)
    errors = []
    manifest_hash = (Path(args.manifest_hash).read_text().split()[0]
                     if getattr(args,"manifest_hash",None) else "UNIT_TEST_NO_MANIFEST")
    mode_contract = (load_json(args.mode_contract)
                     if getattr(args,"mode_contract",None) else None)
    runtime_provenance_value = (load_json(args.runtime_provenance)
                     if getattr(args,"runtime_provenance",None) else None)

    numeric_tables = (explain, modes, subspace, rank_rows, inner_rows,
                      matrix_rows, reconstructed, authoritative, baseline)
    if any(csv_nonfinite(rows) for rows in numeric_tables):
        errors.append("raw CSV diagnostic contains NaN/Inf")
    if json_contains_nonfinite(runtime) or json_contains_nonfinite(pod) or \
            json_contains_nonfinite(reference) or json_contains_nonfinite(completion):
        errors.append("raw JSON diagnostic contains NaN/Inf")

    if not completion.get("complete") or not completion.get("valid") or completion.get("exit_status") != 0:
        errors.append("numerical completion invalid")
    if runtime.get("outer_preconditioning_orientation") != "RIGHT_FGMRES":
        errors.append("outer orientation mismatch")
    iterations = runtime.get("iterations")
    early_joint = (isinstance(iterations, int) and iterations < 60 and
                   runtime.get("joint_target_reached") is True)
    if (not isinstance(iterations, int) or iterations < 1 or iterations > 60 or
            (iterations < 60 and not early_joint)):
        errors.append("neither frozen 60-iteration trajectory nor genuine early joint convergence")
    if not pod.get("valid") or pod.get("selected_modes", 0) < 1:
        errors.append("POD reconstruction invalid")
    if mode_contract:
        if pod.get("selected_modes") != mode_contract.get("authoritative_mode_count"):
            errors.append("authoritative E2 mode count mismatch")
        if pod.get("mode_ids") != mode_contract.get("authoritative_mode_ids"):
            errors.append("authoritative E2 mode ID/order mismatch")
        observed_weights=[x.get("weight") for x in pod.get("mode_checksums",[])]
        expected_weights=mode_contract.get("authoritative_energy_weights")
        if observed_weights != expected_weights:
            errors.append("authoritative E2 energy-weight mismatch")
        if pod.get("authoritative_mode_set_reused") is not True or pod.get("pod_reselection_performed") is not False:
            errors.append("Q_E2 selection contract mismatch")
    if runtime_provenance_value and (not runtime_provenance_value.get("complete") or
            runtime_provenance_value.get("manifest_sha256")!=manifest_hash):
        errors.append("runtime provenance invalid")
    inner_targets_pass = bool(inner_rows) and all(
        row.get("status") == "CONVERGED" and
        finite(row.get("achieved_relative")) and
        finite(row.get("requested_relative_tolerance")) and
        float(row["achieved_relative"]) <=
        float(row["requested_relative_tolerance"]) *
        (1.0 + t["inner_target_relative_roundoff"])
        for row in inner_rows)
    if not inner_targets_pass:
        errors.append("required inner solve target not achieved")
    if abs(float(pod.get("pod_weight_sum", 0.0)) -
           sum(float(r["mode_energy_weight"]) for r in modes
               if int(r["iteration"]) == 0)) > t["pod_weight_sum_absolute_tolerance"]:
        errors.append("POD weight sum inconsistency")

    auth = {(r["case"], int(r["iteration"]), int(r["pod_mode"])):
            float(r["coefficient"]) for r in authoritative}
    reconstruction_max = 0.0
    for row in reconstructed:
        key = (row["case"], int(row["iteration"]), int(row["pod_mode"]))
        if key not in auth:
            errors.append(f"missing authoritative POD coefficient {key}")
            continue
        reconstruction_max = max(reconstruction_max,
            relative_difference(float(row["coefficient"]), auth[key]))
    if len(reconstructed) != len(auth):
        errors.append("POD coefficient trajectory length mismatch")
    if reconstruction_max > t["pod_reconstruction_coefficient_relative_tolerance"]:
        errors.append("POD coefficient reconstruction mismatch")

    # The established reproduction envelope is exactly zero for these fields.
    base = {int(r["iteration"]): r for r in baseline}
    reproduction_max = 0.0
    for row in explain:
        iteration = int(row["iteration"])
        if iteration == 0:
            continue
        if iteration not in base:
            errors.append(f"baseline missing iteration {iteration}")
            continue
        for field in ("continuity_relative", "momentum_relative"):
            reproduction_max = max(reproduction_max,
                relative_difference(float(row[field]), float(base[iteration][field])))
    reproduction_pass = reproduction_max <= t["base_reproduction_allowed_relative_difference"]
    if not reproduction_pass:
        errors.append("BASE reproduction gate failed")

    for row in rank_rows:
        numeric = ("pre_norm", "post_first_pass_norm", "post_second_pass_norm",
                   "first_pass_orthogonality_defect", "relative_remaining_norm")
        if any(not finite(row[x]) for x in numeric):
            errors.append("nonfinite QR diagnostic")
            break
        if int(row["second_pass_performed"]) != 1:
            errors.append("non-deterministic QR pass count")
            break
        expected_reject=float(row["relative_remaining_norm"]) <= t["qr_rank_rejection_relative_threshold"]
        if expected_reject != (int(row["rejected"])==1):
            errors.append("QR rank rejection does not follow canonical threshold")
            break

    # All-history projection errors must not increase beyond a roundoff guard.
    for mode_id in sorted({int(r["mode_id"]) for r in subspace if r["row_type"] == "MODE"}):
        rows = sorted((r for r in subspace if r["row_type"] == "MODE" and
                       int(r["mode_id"]) == mode_id), key=lambda r: int(r["iteration"]))
        for field in ("E_V_all", "E_Z_all", "E_Y_all"):
            values = [float(r[field]) for r in rows]
            if any(b > a + t["projection_monotonicity_absolute_tolerance"]
                   for a, b in zip(values, values[1:])):
                errors.append(f"nonmonotone {field} mode {mode_id}")

    mode_by_iteration = {}
    for row in modes:
        mode_by_iteration.setdefault(int(row["mode_id"]), {})[int(row["iteration"])] = row
    residual_norm = {int(r["iteration"]): float(r["residual_global_pdot_norm"])
                     for r in explain}
    significant = []
    mode_summary = []
    for mode_id, rows in sorted(mode_by_iteration.items()):
        weight = float(next(iter(rows.values()))["mode_energy_weight"])
        required = set(range(20, 30)) | set(range(31, 41))
        if not required.issubset(rows) or not set(range(20, 30)).issubset(residual_norm):
            pre, late, rpre = None, None, None
            status, ratio = "UNAVAILABLE", None
        else:
            pre = median([float(rows[k]["absolute_modal_amplitude"])
                          for k in range(20, 30)])
            late = median([float(rows[k]["absolute_modal_amplitude"])
                           for k in range(31, 41)])
            rpre = median([residual_norm[k] for k in range(20, 30)])
        if pre is not None and pre < t["modal_significance_floor_relative"] * rpre:
            status, ratio = "BELOW_SIGNIFICANCE_FLOOR", None
        elif pre is not None:
            ratio = late / max(pre, 1e-300)
            if ratio >= t["mode_stagnation_ratio_min"]:
                status = "STAGNATING"
            elif ratio <= t["mode_decay_ratio_max"]:
                status = "DECAYING"
            else:
                status = "WEAK_DECAY"
        if weight >= t["authoritative_mode_energy_min"] and status not in \
                ("BELOW_SIGNIFICANCE_FLOOR", "UNAVAILABLE"):
            significant.append((mode_id, status))
        mode_summary.append({"mode_id": mode_id, "energy_weight": weight,
                             "A_pre": pre, "R_pre": rpre,
                             "rho_mode": ratio, "mode_status": status})

    pod_endpoint_iteration = 40 if iterations >= 40 else max(iterations - 1, 0)
    pod_window = [float(r["f_Q"]) for r in explain
                  if 20 <= int(r["iteration"]) <= pod_endpoint_iteration]
    f_late = median(pod_window)
    endpoint_row = next((r for r in explain
                         if int(r["iteration"]) == pod_endpoint_iteration), None)
    f_endpoint = float(endpoint_row["f_Q"]) if endpoint_row else None
    pod_explains = (f_late is not None and f_endpoint is not None and
                    f_late >= t["POD_late_median_explanation_min"] and
                    f_endpoint >= t["POD_iteration40_explanation_min"])

    reachability_iteration = 40 if iterations >= 40 else max(iterations - 1, 0)
    reach_rows = [r for r in subspace
                  if int(r["iteration"]) == reachability_iteration]
    reach_modes = {int(r["mode_id"]): r for r in reach_rows
                   if r["row_type"] == "MODE"}
    reach_aggregate = next((r for r in reach_rows
                            if r["row_type"] == "AGGREGATE"), None)
    important_ids = [m["mode_id"] for m in mode_summary
                     if m["energy_weight"] >= t["authoritative_mode_energy_min"]]
    y_reachable = bool(reach_aggregate) and all(mid in reach_modes and
        float(reach_modes[mid]["E_Y_all"]) <= t["Y_mode_reachability_error_limit"]
        for mid in important_ids) and \
        float(reach_aggregate["E_Y_weighted_all"]) <= t["Y_weighted_reachability_error_limit"]

    try:
        Htrue = matrix_from_rows(matrix_rows, "H_true")
        Hright = matrix_from_rows(matrix_rows, "H_right")
        if Htrue.shape != Hright.shape or Htrue.shape[0] != pod.get("selected_modes"):
            raise ValueError("projected matrix dimension mismatch")
    except (KeyError, TypeError, ValueError) as error:
        errors.append(str(error))
        Htrue = np.zeros((1, 1))
        Hright = np.zeros((1, 1))
    if not np.isfinite(Htrue).all() or not np.isfinite(Hright).all():
        errors.append("nonfinite projected matrix")
    right_metrics = projected_metrics(Hright, t["sigma_zero_relative_threshold"])
    true_metrics = projected_metrics(Htrue, t["sigma_zero_relative_threshold"])
    try:
        leakages = [float(x) for x in reference["leakage_per_mode"]]
        weighted_leakage = float(reference["leakage_weighted"])
        if (len(leakages) != pod.get("selected_modes") or
                not all(math.isfinite(x) for x in leakages) or
                not math.isfinite(weighted_leakage)):
            raise ValueError("invalid projected-subspace leakage")
    except (KeyError, TypeError, ValueError) as error:
        errors.append(str(error))
        leakages = [math.inf] * int(pod.get("selected_modes", 0))
        weighted_leakage = math.inf
    projected_closed = all(leakages[mid - 1] <= t["projected_subspace_mode_leakage_limit"]
                           for mid in important_ids) and \
        weighted_leakage <= t["projected_subspace_weighted_leakage_limit"]
    if not reference.get("S_ref_validation_pass"):
        errors.append("S_ref validation failed")
    if not reference.get("H_true_validation_pass"):
        errors.append("H_true validation failed")
    if not reference.get("H_right_fixed_operator_valid"):
        errors.append("H_right fixed-operator validation failed")
    if reference.get("raw_vector_action_nonfinite") is not False:
        errors.append("raw vector/action nonfinite or unverified")
    if (true_metrics["symmetry_relative_defect"] >
            t["H_true_symmetry_relative_tolerance"]):
        errors.append("H_true symmetry threshold failed")
    raw_minimum = reference.get("H_true_minimum_symmetric_eigenvalue")
    raw_scale = reference.get("reference_action_norm_scale")
    raw_floor = reference.get("numerical_floor")
    if not all(finite(x) for x in (raw_minimum, raw_scale, raw_floor)):
        errors.append("H_true curvature validation diagnostics nonfinite")
    elif float(raw_minimum) < -t["H_true_negative_curvature_relative_tolerance"] * \
            max(abs(float(raw_scale)), float(raw_floor)):
        errors.append("H_true negative-curvature threshold failed")

    plateau = (reproduction_pass and runtime.get("iterations") == 60 and
               not runtime.get("joint_target_reached") and
               float(runtime["best_R_cont"]) > float(runtime["R_cont_target"]))
    adverse = []
    if not projected_closed:
        adverse.append("POD_SUBSPACE_LEAKAGE")
    elif reference.get("H_right_fixed_operator_valid"):
        if right_metrics["kappa_2_is_infinite"] or \
                right_metrics["kappa_2"] >= t["projected_condition_number_adverse"]:
            adverse.append("ILL_CONDITIONED")
        if right_metrics["eta_mix"] >= t["projected_mixing_adverse"]:
            adverse.append("STRONG_MODE_MIXING")
        if right_metrics["eta_non_normal"] >= t["projected_non_normality_adverse"]:
            adverse.append("STRONGLY_NON_NORMAL")
    slow_important = any(status in ("STAGNATING", "WEAK_DECAY")
                         for _, status in significant)

    if errors:
        decision = "INVALID_EXPERIMENT"
        next_task = "REPAIR_AND_REPEAT_STRICT_STAGE_F3"
    elif plateau and not pod_explains:
        decision = "POD_SUBSPACE_NO_LONGER_EXPLAINS_PLATEAU"
        next_task = "STOP_AND_REVIEW_BEFORE_POD_REFRESH"
    elif not y_reachable and slow_important:
        decision = "OUTER_SCHUR_IMAGE_REACHABILITY_LIMITED"
        next_task = "GLOBAL_OR_DEFLATION_REACHABILITY_DESIGN"
    elif y_reachable and slow_important and reference["H_right_fixed_operator_valid"] and adverse:
        decision = "OUTER_SPACE_REACHABLE_BUT_PROJECTED_OPERATOR_ADVERSE"
        next_task = "STOP_AND_SELECT_PROJECTED_OPERATOR_REMEDY"
    else:
        decision = "OUTER_BOTTLENECK_UNRESOLVED"
        next_task = "STOP_AND_REVIEW_BEFORE_NEXT_SOLVER_DESIGN"

    projected = {
        "schema_id": "strict-ala-stage-F3-projected-operator",
        "schema_version":"1","stage_id":STAGE_ID,"manifest_sha256":manifest_hash,
        "outer_preconditioning_orientation": "RIGHT_FGMRES",
        "H_true_definition": "Q^T S_ref Q",
        "H_right_definition": "Q^T S_ref M^-1 Q",
        "H_true": Htrue.tolist(), "H_right": Hright.tolist(),
        "H_true_metrics": true_metrics, "H_right_metrics": right_metrics,
        "projected_subspace_closed": projected_closed,
        "leakage_per_mode": json_safe(leakages),
        "leakage_weighted": json_safe(weighted_leakage),
        "spectral_interpretation_authoritative": bool(projected_closed and
            reference["H_right_fixed_operator_valid"]),
        "Y_space_interpretation":"empirical_actual_production_Schur_image_space_not_fixed_operator_spectrum",
        "complete":True,
    }
    write_json(args.projected_operator_output, projected)
    reference_final = json_safe(dict(reference))
    reference_final.update({
        "schema_id": "strict-ala-stage-F3-reference-operator-validation",
        "schema_version":"1","stage_id":STAGE_ID,"manifest_sha256":manifest_hash,
        "H_true_metrics": true_metrics,
        "H_right_metrics": right_metrics,
        "derived_infinite_condition_number_is_valid": True,
        "raw_vector_action_nonfinite_invalid": True,
        "complete":True,
    })
    write_json(args.reference_validation_output, reference_final)
    result = {
        "schema_id": "strict-ala-stage-F3-decision","schema_version":"1",
        "stage_id":STAGE_ID,"manifest_sha256":manifest_hash,
        "decision": decision, "next_authorized_task": next_task,
        "valid": not errors, "complete": True, "validity_errors": errors,
        "base_reproduction_pass": reproduction_pass,
        "base_reproduction_max_relative_difference": reproduction_max,
        "base_plateau_confirmed": plateau,
        "base_plateau_definition": "reproduction_gate_pass && iterations==60 && !joint_target_reached && best_R_cont>target",
        "pod_reconstruction_max_relative_coefficient_error": reconstruction_max,
        "all_required_inner_solves_satisfy_requested_target": inner_targets_pass,
        "pod_explanation_late_median_20_40": f_late,
        "pod_explanation_endpoint_iteration": pod_endpoint_iteration,
        "pod_explanation_endpoint": f_endpoint,
        "shortened_window_due_to_early_joint_convergence": early_joint,
        "pod_subspace_explains_plateau": pod_explains,
        "mode_summary": mode_summary,
        "reachability_evaluation_iteration": reachability_iteration,
        "reachability_iteration_40_substituted": reachability_iteration != 40,
        "Y_REACHABLE_AT_40": y_reachable if reachability_iteration == 40 else None,
        "Y_REACHABLE_AT_EVALUATION_ITERATION": y_reachable,
        "V_Z_descriptive_geometry_only": True,
        "projected_subspace_closed": projected_closed,
        "adverse_causes": adverse,
        "restart_retention_status": "INCONCLUSIVE_SINGLE_BOUNDARY",
        "restart_boundary_observed": runtime.get("restart_boundary_observed", False),
        "recovery_metric_available": False,
        "production_default_change_authorized": False,
        "thresholds_sha256": sha256(args.thresholds),
    }
    write_json(args.output, result)
    return 0


def audit(args):
    decision = load_json(args.decision)
    pairs = [(args.binary_pre, args.binary_post),
             (args.inputs_pre, args.inputs_post),
             (args.source_pre, args.source_post)]
    provenance = all(Path(a).read_bytes() == Path(b).read_bytes() for a, b in pairs)
    manifest_hash=Path(args.manifest_hash).read_text().split()[0]
    manifest_unchanged=sha256(args.manifest)==manifest_hash
    schema=load_json(args.artifact_schema); index=load_json(args.evidence_index)
    runtime=load_json(args.runtime_provenance); manifest=load_json(args.manifest)
    errors=[]
    if not manifest_unchanged: errors.append("manifest hash changed")
    if sha256(args.mode_contract)!=manifest.get("authoritative_E2_mode_contract_hash"):
        errors.append("E2 mode contract hash mismatch")
    if sha256(args.artifact_schema)!=manifest.get("artifact_schema_contract_hash"):
        errors.append("artifact schema hash mismatch")
    if runtime.get("manifest_sha256")!=manifest_hash or not runtime.get("complete"):
        errors.append("runtime provenance lineage/closure failure")
    if (runtime.get("actual_mpi_world_size"),runtime.get("actual_restart"),runtime.get("actual_outer_budget"))!=(384,50,60):
        errors.append("runtime execution contract mismatch")
    if index.get("manifest_sha256")!=manifest_hash or not index.get("complete"):
        errors.append("evidence index lineage/closure failure")
    indexed={entry["relative_path"]:entry for entry in index.get("artifacts",[])}
    for spec in schema.get("artifacts",[]):
        rel=spec["relative_path"]
        if rel.endswith("strict_ala_stage_F3_evidence_index.json") or rel.endswith("strict_ala_stage_F3_final_audit.json"):
            continue
        path=resolve_relative(args.stage_root,rel); ok,status=validate_artifact(path,spec)
        entry=indexed.get(rel)
        if not ok or not entry or entry.get("sha256")!=sha256(path):
            errors.append("artifact closure failed: "+rel+":"+status)
    allowed={"POD_SUBSPACE_NO_LONGER_EXPLAINS_PLATEAU",
      "OUTER_SCHUR_IMAGE_REACHABILITY_LIMITED",
      "OUTER_SPACE_REACHABLE_BUT_PROJECTED_OPERATOR_ADVERSE",
      "OUTER_BOTTLENECK_UNRESOLVED","INVALID_EXPERIMENT"}
    if decision.get("decision") not in allowed:
        errors.append("primary decision missing or invalid")
    if decision.get("production_default_change_authorized") is not False:
        errors.append("production authorization changed")
    value = {
        "schema_id": "strict-ala-stage-F3-final-audit","schema_version":"1",
        "stage_id":STAGE_ID,"manifest_sha256":manifest_hash,
        "valid": bool(decision.get("valid") and provenance and not errors),
        "complete": bool(decision.get("complete") and provenance and not errors),
        "decision": decision.get("decision"),
        "next_authorized_task": decision.get("next_authorized_task"),
        "provenance_unchanged": provenance,
        "manifest_hash_unchanged":manifest_unchanged,
        "evidence_chain_errors":errors,
        "production_default_change_authorized": False,
    }
    write_json(args.output, value)
    return 0 if provenance and not errors else 1


def parser():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="command", required=True)
    q = sub.add_parser("cfg-contract")
    q.add_argument("--cfg", required=True); q.add_argument("--thresholds", required=True)
    q.add_argument("--output", required=True); q.set_defaults(fn=cfg_contract)
    q = sub.add_parser("analyze")
    for name in ("thresholds", "runtime", "pod-reconstruction",
                 "reference-validation-raw", "completion", "pod-explanation",
                 "inner-solves",
                 "mode-coefficients", "cumulative-subspace", "basis-rank",
                 "projected-matrices", "reconstructed-coefficients",
                 "authoritative-coefficients", "baseline-iterations",
                 "projected-operator-output", "reference-validation-output",
                 "output"):
        q.add_argument("--" + name, required=True)
    q.set_defaults(fn=analyze)
    q.add_argument("--manifest-hash")
    q.add_argument("--mode-contract")
    q.add_argument("--runtime-provenance")
    q = sub.add_parser("runtime-provenance")
    for name in ("manifest","manifest-hash","runtime-raw","pod-reconstruction",
                 "thresholds","threshold-contract","viscosity-rank-records","output"):
        q.add_argument("--"+name,required=True)
    q.set_defaults(fn=runtime_provenance)
    q = sub.add_parser("evidence-index")
    for name in ("stage-root","artifact-schema","manifest-hash","runtime-provenance",
                 "thresholds","mode-contract","output"):
        q.add_argument("--"+name,required=True)
    q.set_defaults(fn=evidence_index)
    q = sub.add_parser("audit")
    for name in ("decision", "binary-pre", "binary-post", "inputs-pre",
                 "inputs-post", "source-pre", "source-post", "output"):
        q.add_argument("--" + name, required=True)
    for name in ("stage-root","manifest","manifest-hash","mode-contract",
                 "artifact-schema","runtime-provenance","evidence-index"):
        q.add_argument("--"+name,required=True)
    q.set_defaults(fn=audit)
    return p


if __name__ == "__main__":
    args = parser().parse_args()
    raise SystemExit(args.fn(args))
