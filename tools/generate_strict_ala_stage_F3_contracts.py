#!/usr/bin/env python3
"""Generate deterministic pre-run contracts for Strict-ALA Stage F3."""
import argparse, csv, hashlib, json, os
from pathlib import Path

STAGE_ID = "STRICT_STAGE_F3_OUTER_SCHUR_KRYLOV_RELOCALIZATION"

def sha256(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def rows(path):
    with Path(path).open(newline="") as stream: return list(csv.DictReader(stream))
def write_json(path, value, atomic=False):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    data=json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+"\n"
    if not atomic: path.write_text(data); return
    tmp=path.with_name(path.name+".tmp")
    with tmp.open("w") as stream:
        stream.write(data); stream.flush(); os.fsync(stream.fileno())
    os.replace(tmp,path)

def derive_modes(pod, coefficients):
    eigen={int(r["mode"]):r for r in rows(pod)
           if r.get("case")=="JOINT" and r.get("basis_type")=="late_direction_case_balanced"}
    if not eigen: raise ValueError("missing authoritative JOINT POD eigensystem")
    ids=[]
    for row in rows(coefficients):
        if row.get("case") in ("E0","E1"):
            mode=int(row["pod_mode"])
            if mode not in ids: ids.append(mode)
    if not ids or ids!=sorted(set(ids)) or any(i not in eigen for i in ids):
        raise ValueError("authoritative E2 mode mapping/order is ambiguous")
    weights=[float(eigen[i]["energy_fraction"]) for i in ids]
    if any(w<0 for w in weights): raise ValueError("negative authoritative POD energy")
    return {"authoritative_mode_count":len(ids),"authoritative_mode_ids":ids,
            "authoritative_mode_order":ids,"authoritative_energy_weights":weights,
            "authoritative_selected_energy_fraction":sum(weights)}

def mode_contract(a):
    pod,coeff=Path(a.pod),Path(a.coefficients)
    write_json(a.output,{"schema_id":"strict-ala-stage-F3-authoritative-E2-modes",
      "schema_version":"1","stage_id":STAGE_ID,
      "source_E2_artifact_paths":[str(pod.resolve()),str(coeff.resolve())],
      "source_E2_artifact_hashes":[sha256(pod),sha256(coeff)],
      **derive_modes(pod,coeff),"coefficient_trajectory_source":str(coeff.resolve()),
      "singular_value_energy_source":str(pod.resolve()),
      "derivation_method":"coefficient_first_appearance_order_mapped_to_JOINT_late_direction_case_balanced_energy",
      "derivation_complete":True,"pod_reselection_performed":False,
      "pod_retraining_performed":False,"complete":True})

NUMERIC={
 "S_ref_inner_relative_tolerance":"S_REF_INNER_RELATIVE_TOLERANCE",
 "S_ref_max_MG_cycles":"S_REF_MAX_MG_CYCLES",
 "repeat_action_relative_tolerance":"REPEAT_ACTION_RELATIVE_TOLERANCE",
 "scaling_scalar_c":"SCALING_SCALAR_C",
 "scaling_relative_tolerance":"SCALING_RELATIVE_TOLERANCE",
 "additivity_relative_tolerance":"ADDITIVITY_RELATIVE_TOLERANCE",
 "modal_significance_floor_relative":"MODAL_SIGNIFICANCE_FLOOR_RELATIVE",
 "qr_first_pass_orthogonality_threshold":"QR_FIRST_PASS_ORTHOGONALITY_THRESHOLD",
 "qr_rank_rejection_relative_threshold":"QR_RANK_REJECTION_RELATIVE_THRESHOLD",
 "projected_subspace_mode_leakage_limit":"PROJECTED_SUBSPACE_MODE_LEAKAGE_LIMIT",
 "projected_subspace_weighted_leakage_limit":"PROJECTED_SUBSPACE_WEIGHTED_LEAKAGE_LIMIT",
 "Y_mode_reachability_error_limit":"Y_MODE_REACHABILITY_ERROR_LIMIT",
 "Y_weighted_reachability_error_limit":"Y_WEIGHTED_REACHABILITY_ERROR_LIMIT",
 "POD_late_median_explanation_min":"POD_LATE_MEDIAN_EXPLANATION_MIN",
 "POD_iteration40_explanation_min":"POD_ITERATION40_EXPLANATION_MIN",
 "mode_stagnation_ratio_min":"MODE_STAGNATION_RATIO_MIN",
 "mode_decay_ratio_max":"MODE_DECAY_RATIO_MAX",
 "authoritative_mode_energy_min":"AUTHORITATIVE_MODE_ENERGY_MIN",
 "projected_condition_number_adverse":"PROJECTED_CONDITION_NUMBER_ADVERSE",
 "projected_mixing_adverse":"PROJECTED_MIXING_ADVERSE",
 "projected_non_normality_adverse":"PROJECTED_NON_NORMALITY_ADVERSE",
 "sigma_zero_relative_threshold":"SIGMA_ZERO_RELATIVE_THRESHOLD",
 "H_true_symmetry_relative_tolerance":"H_TRUE_SYMMETRY_RELATIVE_TOLERANCE",
 "H_true_negative_curvature_relative_tolerance":"H_TRUE_NEGATIVE_CURVATURE_RELATIVE_TOLERANCE",
 "pod_reconstruction_coefficient_relative_tolerance":"POD_RECONSTRUCTION_COEFFICIENT_RELATIVE_TOLERANCE",
 "projection_monotonicity_absolute_tolerance":"PROJECTION_MONOTONICITY_ABSOLUTE_TOLERANCE",
 "pod_weight_sum_absolute_tolerance":"POD_WEIGHT_SUM_ABSOLUTE_TOLERANCE",
 "inner_target_relative_roundoff":"INNER_TARGET_RELATIVE_ROUNDOFF"}

def threshold_contract(a):
    source=Path(a.thresholds); values=json.loads(source.read_text())
    missing=[k for k in NUMERIC if k not in values]
    if missing: raise ValueError("missing threshold keys: "+",".join(missing))
    canonical=sha256(source)
    mirror={"schema_id":"strict-ala-stage-F3-runtime-threshold-contract",
      "schema_version":"1","stage_id":STAGE_ID,"canonical_threshold_sha256":canonical,
      "values":values,"complete":True}
    write_json(a.mirror,mirror)
    lines=["/* Generated from canonical F3 threshold JSON. */",
      "#ifndef STRICT_ALA_STAGE_F3_THRESHOLDS_GENERATED_H",
      "#define STRICT_ALA_STAGE_F3_THRESHOLDS_GENERATED_H",
      '#define STRICT_ALA_STAGE_F3_THRESHOLD_SOURCE_SHA256 "'+canonical+'"']
    for key,macro in NUMERIC.items():
        value=values[key]
        if isinstance(value,bool) or not isinstance(value,(int,float)):
            raise ValueError("threshold is not numeric: "+key)
        rendered=str(value) if isinstance(value,int) else "%.17g"%float(value)
        lines.append("#define STRICT_ALA_STAGE_F3_%s %s"%(macro,rendered))
    lines.append("#endif")
    Path(a.header).parent.mkdir(parents=True,exist_ok=True)
    Path(a.header).write_text("\n".join(lines)+"\n")
    if a.verify:
        if json.loads(Path(a.mirror).read_text())!=mirror or canonical not in Path(a.header).read_text():
            raise ValueError("generated threshold contract mismatch")

CSV={
 "01_outer_krylov_relocalization/strict_ala_stage_F3_iteration_pod_projection.csv":["iteration","restart_cycle","mode_id","mode_energy_weight","q_dot_v","q_dot_z","q_dot_s","v_norm","z_norm","s_norm","requested_inner_tolerance","production_action_sampled_before_arnoldi"],
 "01_outer_krylov_relocalization/strict_ala_stage_F3_basis_rank.csv":["iteration","space","scope","accepted","rejected","rejection_reason","pre_norm","post_first_pass_norm","post_second_pass_norm","first_pass_orthogonality_defect","first_pass_defect_exceeds_threshold","relative_remaining_norm","qr_first_pass_orthogonality_threshold_loaded","qr_rank_rejection_relative_threshold_loaded","second_pass_performed","resulting_rank"],
 "01_outer_krylov_relocalization/strict_ala_stage_F3_mode_coefficients.csv":["iteration","mode_id","mode_energy_weight","coefficient","absolute_modal_amplitude","residual_global_pdot_norm","mode_status"],
 "01_outer_krylov_relocalization/strict_ala_stage_F3_pod_explanation.csv":["iteration","residual_global_pdot_norm","POD_projected_norm","f_Q","continuity_relative","momentum_relative","krylov_recursive_absolute","krylov_explicit_absolute","cumulative_inner_solves","cumulative_MG_cycles","cumulative_schur_actions","elapsed_seconds"],
 "01_outer_krylov_relocalization/strict_ala_stage_F3_cumulative_subspace.csv":["iteration","restart_cycle","row_type","mode_id","mode_energy_weight","E_V_all","E_Z_all","E_Y_all","E_V_cycle","E_Z_cycle","E_Y_cycle","E_V_weighted_all","E_Z_weighted_all","E_Y_weighted_all","E_V_weighted_cycle","E_Z_weighted_cycle","E_Y_weighted_cycle","descriptive_geometry_only_VZ"],
 "01_outer_krylov_relocalization/strict_ala_stage_F3_restart_retention.csv":["restart_boundary_iteration","restart_cycle","mode_id","pre_restart_Y_reachability","post_restart_iteration","post_restart_cycle_local_Y_reachability","recovery_metric_available","recovery_iterations"]}
JSON={
 "00_manifest/strict_ala_stage_F3_manifest.json":["stage_id","schema_version","expected_mpi_world_size","required_runtime_evidence_relative_paths","required_final_artifact_relative_paths","production_default_change_authorized"],
 "00_manifest/strict_ala_stage_F3_authoritative_E2_modes.json":["schema_id","schema_version","stage_id","authoritative_mode_count","authoritative_mode_ids","authoritative_mode_order","authoritative_energy_weights","derivation_complete","complete"],
 "00_manifest/strict_ala_stage_F3_artifact_schema.json":["schema_id","schema_version","stage_id","artifacts","complete"],
 "01_outer_krylov_relocalization/strict_ala_stage_F3_pod_reconstruction.json":["schema_id","schema_version","stage_id","manifest_sha256","authoritative_mode_set_reused","Q_E2_validation_complete","complete"],
 "01_outer_krylov_relocalization/strict_ala_stage_F3_runtime_provenance.json":["schema_id","schema_version","stage_id","manifest_sha256","actual_mpi_world_size","actual_restart","actual_outer_budget","runtime_threshold_values","viscosity_provenance","complete"],
 "01_outer_krylov_relocalization/strict_ala_stage_F3_numerical_completion.json":["schema_id","schema_version","stage_id","manifest_sha256","complete","valid","exit_status"],
 "02_analysis/strict_ala_stage_F3_reference_operator_validation.json":["schema_id","schema_version","stage_id","manifest_sha256","H_right_fixed_operator_valid","complete"],
 "02_analysis/strict_ala_stage_F3_projected_operator.json":["schema_id","schema_version","stage_id","manifest_sha256","H_true_definition","H_right_definition","complete"],
 "02_analysis/strict_ala_stage_F3_decision.json":["schema_id","schema_version","stage_id","manifest_sha256","decision","production_default_change_authorized","complete"],
 "03_final/strict_ala_stage_F3_evidence_index.json":["schema_id","schema_version","stage_id","manifest_sha256","artifacts","errors","complete"],
 "03_final/strict_ala_stage_F3_final_audit.json":["schema_id","schema_version","stage_id","manifest_sha256","decision","evidence_chain_errors","production_default_change_authorized","complete"]}

def schema_contract(a):
    artifacts=[]
    for path,keys in JSON.items(): artifacts.append({"relative_path":path,"artifact_type":"json","schema_id":Path(path).stem,"schema_version":"1","required_fields":{k:"present" for k in keys},"completeness_rule":"exists_nonempty_and_valid_json"})
    artifacts.append({"relative_path":"00_manifest/strict_ala_stage_F3_manifest.sha256","artifact_type":"sha256","schema_id":"sha256sum-line","schema_version":"1","required_fields":{},"completeness_rule":"exists_nonempty_and_matches_manifest"})
    for path,header in CSV.items():
        minimum=0 if path.endswith("strict_ala_stage_F3_restart_retention.csv") else 1
        artifacts.append({"relative_path":path,"artifact_type":"csv","schema_id":Path(path).stem,"schema_version":"1","exact_ordered_header":header,"minimum_record_count":minimum,"required_finite_non_null_columns":[],"completeness_rule":("header_exact; records_optional_for_early_joint_convergence_before_restart" if minimum==0 else "header_exact_and_at_least_one_record")})
    write_json(a.output,{"schema_id":"strict-ala-stage-F3-artifact-schema","schema_version":"1","stage_id":STAGE_ID,"artifacts":artifacts,"complete":True})

def main():
    p=argparse.ArgumentParser(); sub=p.add_subparsers(dest="command",required=True)
    q=sub.add_parser("mode-contract"); q.add_argument("--pod",required=True); q.add_argument("--coefficients",required=True); q.add_argument("--output",required=True); q.set_defaults(fn=mode_contract)
    q=sub.add_parser("threshold-contract"); q.add_argument("--thresholds",required=True); q.add_argument("--header",required=True); q.add_argument("--mirror",required=True); q.add_argument("--verify",action="store_true"); q.set_defaults(fn=threshold_contract)
    q=sub.add_parser("artifact-schema"); q.add_argument("--output",required=True); q.set_defaults(fn=schema_contract)
    a=p.parse_args(); a.fn(a)
if __name__=="__main__": main()
