import csv
import importlib.util
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "f3", ROOT / "tools/analyze_strict_ala_stage_F3.py")
F3 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(F3)
GEN_SPEC = importlib.util.spec_from_file_location(
    "f3_contracts", ROOT / "tools/generate_strict_ala_stage_F3_contracts.py")
GEN = importlib.util.module_from_spec(GEN_SPEC)
GEN_SPEC.loader.exec_module(GEN)
RUNS = ROOT.parents[1] / "runs"


def write_csv(path, fields, rows):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


class StageF3Tests(unittest.TestCase):
    def test_source_hooks_real_right_fgmres_before_arnoldi(self):
        source = (ROOT / "lib/Stokes_flow_Incomp.c").read_text()
        stage = (ROOT / "lib/Strict_ala_stage_f3.inc").read_text()
        call = source.index("strict_ala_stage_f3_direction")
        orth = source.index("for(i=0;i<=j;i++)", call)
        self.assertLess(call, orth)
        self.assertIn("RIGHT_FGMRES", stage)
        self.assertIn("Deterministic two-pass distributed MGS", stage)
        self.assertIn("global_pdot", stage)
        self.assertIn("base_plateau_confirmed", stage)

    def test_reference_action_uses_cumulative_bounded_residual_correction(self):
        stage = (ROOT / "lib/Strict_ala_stage_f3.inc").read_text()
        action = stage[stage.index("static int ala_f3_reference_action"):
                       stage.index("static void strict_ala_stage_f3_finalize")]
        marker_set = action.index(
            'setenv("STRICT_ALA_STAGE_F3_REFERENCE_ACTION","1",1)')
        first_solve = action.index("valid=solve_del2_u_bounded", marker_set)
        correction_solve = action.index(
            "correction_valid=solve_del2_u_bounded", first_solve)
        marker_unset = action.index(
            'unsetenv("STRICT_ALA_STAGE_F3_REFERENCE_ACTION")',
            correction_solve)
        self.assertLess(marker_set, first_solve)
        self.assertLess(first_solve, correction_solve)
        self.assertLess(correction_solve, marker_unset)
        self.assertIn("total_cycles=0;", action)
        self.assertIn(
            "remaining_cycles=STRICT_ALA_STAGE_F3_S_REF_MAX_MG_CYCLES-",
            action)
        self.assertIn("correction_rhs2>0.0", action)
        self.assertIn("STRICT_ALA_STAGE_F3_S_REF_INNER_RELATIVE_TOLERANCE*",
                      action)
        self.assertIn("total_cycles += E->control.ala_stage_f3_last_inner_cycles",
                      action)
        rhs_restore = action.index(
            "rhs[m][e]=saved[m][e];\n"
            "        (void)ala_schur_velocity_residual",
            first_solve)
        self.assertLess(first_solve, rhs_restore)
        self.assertLess(rhs_restore, correction_solve)
        self.assertNotIn(
            "0.1*STRICT_ALA_STAGE_F3_S_REF_INNER_RELATIVE_TOLERANCE",
            action)

    def test_threshold_contract_contains_four_audit_corrections(self):
        thresholds = json.loads((RUNS /
            "cmbhf_ALA_strict_stage_F3_thresholds.json").read_text())
        F3.validate_thresholds(thresholds)
        self.assertTrue(thresholds["derived_kappa_infinity_valid"])
        self.assertEqual(thresholds["qr_policy"],
            "deterministic_always_two_pass_batched_MGS_global_pdot")
        self.assertEqual(thresholds["H_true_symmetry_relative_tolerance"], 1e-8)
        self.assertIn("iterations==60", thresholds["base_plateau_definition"])

    def test_derived_infinite_condition_number_is_json_safe(self):
        metrics = F3.projected_metrics(
            F3.np.array([[1.0, 0.0], [0.0, 0.0]]), 1e-14)
        self.assertIsNone(metrics["kappa_2"])
        self.assertTrue(metrics["kappa_2_is_infinite"])
        json.dumps(metrics, allow_nan=False)

    def make_case(self, directory, htrue_valid=True):
        root = Path(directory)
        thresholds = RUNS / "cmbhf_ALA_strict_stage_F3_thresholds.json"
        weights = [0.85, 0.095, 0.016]
        runtime = root / "runtime.json"
        runtime.write_text(json.dumps({
            "outer_preconditioning_orientation": "RIGHT_FGMRES",
            "iterations": 60, "joint_target_reached": False,
            "best_R_cont": .05, "final_R_cont": .06,
            "R_cont_target": .01, "restart_boundary_observed": True}))
        pod = root / "pod.json"
        pod.write_text(json.dumps({"valid": True, "selected_modes": 3,
                                   "pod_weight_sum": sum(weights)}))
        raw = root / "reference.json"
        raw.write_text(json.dumps({
            "S_ref_validation_pass": True,
            "H_true_validation_pass": htrue_valid,
            "H_right_fixed_operator_valid": True,
            "H_true_minimum_symmetric_eigenvalue": 1.0,
            "reference_action_norm_scale": 1.0,
            "numerical_floor": 1e-14,
            "raw_vector_action_nonfinite": False,
            "leakage_per_mode": [.01, .01, .01],
            "leakage_weighted": .01}))
        completion = root / "completion.json"
        completion.write_text(json.dumps({"complete": True, "valid": True,
                                          "exit_status": 0}))
        inner = root / "inner.csv"
        write_csv(inner, ("status", "requested_relative_tolerance",
                          "achieved_relative"), [{"status": "CONVERGED",
                          "requested_relative_tolerance": 1e-2,
                          "achieved_relative": 9e-3}])
        explain_rows = [{
            "iteration": k, "residual_global_pdot_norm": 1.0,
            "POD_projected_norm": .95, "f_Q": .9025,
            "continuity_relative": .05, "momentum_relative": .0005,
            "krylov_recursive_absolute": 1.0,
            "krylov_explicit_absolute": 1.0,
            "cumulative_inner_solves": k, "cumulative_MG_cycles": 10*k,
            "cumulative_schur_actions": k, "elapsed_seconds": k}
            for k in range(0, 61)]
        explain = root / "explain.csv"
        write_csv(explain, explain_rows[0].keys(), explain_rows)
        mode_rows = []
        for k in range(0, 61):
            for mode, weight in enumerate(weights, 1):
                amplitude = .9 if mode == 1 else .02 / mode
                mode_rows.append({"iteration": k, "mode_id": mode,
                    "mode_energy_weight": weight, "coefficient": amplitude,
                    "absolute_modal_amplitude": amplitude,
                    "residual_global_pdot_norm": 1.0,
                    "mode_status": "PENDING_ANALYSIS"})
        modes = root / "modes.csv"
        write_csv(modes, mode_rows[0].keys(), mode_rows)
        sub_rows = []
        for mode, weight in enumerate(weights, 1):
            sub_rows.append({"iteration": 40, "restart_cycle": 1,
                "row_type": "MODE", "mode_id": mode,
                "mode_energy_weight": weight, "E_V_all": .5,
                "E_Z_all": .5, "E_Y_all": .1, "E_V_cycle": .5,
                "E_Z_cycle": .5, "E_Y_cycle": .1,
                "E_V_weighted_all": "", "E_Z_weighted_all": "",
                "E_Y_weighted_all": "", "E_V_weighted_cycle": "",
                "E_Z_weighted_cycle": "", "E_Y_weighted_cycle": "",
                "descriptive_geometry_only_VZ": 1})
        sub_rows.append({"iteration": 40, "restart_cycle": 1,
            "row_type": "AGGREGATE", "mode_id": 0,
            "mode_energy_weight": sum(weights), "E_V_all": "",
            "E_Z_all": "", "E_Y_all": "", "E_V_cycle": "",
            "E_Z_cycle": "", "E_Y_cycle": "", "E_V_weighted_all": .5,
            "E_Z_weighted_all": .5, "E_Y_weighted_all": .1,
            "E_V_weighted_cycle": .5, "E_Z_weighted_cycle": .5,
            "E_Y_weighted_cycle": .1, "descriptive_geometry_only_VZ": 1})
        subspace = root / "subspace.csv"
        write_csv(subspace, sub_rows[0].keys(), sub_rows)
        rank = root / "rank.csv"
        write_csv(rank, ("iteration", "space", "scope", "accepted",
            "rejected", "rejection_reason", "pre_norm",
            "post_first_pass_norm", "post_second_pass_norm",
            "first_pass_orthogonality_defect",
            "first_pass_defect_exceeds_threshold", "relative_remaining_norm",
            "qr_first_pass_orthogonality_threshold_loaded",
            "qr_rank_rejection_relative_threshold_loaded", "second_pass_performed",
            "resulting_rank"), [{"iteration": 1, "space": "Y",
            "scope": "ALL", "accepted": 1, "rejected": 0,
            "rejection_reason": "NONE", "pre_norm": 1,
            "post_first_pass_norm": 1, "post_second_pass_norm": 1,
            "first_pass_orthogonality_defect": 0,
            "first_pass_defect_exceeds_threshold": 0,
            "relative_remaining_norm": 1,
            "qr_first_pass_orthogonality_threshold_loaded": 1e-10,
            "qr_rank_rejection_relative_threshold_loaded": 1e-12,
            "second_pass_performed": 1, "resulting_rank": 1}])
        matrices = root / "matrices.csv"
        matrix_rows = []
        for name, values in (("H_true", [1, 1, 1]),
                             ("H_right", [1, 1e-4, 1])):
            for i in range(3):
                for j in range(3):
                    matrix_rows.append({"matrix": name, "row_mode": i+1,
                        "column_mode": j+1,
                        "value": values[i] if i == j else 0})
        write_csv(matrices, matrix_rows[0].keys(), matrix_rows)
        reconstructed = root / "reconstructed.csv"
        authoritative = root / "authoritative.csv"
        coeff_rows = [{"case": case, "iteration": iteration,
                       "pod_mode": mode, "coefficient": mode*.01+iteration}
                      for case in ("E0", "E1") for iteration in range(19)
                      for mode in range(1, 4)]
        write_csv(reconstructed, coeff_rows[0].keys(), coeff_rows)
        auth_rows = [dict(r, contraction=1, sign_change=0) for r in coeff_rows]
        write_csv(authoritative, auth_rows[0].keys(), auth_rows)
        baseline = root / "baseline.csv"
        base_rows = [{"case": "INNER_1E2", "iteration": k,
                      "continuity_relative": .05,
                      "momentum_relative": .0005} for k in range(1, 61)]
        write_csv(baseline, base_rows[0].keys(), base_rows)
        return SimpleNamespace(thresholds=thresholds, runtime=runtime,
            pod_reconstruction=pod, reference_validation_raw=raw,
            completion=completion, pod_explanation=explain,
            inner_solves=inner,
            mode_coefficients=modes, cumulative_subspace=subspace,
            basis_rank=rank, projected_matrices=matrices,
            reconstructed_coefficients=reconstructed,
            authoritative_coefficients=authoritative,
            baseline_iterations=baseline,
            projected_operator_output=root/"projected.json",
            reference_validation_output=root/"reference_final.json",
            output=root/"decision.json")

    def test_machine_classification_projected_operator_adverse(self):
        with tempfile.TemporaryDirectory() as directory:
            args = self.make_case(directory)
            self.assertEqual(F3.analyze(args), 0)
            decision = json.loads(args.output.read_text())
            self.assertEqual(decision["decision"],
                             "OUTER_SPACE_REACHABLE_BUT_PROJECTED_OPERATOR_ADVERSE")
            self.assertIn("ILL_CONDITIONED", decision["adverse_causes"])
            self.assertTrue(decision["base_plateau_confirmed"])

    def test_htrue_validation_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            args = self.make_case(directory, htrue_valid=False)
            F3.analyze(args)
            decision = json.loads(args.output.read_text())
            self.assertEqual(decision["decision"], "INVALID_EXPERIMENT")

    def test_early_joint_convergence_is_valid_shortened_window(self):
        with tempfile.TemporaryDirectory() as directory:
            args = self.make_case(directory)
            runtime = json.loads(args.runtime.read_text())
            runtime.update({"iterations": 35, "joint_target_reached": True,
                            "best_R_cont": .005})
            args.runtime.write_text(json.dumps(runtime))
            rows = F3.load_csv(args.pod_explanation)
            write_csv(args.pod_explanation, rows[0].keys(),
                      [row for row in rows if int(row["iteration"]) <= 35])
            modes = F3.load_csv(args.mode_coefficients)
            write_csv(args.mode_coefficients, modes[0].keys(),
                      [row for row in modes if int(row["iteration"]) <= 35])
            subspace = F3.load_csv(args.cumulative_subspace)
            for row in subspace:
                row["iteration"] = 34
            write_csv(args.cumulative_subspace, subspace[0].keys(), subspace)
            baseline = F3.load_csv(args.baseline_iterations)
            write_csv(args.baseline_iterations, baseline[0].keys(),
                      [row for row in baseline if int(row["iteration"]) <= 35])
            self.assertEqual(F3.analyze(args), 0)
            decision = json.loads(args.output.read_text())
            self.assertTrue(decision["valid"])
            self.assertFalse(decision["base_plateau_confirmed"])
            self.assertTrue(decision["shortened_window_due_to_early_joint_convergence"])
            self.assertEqual(decision["reachability_evaluation_iteration"], 34)

    def test_missing_projected_matrix_fails_closed_without_crash(self):
        with tempfile.TemporaryDirectory() as directory:
            args = self.make_case(directory)
            rows = [row for row in F3.load_csv(args.projected_matrices)
                    if row["matrix"] != "H_right"]
            write_csv(args.projected_matrices, rows[0].keys(), rows)
            self.assertEqual(F3.analyze(args), 0)
            decision = json.loads(args.output.read_text())
            self.assertEqual(decision["decision"], "INVALID_EXPERIMENT")

    def test_raw_nonfinite_csv_fails_closed_and_outputs_strict_json(self):
        with tempfile.TemporaryDirectory() as directory:
            args = self.make_case(directory)
            rows = F3.load_csv(args.mode_coefficients)
            rows[0]["coefficient"] = "nan"
            write_csv(args.mode_coefficients, rows[0].keys(), rows)
            self.assertEqual(F3.analyze(args), 0)
            decision = json.loads(args.output.read_text())
            self.assertEqual(decision["decision"], "INVALID_EXPERIMENT")

    def test_lsf_single_384_rank_production_process(self):
        text = (RUNS / "cmbhf_ALA_strict_stage_F3.lsf").read_text()
        self.assertIn("#BSUB -n 400", text)
        self.assertIn("--nodes=384 --macros.nodes=384", text)
        self.assertEqual(text.count("pythia pyre.schedulers:jobstart"), 1)
        self.assertIn("-c 'import numpy'", text)
        self.assertIn("LINEAGE_E2_ROOT", text)
        self.assertIn("reproduction_envelope.json", text)
        self.assertIn("E2_snapshots.sha256", text)
        self.assertIn("MPI_cap_decomposition", text)
        self.assertIn("STRICT_STAGE_F3_${LSB_JOBID}", text)
        self.assertNotIn("F2H", text)

    def test_lsf_reproduction_reference_resolves_e2_root(self):
        text = (RUNS / "cmbhf_ALA_strict_stage_F3.lsf").read_text()
        self.assertIn("print(path.parents[3])", text)
        self.assertNotIn("print(path.parents[2])", text)
        self.assertIn("restored E2 evidence hash mismatch", text)
        self.assertIn("authoritative_E2_lineage_root", text)
        self.assertIn("authoritative_E2_runtime_root", text)
        self.assertNotIn(
            "pathlib.Path(sys.argv[1]).resolve())' \"${E2_ROOT}\")\" =",
            text)

    def test_no_pod_reselection_and_single_q_e2_source(self):
        stage=(ROOT/"lib/Strict_ala_stage_f3.inc").read_text()
        self.assertNotIn("ala_f1a_selected_count",stage)
        self.assertIn("authoritative_mode_set_reused",stage)
        self.assertIn("s->mode[i]",stage)

    def test_qr_always_two_pass_and_generated_rank_gate(self):
        stage=(ROOT/"lib/Strict_ala_stage_f3.inc").read_text()
        body=stage[stage.index("static void ala_f3_qr_insert"):stage.index("static double ala_f3_projection_error")]
        self.assertGreaterEqual(body.count("for(i=0;i<basis->rank;i++)"),3)
        self.assertIn("QR_RANK_REJECTION_RELATIVE_THRESHOLD",body)
        self.assertIn("first_pass_defect_exceeds_threshold",stage)

    def test_sref_residual_nonzero_ratio(self):
        row=F3.sref_residual_record(4.0,1.0,2,0.6,3,10)
        self.assertEqual(row["relative_residual"],0.5)
        self.assertTrue(row["residual_gate_pass"])

    def test_sref_zero_rhs_exact_is_valid_action(self):
        row=F3.sref_residual_record(0.0,0.0,10,1e-10,0,1024)
        self.assertEqual(row["residual_status"],"ZERO_RHS_EXACT")
        self.assertEqual(row["relative_residual"],0.0)
        self.assertTrue(row["residual_gate_pass"])

    def test_sref_zero_rhs_nonzero_residual_fails(self):
        row=F3.sref_residual_record(0.0,1.0,10,1e-10,0,1024)
        self.assertIsNone(row["relative_residual"])
        self.assertFalse(row["residual_gate_pass"])

    def test_ordered_viscosity_checksum_determinism_and_sensitivity(self):
        records=[{"rank_id":0,"local_entry_count":2,"local_viscosity_hash":"aa"},
                 {"rank_id":1,"local_entry_count":2,"local_viscosity_hash":"bb"}]
        baseline=F3.ordered_viscosity_checksum(records,2)
        self.assertEqual(baseline,F3.ordered_viscosity_checksum(list(records),2))
        changed=[dict(x) for x in records]; changed[1]["local_viscosity_hash"]="bc"
        self.assertNotEqual(baseline,F3.ordered_viscosity_checksum(changed,2))

    def test_viscosity_checksum_rejects_rank_count_and_permutation_identity(self):
        records=[{"rank_id":1,"local_entry_count":2,"local_viscosity_hash":"aa"},
                 {"rank_id":0,"local_entry_count":2,"local_viscosity_hash":"bb"}]
        digest=F3.ordered_viscosity_checksum(records,2)
        swapped=[dict(records[0],rank_id=0),dict(records[1],rank_id=1)]
        self.assertNotEqual(digest,F3.ordered_viscosity_checksum(swapped,2))
        with self.assertRaises(ValueError): F3.ordered_viscosity_checksum(records,400)

    def test_root_relative_resolution_rejects_escape(self):
        self.assertEqual(F3.resolve_relative("/a","b/c"),Path("/a/b/c"))
        with self.assertRaises(ValueError): F3.resolve_relative("/a","../c")
        with self.assertRaises(ValueError): F3.resolve_relative("/a","/c")

    def test_json_and_csv_schema_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); bad=root/"bad.json"; bad.write_text("{}")
            ok,status=F3.validate_artifact(bad,{"artifact_type":"json","required_fields":{"x":"present"}})
            self.assertFalse(ok); self.assertIn("missing_json_fields",status)
            csv_path=root/"bad.csv"; csv_path.write_text("b,a\n1,2\n")
            ok,status=F3.validate_artifact(csv_path,{"artifact_type":"csv","exact_ordered_header":["a","b"],"minimum_record_count":1})
            self.assertFalse(ok); self.assertEqual(status,"csv_header_mismatch")

    def test_missing_and_empty_artifacts_fail_schema(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"x"; spec={"artifact_type":"json","required_fields":{}}
            self.assertFalse(F3.validate_artifact(path,spec)[0])
            path.write_text(""); self.assertFalse(F3.validate_artifact(path,spec)[0])

    def test_threshold_generator_is_deterministic_and_hash_grounded(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); outputs=[]
            for index in range(2):
                header=root/f"h{index}.h"; mirror=root/f"m{index}.json"
                args=SimpleNamespace(thresholds=RUNS/"cmbhf_ALA_strict_stage_F3_thresholds.json",header=header,mirror=mirror,verify=True)
                GEN.threshold_contract(args); outputs.append((header.read_bytes(),mirror.read_bytes()))
            self.assertEqual(outputs[0],outputs[1])

    def test_threshold_identity_mismatch_fails(self):
        canonical=json.loads((RUNS/"cmbhf_ALA_strict_stage_F3_thresholds.json").read_text())
        digest=hashlib.sha256((RUNS/"cmbhf_ALA_strict_stage_F3_thresholds.json").read_bytes()).hexdigest()
        mirror={"canonical_threshold_sha256":digest,"values":canonical}
        runtime={"canonical_threshold_sha256":digest,"runtime_threshold_values":{"S_ref_max_MG_cycles":1}}
        self.assertTrue(F3.validate_threshold_identity(canonical,digest,mirror,runtime))

    def test_mode_contract_unique_ids_order_and_weights(self):
        evidence=Path("/Volumes/Pōwehi/CitcomS_DATA/STRICT_STAGE_E2_12121889/03_subspace")
        if not evidence.is_dir(): self.skipTest("local E2 evidence unavailable")
        value=GEN.derive_modes(evidence/"strict_ala_stage_E2_pod.csv",evidence/"strict_ala_stage_E2_joint_pod_coefficients.csv")
        self.assertEqual(value["authoritative_mode_count"],2)
        self.assertEqual(value["authoritative_mode_ids"],[1,2])
        self.assertEqual(value["authoritative_mode_order"],[1,2])
        self.assertAlmostEqual(sum(value["authoritative_energy_weights"]),0.9438460932970015)

    def test_source_operator_viscosity_is_finest_evi(self):
        element=(ROOT/"lib/Element_calculations.c").read_text()
        stage=(ROOT/"lib/Strict_ala_stage_f3.inc").read_text()
        self.assertIn("E->EVI[lev][m][(el-1)*vpts+k]",element)
        self.assertIn("E->EVI[E->mesh.levmax]",stage)
        self.assertIn("MPI_Gather",stage)

    def test_manifest_hash_modification_is_detectable(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"manifest.json"; path.write_text("{}\n")
            frozen=F3.sha256(path); path.write_text('{"changed":true}\n')
            self.assertNotEqual(frozen,F3.sha256(path))

    def test_artifact_schema_declares_final_outputs_without_self_hash_cycle(self):
        with tempfile.TemporaryDirectory() as directory:
            output=Path(directory)/"schema.json"
            GEN.schema_contract(SimpleNamespace(output=output))
            artifacts={x["relative_path"]:x for x in
                       json.loads(output.read_text())["artifacts"]}
            self.assertIn("03_final/strict_ala_stage_F3_evidence_index.json",artifacts)
            self.assertIn("03_final/strict_ala_stage_F3_final_audit.json",artifacts)
            restart=artifacts["01_outer_krylov_relocalization/strict_ala_stage_F3_restart_retention.csv"]
            self.assertEqual(restart["minimum_record_count"],0)
            analyzer=(ROOT/"tools/analyze_strict_ala_stage_F3.py").read_text()
            self.assertIn('rel.endswith("strict_ala_stage_F3_evidence_index.json")',analyzer)
            self.assertIn('rel.endswith("strict_ala_stage_F3_final_audit.json")',analyzer)

    def test_fq_formula_uses_global_pressure_norm_and_exact_modes(self):
        source=(ROOT/"lib/Strict_ala_stage_f3.inc").read_text()
        self.assertIn("projected2+=coefficient*coefficient",source)
        self.assertIn("projected2/max(norm*norm",source)


if __name__ == "__main__":
    unittest.main()
