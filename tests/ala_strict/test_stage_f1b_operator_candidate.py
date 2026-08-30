#!/usr/bin/env python3
import csv
import importlib.util
import json
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT.parents[1] / "runs"
SPEC = importlib.util.spec_from_file_location(
    "f1b", ROOT / "tools/analyze_strict_ala_stage_F1b.py")
F1B = importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(F1B)


class Args:
    def __init__(self, **values): self.__dict__.update(values)


def write_cfg(path, selector):
    values = {
        "ala_shallow_patch_local_operator": selector,
        "ala_shallow_patch_preconditioner": "on",
        "ala_shallow_patch_velocity_solver": "element_vanka",
        "ala_augmented_lagrangian_gamma": "10.0",
        "ala_outer_solver": "fgmres", "ala_pcg_restart_interval": "50",
        "ala_stage_abc_production_logging": "on",
        "ala_stage_e_diagnostic": "on", "piterations": "30",
        "nprocx": "4", "nprocy": "4", "nprocz": "2",
        "nodex": "129", "nodey": "129", "nodez": "65", "steps": "1"}
    mesher = {"nprocx", "nprocy", "nprocz", "nodex", "nodey", "nodez"}
    lines = ["[CitcomS]\n", "steps = 1\n", "[CitcomS.solver.mesher]\n"]
    lines.extend(f"{key} = {values[key]}\n" for key in sorted(mesher))
    lines.append("[CitcomS.solver.vsolver]\n")
    lines.extend(f"{key} = {value}\n" for key, value in values.items()
                 if key not in mesher and key != "steps")
    Path(path).write_text("".join(lines))


def write_stage(path, case, terminal_cont=0.5, terminal_mom=1e-4):
    fields = ("case", "iteration", "true_continuity_relative",
              "true_momentum_relative", "recursive_residual",
              "explicit_residual", "velocity_MG_cycles",
              "K_gamma_operator_applications", "elapsed_seconds")
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader()
        for iteration in range(31):
            fraction = iteration / 30.0
            writer.writerow({"case": case, "iteration": iteration,
                "true_continuity_relative": 1.0 + fraction * (terminal_cont - 1.0),
                "true_momentum_relative": terminal_mom,
                "recursive_residual": 1.0 - 0.5 * fraction,
                "explicit_residual": 1.0 - 0.5 * fraction,
                "velocity_MG_cycles": 100 + 10 * iteration,
                "K_gamma_operator_applications": iteration,
                "elapsed_seconds": 1 + iteration})


def write_memory(path, case, rss, retained=None):
    retained = rss * 1024 if retained is None else retained
    Path(path).write_text(
        "case,rank_rss_max_kib,rank_rss_sum_kib,ranks,"
        "retained_cache_bytes_per_rank,temporary_peak_bytes_per_rank,source\n"
        f"{case},{rss},{rss*384},384,{retained},0,"
        "getrusage_RUSAGE_SELF_MPI_reduce\n")


def write_dict_csv(path, rows):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(rows[0]))
        writer.writeheader(); writer.writerows(rows)


class F1bContractTests(unittest.TestCase):
    @staticmethod
    def audit_args(td, offline, output):
        evidence = {}
        for name in ("binary", "inputs", "source"):
            before, after = td/f"{name}.pre", td/f"{name}.post"
            before.write_text(f"frozen-{name}\n"); after.write_bytes(before.read_bytes())
            evidence[f"{name}_pre"] = before
            evidence[f"{name}_post"] = after
        return Args(offline=offline, short=None, output=output, **evidence)

    def test_source_keeps_legacy_default_and_outer_projection(self):
        instructions = (ROOT / "lib/Instructions.c").read_text()
        bridge = (ROOT / "module/setProperties.c").read_text()
        source = (ROOT / "lib/Stokes_flow_Incomp.c").read_text()
        candidate = (ROOT / "lib/Strict_ala_stage_f1b.inc").read_text()
        self.assertIn('ala_shallow_patch_local_operator,"legacy"', instructions)
        self.assertIn('"operator_consistent"', instructions)
        self.assertIn(
            'getStringProperty(properties, "ala_shallow_patch_local_operator"',
            bridge)
        self.assertIn(
            'strcmp(E->control.ala_shallow_patch_local_operator,', bridge)
        self.assertIn("ala_f1b_build_operator_consistent_patch", source)
        self.assertIn("STRICT_ALA_STAGE_F1B_EXPECTED_OPERATOR", source)
        self.assertIn("STRICT_ALA_STAGE_F1B_CONFIG_ACTIVE", source)
        self.assertIn("STRICT_ALA_STAGE_F1B_CONFIG_MISMATCH", source)
        self.assertIn("operator_consistent(GKgamma^-1G^T)", source)
        self.assertIn("direct-face plus two-hop corner element records", candidate)
        self.assertIn("Kgamma_replay_relative_defect", candidate)
        self.assertIn("ALA_F1B_RECORD_GLOBAL_DIAG", candidate)
        self.assertIn("1.0/E->BI[lev][m][eq]", candidate)
        self.assertNotIn("K[i*nv+i]=vdiag[i]", candidate)
        self.assertIn("ALA_F1B_NODE_RECORD_COEFFICIENT", candidate)
        self.assertIn("E->Node_map[lev][m]", candidate)
        self.assertIn("E->Eqn_k1[lev][m]", candidate)
        self.assertIn("E->Eqn_k2[lev][m]", candidate)
        self.assertIn("E->Eqn_k3[lev][m]", candidate)
        self.assertIn("pool->local_node_records", candidate)
        self.assertIn("pool->remote_node_records", candidate)
        self.assertIn("snapshot->element_matrix", candidate)
        self.assertIn("diagonal_completion_relative_max", candidate)
        self.assertIn("production_nodal_Eqn_k", candidate)
        self.assertIn("MPI_Abort(E->parallel.world,8)", candidate)
        self.assertIn("pivot_index=%d", candidate)
        self.assertIn("production_output_replica_relative_difference", candidate)
        self.assertIn("B_restriction_relative_defect", candidate)
        self.assertIn("Bt_restriction_relative_defect", candidate)
        self.assertIn("B_Bt_bilinear_relative_defect", candidate)
        self.assertIn("STRICT_ALA_STAGE_F1B_K_REPLAY_FORENSIC", candidate)
        self.assertIn("element_sum_relative_defect", candidate)
        self.assertIn("diagonal_completion_action_relative", candidate)
        self.assertIn("maximum_error_relative_to_action_rms", candidate)
        self.assertIn("snapshot->element_diagonal", candidate)
        self.assertIn("snapshot->incident_elements", candidate)
        self.assertIn("MPI_Bcast(forensic,sizeof(*forensic),MPI_BYTE", candidate)
        self.assertIn("Kgamma_replay_forensics", candidate)
        self.assertIn("assemble_div_rho_u", candidate)
        self.assertIn("assemble_grad_rho_p", candidate)
        self.assertIn("STRICT_ALA_STAGE_F1B_MEMORY_PREFLIGHT", source)
        self.assertIn("ala_f1b_replay_global_k(E,lev,0", source)
        self.assertIn("ala_f1b_replay_global_k(E,lev,1", source)
        # Existing pre/post Q_i and multiplicity path remains in the common
        # apply routine; the candidate changes only cached A_i factors.
        self.assertGreaterEqual(source.count("rhs[i] -= sum*constant_mode[i]"), 2)
        self.assertGreaterEqual(source.count("solution[i] -= sum*constant_mode[i]"), 2)

    def test_cfg_diff_is_exactly_one_selector(self):
        with tempfile.TemporaryDirectory() as td:
            base, cand, diff = (Path(td) / name for name in ("base", "cand", "diff"))
            write_cfg(base, "legacy"); write_cfg(cand, "operator_consistent")
            self.assertEqual(F1B.check_cfg(Args(base=base, candidate=cand, diff=diff)), 0)
            self.assertEqual(diff.read_text().strip(),
                             "ala_shallow_patch_local_operator: legacy -> operator_consistent")

    def test_cfg_editor_inserts_new_keys_in_vsolver_section(self):
        with tempfile.TemporaryDirectory() as td:
            cfg = Path(td)/"case.cfg"
            shutil.copy(RUNS/"cmbhf_ALA_strict.cfg", cfg)
            for key, value in (("ala_stage_abc_production_logging", "on"),
                               ("ala_stage_e_diagnostic", "on"),
                               ("ala_shallow_patch_local_operator",
                                "operator_consistent")):
                F1B.set_cfg_value(Args(cfg=cfg, section="CitcomS.solver.vsolver",
                                       key=key, value=value))
            entries = F1B.cfg_entries(cfg)
            for key in ("ala_stage_abc_production_logging",
                        "ala_stage_e_diagnostic",
                        "ala_shallow_patch_local_operator"):
                self.assertEqual(entries[key][1], "CitcomS.solver.vsolver")
            self.assertNotIn("full.phase.ala_", cfg.read_text())

    def test_lsf_times_executable_and_preserves_manifest_templates(self):
        lsf = (RUNS / "cmbhf_ALA_strict_stage_F1b.lsf").read_text()
        self.assertNotIn('/usr/bin/time -v -o "${B}/case_time.txt" launch_case', lsf)
        self.assertIn('/usr/bin/time -v -o "${time_file}"', lsf)
        self.assertIn('/python/pythia-0.8.1.15-py2.6.egg:', lsf)
        self.assertIn('require_path "${CODE_DIR}/python/pythia-0.8.1.15-py2.6.egg"', lsf)
        self.assertNotIn('test -f "${CODE_DIR}/python/pythia-0.8.1.15-py2.6.egg"', lsf)
        self.assertIn('STRICT_ALA_STAGE_F1B_MISSING_PATH', lsf)
        self.assertIn('if [ "${rc}" -ne 0 ]; then return "${rc}"; fi', lsf)
        self.assertIn('stage_inputs "${O}" "${OFFLINE_CFG}"', lsf)
        self.assertNotIn('stage_inputs "${O}" "${CAND_CFG}"', lsf)
        self.assertIn('completion_valid "${B}/phase3_base_completion.json"', lsf)
        self.assertIn('--binary-pre "${M}/binary.pre.sha256"', lsf)
        self.assertIn("require_operator_handshake", lsf)
        self.assertIn(
            "STRICT_ALA_STAGE_F1B_EXPECTED_OPERATOR=operator_consistent", lsf)
        self.assertIn("STRICT_ALA_STAGE_F1B_EXPECTED_OPERATOR=legacy", lsf)
        self.assertIn("binary_F1b_diagonal_completion_marker.txt", lsf)
        self.assertIn("binary_F1b_nodal_operator_marker.txt", lsf)
        self.assertIn("production_nodal_Eqn_k", lsf)
        self.assertIn("binary_F1b_replica_conformity_marker.txt", lsf)
        self.assertIn("binary_F1b_B_restriction_marker.txt", lsf)
        self.assertIn("binary_F1b_Bt_restriction_marker.txt", lsf)
        self.assertIn("binary_F1b_bilinear_marker.txt", lsf)
        self.assertIn("strict_ala_stage_F1b_factor_forensic.json", lsf)
        self.assertIn("offline_numerical_evidence_valid", lsf)
        self.assertIn("STRICT_ALA_STAGE_F1B_RESUME_NUMERICAL_EVIDENCE", lsf)
        self.assertIn("--allow-configured-replay-change", lsf)
        self.assertIn("numerical_completion_valid", lsf)
        self.assertIn("STRICT_ALA_STAGE_F1B_IDENTITY_PREFLIGHT=1", lsf)
        self.assertIn("run_f1b identity", lsf)
        self.assertIn("reason=identity", lsf)
        self.assertIn("BASE_MEMORY_PREFLIGHT", lsf)
        self.assertIn("CAND_MEMORY_PREFLIGHT", lsf)
        self.assertIn('F1B_ALLOCATED_TASKS=400', lsf)
        self.assertIn('F1B_SOLVER_WORLD_SIZE=384', lsf)
        self.assertIn('--thresholds "${M}/strict_ala_stage_F1b_thresholds.json"', lsf)
        self.assertIn('"${M}/strict_ala_stage_F1b_layout.json"', lsf)

    def test_threshold_and_mpi_layout_are_frozen_before_launch(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td); cfg = td/"case.cfg"; write_cfg(cfg, "legacy")
            args = Args(output_dir=td, allocated_tasks="400",
                        solver_world_size="384", nprocx="4", nprocy="4",
                        nprocz="2", caps="12", ranks_per_node="40",
                        node_memory_kib=str(256 * 1024 * 1024),
                        allocated_wall_seconds="604800",
                        startup_budget_seconds="7200")
            F1B.freeze_contract(args)
            threshold_path = td/"strict_ala_stage_F1b_thresholds.json"
            frozen = json.loads(threshold_path.read_text())
            self.assertEqual(frozen["threshold_status"], "FROZEN_PRELAUNCH")
            self.assertEqual(frozen["Kgamma_replay_relative_tolerance"], 1e-10)
            self.assertEqual(
                frozen["Kgamma_replay_stop_review_upper_bound"], 1e-8)
            self.assertTrue(frozen["production_freeze_contract"]
                            ["production_representative_smoke_soak_required"])
            self.assertEqual(
                frozen["declared_project_startup_budget_seconds"], 7200)
            schema = json.loads(
                (td/"strict_ala_stage_F1b_schema.json").read_text())
            self.assertEqual(
                [entry["patch_category"] for entry in schema["collective_manifest"]],
                ["LOCAL", "INTERFACE"])
            F1B.verify_layout(Args(cfg=cfg, thresholds=threshold_path,
                allocated_tasks="400", solver_world_size="384", caps="12",
                ranks_per_node="40", output=td/"layout.json"))
            layout = json.loads((td/"layout.json").read_text())
            self.assertTrue(layout["valid"])
            self.assertEqual(layout["unused_allocation_tasks"], 16)

    def test_identity_gate_precedes_authorized_candidate(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            freeze = Args(output_dir=td, allocated_tasks="400",
                solver_world_size="384", nprocx="4", nprocy="4",
                nprocz="2", caps="12", ranks_per_node="40",
                node_memory_kib=str(256*1024*1024),
                allocated_wall_seconds="604800", startup_budget_seconds="7200")
            F1B.freeze_contract(freeze)
            schema = td/"strict_ala_stage_F1b_schema.json"
            feasibility = {
                "collective_manifest_sha256": F1B.sha256(schema),
                "factorization_phase": "diagnostic_factor_attempt",
                "Kgamma_replay_relative_defect": 1e-12,
                "B_restriction_relative_defect": 1e-12,
                "Bt_restriction_relative_defect": 1e-12,
                "B_Bt_bilinear_relative_defect": 1e-12,
                "production_output_replica_relative_difference": 1e-12}
            (td/"feasibility.json").write_text(json.dumps(feasibility))
            rows = []
            for category in ("LOCAL", "INTERFACE"):
                rows.append(dict(zip(F1B.ASSEMBLY_COLUMNS,
                    ("AGGREGATE", category, "-1", "-1", "0", "0", "0", "0",
                     "PASS", "1", "0", "0", "0", "0", "1e-12", "1e-12",
                     "1e-12", "1e-12", "1e-12"))))
            write_dict_csv(td/"assembly.csv", rows)
            F1B.identity(Args(feasibility=td/"feasibility.json",
                assembly=td/"assembly.csv",
                thresholds=td/"strict_ala_stage_F1b_thresholds.json",
                schema=schema, output=td/"identity.json"))
            result = json.loads((td/"identity.json").read_text())
            self.assertTrue(result["identity_gate_pass"])
            self.assertEqual(result["factorization_phase"],
                             "diagnostic_factor_attempt")
            self.assertFalse(result["authorized_candidate_factorization"])

    def test_restriction_review_never_accepts_candidate(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            offline = {"experiment_evidence_valid": True,
                       "offline_gate_pass": False,
                       "restriction_tolerance_review_required": True,
                       "gates": {"candidate_isolatable": False,
                                 "memory": True, "assembly": False}}
            (td/"offline.json").write_text(json.dumps(offline))
            F1B.final_audit(self.audit_args(td, td/"offline.json",
                                            td/"audit.json"))
            result = json.loads((td/"audit.json").read_text())
            self.assertEqual(result["decision"],
                             "F1B_RESTRICTION_TOLERANCE_REVIEW")
            self.assertFalse(result["production_freeze_authorized"])

    def test_offline_gate_consumes_frozen_layout_memory_and_manifest(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td); out = td/"out"
            freeze = Args(output_dir=td, allocated_tasks="400",
                solver_world_size="384", nprocx="4", nprocy="4",
                nprocz="2", caps="12", ranks_per_node="40",
                node_memory_kib=str(256*1024*1024),
                allocated_wall_seconds="604800", startup_budget_seconds="7200")
            F1B.freeze_contract(freeze)
            schema = td/"strict_ala_stage_F1b_schema.json"
            thresholds = td/"strict_ala_stage_F1b_thresholds.json"
            (td/"lineage.json").write_text(json.dumps(
                {"valid": True, "H_RAW_legacy": 1.0}))
            (td/"identity.json").write_text(json.dumps(
                {"identity_gate_pass": True}))
            feasibility = {
                "collective_manifest_sha256": F1B.sha256(schema),
                "factorization_phase": "authorized_candidate_factorization",
                "candidate_isolatable": True, "factorization_count": 2,
                "RHS_solve_count": 2, "estimated_dense_factor_and_solve_work": 1,
                "temporary_peak_bytes_per_rank": 100,
                "retained_cache_bytes_per_rank": 100,
                "legacy_retained_cache_bytes_per_rank": 100,
                "retained_cache_ratio_vs_legacy": 1.0,
                "projected_peak_ratio_vs_legacy_cache": 2.0,
                "new_mpi_payload_bytes_max_per_rank": 1,
                "diagonal_completion_relative_max": 0,
                "Kgamma_replay_relative_defect": 1e-12,
                "production_output_replica_relative_difference": 1e-12,
                "B_restriction_relative_defect": 1e-12,
                "Bt_restriction_relative_defect": 1e-12,
                "B_Bt_bilinear_relative_defect": 1e-12,
                "startup_seconds_max": 60, "memory_gate_pass": True}
            (td/"feasibility.json").write_text(json.dumps(feasibility))
            assembly = []
            for category in ("LOCAL", "INTERFACE"):
                assembly.append(dict(zip(F1B.ASSEMBLY_COLUMNS,
                    ("AGGREGATE", category, "-1", "-1", "1", "1", "0", "0",
                     "PASS", "1", "0", "0", "8", "8", "1e-12", "1e-12",
                     "1e-12", "1e-12", "1e-12"))))
            write_dict_csv(td/"assembly.csv", assembly)
            local = [{"patch_or_bin_ID":"0", "POD_mode":"0",
                      "local_action_error":"1", "contribution_weight":"1"}]
            cand_local = [dict(local[0], local_action_error="0.5")]
            write_dict_csv(td/"legacy_local.csv", local)
            write_dict_csv(td/"cand_local.csv", cand_local)
            tel = [{"POD_mode":"0", "stage":"CONFIGURED",
                    "POD_energy_weight":"1", "E_P":"1", "cosine":"1",
                    "qTPq":"1", "qTSPq":"1", "tight_solve_achieved":"1"}]
            cand_tel = [dict(tel[0], E_P="0.5")]
            write_dict_csv(td/"legacy_tel.csv", tel)
            write_dict_csv(td/"cand_tel.csv", cand_tel)
            (td/"legacy_decision.json").write_text(json.dumps(
                {"weighted_telescoping_damage":{"RAW_LOCAL":1.0}}))
            (td/"cand_decision.json").write_text(json.dumps(
                {"weighted_telescoping_damage":{"RAW_LOCAL":0.4}}))
            write_dict_csv(td/"tight.csv", [{"converged":"1"}])
            (td/"snapshots.csv").write_text("snapshot\n0\n")
            write_memory(td/"base_memory.csv", "BASE_MEMORY_PREFLIGHT", 1000)
            write_memory(td/"cand_memory.csv", "CAND_MEMORY_PREFLIGHT", 1100)
            F1B.offline(Args(lineage=td/"lineage.json",
                identity_gate=td/"identity.json",
                feasibility=td/"feasibility.json", assembly=td/"assembly.csv",
                legacy_local=td/"legacy_local.csv", candidate_local=td/"cand_local.csv",
                legacy_telescoping=td/"legacy_tel.csv",
                candidate_telescoping=td/"cand_tel.csv",
                legacy_decision=td/"legacy_decision.json",
                candidate_decision=td/"cand_decision.json",
                candidate_tight=td/"tight.csv", snapshot_manifest=td/"snapshots.csv",
                thresholds=thresholds, schema=schema,
                base_memory=td/"base_memory.csv", cand_memory=td/"cand_memory.csv",
                output_dir=out))
            result = json.loads((out/"strict_ala_stage_F1b_offline_gate.json").read_text())
            self.assertTrue(result["offline_gate_pass"])
            self.assertTrue(all(result["memory_gates"][key] for key in
                ("M1_retained_cache", "M2_total_process_peak", "M3_node_safety")))

            # M1 is a real matched BASE/CAND retained-cache comparison, not
            # the former candidate/self circular ratio.
            write_memory(td/"cand_memory.csv", "CAND_MEMORY_PREFLIGHT", 1100,
                         retained=1.30 * 1000 * 1024)
            args = Args(lineage=td/"lineage.json", identity_gate=td/"identity.json",
                feasibility=td/"feasibility.json", assembly=td/"assembly.csv",
                legacy_local=td/"legacy_local.csv", candidate_local=td/"cand_local.csv",
                legacy_telescoping=td/"legacy_tel.csv",
                candidate_telescoping=td/"cand_tel.csv",
                legacy_decision=td/"legacy_decision.json",
                candidate_decision=td/"cand_decision.json",
                candidate_tight=td/"tight.csv", snapshot_manifest=td/"snapshots.csv",
                thresholds=thresholds, schema=schema,
                base_memory=td/"base_memory.csv", cand_memory=td/"cand_memory.csv",
                output_dir=td/"out_m1_fail")
            F1B.offline(args)
            failed = json.loads((td/"out_m1_fail/strict_ala_stage_F1b_offline_gate.json").read_text())
            self.assertFalse(failed["memory_gates"]["M1_retained_cache"])
            self.assertFalse(failed["offline_gate_pass"])

    def test_exact_global_diagonal_closes_restricted_spd_operator(self):
        # Exterior elements contribute to the assembled diagonal even when
        # their pressure rows are outside this patch.  Completing that exact
        # diagonal turns the incident-only singular matrix into the true SPD
        # principal restriction, without adding regularization.
        incomplete = [[1.0, -1.0], [-1.0, 1.0]]
        production_diagonal = [2.0, 3.0]
        restricted = [row[:] for row in incomplete]
        for index, value in enumerate(production_diagonal):
            restricted[index][index] = value
        self.assertEqual(restricted, [[2.0, -1.0], [-1.0, 3.0]])
        self.assertEqual(incomplete[0][0] * incomplete[1][1]
                         - incomplete[0][1] ** 2, 0.0)
        self.assertGreater(restricted[0][0] * restricted[1][1]
                           - restricted[0][1] ** 2, 0.0)

    def test_f1b_cholesky_preserves_each_diagonal_factor(self):
        candidate = (ROOT / "lib/Strict_ala_stage_f1b.inc").read_text()
        factor = candidate[candidate.index("static int ala_f1b_factor"):
                           candidate.index("static void ala_f1b_solve")]
        self.assertIn("if(j<i)\n                matrix[j*n+i]=0.0;", factor)
        self.assertNotIn("\n            matrix[j*n+i]=0.0;", factor)

        # Mirror the two-by-two recurrence that exposed the defect on HPC.
        # L00 must survive long enough to form L10 and the second pivot.
        matrix = [[4.0, 2.0], [2.0, 3.0]]
        l00 = matrix[0][0] ** 0.5
        l10 = matrix[1][0] / l00
        l11 = (matrix[1][1] - l10 * l10) ** 0.5
        self.assertEqual(l00, 2.0)
        self.assertEqual(l10, 1.0)
        self.assertEqual(l11, 2.0 ** 0.5)

    def test_interface_velocity_identity_keeps_the_441_dof_contract(self):
        candidate = (ROOT / "lib/Strict_ala_stage_f1b.inc").read_text()
        self.assertIn("#define ALA_F1B_MAX_VELOCITY_DOF", candidate)
        self.assertIn("2*(ALA_PATCH_MAX_MPI_OVERLAP+1)", candidate)
        self.assertIn("#define ALA_F1B_RECORD_NODE_KEY", candidate)
        self.assertIn("E->sphere.capid[m]", candidate)
        self.assertIn("E->lmesh.EXS[lev]+ix", candidate)
        self.assertIn("E->lmesh.EYS[lev]+iy", candidate)
        self.assertIn("E->lmesh.EZS[lev]+iz", candidate)
        self.assertIn("ala_f1b_same_node_key(node_key,vkey+4*i)", candidate)
        self.assertIn("ala_f1b_same_point(coord,vcoord+3*i)", candidate)
        self.assertIn("nearest_same_component_distance", candidate)
        self.assertIn("ala_f1b_measure_patch_topology", candidate)
        self.assertIn("ALA_F1B_FAILURE_INTERFACE_TOPOLOGY", candidate)
        self.assertIn("topology_expected_velocity_dofs", candidate)
        self.assertIn("STRICT_ALA_STAGE_F1B_TOPOLOGY_RECORD", candidate)
        self.assertIn("STRICT_ALA_STAGE_F1B_TOPOLOGY_NODE", candidate)
        self.assertIn("nearest_same_distance=%.17e", candidate)
        self.assertIn("nearest_other_distance=%.17e", candidate)
        self.assertIn("local_occurrences[index]++", candidate)
        self.assertIn("remote_occurrences[index]++", candidate)
        self.assertIn("ala_f1b_find_outward_remote_record", candidate)
        self.assertIn("shared!=4 || (incoming!=NULL && overlap!=0)",
                      candidate)
        self.assertIn("curved-grid nearest-centre shear", candidate)
        source = (ROOT / "lib/Stokes_flow_Incomp.c").read_text()
        self.assertIn("face+1,tangent,ez", source)
        self.assertNotIn("f1b_failure.interface_face=face+1", source)
        self.assertIn("ala_find_halo_record_exact", source)
        self.assertIn("if(f1b_candidate) {", source)
        self.assertIn("interface_halo_topology", source)

        # The legacy nearest-centre selector remains the BASE path; only the
        # operator-consistent candidate uses exact topology and exact centre
        # lookup, preserving the single-variable A/B contract.
        selector = source.index("f1b_remote_record=")
        candidate_branch = source[selector-1800:
                                  source.index("if(ghost_index<0)", selector)]
        self.assertIn("ala_f1b_find_outward_remote_record", candidate_branch)
        self.assertIn("ala_find_halo_record_exact", candidate_branch)
        self.assertIn("else {", candidate_branch)
        self.assertIn("ala_find_halo_record(", candidate_branch)

        # Two three-element-wide strips share one face.  Global structured
        # node keys collapse that face to 7*7*3 physical nodes, exactly the
        # frozen 441-component capacity, independently of rank-local indices.
        left = {(12, gx, gy, gz)
                for gx in range(4) for gy in range(7) for gz in range(3)}
        right = {(12, gx, gy, gz)
                 for gx in range(3, 7) for gy in range(7) for gz in range(3)}
        self.assertEqual(len(left | right), 147)
        self.assertEqual(3 * len(left | right), 441)

    def test_interface_topology_chain_rejects_tangential_shear(self):
        def element(x, y, z=0):
            return {(x+dx, y+dy, z+dz)
                    for dx in (0, 1) for dy in (0, 1)
                    for dz in (0, 1)}

        local = element(3, 0)
        q0 = element(2, 0)
        outward = element(1, 0)
        tangential_shear = element(2, 1)
        radial_shear = element(2, 0, 1)
        incoming_face = local & q0

        self.assertEqual(len(incoming_face), 4)
        self.assertEqual(len(q0 & outward), 4)
        self.assertEqual(len((q0 & outward) & incoming_face), 0)
        self.assertEqual(len(q0 & tangential_shear), 4)
        self.assertEqual(len((q0 & tangential_shear) & incoming_face), 2)
        self.assertEqual(len(q0 & radial_shear), 4)
        self.assertEqual(len((q0 & radial_shear) & incoming_face), 2)

    def test_lsf_accepts_unpacked_pythia_egg_directory(self):
        lsf = (RUNS / "cmbhf_ALA_strict_stage_F1b.lsf").read_text()
        function = re.search(r"(require_path\(\) \{.*?\n\})", lsf,
                             flags=re.DOTALL).group(1)
        with tempfile.TemporaryDirectory() as td:
            egg = Path(td)/"pythia-0.8.1.15-py2.6.egg"; egg.mkdir()
            accepted = subprocess.run(
                ["bash", "-c", function + '\nrequire_path "$1"', "bash", str(egg)],
                text=True, capture_output=True)
            self.assertEqual(accepted.returncode, 0, accepted.stderr)
            missing = subprocess.run(
                ["bash", "-c", function + '\nrequire_path "$1"', "bash",
                 str(egg/"missing")], text=True, capture_output=True)
            self.assertNotEqual(missing.returncode, 0)
            self.assertIn("STRICT_ALA_STAGE_F1B_MISSING_PATH", missing.stderr)

    def test_hpc_evidence_root_is_valid_lineage(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td); evidence = td/"f1a"
            (evidence/"01_operator_replay").mkdir(parents=True)
            (evidence/"01_operator_replay/strict_ala_stage_F1a_projected_matrices.csv").write_text("x\n")
            final = {"experiment_valid": True,
                     "primary_defect_class": "LOCAL_OPERATOR_MISMATCH",
                     "dominant_damage_stage": "RAW_LOCAL", "F1b_authorized": True,
                     "next_authorized_task": "F1B_LOCAL_CANDIDATE",
                     "numerical_evidence_root": str(evidence),
                     "weighted_telescoping_damage": {"RAW_LOCAL": 1.0}}
            (td/"final.json").write_text(json.dumps(final))
            F1B.verify_lineage(Args(final_audit=td/"final.json",
                evidence_root=evidence, expected_hpc_root=str(evidence),
                output=td/"lineage.json"))
            self.assertTrue(json.loads((td/"lineage.json").read_text())["valid"])

    def test_short_routes_use_discrete_same_state_rows(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td); write_stage(td/"base.csv", "BASE", 0.5)
            write_stage(td/"cand.csv", "CAND", 0.35)
            write_memory(td/"bm.csv", "BASE", 1000)
            write_memory(td/"cm.csv", "CAND", 1100)
            (td/"bw").write_text("40\n"); (td/"cw").write_text("42\n")
            F1B.short(Args(base_iterations=td/"base.csv", cand_iterations=td/"cand.csv",
                base_memory=td/"bm.csv", cand_memory=td/"cm.csv",
                base_wall=td/"bw", cand_wall=td/"cw", output_dir=td/"out"))
            result = json.loads((td/"out/strict_ala_stage_F1b_short_decision.json").read_text())
            self.assertTrue(result["route_A_pass"])
            self.assertTrue(result["memory_gate_pass"])
            self.assertTrue(result["short_AB_pass"])
            with (td/"out/strict_ala_stage_F1b_short_iterations.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(rows[0]["fgmres_candidate_iteration"], "0")
            self.assertEqual(rows[-1]["fgmres_candidate_iteration"], "30")

    def test_scientific_rejection_is_not_invalid_experiment(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            offline = {"experiment_evidence_valid": True, "offline_gate_pass": False,
                       "gates": {"candidate_isolatable": True, "memory": True,
                                 "assembly": True, "local_RMS": True,
                                 "local_heavy_mode": True, "H_RAW": False}}
            (td/"offline.json").write_text(json.dumps(offline))
            F1B.final_audit(self.audit_args(td, td/"offline.json", td/"audit.json"))
            result = json.loads((td/"audit.json").read_text())
            self.assertTrue(result["experiment_valid"])
            self.assertFalse(result["candidate_valid"])
            self.assertEqual(result["decision"], "CANDIDATE_REJECTED_RAW_LOCAL")

    def test_provenance_drift_fails_closed(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            offline = {"experiment_evidence_valid": True, "offline_gate_pass": False,
                       "gates": {"candidate_isolatable": True, "memory": True,
                                 "assembly": True, "local_RMS": True,
                                 "local_heavy_mode": True, "H_RAW": False}}
            (td/"offline.json").write_text(json.dumps(offline))
            args = self.audit_args(td, td/"offline.json", td/"audit.json")
            Path(args.binary_post).write_text("drift\n")
            F1B.final_audit(args)
            result = json.loads((td/"audit.json").read_text())
            self.assertFalse(result["experiment_valid"])
            self.assertEqual(result["decision"], "INVALID_EXPERIMENT")
            self.assertEqual(result["next_authorized_task"],
                             "REPEAT_F1B_VALID_EXPERIMENT")


if __name__ == "__main__": unittest.main()
