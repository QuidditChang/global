import importlib.util
import json
import pathlib
import tempfile
import types
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[2]
TOOL_PATH = ROOT / "tools" / "analyze_strict_ala_stage_VC1.py"
DRIVER = ROOT / "lib" / "Drive_solvers.c"
STOKES = ROOT / "lib" / "Stokes_flow_Incomp.c"
LSF = ROOT.parents[1] / "runs" / "cmbhf_ALA_strict_stage_VC1.lsf"
spec = importlib.util.spec_from_file_location("vc1", TOOL_PATH)
vc1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(vc1)


CFG = """[CitcomS]
steps = {steps}
[CitcomS.solver]
datadir = {datadir}
datadir_old = {old}
[CitcomS.solver.vsolver]
piterations = 220
ala_outer_solver = coupled_fgmres
ala_element_vanka_smoother = on
ala_coupled_element_vanka = on
[CitcomS.solver.phase]
phase_delta_rho = 180.0, 58.0, 459.0
phase_delta_s = -0.032902012641276054, -0.022663818742040032, 0.017520140281939885
[CitcomS.solver.visc]
SDEPV = off
PDEPV = off
visc_max = {vmax}
"""


class VC1ContinuationTests(unittest.TestCase):
    def test_runtime_cfg_adds_only_step_zero_phase_parser_compatibility(self):
        upstream_text = CFG.format(steps=30000, datadir="old", old="old", vmax=100)
        upstream_text = upstream_text.replace(
            "phase_delta_s = -0.032902012641276054, -0.022663818742040032, "
            "0.017520140281939885\n", "")
        with tempfile.TemporaryDirectory() as directory:
            directory = pathlib.Path(directory)
            upstream = directory / "upstream.cfg"
            runtime = directory / "runtime.cfg"
            audit = directory / "audit.json"
            upstream.write_text(upstream_text)
            vc1.prepare_runtime_cfg(types.SimpleNamespace(
                upstream=str(upstream), output=str(runtime), audit=str(audit),
                expected_upstream_sha256=vc1.digest(upstream)))
            value = json.loads(audit.read_text())
            self.assertTrue(value["valid"])
            self.assertEqual(value["changed_fields"],
                             ["CitcomS.solver.phase.phase_delta_s"])
            self.assertFalse(value["stokes_operator_change"])
            self.assertEqual(vc1.parse_cfg(runtime)[
                ("CitcomS.solver.phase", "phase_delta_s")],
                ", ".join(format(item, ".17g")
                          for item in vc1.PHASE_DELTA_S_COMPAT))

    def test_runtime_cfg_rejects_wrong_upstream_hash(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = pathlib.Path(directory)
            upstream = directory / "upstream.cfg"
            upstream.write_text("[CitcomS.solver.phase]\nphase_delta_rho = 1, 2, 3\n")
            with self.assertRaises(SystemExit):
                vc1.prepare_runtime_cfg(types.SimpleNamespace(
                    upstream=str(upstream), output=str(directory / "runtime.cfg"),
                    audit=str(directory / "audit.json"),
                    expected_upstream_sha256="0" * 64))

    def test_momentum_only_transfer_and_fresh_process_contract(self):
        source = DRIVER.read_text()
        self.assertIn('STRICT_ALA_VC1_WARM_INPUT', source)
        self.assertIn('STRICT_ALA_VC1_WARM_OUTPUT', source)
        self.assertIn('pressure_gauge=exact_no_regauge', source)
        transfer = source[source.index("static void strict_ala_vc1_transfer_warm_state"):
                          source.index("static void strict_ala_vc1_viscosity_audit")]
        self.assertIn("E->P[m]", transfer)
        self.assertIn("E->U[m]", transfer)
        self.assertNotIn("E->T[m]", transfer)
        self.assertNotIn("E->trace", transfer)
        self.assertLess(source.index("strict_ala_vc1_transfer_warm_state(E,vc1_warm_input,0)"),
                        source.index("assemble_forces(E,0)"))
        lsf = LSF.read_text()
        self.assertIn('CitcomS.SimpleApp:SimpleApp "${case_dir}/case.cfg"', lsf)
        self.assertNotIn("restart              = on", lsf)

    def test_frozen_state_and_cache_rebuild_order_are_guarded(self):
        source = DRIVER.read_text()
        self.assertIn("strict_ala_vc1_physical_state_hash", source)
        self.assertIn("VC1 frozen T/C/tracer/time state changed", source)
        warm = source.index("strict_ala_vc1_transfer_warm_state(E,vc1_warm_input,0)")
        viscosity = source.index("get_system_viscosity(E,1", warm)
        stiffness = source.index("construct_stiffness_B_matrix(E);", viscosity)
        solve = source.index("solve_constrained_flow_iterative(E);", stiffness)
        self.assertLess(warm, viscosity)
        self.assertLess(viscosity, stiffness)
        self.assertLess(stiffness, solve)

    def test_cfg_contract_changes_only_stage_fields(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = pathlib.Path(directory)
            canonical = directory / "canonical.cfg"
            candidate = directory / "candidate.cfg"
            output = directory / "contract.json"
            canonical.write_text(CFG.format(steps=30000, datadir="old", old="old", vmax=100))
            candidate.write_text(CFG.format(steps=1, datadir="new", old="restart", vmax=10))
            vc1.cfg_contract(types.SimpleNamespace(canonical=str(canonical), cfg=str(candidate),
                                                   visc_max=10.0, output=str(output)))
            value = json.loads(output.read_text())
            self.assertTrue(value["valid"])
            self.assertEqual(value["visc_max"], 10.0)

    def test_stage_parser_accepts_numerical_nonconvergence_without_warm_output(self):
        with tempfile.TemporaryDirectory() as directory:
            case = pathlib.Path(directory)
            rank0 = case / "DATA" / "0"
            rank0.mkdir(parents=True)
            (rank0 / "raw.log").write_text(
                "STRICT_ALA_VC1_VISCOSITY configured_visc_max=1.00000000000000000e+02 "
                "eta_min=1.00000000000000000e-01 eta_max=1.00000000000000000e+02 "
                "eta_ratio=1.00000000000000000e+03 lower_clamp_fraction=1.0e-01 "
                "upper_clamp_fraction=2.0e-01 sample_count=8 global_xor_checksum=abc\n"
                "STRICT_ALA_VC1_FROZEN_STATE_GUARD before=def solver_mutable_scope=U_P_only\n"
                "STRICT_ALA_VC1_TIMING operator_rebuild_seconds=1.0e+00 "
                "fgmres_and_preconditioner_seconds=0.0e+00\n"
                "ALA COUPLED FGMRES startup restart=50 block_norm=1 cancellation=1.0e+00 "
                "raw_momentum_relative=2.0e-04\n"
                "STRICT_ALA_VC1_COUPLED_COST K_gamma_rhs_solves=221 "
                "K_gamma_operator_applications=221 velocity_MG_cycles=1234 "
                "preconditioner_applications=220\n"
                "STRICT_ALA_STAGE_C_CASE_COMPLETE case=COLD status=NUMERICAL_NONCONVERGENCE "
                "iterations=220 continuity=1.938367e-2 momentum=1.303084e-5\n")
            output = case / "completion.json"
            vc1.stage_summary(types.SimpleNamespace(path_id="COLD_100", stage_index=0,
                visc_max=100.0, warm_source="COLD", case_dir=str(case), exit_status=1,
                wall_seconds=10.0, output=str(output)))
            value = json.loads(output.read_text())
            self.assertTrue(value["valid"])
            self.assertEqual(value["convergence_status"], "NUMERICAL_NONCONVERGENCE")

    def test_stage_parser_fails_closed_after_writing_infrastructure_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            case = pathlib.Path(directory)
            rank0 = case / "DATA" / "0"
            rank0.mkdir(parents=True)
            (rank0 / "raw.log").write_text("phase parser stopped before solve\n")
            output = case / "completion.json"
            with self.assertRaises(SystemExit):
                vc1.stage_summary(types.SimpleNamespace(
                    path_id="COLD_100", stage_index=0, visc_max=100.0,
                    warm_source="COLD", case_dir=str(case), exit_status=1,
                    wall_seconds=2.0, output=str(output)))
            value = json.loads(output.read_text())
            self.assertFalse(value["complete"])
            self.assertFalse(value["valid"])
            self.assertEqual(value["convergence_status"], "INFRASTRUCTURE_FAILURE")

    def test_path_sequences_and_independent_stop_logic(self):
        text = LSF.read_text()
        self.assertIn('run_path_stage DIRECT_WARM 1 100', text)
        self.assertIn('run_path_stage COARSE 1 10', text)
        self.assertIn('run_path_stage COARSE 2 100', text)
        self.assertIn('for fine_vmax in 2 3 5 10 30 100', text)
        self.assertIn("break", text)
        self.assertIn("|| true", text)
        self.assertIn("prepare-runtime-cfg", text)
        self.assertIn("EXPECTED_UPSTREAM_CFG_SHA256", text)
        self.assertIn('"${UPSTREAM_CANONICAL}" "${CANONICAL}"', text)

    def test_continuation_path_completion_accepts_early_numerical_stop(self):
        rows = [{"path_id": "FINE", "continuation_stage_index": 1,
                 "visc_max": 2.0,
                 "convergence_status": "NUMERICAL_NONCONVERGENCE"}]
        self.assertTrue(vc1.continuation_path_complete(
            rows, "FINE", (2, 3, 5, 10, 30, 100)))
        rows[0]["convergence_status"] = "CONVERGED"
        self.assertFalse(vc1.continuation_path_complete(
            rows, "FINE", (2, 3, 5, 10, 30, 100)))

    def test_aggregate_marks_incomplete_invalid_run_nonzero(self):
        def failed_stage(path_id, visc_max):
            return {
                "path_id": path_id, "continuation_stage_index": 0,
                "visc_max": visc_max, "warm_start_source": "COLD",
                "convergence_status": "INFRASTRUCTURE_FAILURE",
                "FGMRES_iterations": 0, "final_R_cont": None,
                "final_R_mom": None, "K_gamma_solve_count": 0,
                "total_MG_cycles": 0, "operator_rebuild_time_seconds": None,
                "preconditioner_cache_rebuild_time_seconds": 0.0,
                "FGMRES_solve_time_seconds": 0.0,
                "total_stage_wall_time_seconds": 1.0,
                "R_cont_trajectory": [], "R_mom_trajectory": [],
                "physical_state_fingerprint": None, "valid": False,
                "complete": False,
            }

        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            for name, viscosity in (("COLD_100", 100.0), ("SEED_1", 1.0)):
                stage = root / name
                stage.mkdir()
                (stage / "completion.json").write_text(
                    json.dumps(failed_stage(name, viscosity)))
            final = root / "FINAL"
            with self.assertRaises(SystemExit):
                vc1.aggregate(types.SimpleNamespace(root=str(root),
                                                   output_dir=str(final)))
            audit = json.loads((final / "strict_ala_stage_VC1_final_audit.json").read_text())
            self.assertFalse(audit["valid"])
            self.assertFalse(audit["complete"])
            self.assertEqual(audit["exit_status"], 1)

    def test_coupled_solver_exports_required_cost_counters(self):
        text = STOKES.read_text()
        self.assertIn("STRICT_ALA_VC1_COUPLED_COST", text)
        self.assertIn("K_gamma_rhs_solves=%lld", text)
        self.assertIn("velocity_MG_cycles=%lld", text)


if __name__ == "__main__":
    unittest.main()
