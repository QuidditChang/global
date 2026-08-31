import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "f2", ROOT / "tools/analyze_strict_ala_stage_F2.py")
F2 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(F2)


def write_csv(path, fields, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


class F2ContractTests(unittest.TestCase):
    def test_runtime_gate_and_single_changed_object(self):
        source = (ROOT / "lib/Strict_ala_stage_f2.inc").read_text()
        driver = (ROOT / "lib/Stokes_flow_Incomp.c").read_text()
        self.assertIn('getenv("STRICT_ALA_STAGE_F2_REQUIRED")', source)
        self.assertIn("CONFIGURED_VANKA", source)
        self.assertIn("DIAGONAL_ROLLBACK", source)
        self.assertIn("ala_element_vanka_smoother=variant==0", source)
        self.assertIn(
            "E->control.ala_element_vanka_smoother=original_vanka;\n"
            "        q=ala_e2_build_mode", source)
        self.assertIn("E->control.ala_element_vanka_smoother=original_vanka", source)
        self.assertIn("strict_ala_stage_f2_run(E,lev)", driver)
        self.assertIn('#include "Strict_ala_stage_f2.inc"', driver)

    def test_submission_contract_is_dynamic_and_resumable(self):
        lsf = ROOT.parents[1] / "runs/cmbhf_ALA_strict_stage_F2.lsf"
        text = lsf.read_text()
        self.assertIn("STRICT_STAGE_F2_${LSB_JOBID}", text)
        self.assertIn("STRICT_ALA_STAGE_F2_RESUME_ROOT", text)
        self.assertIn("STRICT_ALA_STAGE_F2_F1D_ROOT", text)
        self.assertIn("no valid F1d rejection authorizing F2", text)
        self.assertNotIn("STRICT_STAGE_F1D_12133673", text)
        self.assertIn("--nodes=384", text)
        self.assertIn("rank_zero_log_closure", text)
        self.assertIn(
            "CitcomS.solver.vsolver ala_shallow_patch_preconditioner off",
            text)

    def make_inputs(self, directory, configured=0.10, rollback=0.20):
        directory = Path(directory)
        f1d = {
            "experiment_evidence_valid": True,
            "decision": "F1D_RESIDUAL_CLOSED_COARSE_REJECTED",
            "next_authorized_task": "F2_AVV_PRECONDITIONER_REVIEW",
            "production_default_change_authorized": False,
        }
        (directory / "f1d.json").write_text(json.dumps(f1d))
        runtime = {
            "schema": "strict-ala-stage-F2-runtime-v1",
            "selected_modes": 3,
            "selected_energy": 0.96,
            "variants": ["CONFIGURED_VANKA", "DIAGONAL_ROLLBACK"],
            "maximum_MG_cycles": 64,
            "tight_solve_count": 3,
            "operator": "frozen_current_K_gamma",
            "changed_object": "velocity_multigrid_smoother_only",
            "production_default_change_authorized": False,
        }
        (directory / "runtime.json").write_text(json.dumps(runtime))
        fields = (
            "variant", "POD_mode", "POD_energy_weight", "MG_cycles",
            "momentum_residual_relative", "velocity_error_relative",
            "K_energy_error_relative", "Schur_action_error_relative",
            "velocity_cosine", "exact_tight_achieved", "valid")
        rows = []
        weights = (0.85, 0.09, 0.02)
        for variant, final in (("CONFIGURED_VANKA", configured),
                               ("DIAGONAL_ROLLBACK", rollback)):
            for mode, weight in enumerate(weights, 1):
                for cycle in F2.CHECKPOINTS:
                    value = final ** (cycle / 64.0)
                    rows.append({
                        "variant": variant, "POD_mode": mode,
                        "POD_energy_weight": weight, "MG_cycles": cycle,
                        "momentum_residual_relative": value,
                        "velocity_error_relative": value,
                        "K_energy_error_relative": value,
                        "Schur_action_error_relative": value,
                        "velocity_cosine": 1.0 - 0.01 * value,
                        "exact_tight_achieved": 1e-11, "valid": 1})
        write_csv(directory / "trajectory.csv", fields, rows)
        tight_fields = (
            "call_id", "POD_mode", "rhs_norm",
            "requested_relative_tolerance", "target_residual",
            "achieved_relative_residual", "cycles", "max_cycles",
            "converged")
        write_csv(directory / "tight.csv", tight_fields, [{
            "call_id": mode, "POD_mode": mode, "rhs_norm": 1.0,
            "requested_relative_tolerance": 1e-10,
            "target_residual": 1e-12,
            "achieved_relative_residual": 9e-11,
            "cycles": 400, "max_cycles": 2000, "converged": 1}
            for mode in range(1, 4)])
        return SimpleNamespace(
            f1d_decision=directory / "f1d.json",
            trajectory=directory / "trajectory.csv",
            tight=directory / "tight.csv", runtime=directory / "runtime.json",
            output=directory / "decision.json")

    def test_analyzer_classifies_effective_vanka(self):
        with tempfile.TemporaryDirectory() as name:
            args = self.make_inputs(name, configured=0.10, rollback=0.20)
            self.assertEqual(F2.analyze(args), 0)
            result = json.loads(args.output.read_text())
            self.assertTrue(result["experiment_evidence_valid"])
            self.assertEqual(
                result["decision"],
                "F2_VANKA_EFFECTIVE_BUT_AVV_CONVERGENCE_LIMITED")
            self.assertEqual(
                result["next_authorized_task"],
                "F2B_COARSE_GRID_TRANSFER_LOCALIZATION")

    def test_analyzer_classifies_vanka_regression(self):
        with tempfile.TemporaryDirectory() as name:
            args = self.make_inputs(name, configured=0.30, rollback=0.10)
            F2.analyze(args)
            result = json.loads(args.output.read_text())
            self.assertEqual(result["decision"], "F2_VANKA_SMOOTHER_REGRESSION")
            self.assertEqual(result["next_authorized_task"],
                             "F2B_VANKA_SMOOTHER_REDESIGN")

    def test_analyzer_rejects_missing_checkpoint(self):
        with tempfile.TemporaryDirectory() as name:
            args = self.make_inputs(name)
            with args.trajectory.open() as stream:
                rows = list(csv.DictReader(stream))
                fields = rows[0].keys()
            write_csv(args.trajectory, fields, rows[:-1])
            with self.assertRaisesRegex(ValueError, "incomplete or duplicate"):
                F2.analyze(args)


if __name__ == "__main__":
    unittest.main()
