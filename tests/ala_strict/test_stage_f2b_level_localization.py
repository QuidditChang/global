import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "f2b", ROOT / "tools/analyze_strict_ala_stage_F2b.py")
F2B = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(F2B)


def write_csv(path, fields, rows):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)


class F2bContractTests(unittest.TestCase):
    def test_observer_and_runtime_gate(self):
        mg = (ROOT / "lib/General_matrix_functions.c").read_text()
        stage = (ROOT / "lib/Strict_ala_stage_f2b.inc").read_text()
        driver = (ROOT / "lib/Stokes_flow_Incomp.c").read_text()
        for marker in ("DOWN_SMOOTH", "RESTRICTION", "BOTTOM_SOLVE",
                       "UP_CORRECTION"):
            self.assertIn(marker, mg)
        self.assertIn('getenv("STRICT_ALA_STAGE_F2B_REQUIRED")', stage)
        self.assertIn("STRICT_ALA_STAGE_F2B_OBSERVE", stage)
        self.assertIn("strict_ala_stage_f2b_run(E,lev)", driver)
        self.assertIn('#include "Strict_ala_stage_f2b.inc"', driver)

    def test_lsf_is_dynamic_resumable_and_skips_pressure_cache(self):
        text = (ROOT.parents[1] / "runs/cmbhf_ALA_strict_stage_F2b.lsf").read_text()
        self.assertIn("STRICT_STAGE_F2B_${LSB_JOBID}", text)
        self.assertIn("STRICT_ALA_STAGE_F2B_RESUME_ROOT", text)
        self.assertIn("STRICT_ALA_STAGE_F2B_F2_ROOT", text)
        self.assertIn("no valid F2 experiment authorizing F2b", text)
        self.assertNotIn("STRICT_STAGE_F2_12134701", text)
        self.assertIn("ala_shallow_patch_preconditioner off", text)
        self.assertIn("--nodes=384", text)

    def make_inputs(self, directory, worst_phase="UP_CORRECTION"):
        td = Path(directory)
        (td / "f2.json").write_text(json.dumps({
            "experiment_evidence_valid": True,
            "decision": "F2_VANKA_EFFECTIVE_BUT_AVV_CONVERGENCE_LIMITED",
            "next_authorized_task": "F2B_COARSE_GRID_TRANSFER_LOCALIZATION",
            "production_default_change_authorized": False}))
        (td / "runtime.json").write_text(json.dumps({
            "schema": "strict-ala-stage-F2b-runtime-v1", "selected_modes": 3,
            "selected_energy": 0.96, "maximum_MG_cycles": 64,
            "checkpoints": list(F2B.CHECKPOINTS),
            "operator": "frozen_current_K_gamma",
            "velocity_smoother": "configured_element_vanka",
            "production_default_change_authorized": False}))
        weights = (0.85, 0.09, 0.02)
        trajectory_fields = ("POD_mode", "POD_energy_weight", "MG_cycles",
                             "momentum_residual_relative", "valid")
        trajectory = []
        f2_fields = ("variant", "POD_mode", "POD_energy_weight", "MG_cycles",
                     "momentum_residual_relative", "velocity_error_relative",
                     "K_energy_error_relative", "Schur_action_error_relative",
                     "velocity_cosine", "exact_tight_achieved", "valid")
        f2_rows = []
        for mode, weight in enumerate(weights, 1):
            for cycle in F2B.CHECKPOINTS:
                value = 0.9 ** cycle
                trajectory.append({"POD_mode": mode, "POD_energy_weight": weight,
                    "MG_cycles": cycle, "momentum_residual_relative": value,
                    "valid": 1})
                f2_rows.append({"variant": "CONFIGURED_VANKA", "POD_mode": mode,
                    "POD_energy_weight": weight, "MG_cycles": cycle,
                    "momentum_residual_relative": value,
                    "velocity_error_relative": value,
                    "K_energy_error_relative": value,
                    "Schur_action_error_relative": value,
                    "velocity_cosine": 1.0, "exact_tight_achieved": 1e-10,
                    "valid": 1})
        write_csv(td / "trajectory.csv", trajectory_fields, trajectory)
        write_csv(td / "f2trajectory.csv", f2_fields, f2_rows)
        rhs_fields = ("POD_mode", "POD_energy_weight", "rhs_norm", "valid")
        write_csv(td / "rhs.csv", rhs_fields, [{"POD_mode": i,
            "POD_energy_weight": weights[i-1], "rhs_norm": i, "valid": 1}
            for i in range(1, 4)])
        tight_fields = ("call_id", "POD_mode", "rhs_norm",
            "requested_relative_tolerance", "target_residual",
            "achieved_relative_residual", "cycles", "max_cycles", "converged")
        write_csv(td / "tight.csv", tight_fields, [{"call_id": i,
            "POD_mode": i, "rhs_norm": i, "requested_relative_tolerance": 1e-10,
            "target_residual": 1e-12, "achieved_relative_residual": 9e-11,
            "cycles": 400, "max_cycles": 2000, "converged": 1}
            for i in range(1, 4)])
        level_rows = []
        ratios = {"DOWN_SMOOTH": 0.4, "RESTRICTION": 0.7,
                  "BOTTOM_SOLVE": 0.1, "UP_CORRECTION": 0.3}
        ratios[worst_phase] = 0.95
        for mode, weight in enumerate(weights, 1):
            for cycle in F2B.CHECKPOINTS:
                for phase, level in (("DOWN_SMOOTH", 4), ("RESTRICTION", 4),
                                     ("BOTTOM_SOLVE", 0),
                                     ("UP_CORRECTION", 4)):
                    level_rows.append({"POD_mode": mode,
                        "POD_energy_weight": weight, "MG_cycle": cycle,
                        "top_level": 4, "level": level, "phase": phase,
                        "input_rms": 1.0, "output_rms": ratios[phase],
                        "reduction_ratio": ratios[phase], "alpha": 1.0,
                        "valid": 1})
        write_csv(td / "levels.csv", F2B.LEVEL_COLUMNS, level_rows)
        return SimpleNamespace(f2_decision=td/"f2.json",
            f2_trajectory=td/"f2trajectory.csv", f2_tight=td/"tight.csv",
            level_transfer=td/"levels.csv", trajectory=td/"trajectory.csv",
            mode_rhs=td/"rhs.csv", runtime=td/"runtime.json",
            output=td/"decision.json")

    def test_analyzer_localizes_upward_transfer(self):
        with tempfile.TemporaryDirectory() as name:
            args = self.make_inputs(name, "UP_CORRECTION")
            self.assertEqual(F2B.analyze(args), 0)
            result = json.loads(args.output.read_text())
            self.assertTrue(result["experiment_evidence_valid"])
            self.assertEqual(result["decision"],
                             "F2B_UPWARD_TRANSFER_BOTTLENECK")
            self.assertEqual(result["next_authorized_task"],
                             "F2C_PROLONGATION_GALERKIN_REVIEW")

    def test_analyzer_rejects_changed_F2_trajectory(self):
        with tempfile.TemporaryDirectory() as name:
            args = self.make_inputs(name)
            with args.trajectory.open() as stream:
                rows = list(csv.DictReader(stream))
            rows[0]["momentum_residual_relative"] = "0.5"
            write_csv(args.trajectory, rows[0].keys(), rows)
            F2B.analyze(args)
            result = json.loads(args.output.read_text())
            self.assertFalse(result["experiment_evidence_valid"])
            self.assertEqual(result["decision"],
                             "F2B_INVALID_STRUCTURAL_EVIDENCE")


if __name__ == "__main__":
    unittest.main()
