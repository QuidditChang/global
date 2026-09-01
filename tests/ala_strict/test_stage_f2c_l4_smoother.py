import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "f2c", ROOT / "tools/analyze_strict_ala_stage_F2c.py")
F2C = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(F2C)


def write_csv(path, fields, rows):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)


class F2cContractTests(unittest.TestCase):
    def test_runtime_override_is_finest_only_and_production_gated(self):
        mg = (ROOT / "lib/General_matrix_functions.c").read_text()
        stage = (ROOT / "lib/Strict_ala_stage_f2c.inc").read_text()
        driver = (ROOT / "lib/Stokes_flow_Incomp.c").read_text()
        self.assertIn("level!=E->mesh.levmax", mg)
        self.assertIn("STRICT_ALA_STAGE_F2C_FINEST_SWEEPS", mg)
        self.assertIn("STRICT_ALA_STAGE_F2C_FINEST_DAMPING", mg)
        self.assertIn('getenv("STRICT_ALA_STAGE_F2C_REQUIRED")', stage)
        self.assertIn("strict_ala_stage_f2c_run(E,lev)", driver)
        self.assertIn('#include "Strict_ala_stage_f2c.inc"', driver)

    def test_lsf_contract_is_dynamic_and_resumable(self):
        text = (ROOT.parents[1] / "runs/cmbhf_ALA_strict_stage_F2c.lsf").read_text()
        self.assertIn("STRICT_STAGE_F2C_${LSB_JOBID}", text)
        self.assertIn("STRICT_ALA_STAGE_F2C_RESUME_ROOT", text)
        self.assertIn("STRICT_ALA_STAGE_F2C_F2B_ROOT", text)
        self.assertIn("no valid F2b experiment authorizing F2c", text)
        self.assertNotIn("STRICT_STAGE_F2B_12137549", text)
        self.assertIn("--nodes=384", text)

    def make_inputs(self, directory, improve=True):
        td = Path(directory)
        (td / "f2b.json").write_text(json.dumps({
            "experiment_evidence_valid": True,
            "decision": "F2B_LEVEL_SMOOTHER_BOTTLENECK",
            "next_authorized_task": "F2C_LEVEL_SMOOTHER_REDESIGN",
            "production_default_change_authorized": False}))
        (td / "runtime.json").write_text(json.dumps({
            "schema": "strict-ala-stage-F2c-runtime-v1", "selected_modes": 3,
            "selected_energy": 0.97, "candidate_count": 3,
            "scope": "finest_level_only", "viscosity": "frozen_current",
            "operator": "frozen_current_K_gamma_D_plus_C",
            "production_default_change_authorized": False}))
        weights = (0.85, 0.10, 0.02)
        trajectory_fields = ("candidate", "finest_sweeps", "finest_damping",
            "POD_mode", "POD_energy_weight", "MG_cycles",
            "momentum_residual_relative", "valid")
        prior_fields = ("POD_mode", "POD_energy_weight", "MG_cycles",
                        "momentum_residual_relative", "valid")
        trajectory, prior = [], []
        factors = {"CONFIGURED": 1.0, "L4_DOUBLE_SWEEPS": 0.80 if improve else 0.95,
                   "L4_DOUBLE_SWEEPS_DAMPING_1P5": 0.70 if improve else 1.1}
        for mode, weight in enumerate(weights, 1):
            for cycle in F2C.CHECKPOINTS:
                base = 0.9 ** cycle
                prior.append({"POD_mode": mode, "POD_energy_weight": weight,
                    "MG_cycles": cycle, "momentum_residual_relative": base,
                    "valid": 1})
                for candidate in F2C.CANDIDATES:
                    trajectory.append({"candidate": candidate,
                        "finest_sweeps": 1 if candidate == "CONFIGURED" else 2,
                        "finest_damping": .2 if candidate !=
                            "L4_DOUBLE_SWEEPS_DAMPING_1P5" else .3,
                        "POD_mode": mode, "POD_energy_weight": weight,
                        "MG_cycles": cycle,
                        "momentum_residual_relative": base*factors[candidate],
                        "valid": 1})
        write_csv(td/"trajectory.csv", trajectory_fields, trajectory)
        write_csv(td/"prior.csv", prior_fields, prior)
        rhs_fields = ("POD_mode", "POD_energy_weight", "rhs_norm", "valid")
        rhs = [{"POD_mode": mode, "POD_energy_weight": weight,
                "rhs_norm": mode, "valid": 1}
               for mode, weight in enumerate(weights, 1)]
        write_csv(td/"rhs.csv", rhs_fields, rhs)
        write_csv(td/"prior_rhs.csv", rhs_fields, rhs)
        level_fields = ("candidate", "POD_mode", "POD_energy_weight",
            "MG_cycle", "top_level", "level", "phase", "input_rms",
            "output_rms", "reduction_ratio", "alpha", "valid")
        old_level_fields = level_fields[1:]
        levels, old_levels = [], []
        for mode, weight in enumerate(weights, 1):
            for cycle in F2C.CHECKPOINTS:
                old_levels.append({"POD_mode": mode, "POD_energy_weight": weight,
                    "MG_cycle": cycle, "top_level": 4, "level": 4,
                    "phase": "DOWN_SMOOTH", "input_rms": 1,
                    "output_rms": .98, "reduction_ratio": .98,
                    "alpha": 1, "valid": 1})
                ratios = {"CONFIGURED": .98, "L4_DOUBLE_SWEEPS": .90,
                          "L4_DOUBLE_SWEEPS_DAMPING_1P5": .82}
                for candidate in F2C.CANDIDATES:
                    levels.append({"candidate": candidate, "POD_mode": mode,
                        "POD_energy_weight": weight, "MG_cycle": cycle,
                        "top_level": 4, "level": 4, "phase": "DOWN_SMOOTH",
                        "input_rms": 1, "output_rms": ratios[candidate],
                        "reduction_ratio": ratios[candidate], "alpha": 1,
                        "valid": 1})
        write_csv(td/"levels.csv", level_fields, levels)
        write_csv(td/"prior_levels.csv", old_level_fields, old_levels)
        return SimpleNamespace(f2b_decision=td/"f2b.json",
            f2b_trajectory=td/"prior.csv", f2b_levels=td/"prior_levels.csv",
            f2b_rhs=td/"prior_rhs.csv", level_transfer=td/"levels.csv",
            trajectory=td/"trajectory.csv", mode_rhs=td/"rhs.csv",
            runtime=td/"runtime.json", output=td/"decision.json")

    def test_analyzer_selects_improved_candidate(self):
        with tempfile.TemporaryDirectory() as name:
            args = self.make_inputs(name)
            self.assertEqual(F2C.analyze(args), 0)
            result = json.loads(args.output.read_text())
            self.assertTrue(result["experiment_evidence_valid"])
            self.assertEqual(result["decision"],
                             "F2C_L4_SMOOTHER_CANDIDATE_IDENTIFIED")
            self.assertEqual(result["selected_candidate"],
                             "L4_DOUBLE_SWEEPS_DAMPING_1P5")
            self.assertFalse(result["production_default_change_authorized"])

    def test_analyzer_rejects_baseline_lineage_change(self):
        with tempfile.TemporaryDirectory() as name:
            args = self.make_inputs(name)
            with args.trajectory.open() as stream:
                rows = list(csv.DictReader(stream))
            rows[0]["momentum_residual_relative"] = "0.123"
            write_csv(args.trajectory, rows[0].keys(), rows)
            F2C.analyze(args)
            result = json.loads(args.output.read_text())
            self.assertFalse(result["experiment_evidence_valid"])
            self.assertEqual(result["decision"],
                             "F2C_INVALID_STRUCTURAL_EVIDENCE")


if __name__ == "__main__":
    unittest.main()
