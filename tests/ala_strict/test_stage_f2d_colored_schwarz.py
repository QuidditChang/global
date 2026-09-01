import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "f2d", ROOT / "tools/analyze_strict_ala_stage_F2d.py")
F2D = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(F2D)


def write_csv(path, fields, rows):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)


class F2dContractTests(unittest.TestCase):
    def test_colored_composition_is_finest_only_and_production_gated(self):
        mg = (ROOT / "lib/General_matrix_functions.c").read_text()
        stage = (ROOT / "lib/Strict_ala_stage_f2d.inc").read_text()
        driver = (ROOT / "lib/Stokes_flow_Incomp.c").read_text()
        self.assertIn("level!=E->mesh.levmax", mg)
        self.assertIn("STRICT_ALA_STAGE_F2D_L4_MODE", mg)
        self.assertIn("colored_alternating", mg)
        self.assertIn("n_assemble_del2_u(E,d0,Ad,level,1)", mg)
        self.assertIn('getenv("STRICT_ALA_STAGE_F2D_REQUIRED")', stage)
        self.assertIn("strict_ala_stage_f2d_run(E,lev)", driver)
        self.assertIn('#include "Strict_ala_stage_f2d.inc"', driver)

    def test_lsf_contract_is_dynamic_and_resumable(self):
        text = (ROOT.parents[1] / "runs/cmbhf_ALA_strict_stage_F2d.lsf").read_text()
        self.assertIn("STRICT_STAGE_F2D_${LSB_JOBID}", text)
        self.assertIn("STRICT_ALA_STAGE_F2D_RESUME_ROOT", text)
        self.assertIn("STRICT_ALA_STAGE_F2D_F2C_ROOT", text)
        self.assertIn("no valid F2c experiment authorizing F2d", text)
        self.assertNotIn("STRICT_STAGE_F2C_12141823", text)
        self.assertIn("--nodes=384", text)

    def make_inputs(self, directory, improve=True):
        td = Path(directory)
        (td / "f2c.json").write_text(json.dumps({
            "experiment_evidence_valid": True,
            "decision": "F2C_NO_ROBUST_L4_SMOOTHER_CANDIDATE",
            "next_authorized_task": "F2C_BLOCK_STRUCTURE_REDESIGN",
            "production_default_change_authorized": False}))
        (td / "runtime.json").write_text(json.dumps({
            "schema": "strict-ala-stage-F2d-runtime-v1", "selected_modes": 3,
            "selected_energy": 0.97, "candidate_count": 3,
            "configured_finest_sweeps": 3,
            "configured_finest_damping": .2,
            "scope": "finest_level_only", "viscosity": "frozen_current",
            "operator": "frozen_current_K_gamma_D_plus_C",
            "candidate_structure": "eight_color_multiplicative_residual_update",
            "cached_blocks": "unchanged_24x24_element_K_gamma", "colors": 8,
            "production_default_change_authorized": False}))
        weights = (0.85, 0.10, 0.02)
        trajectory_fields = ("candidate", "smoother_composition",
            "finest_sweeps", "finest_damping", "POD_mode",
            "POD_energy_weight", "MG_cycles",
            "momentum_residual_relative", "valid")
        prior_fields = ("candidate", "finest_sweeps", "finest_damping",
            "POD_mode", "POD_energy_weight", "MG_cycles",
            "momentum_residual_relative", "valid")
        trajectory, prior = [], []
        contractions = {"CONFIGURED": .90,
                        "L4_COLORED_FORWARD": .86 if improve else .91,
                        "L4_COLORED_ALTERNATING": .80 if improve else .92}
        compositions = {"CONFIGURED": "configured",
                        "L4_COLORED_FORWARD": "colored_forward",
                        "L4_COLORED_ALTERNATING": "colored_alternating"}
        for mode, weight in enumerate(weights, 1):
            for cycle in F2D.CHECKPOINTS:
                base = 0.9 ** cycle
                prior.append({"candidate": "CONFIGURED", "finest_sweeps": 3,
                    "finest_damping": .2, "POD_mode": mode,
                    "POD_energy_weight": weight,
                    "MG_cycles": cycle, "momentum_residual_relative": base,
                    "valid": 1})
                for candidate in F2D.CANDIDATES:
                    trajectory.append({"candidate": candidate,
                        "smoother_composition": compositions[candidate],
                        "finest_sweeps": 3, "finest_damping": .2,
                        "POD_mode": mode, "POD_energy_weight": weight,
                        "MG_cycles": cycle,
                        "momentum_residual_relative": contractions[candidate]**cycle,
                        "valid": 1})
        write_csv(td/"trajectory.csv", trajectory_fields, trajectory)
        write_csv(td/"prior.csv", prior_fields, prior)
        rhs_fields = ("POD_mode", "POD_energy_weight", "rhs_norm", "valid")
        rhs = [{"POD_mode": mode, "POD_energy_weight": weight,
                "rhs_norm": mode, "valid": 1}
               for mode, weight in enumerate(weights, 1)]
        write_csv(td/"rhs.csv", rhs_fields, rhs)
        write_csv(td/"prior_rhs.csv", rhs_fields, rhs)
        performance_fields = ("candidate", "smoother_composition",
            "finest_sweeps", "finest_damping", "colors_per_sweep",
            "seconds_max", "valid")
        write_csv(td/"performance.csv", performance_fields, [{
            "candidate": candidate,
            "smoother_composition": compositions[candidate],
            "finest_sweeps": 3, "finest_damping": .2,
            "colors_per_sweep": 1 if candidate == "CONFIGURED" else 8,
            "seconds_max": 10 if candidate == "CONFIGURED" else 60,
            "valid": 1} for candidate in F2D.CANDIDATES])
        level_fields = ("candidate", "POD_mode", "POD_energy_weight",
            "MG_cycle", "top_level", "level", "phase", "input_rms",
            "output_rms", "reduction_ratio", "alpha", "valid")
        old_level_fields = level_fields
        levels, old_levels = [], []
        for mode, weight in enumerate(weights, 1):
            for cycle in F2D.CHECKPOINTS:
                old_levels.append({"candidate": "CONFIGURED",
                    "POD_mode": mode, "POD_energy_weight": weight,
                    "MG_cycle": cycle, "top_level": 4, "level": 4,
                    "phase": "DOWN_SMOOTH", "input_rms": 1,
                    "output_rms": .98, "reduction_ratio": .98,
                    "alpha": 1, "valid": 1})
                ratios = {"CONFIGURED": .98, "L4_COLORED_FORWARD": .90,
                          "L4_COLORED_ALTERNATING": .82}
                for candidate in F2D.CANDIDATES:
                    levels.append({"candidate": candidate, "POD_mode": mode,
                        "POD_energy_weight": weight, "MG_cycle": cycle,
                        "top_level": 4, "level": 4, "phase": "DOWN_SMOOTH",
                        "input_rms": 1, "output_rms": ratios[candidate],
                        "reduction_ratio": ratios[candidate], "alpha": 1,
                        "valid": 1})
        write_csv(td/"levels.csv", level_fields, levels)
        write_csv(td/"prior_levels.csv", old_level_fields, old_levels)
        return SimpleNamespace(f2c_decision=td/"f2c.json",
            f2c_trajectory=td/"prior.csv", f2c_levels=td/"prior_levels.csv",
            f2c_rhs=td/"prior_rhs.csv", level_transfer=td/"levels.csv",
            trajectory=td/"trajectory.csv", mode_rhs=td/"rhs.csv",
            performance=td/"performance.csv",
            runtime=td/"runtime.json", output=td/"decision.json")

    def test_analyzer_selects_improved_candidate(self):
        with tempfile.TemporaryDirectory() as name:
            args = self.make_inputs(name)
            self.assertEqual(F2D.analyze(args), 0)
            result = json.loads(args.output.read_text())
            self.assertTrue(result["experiment_evidence_valid"])
            self.assertEqual(result["decision"],
                             "F2D_L4_SMOOTHER_CANDIDATE_IDENTIFIED")
            self.assertEqual(result["selected_candidate"],
                             "L4_COLORED_ALTERNATING")
            self.assertFalse(result["production_default_change_authorized"])

    def test_analyzer_rejects_baseline_lineage_change(self):
        with tempfile.TemporaryDirectory() as name:
            args = self.make_inputs(name)
            with args.trajectory.open() as stream:
                rows = list(csv.DictReader(stream))
            rows[0]["momentum_residual_relative"] = "0.123"
            write_csv(args.trajectory, rows[0].keys(), rows)
            F2D.analyze(args)
            result = json.loads(args.output.read_text())
            self.assertFalse(result["experiment_evidence_valid"])
            self.assertEqual(result["decision"],
                             "F2D_INVALID_STRUCTURAL_EVIDENCE")


if __name__ == "__main__":
    unittest.main()
