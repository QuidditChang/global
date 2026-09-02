import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "f2e", ROOT / "tools/analyze_strict_ala_stage_F2e.py")
F2E = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(F2E)


def write_csv(path, fields, rows):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


class F2eContractTests(unittest.TestCase):
    def test_face_patch_is_bounded_finest_only_and_gated(self):
        matrix = (ROOT / "lib/General_matrix_functions.c").read_text()
        stage = (ROOT / "lib/Strict_ala_stage_f2e.inc").read_text()
        driver = (ROOT / "lib/Stokes_flow_Incomp.c").read_text()
        self.assertIn("#define ALA_F2E_FACE_DOF 36", matrix)
        self.assertIn("#define ALA_F2E_FACE_CHOL_SIZE 666", matrix)
        self.assertIn("level!=E->mesh.levmax", matrix)
        self.assertIn("ala_f2e_free_face_cache()", matrix)
        self.assertIn("one temporary cache exists at a time", stage)
        self.assertIn('getenv("STRICT_ALA_STAGE_F2E_REQUIRED")', stage)
        self.assertIn('#include "Strict_ala_stage_f2e.inc"', driver)
        self.assertIn("strict_ala_stage_f2e_run(E,lev)", driver)

    def test_lsf_uses_dynamic_f2d_lineage(self):
        text = (ROOT.parents[1] /
                "runs/cmbhf_ALA_strict_stage_F2e.lsf").read_text()
        self.assertIn("STRICT_STAGE_F2E_${LSB_JOBID}", text)
        self.assertIn("STRICT_ALA_STAGE_F2E_F2D_ROOT", text)
        self.assertIn("no valid F2d experiment authorizing F2e", text)
        self.assertNotIn("STRICT_STAGE_F2D_12145805", text)
        self.assertIn("--nodes=384", text)

    def make_inputs(self, directory):
        root = Path(directory)
        (root / "f2d.json").write_text(json.dumps({
            "experiment_evidence_valid": True,
            "decision": "F2D_NO_ROBUST_L4_SMOOTHER_CANDIDATE",
            "next_authorized_task": "F2E_FACE_PATCH_BLOCK_REDESIGN",
            "production_default_change_authorized": False}))
        (root / "runtime.json").write_text(json.dumps({
            "schema": "strict-ala-stage-F2e-runtime-v1",
            "selected_modes": 1, "selected_energy": .97,
            "candidate_count": 4, "configured_finest_sweeps": 3,
            "configured_finest_damping": .2, "scope": "finest_level_only",
            "candidate_structure": "nonoverlapping_two_element_face_pairs",
            "cached_blocks": "one_orientation_at_a_time_36x36_K_gamma",
            "face_patch_dof": 36, "viscosity": "frozen_current",
            "operator": "frozen_current_K_gamma_D_plus_C",
            "production_default_change_authorized": False}))
        tf = ("candidate", "smoother_composition", "finest_sweeps",
              "finest_damping", "POD_mode", "POD_energy_weight",
              "MG_cycles", "momentum_residual_relative", "valid")
        contraction = {"CONFIGURED": .90, "L4_FACE_Z": .88,
                       "L4_FACE_X": .87, "L4_FACE_Y": .80}
        trajectory = [{"candidate": candidate,
            "smoother_composition": F2E.BLOCKS[candidate],
            "finest_sweeps": 3, "finest_damping": .2, "POD_mode": 1,
            "POD_energy_weight": .97, "MG_cycles": cycle,
            "momentum_residual_relative": contraction[candidate] ** cycle,
            "valid": 1} for candidate in F2E.CANDIDATES
            for cycle in F2E.CHECKPOINTS]
        write_csv(root / "trajectory.csv", tf, trajectory)
        write_csv(root / "prior_trajectory.csv", tf,
                  [row for row in trajectory
                   if row["candidate"] == "CONFIGURED"])
        rf = ("POD_mode", "POD_energy_weight", "rhs_norm", "valid")
        rhs = [{"POD_mode": 1, "POD_energy_weight": .97,
                "rhs_norm": 1, "valid": 1}]
        write_csv(root / "rhs.csv", rf, rhs)
        write_csv(root / "prior_rhs.csv", rf, rhs)
        lf = ("candidate", "POD_mode", "POD_energy_weight", "MG_cycle",
              "top_level", "level", "phase", "input_rms", "output_rms",
              "reduction_ratio", "alpha", "valid")
        levels = [{"candidate": candidate, "POD_mode": 1,
            "POD_energy_weight": .97, "MG_cycle": cycle,
            "top_level": 4, "level": 4, "phase": "DOWN_SMOOTH",
            "input_rms": 1, "output_rms": .98, "reduction_ratio": .98,
            "alpha": 1, "valid": 1} for candidate in F2E.CANDIDATES
            for cycle in F2E.CHECKPOINTS]
        write_csv(root / "levels.csv", lf, levels)
        write_csv(root / "prior_levels.csv", lf,
                  [row for row in levels
                   if row["candidate"] == "CONFIGURED"])
        pf = ("candidate", "smoother_block", "finest_sweeps",
              "finest_damping", "block_dof", "seconds_max", "valid")
        write_csv(root / "performance.csv", pf, [{
            "candidate": candidate, "smoother_block": F2E.BLOCKS[candidate],
            "finest_sweeps": 3, "finest_damping": .2,
            "block_dof": 24 if candidate == "CONFIGURED" else 36,
            "seconds_max": 10, "valid": 1} for candidate in F2E.CANDIDATES])
        cf = ("orientation", "patches_per_cap", "block_dof",
              "cache_mb_max", "build_seconds_max", "min_pivot_ratio", "valid")
        write_csv(root / "cache.csv", cf, [{"orientation": orientation,
            "patches_per_cap": 256, "block_dof": 36, "cache_mb_max": 2,
            "build_seconds_max": 4, "min_pivot_ratio": 1e-8, "valid": 1}
            for orientation in ("face_z", "face_x", "face_y")])
        return SimpleNamespace(
            f2d_decision=root / "f2d.json",
            f2d_trajectory=root / "prior_trajectory.csv",
            f2d_levels=root / "prior_levels.csv", f2d_rhs=root / "prior_rhs.csv",
            level_transfer=root / "levels.csv", trajectory=root / "trajectory.csv",
            mode_rhs=root / "rhs.csv", performance=root / "performance.csv",
            cache=root / "cache.csv", runtime=root / "runtime.json",
            output=root / "decision.json")

    def test_analyzer_selects_improved_face_patch(self):
        with tempfile.TemporaryDirectory() as directory:
            args = self.make_inputs(directory)
            self.assertEqual(F2E.analyze(args), 0)
            result = json.loads(args.output.read_text())
            self.assertTrue(result["experiment_evidence_valid"])
            self.assertEqual(result["selected_candidate"], "L4_FACE_Y")
            self.assertFalse(result["production_default_change_authorized"])

    def test_analyzer_rejects_cache_over_64_mb(self):
        with tempfile.TemporaryDirectory() as directory:
            args = self.make_inputs(directory)
            with args.cache.open() as stream:
                rows = list(csv.DictReader(stream))
            rows[0]["cache_mb_max"] = "65"
            write_csv(args.cache, rows[0].keys(), rows)
            F2E.analyze(args)
            result = json.loads(args.output.read_text())
            self.assertFalse(result["experiment_evidence_valid"])


if __name__ == "__main__":
    unittest.main()
