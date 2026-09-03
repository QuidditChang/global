import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "f2g", ROOT / "tools/analyze_strict_ala_stage_F2g.py")
F2G = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(F2G)


def write_csv(path, fields, rows):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


class F2gContractTests(unittest.TestCase):
    def test_runtime_is_gated_and_collective_action_is_outside_root_guard(self):
        stage = (ROOT / "lib/Strict_ala_stage_f2g.inc").read_text()
        matrix = (ROOT / "lib/General_matrix_functions.c").read_text()
        driver = (ROOT / "lib/Stokes_flow_Incomp.c").read_text()
        self.assertIn('getenv("STRICT_ALA_STAGE_F2G_REQUIRED")', stage)
        self.assertIn('getenv("STRICT_ALA_STAGE_F2G_ACTIVE")', matrix)
        self.assertIn('strcmp(mode,"overlap_face_z")', matrix)
        self.assertIn("rhs_dot_Ku=global_vdot(E,rhs_saved,Ku,lev);", stage)
        self.assertIn("correction2=global_vdot(E,u,u,lev);", stage)
        self.assertIn('#include "Strict_ala_stage_f2g.inc"', driver)
        self.assertIn("strict_ala_stage_f2g_run(E,lev)", driver)
        collective = stage.index("correction2=global_vdot")
        root_write = stage.index("if(E->parallel.me==0) {", collective)
        self.assertLess(collective, root_write)

    def test_lsf_uses_dynamic_f2f_lineage(self):
        text = (ROOT.parents[1] /
                "runs/cmbhf_ALA_strict_stage_F2g.lsf").read_text()
        self.assertIn("STRICT_STAGE_F2G_${LSB_JOBID}", text)
        self.assertIn("STRICT_ALA_STAGE_F2G_F2F_ROOT", text)
        self.assertIn("no valid F2f experiment authorizing F2g", text)
        self.assertNotIn("STRICT_STAGE_F2F_12146335", text)
        self.assertIn("--nodes=384", text)

    def make_inputs(self, directory):
        root = Path(directory)
        (root / "f2f.json").write_text(json.dumps({
            "experiment_evidence_valid": True,
            "decision": "F2F_NO_ROBUST_OVERLAPPING_FACE_PATCH_CANDIDATE",
            "next_authorized_task": "F2G_CROSS_RANK_FACE_PATCH_REDESIGN",
            "production_default_change_authorized": False}))
        (root / "runtime.json").write_text(json.dumps({
            "schema": "strict-ala-stage-F2g-runtime-v1",
            "selected_modes": 2, "selected_energy": .97,
            "candidate_count": 4, "configured_finest_sweeps": 3,
            "configured_finest_damping": .2, "scope": "finest_level_only",
            "candidate_structure": "F2f_best_z_overlap_damping_calibration",
            "action_diagnostic": "first_cycle_M_inverse_r_and_K_gamma_action",
            "cross_rank_faces": False,
            "cross_rank_deferred_reason":
                "calibrate_rank_local_action_before_MPI_redesign",
            "viscosity": "frozen_current",
            "operator": "frozen_current_K_gamma_D_plus_C",
            "production_default_change_authorized": False}))
        trajectory_fields = (
            "candidate", "smoother_composition", "finest_sweeps",
            "finest_damping", "POD_mode", "POD_energy_weight",
            "MG_cycles", "momentum_residual_relative", "valid")
        factors = {"CONFIGURED": .90, "L4_OVERLAP_FACE_Z_D0P2": .98,
                   "L4_OVERLAP_FACE_Z_D0P5": .82,
                   "L4_OVERLAP_FACE_Z_D1P0": .85}
        damping = {"CONFIGURED": .2, **F2G.DAMPING}
        trajectory = [{
            "candidate": candidate,
            "smoother_composition": F2G.BLOCKS[candidate],
            "finest_sweeps": 3, "finest_damping": damping[candidate],
            "POD_mode": mode, "POD_energy_weight": .3,
            "MG_cycles": cycle,
            "momentum_residual_relative": factors[candidate] ** cycle,
            "valid": 1}
            for candidate in F2G.CANDIDATES for mode in (1, 2)
            for cycle in F2G.CHECKPOINTS]
        write_csv(root / "trajectory.csv", trajectory_fields, trajectory)
        prior = []
        for row in trajectory:
            if row["candidate"] == "CONFIGURED":
                prior.append(dict(row))
            elif row["candidate"] == "L4_OVERLAP_FACE_Z_D0P2":
                item = dict(row); item["candidate"] = "L4_OVERLAP_FACE_Z"
                prior.append(item)
        write_csv(root / "prior_trajectory.csv", trajectory_fields, prior)
        rhs_fields = ("POD_mode", "POD_energy_weight", "rhs_norm", "valid")
        rhs = [{"POD_mode": mode, "POD_energy_weight": .3,
                "rhs_norm": 1, "valid": 1} for mode in (1, 2)]
        write_csv(root / "rhs.csv", rhs_fields, rhs)
        write_csv(root / "prior_rhs.csv", rhs_fields, rhs)
        level_fields = (
            "candidate", "POD_mode", "POD_energy_weight", "MG_cycle",
            "top_level", "level", "phase", "input_rms", "output_rms",
            "reduction_ratio", "alpha", "valid")
        levels = [{
            "candidate": candidate, "POD_mode": mode,
            "POD_energy_weight": .3, "MG_cycle": cycle,
            "top_level": 4, "level": 4, "phase": "DOWN_SMOOTH",
            "input_rms": 1, "output_rms": .9,
            "reduction_ratio": .9, "alpha": 1, "valid": 1}
            for candidate in F2G.CANDIDATES for mode in (1, 2)
            for cycle in F2G.CHECKPOINTS]
        write_csv(root / "levels.csv", level_fields, levels)
        prior_levels = []
        for row in levels:
            if row["candidate"] == "CONFIGURED": prior_levels.append(dict(row))
            elif row["candidate"] == "L4_OVERLAP_FACE_Z_D0P2":
                item = dict(row); item["candidate"] = "L4_OVERLAP_FACE_Z"
                prior_levels.append(item)
        write_csv(root / "prior_levels.csv", level_fields, prior_levels)
        action_fields = (
            "candidate", "finest_damping", "POD_mode", "POD_energy_weight",
            "rhs_norm", "correction_norm", "rhs_dot_correction",
            "rhs_dot_Kcorrection", "Kcorrection_norm", "alpha_opt",
            "residual_ratio_alpha_1", "residual_ratio_alpha_opt",
            "positive_action", "valid")
        write_csv(root / "action.csv", action_fields, [{
            "candidate": candidate, "finest_damping": damping[candidate],
            "POD_mode": mode, "POD_energy_weight": .3, "rhs_norm": 1,
            "correction_norm": 1, "rhs_dot_correction": .4,
            "rhs_dot_Kcorrection": .5, "Kcorrection_norm": 1,
            "alpha_opt": .5, "residual_ratio_alpha_1": .9,
            "residual_ratio_alpha_opt": .7, "positive_action": 1, "valid": 1}
            for candidate in F2G.CANDIDATES for mode in (1, 2)])
        performance_fields = (
            "candidate", "smoother_block", "finest_sweeps",
            "finest_damping", "block_dof", "seconds_max", "valid")
        write_csv(root / "performance.csv", performance_fields, [{
            "candidate": candidate,
            "smoother_block": F2G.BLOCKS[candidate], "finest_sweeps": 3,
            "finest_damping": damping[candidate],
            "block_dof": 24 if candidate == "CONFIGURED" else 36,
            "seconds_max": 10, "valid": 1}
            for candidate in F2G.CANDIDATES])
        cache_fields = ("orientation", "patches_per_cap", "block_dof",
                        "cache_mb_max", "build_seconds_max",
                        "min_pivot_ratio", "valid")
        write_csv(root / "cache.csv", cache_fields, [{
            "orientation": "overlap_face_z", "patches_per_cap": 31744,
            "block_dof": 36, "cache_mb_max": 86,
            "build_seconds_max": 2, "min_pivot_ratio": 1e-4, "valid": 1}])
        return SimpleNamespace(
            f2f_decision=root / "f2f.json",
            f2f_trajectory=root / "prior_trajectory.csv",
            f2f_levels=root / "prior_levels.csv", f2f_rhs=root / "prior_rhs.csv",
            level_transfer=root / "levels.csv", trajectory=root / "trajectory.csv",
            mode_rhs=root / "rhs.csv", performance=root / "performance.csv",
            cache=root / "cache.csv", action=root / "action.csv",
            runtime=root / "runtime.json", output=root / "decision.json")

    def test_analyzer_selects_calibrated_candidate(self):
        with tempfile.TemporaryDirectory() as directory:
            args = self.make_inputs(directory)
            self.assertEqual(F2G.analyze(args), 0)
            result = json.loads(args.output.read_text())
            self.assertTrue(result["experiment_evidence_valid"])
            self.assertEqual(result["selected_candidate"],
                             "L4_OVERLAP_FACE_Z_D0P5")
            self.assertAlmostEqual(result["pod_weight_sum"], 0.6)
            self.assertTrue(result["normalized_weighted_rms"])
            self.assertTrue(result["normalized_weighted_mean"])
            self.assertAlmostEqual(result["metrics"]["CONFIGURED"]
                                   ["momentum_residual_RMS_64"], .9 ** 64)
            self.assertAlmostEqual(result["metrics"]["CONFIGURED"]
                                   ["alpha_opt_weighted_mean"], .5)
            self.assertFalse(result["production_default_change_authorized"])

    def test_nonpositive_candidate_is_scientifically_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            args = self.make_inputs(directory)
            with args.action.open() as stream:
                rows = list(csv.DictReader(stream))
            for row in rows:
                if row["candidate"] == "L4_OVERLAP_FACE_Z_D0P5":
                    row["positive_action"] = "0"
                    row["rhs_dot_correction"] = "-0.4"
            write_csv(args.action, rows[0].keys(), rows)
            F2G.analyze(args)
            result = json.loads(args.output.read_text())
            self.assertTrue(result["experiment_evidence_valid"])
            self.assertFalse(result["metrics"][
                "L4_OVERLAP_FACE_Z_D0P5"]["candidate_gate_pass"])


if __name__ == "__main__":
    unittest.main()
