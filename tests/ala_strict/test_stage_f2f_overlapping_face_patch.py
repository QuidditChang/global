import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "f2f", ROOT / "tools/analyze_strict_ala_stage_F2f.py")
F2F = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(F2F)


def write_csv(path, fields, rows):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


class F2fContractTests(unittest.TestCase):
    def test_each_rank_local_internal_face_is_enumerated_once(self):
        nx, ny, nz = 4, 3, 5
        dimensions = {1: nz, 2: nx, 3: ny}
        for orientation, paired in dimensions.items():
            patches = nx * ny * nz * (paired - 1) // paired
            pairs = []
            for patch in range(patches):
                if orientation == 1:
                    z = patch % (nz - 1) + 1
                    q = patch // (nz - 1)
                    x, y = q % nx + 1, q // nx + 1
                    second = (x, y, z + 1)
                elif orientation == 2:
                    z = patch % nz + 1
                    q = patch // nz
                    x, y = q % (nx - 1) + 1, q // (nx - 1) + 1
                    second = (x + 1, y, z)
                else:
                    z = patch % nz + 1
                    q = patch // nz
                    x, y = q % nx + 1, q // nx + 1
                    second = (x, y + 1, z)
                pairs.append(((x, y, z), second))
            self.assertEqual(len(pairs), len(set(pairs)))
            self.assertEqual(len(pairs), patches)
            self.assertTrue(all(
                1 <= x <= nx and 1 <= y <= ny and 1 <= z <= nz
                for pair in pairs for x, y, z in pair))

    def test_overlapping_patch_is_finest_only_bounded_and_gated(self):
        matrix = (ROOT / "lib/General_matrix_functions.c").read_text()
        stage = (ROOT / "lib/Strict_ala_stage_f2f.inc").read_text()
        driver = (ROOT / "lib/Stokes_flow_Incomp.c").read_text()
        self.assertIn('getenv("STRICT_ALA_STAGE_F2F_ACTIVE")', matrix)
        self.assertIn('strcmp(mode,"overlap_face_z")', matrix)
        self.assertIn('strcmp(mode,"overlap_face_x")', matrix)
        self.assertIn('strcmp(mode,"overlap_face_y")', matrix)
        self.assertIn("*(paired_dimension-1)/paired_dimension", matrix)
        self.assertIn("one cache exists at a time", stage)
        self.assertIn(r'\"cross_rank_faces\": false', stage)
        self.assertIn('#include "Strict_ala_stage_f2f.inc"', driver)
        self.assertIn("strict_ala_stage_f2f_run(E,lev)", driver)

    def test_lsf_uses_dynamic_f2e_lineage(self):
        text = (ROOT.parents[1] /
                "runs/cmbhf_ALA_strict_stage_F2f.lsf").read_text()
        self.assertIn("STRICT_STAGE_F2F_${LSB_JOBID}", text)
        self.assertIn("STRICT_ALA_STAGE_F2F_F2E_ROOT", text)
        self.assertIn("no valid F2e experiment authorizing F2f", text)
        self.assertNotIn("STRICT_STAGE_F2E_12146055", text)
        self.assertIn("--nodes=384", text)

    def make_inputs(self, directory):
        root = Path(directory)
        (root / "f2e.json").write_text(json.dumps({
            "experiment_evidence_valid": True,
            "decision": "F2E_NO_ROBUST_FACE_PATCH_CANDIDATE",
            "next_authorized_task":
                "F2F_OVERLAPPING_OR_CROSS_RANK_BLOCK_REDESIGN",
            "production_default_change_authorized": False}))
        (root / "runtime.json").write_text(json.dumps({
            "schema": "strict-ala-stage-F2f-runtime-v1",
            "selected_modes": 2, "selected_energy": .97,
            "candidate_count": 4, "configured_finest_sweeps": 3,
            "configured_finest_damping": .2, "scope": "finest_level_only",
            "candidate_structure":
                "overlapping_two_element_face_pairs_all_local_internal_faces",
            "cached_blocks": "one_orientation_at_a_time_36x36_K_gamma",
            "face_patch_dof": 36, "cross_rank_faces": False,
            "viscosity": "frozen_current",
            "operator": "frozen_current_K_gamma_D_plus_C",
            "production_default_change_authorized": False}))
        trajectory_fields = (
            "candidate", "smoother_composition", "finest_sweeps",
            "finest_damping", "POD_mode", "POD_energy_weight",
            "MG_cycles", "momentum_residual_relative", "valid")
        factors = {"CONFIGURED": .90, "L4_OVERLAP_FACE_Z": .86,
                   "L4_OVERLAP_FACE_X": .84,
                   "L4_OVERLAP_FACE_Y": .80}
        trajectory = [{
            "candidate": candidate,
            "smoother_composition": F2F.BLOCKS[candidate],
            "finest_sweeps": 3, "finest_damping": .2,
            "POD_mode": mode, "POD_energy_weight": .5,
            "MG_cycles": cycle,
            "momentum_residual_relative": factors[candidate] ** cycle,
            "valid": 1}
            for candidate in F2F.CANDIDATES for mode in (1, 2)
            for cycle in F2F.CHECKPOINTS]
        write_csv(root / "trajectory.csv", trajectory_fields, trajectory)
        write_csv(root / "prior_trajectory.csv", trajectory_fields,
                  [row for row in trajectory
                   if row["candidate"] == "CONFIGURED"])
        rhs_fields = ("POD_mode", "POD_energy_weight", "rhs_norm", "valid")
        rhs = [{"POD_mode": mode, "POD_energy_weight": .5,
                "rhs_norm": 1, "valid": 1} for mode in (1, 2)]
        write_csv(root / "rhs.csv", rhs_fields, rhs)
        write_csv(root / "prior_rhs.csv", rhs_fields, rhs)
        level_fields = (
            "candidate", "POD_mode", "POD_energy_weight", "MG_cycle",
            "top_level", "level", "phase", "input_rms", "output_rms",
            "reduction_ratio", "alpha", "valid")
        levels = [{
            "candidate": candidate, "POD_mode": mode,
            "POD_energy_weight": .5, "MG_cycle": cycle,
            "top_level": 4, "level": 4, "phase": "DOWN_SMOOTH",
            "input_rms": 1, "output_rms": .95,
            "reduction_ratio": .95, "alpha": 1, "valid": 1}
            for candidate in F2F.CANDIDATES for mode in (1, 2)
            for cycle in F2F.CHECKPOINTS]
        write_csv(root / "levels.csv", level_fields, levels)
        write_csv(root / "prior_levels.csv", level_fields,
                  [row for row in levels
                   if row["candidate"] == "CONFIGURED"])
        performance_fields = (
            "candidate", "smoother_block", "finest_sweeps",
            "finest_damping", "block_dof", "seconds_max", "valid")
        write_csv(root / "performance.csv", performance_fields, [{
            "candidate": candidate,
            "smoother_block": F2F.BLOCKS[candidate],
            "finest_sweeps": 3, "finest_damping": .2,
            "block_dof": 24 if candidate == "CONFIGURED" else 36,
            "seconds_max": 10, "valid": 1}
            for candidate in F2F.CANDIDATES])
        cache_fields = (
            "orientation", "patches_per_cap", "block_dof", "cache_mb_max",
            "build_seconds_max", "min_pivot_ratio", "valid")
        write_csv(root / "cache.csv", cache_fields, [{
            "orientation": orientation, "patches_per_cap": 31744,
            "block_dof": 36, "cache_mb_max": 87,
            "build_seconds_max": 2, "min_pivot_ratio": 1e-4, "valid": 1}
            for orientation in
            ("overlap_face_z", "overlap_face_x", "overlap_face_y")])
        return SimpleNamespace(
            f2e_decision=root / "f2e.json",
            f2e_trajectory=root / "prior_trajectory.csv",
            f2e_levels=root / "prior_levels.csv",
            f2e_rhs=root / "prior_rhs.csv",
            level_transfer=root / "levels.csv",
            trajectory=root / "trajectory.csv", mode_rhs=root / "rhs.csv",
            performance=root / "performance.csv", cache=root / "cache.csv",
            runtime=root / "runtime.json", output=root / "decision.json")

    def test_analyzer_selects_improved_overlapping_patch(self):
        with tempfile.TemporaryDirectory() as directory:
            args = self.make_inputs(directory)
            self.assertEqual(F2F.analyze(args), 0)
            result = json.loads(args.output.read_text())
            self.assertTrue(result["experiment_evidence_valid"])
            self.assertEqual(result["selected_candidate"],
                             "L4_OVERLAP_FACE_Y")
            self.assertFalse(result["production_default_change_authorized"])

    def test_analyzer_rejects_cache_over_128_mb(self):
        with tempfile.TemporaryDirectory() as directory:
            args = self.make_inputs(directory)
            with args.cache.open() as stream:
                rows = list(csv.DictReader(stream))
            rows[0]["cache_mb_max"] = "129"
            write_csv(args.cache, rows[0].keys(), rows)
            F2F.analyze(args)
            result = json.loads(args.output.read_text())
            self.assertFalse(result["experiment_evidence_valid"])


if __name__ == "__main__":
    unittest.main()
