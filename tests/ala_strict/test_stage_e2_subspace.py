import importlib.util
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
TOOL = ROOT / "tools" / "analyze_strict_ala_stage_E2.py"
SPEC = importlib.util.spec_from_file_location("stage_e2", TOOL)
stage_e2 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(stage_e2)


class PodTests(unittest.TestCase):
    def test_repeated_snapshots_rank_one(self):
        result = stage_e2.pod_metrics(np.ones((8, 5)), normalized=True)
        self.assertEqual(result["rank"], 1)
        self.assertAlmostEqual(result["fractions"][0], 1.0)

    def test_orthogonal_snapshots(self):
        result = stage_e2.pod_metrics(np.eye(5))
        np.testing.assert_allclose(result["gram"], np.eye(5))
        np.testing.assert_allclose(result["modes"].T @ result["modes"], np.eye(5))

    def test_rank_deficient_reconstruction(self):
        matrix = np.array([[1., 2., 3.], [2., 4., 6.], [0., 1., 1.]])
        result = stage_e2.pod_metrics(matrix)
        reconstructed = result["modes"] @ np.diag(np.sqrt(
            result["values"][:result["rank"]])) @ result["vectors"][:, :result["rank"]].T
        np.testing.assert_allclose(reconstructed, matrix, atol=1e-12)

    def test_principal_angles(self):
        a = np.eye(4)[:, :2]
        b = np.column_stack((a[:, 0], np.eye(4)[:, 2]))
        np.testing.assert_allclose(stage_e2.principal_angles(a, b), [0., 90.], atol=1e-7)


class SpatialTests(unittest.TestCase):
    def test_nested_scale_discrimination(self):
        constant = np.ones((8, 8))
        checker = np.fromfunction(lambda y, x: (-1.) ** (x + y), (8, 8))
        self.assertEqual(int(np.argmax(stage_e2.nested_projection_bands(constant))), 0)
        self.assertEqual(int(np.argmax(stage_e2.nested_projection_bands(checker))), 3)

    def test_patch_scale(self):
        field = np.zeros((8, 8)); field[:2, :2] = 1.; field[:2, 2:4] = -1.
        self.assertIn(int(np.argmax(stage_e2.nested_projection_bands(field))), (2, 3))

    def test_mpi_normalization(self):
        self.assertAlmostEqual(stage_e2.mpi_concentration([20, 0, 0, 5], [10, 1, 1, 10]), 4.)

    def test_scale_energy_reconciles(self):
        field = np.arange(64., dtype=float).reshape(8, 8)
        self.assertAlmostEqual(stage_e2.nested_projection_bands(field).sum(),
                               np.sum(field * field))


class IntegrityAndOperatorTests(unittest.TestCase):
    def test_snapshot_checksum_is_deterministic_and_sensitive(self):
        self.assertEqual(stage_e2.fnv1a64(b"actual residual"),
                         stage_e2.fnv1a64(b"actual residual"))
        self.assertNotEqual(stage_e2.fnv1a64(b"actual residual"),
                            stage_e2.fnv1a64(b"recursive residual"))

    def test_exact_spd_preconditioner(self):
        operator = np.diag([2., 4., 8.])
        metrics = stage_e2.preconditioner_metrics(
            operator, np.linalg.inv(operator), [1., 2., 3.])
        self.assertAlmostEqual(metrics["E_P"], 0.)
        self.assertAlmostEqual(metrics["cosine"], 1.)
        self.assertGreater(metrics["qTPq"], 0.)
        self.assertGreater(metrics["qTSPq"], 0.)


class ContractTests(unittest.TestCase):
    def test_high_frequency_authorizes_only_local_forensics(self):
        self.assertIn('scale_verdict in ("PATCH_LOCAL", "HIGH_FREQUENCY")',
                      TOOL.read_text())
        self.assertIn('next_path = "LOCAL_SCHWARZ_PATH"', TOOL.read_text())
        self.assertIn('"forensic_path_authorized": next_path == "LOCAL_SCHWARZ_PATH"',
                      TOOL.read_text())
        self.assertIn('"production_schwarz_modification_authorized": False',
                      TOOL.read_text())

    def test_next_action_enums_are_closed(self):
        source = TOOL.read_text()
        expected = {"LOCAL_SCHWARZ_PATH", "MPI_OVERLAP_PATH",
                    "GLOBAL_COARSE_PATH", "TARGETED_LOW_RANK_PATH",
                    "MULTILEVEL_PRESSURE_PATH", "PRECONDITIONER_REDESIGN_PATH",
                    "UNRESOLVED"}
        for value in expected:
            self.assertIn(f'"{value}"', source)
        for obsolete in ("LOCAL_SCHWARZ_REDESIGN_REVIEW",
                         "ALTERNATE_LAYOUT_MPI_OVERLAP_CONFIRMATION",
                         "BROADER_MULTILEVEL_TREATMENT_REVIEW"):
            self.assertNotIn(obsolete, source)

    def test_schedule_and_case_balancing_are_frozen(self):
        self.assertEqual(len(stage_e2.SCHEDULE), 19)
        source = (ROOT / "lib" / "Strict_ala_stage_e2.inc").read_text()
        self.assertIn("late_direction_case_balanced", source)
        self.assertIn("online_collectives=0", source)
        self.assertIn("map==1 ? 3 : 1", source)

    def test_truncated_csv_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "empty.csv"
            path.write_text("case,iteration\n")
            with self.assertRaises(ValueError):
                stage_e2.read_csv(path, ("case", "iteration"))

    def test_complete_synthetic_contract_reaches_classification(self):
        labels = [(case, iteration) for case in ("E0", "E1")
                  for iteration in stage_e2.SCHEDULE]
        data = {key: [] for key in stage_e2.SCHEMAS}
        data["snapshots_manifest"] = [
            {"case": case, "iteration": str(iteration), "residual_norm": "1",
             "global_checksum": "0123456789abcdef",
             "storage_identifier": "snap/E.rank%06d.bin",
             "ranks": "384", "local_npno": "8"} for case, iteration in labels]
        for ci, ii in labels:
            for cj, ij in labels:
                data["snapshot_gram"].append({"case_i": ci, "iteration_i": str(ii),
                    "case_j": cj, "iteration_j": str(ij),
                    "inner_product": "1" if (ci, ii) == (cj, ij) else "0"})
        groups = (("E0", "full_physical", 19), ("E0", "late_physical", 13),
                  ("E0", "late_direction", 13), ("E1", "full_physical", 19),
                  ("E1", "late_physical", 13), ("E1", "late_direction", 13),
                  ("JOINT", "late_physical", 26),
                  ("JOINT", "late_direction_case_balanced", 26))
        for case, basis, count in groups:
            for mode in range(1, count + 1):
                data["pod"].append({"case": case, "basis_type": basis,
                    "mode": str(mode), "singular_value": "1",
                    "energy_fraction": str(1 / count),
                    "cumulative_fraction": str(mode / count), "N50": str((count + 1)//2),
                    "N80": str(math.ceil(.8 * count)), "N90": str(math.ceil(.9 * count)),
                    "N95": str(math.ceil(.95 * count)), "weighting": "test"})
        for label_index, (case, iteration) in enumerate(labels):
            for mode in (1, 2):
                data["joint_pod_coefficients"].append({"case": case,
                    "iteration": str(iteration), "pod_mode": str(mode),
                    "coefficient": "1" if label_index == mode - 1 else "0",
                    "contraction": "1", "sign_change": "0"})
        data["subspace_angles"] = [{"energy_level": "80", "subspace_dimension": "1",
            "angle_index": "1", "angle_degrees": "5", "cosine": ".996",
            "threshold_status": stage_e2.PROVISIONAL}]
        objects = [(case, "pod_mode", str(mode))
                   for case in ("E0", "E1", "JOINT") for mode in (1, 2)]
        for case, kind, identifier in objects:
            for depth in ("0_200", "200_410", "410_660", "660_1000", "1000_cmb"):
                for band in ("GLOBAL_VERY_LONG", "LONG_INTERMEDIATE",
                             "PATCH_LOCAL", "HIGH_FREQUENCY"):
                    data["depth_scale_spectrum"].append({"case": case,
                        "object_type": kind, "object_id": identifier,
                        "depth_band": depth, "horizontal_scale_band": band,
                        "energy": "1", "normalized_fraction": ".05",
                        "method": "geometry_native_nested_v1"})
            for distance in ("0", "1", "2", "3+"):
                data["mpi_distance"].append({"case": case, "object_type": kind,
                    "object_id": identifier, "mpi_distance_band": distance,
                    "element_count": "10", "energy": "10", "energy_per_element": "1"})
        for mode in (1, 2):
            for mapping in ("BPI", "CONFIGURED", "PURE_SCHWARZ"):
                data["true_mode_preconditioner"].append({"mode": str(mode),
                    "map": mapping, "mode_energy_fraction": ".5",
                    "cumulative_energy": str(.5 * mode), "E_P": ".1",
                    "cosine": "1", "alpha_opt": "1", "Pq_norm": "1",
                    "SPq_norm": "1", "qTPq": "1", "qTSPq": "1",
                    "positive_qTPq": "1", "positive_qTSPq": "1",
                    "tight_solve_achieved": "1e-10",
                    "map_semantics": ("configured_minus_BPI_exact_component"
                        if mapping == "PURE_SCHWARZ" else "complete_map")})
        modes = stage_e2.validate(data)
        validation = stage_e2.numerical_validation(data, modes)
        self.assertEqual(validation["status"], "PASS")
        prior = {"case_E0": {"restart_verdict": "RESIDUAL_TAIL_PREEXISTS_RESTART",
                              "density_gauge_fraction": 1e-6},
                 "case_E1": {"restart_verdict": "RESIDUAL_TAIL_PREEXISTS_RESTART",
                              "density_gauge_fraction": 2e-6}}
        result = stage_e2.classify(data, modes, prior)
        self.assertFalse(result["automatic_solver_change_authorized"])
        self.assertTrue(result["density_gauge_remains_non_dominant"])
        self.assertIn("case_E0", result)
        self.assertIn("case_E1", result)
        self.assertIn("thresholds", result)


if __name__ == "__main__":
    unittest.main()
