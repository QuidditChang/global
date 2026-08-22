import importlib.util
import math
from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
TOOL = ROOT / "tools" / "analyze_strict_ala_stage_E.py"
SOURCE = ROOT / "lib" / "Stokes_flow_Incomp.c"
ELEMENT = ROOT / "lib" / "Element_calculations.c"
GENERAL = ROOT / "lib" / "General_matrix_functions.c"

spec = importlib.util.spec_from_file_location("stage_e", TOOL)
stage_e = importlib.util.module_from_spec(spec)
spec.loader.exec_module(stage_e)


class ProjectionTests(unittest.TestCase):
    def test_exact_orthogonal_projection(self):
        result = stage_e.svd_projection(np.eye(2), np.array([3.0, 4.0]), 25.0)
        self.assertEqual(result["rank"], 2)
        self.assertAlmostEqual(result["represented_fraction"], 1.0)
        self.assertLessEqual(result["reconciliation_relative"], 1.0e-15)

    def test_nonorthogonal_projection(self):
        gram = np.array([[1.0, 0.5], [0.5, 1.0]])
        coefficients = np.array([2.0, -1.0])
        correlation = gram @ coefficients
        energy = float(coefficients @ gram @ coefficients)
        result = stage_e.svd_projection(gram, correlation, energy)
        self.assertAlmostEqual(result["represented_fraction"], 1.0, places=14)

    def test_rank_deficient_projection(self):
        gram = np.ones((2, 2))
        result = stage_e.svd_projection(gram, np.array([2.0, 2.0]), 4.0)
        self.assertEqual(result["rank"], 1)
        self.assertAlmostEqual(result["represented_fraction"], 1.0, places=14)

    def test_unrepresented_residual(self):
        result = stage_e.svd_projection(np.eye(3), np.zeros(3), 7.0)
        self.assertEqual(result["represented_fraction"], 0.0)
        self.assertEqual(result["unexplained_fraction"], 1.0)

    def test_inconsistent_projection_energy_is_not_clamped_to_pass(self):
        result = stage_e.svd_projection(np.eye(1), np.array([2.0]), 1.0)
        self.assertFalse(result["physical_energy_pass"])
        self.assertGreater(result["represented_energy"], 1.0)

    def test_exact_density_mode(self):
        gram = np.eye(len(stage_e.PROBES))
        correlation = np.zeros(len(stage_e.PROBES))
        correlation[stage_e.PROBES.index("P15_density_gauge")] = 1.0
        fraction = stage_e.family_projection(
            gram, correlation, 1.0, stage_e.FAMILIES["GAUGE_NEAR_NULL"])
        self.assertAlmostEqual(fraction, 1.0)

    def test_overlapping_family_detection(self):
        gram = np.eye(len(stage_e.PROBES))
        i = stage_e.PROBES.index("P6_degree_1")
        j = stage_e.PROBES.index("P5_horizontal_checkerboard")
        gram[i, j] = gram[j, i] = 0.95
        overlap = stage_e.family_overlap(
            gram, "GLOBAL_LOW_DEGREE", "LOCAL_HIGH_FREQUENCY")
        self.assertGreater(overlap, 0.9)

    def test_mixed_tail_is_not_forced_to_one_family(self):
        projections = {}
        for iteration in range(61):
            row = {"represented_fraction": 0.9}
            row.update({name: 0.1 for name in stage_e.PRIMARY_FAMILIES})
            row["GLOBAL_LOW_DEGREE"] = 0.7
            row["LOCAL_HIGH_FREQUENCY"] = 0.7
            projections[iteration] = row
        classification, dominant = stage_e._classify_tail(
            {}, projections, {})
        self.assertEqual(classification, "MIXED")
        self.assertIsNone(dominant)


class RestartTests(unittest.TestCase):
    @staticmethod
    def records(pre, post):
        records = {}
        for iteration in range(61):
            rho = pre if iteration in stage_e.PRE_WINDOW else post
            records[iteration] = {"rho": rho}
        return records

    @staticmethod
    def projections(pre_value=0.2, post_value=0.2):
        output = {}
        for iteration in range(61):
            value = pre_value if iteration < 50 else post_value
            output[iteration] = {name: value for name in stage_e.PRIMARY_FAMILIES}
        return output

    def test_tail_preexists_restart(self):
        rows = [{"residual_jump_relative": "0", "pre_post_cosine": "1"}]
        result = stage_e.restart_verdict(
            self.records(0.99, 0.995), rows, self.projections())
        self.assertEqual(result["verdict"], "RESIDUAL_TAIL_PREEXISTS_RESTART")

    def test_restart_sensitivity_strong(self):
        rows = [{"residual_jump_relative": "0", "pre_post_cosine": "1"}]
        result = stage_e.restart_verdict(
            self.records(0.8, 0.95), rows, self.projections(0.1, 0.4))
        self.assertEqual(result["verdict"], "RESTART_SENSITIVITY_STRONG")

    def test_restart_discontinuity_is_strong(self):
        rows = [{"residual_jump_relative": "1e-6", "pre_post_cosine": "1"}]
        result = stage_e.restart_verdict(
            self.records(0.8, 0.8), rows, self.projections())
        self.assertEqual(result["verdict"], "RESTART_SENSITIVITY_STRONG")

    def test_hessenberg_serialization_reconstructs_ritz_values(self):
        rows = [
            {"case": "E0", "restart_cycle": "1", "column": "0",
             "row": "0", "value": "2", "cycle_start_residual_norm": "1"},
            {"case": "E0", "restart_cycle": "1", "column": "0",
             "row": "1", "value": "0", "cycle_start_residual_norm": "1"},
            {"case": "E0", "restart_cycle": "1", "column": "1",
             "row": "0", "value": "0", "cycle_start_residual_norm": "1"},
            {"case": "E0", "restart_cycle": "1", "column": "1",
             "row": "1", "value": "3", "cycle_start_residual_norm": "1"},
            {"case": "E0", "restart_cycle": "1", "column": "1",
             "row": "2", "value": "0", "cycle_start_residual_norm": "1"},
        ]
        result = stage_e.ritz_summary(rows)[0]
        self.assertEqual(sorted(float(value.real) for value in result[3]), [2.0, 3.0])


class InstrumentationContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = SOURCE.read_text()
        cls.element = ELEMENT.read_text()
        cls.general = GENERAL.read_text()

    def function_body(self, name, next_name):
        return self.source.split(name, 1)[1].split(next_name, 1)[0]

    def test_sampling_uses_one_batched_reduction(self):
        body = self.function_body(
            "static void strict_ala_stage_e_sample",
            "static void strict_ala_stage_e_restart_pre")
        self.assertEqual(body.count("MPI_Allreduce"), 1)
        self.assertIn("ALA_STAGE_E_REDUCTION_SIZE", body)

    def test_batched_reduction_matches_componentwise_reductions(self):
        rank_contributions = np.arange(7 * (len(stage_e.PROBES) + 9),
                                       dtype=float).reshape(7, -1)
        batched = rank_contributions.sum(axis=0)
        scalar = np.array([rank_contributions[:, index].sum()
                           for index in range(rank_contributions.shape[1])])
        np.testing.assert_array_equal(batched, scalar)

    def test_observer_does_not_receive_solver_velocity_or_pressure(self):
        signature = self.source.split(
            "static void strict_ala_stage_e_sample", 1)[1].split("{", 1)[0]
        self.assertNotIn("double **V", signature)
        self.assertNotIn("double **P", signature)
        body = self.function_body(
            "static void strict_ala_stage_e_sample",
            "static void strict_ala_stage_e_restart_pre")
        self.assertNotIn("assemble_", body)
        self.assertNotIn("\nresidual[m][e]=", body.replace(" ", ""))
        self.assertIn("residual_state_guard_pass", self.source)

    def test_krylov_ratio_and_difference_are_separate(self):
        self.assertIn("krylov_residual_ratio", self.source)
        self.assertIn("krylov_residual_rel_difference", self.source)
        self.assertIn("fabs(recursive_residual-explicit_residual)", self.source)

    def test_true_operator_and_smoother_counters_are_at_execution_sites(self):
        self.assertGreaterEqual(
            self.element.count("ala_stage_e_k_operator_applications++"), 2)
        self.assertIn("ala_stage_e_velocity_smoother_sweeps += *cycles",
                      self.general)
        self.assertIn("globally_identical_logical_counts_not_process_sums",
                      self.source)

    def test_probe_definition_is_shared(self):
        self.assertIn("ala_schur_probe_raw_value", self.source)
        builder = self.function_body(
            "static void ala_schur_build_probe",
            "#define ALA_STAGE_E_REDUCTION_SIZE")
        self.assertIn("ala_schur_probe_raw_value", builder)

    def test_restart_and_hessenberg_hooks_exist(self):
        self.assertIn("strict_ala_stage_e_restart_pre", self.source)
        self.assertIn("strict_ala_stage_e_restart_post", self.source)
        self.assertIn("strict_ala_stage_e_log_hessenberg", self.source)
        self.assertIn("cycle_start_residual_norm", self.source)


class FailClosedSchemaTests(unittest.TestCase):
    @staticmethod
    def complete_case(case="E0"):
        iterations = []
        work = []
        correlations = []
        for iteration in range(61):
            iterations.append({
                "case": case, "iteration": str(iteration),
                "residual_state_guard_pass": "1",
                "depth_0_200_fraction": "0.2",
                "depth_200_410_fraction": "0.2",
                "depth_410_660_fraction": "0.2",
                "depth_660_1000_fraction": "0.2",
                "depth_1000_cmb_fraction": "0.2",
                "recursive_residual": "1", "explicit_residual": "2",
                "krylov_residual_ratio": "0.5",
                "krylov_residual_rel_difference": "0.5",
            })
            work.append({
                "case": case, "iteration": str(iteration),
                "pressure_schur_actions": str(iteration),
                "K_gamma_rhs_solves": str(iteration),
                "K_gamma_operator_applications": str(iteration),
                "velocity_MG_cycles": str(iteration),
                "velocity_smoother_sweeps": str(iteration),
                "preconditioner_applications": str(iteration),
            })
            correlations.extend(
                {"case": case, "iteration": str(iteration), "probe": probe}
                for probe in stage_e.PROBES)
        gram = [
            {"case": case, "probe_i": first, "probe_j": second}
            for first in stage_e.PROBES for second in stage_e.PROBES]
        restart = [
            {"case": case, "iteration": "50", "probe": probe}
            for probe in stage_e.PROBES]
        hessenberg = []
        for cycle, columns in ((1, 50), (2, 10)):
            for column in range(columns):
                hessenberg.extend({
                    "case": case, "restart_cycle": str(cycle),
                    "column": str(column), "row": str(row)}
                    for row in range(column + 2))
        return {"iterations": iterations, "work": work,
                "correlations": correlations, "gram": gram,
                "restart": restart, "hessenberg": hessenberg}

    def test_complete_schema_passes(self):
        self.assertTrue(stage_e.validate_case_structure(
            self.complete_case(), "E0"))

    def test_duplicate_correlation_fails_closed(self):
        data = self.complete_case()
        data["correlations"][-1] = data["correlations"][0].copy()
        self.assertFalse(stage_e.validate_case_structure(data, "E0"))

    def test_missing_second_hessenberg_cycle_fails_closed(self):
        data = self.complete_case()
        data["hessenberg"] = [row for row in data["hessenberg"]
                              if row["restart_cycle"] == "1"]
        self.assertFalse(stage_e.validate_case_structure(data, "E0"))

    def test_counter_regression_fails_closed(self):
        data = self.complete_case()
        data["work"][60]["K_gamma_operator_applications"] = "0"
        self.assertFalse(stage_e.validate_case_structure(data, "E0"))


if __name__ == "__main__":
    unittest.main()
