import importlib.util
import csv
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOL = ROOT / "tools" / "analyze_strict_ala_stage_F1a.py"
SPEC = importlib.util.spec_from_file_location("stage_f1a", TOOL)
f1a = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(f1a)


def rows_from_energies(mode_energies, weights=None, standalone=None, interaction=None):
    weights = weights or [1 / len(mode_energies)] * len(mode_energies)
    standalone = standalone or [values[-1] - 1 for values in mode_energies]
    interaction = interaction or [0.0] * len(mode_energies)
    rows = []
    for mode, (energies, weight, stand, inter) in enumerate(
            zip(mode_energies, weights, standalone, interaction), 1):
        for stage, energy in zip(f1a.STAGES, energies):
            row = {column: "0" for column in f1a.TELESCOPING_COLUMNS}
            row.update({"POD_mode": str(mode), "stage": stage,
                        "POD_energy_weight": str(weight), "residual_energy": str(energy),
                        "standalone_residual_energy": str(1 + stand),
                        "yB_yS": str(inter / 2), "valid": "1"})
            rows.append(row)
    return rows


class TelescopingTests(unittest.TestCase):
    def test_increments_sum_to_full_schwarz_harm(self):
        result = f1a.classify_telescoping(rows_from_energies(
            [[1.0, 1.1, 1.3, 1.25, 1.4, 1.4]]))
        increments = result["weighted_telescoping_damage"]
        self.assertAlmostEqual(sum(increments.values()), 0.4)

    def test_constructive_and_destructive_stages(self):
        result = f1a.classify_telescoping(rows_from_energies(
            [[1.0, .8, .7, .9, 1.4, 1.4]]))
        self.assertLess(result["weighted_telescoping_damage"]["RAW_LOCAL"], 0)
        self.assertGreater(result["weighted_telescoping_damage"]["LEFT_SCALING"], 0)

    def test_unique_primary_requires_sixty_percent(self):
        result = f1a.classify_telescoping(rows_from_energies(
            [[1.0, 1.7, 1.8, 1.9, 2.0, 2.0]]))
        self.assertEqual(result["dominant_damage_stage"], "RAW_LOCAL")
        self.assertGreaterEqual(result["primary_explained_harm_fraction"], .6)

    def test_mixed_25_75_rule(self):
        result = f1a.classify_telescoping(rows_from_energies(
            [[1.0, 1.4, 1.7, 1.7, 2.0, 2.0]]))
        self.assertEqual(result["primary_defect_class"], "MIXED")

    def test_unresolved_plurality(self):
        result = f1a.classify_telescoping(rows_from_energies(
            [[1.0, 1.5, 1.6, 1.7, 1.8, 1.8]]))
        self.assertEqual(result["primary_defect_class"], "UNRESOLVED")

    def test_mode_consistency_needs_eighty_percent(self):
        result = f1a.classify_telescoping(rows_from_energies(
            [[1.0, 1.7, 1.7, 1.7, 1.7, 1.7],
             [1.0, .3, .3, .3, .3, .3]], weights=[.79, .21]))
        self.assertEqual(result["dominant_damage_stage"], "UNRESOLVED")


class AdditiveTests(unittest.TestCase):
    def test_addition_is_independent_not_telescoping_stage(self):
        result = f1a.classify_telescoping(rows_from_energies(
            [[1.0, .9, .8, .7, .6, 1.3]], standalone=[-.4], interaction=[.7]))
        self.assertEqual(result["primary_defect_class"], "ADDITIVE_INTERFERENCE")
        self.assertEqual(result["dominant_damage_stage"], "ADDITION_TO_BPI")
        self.assertAlmostEqual(result["additive_interaction"], .7)

    def test_configured_harm_alone_is_not_additive_causality(self):
        result = f1a.classify_telescoping(rows_from_energies(
            [[1.0, 1.1, 1.1, 1.1, 1.1, 1.3]], standalone=[.1], interaction=[.2]))
        self.assertNotEqual(result["primary_defect_class"], "ADDITIVE_INTERFERENCE")


class StaticContractTests(unittest.TestCase):
    def test_e2_observer_has_exactly_one_owner_free(self):
        source = (ROOT / "lib" / "Strict_ala_stage_e2.inc").read_text()
        body = source.split("static void ala_e2_free_observer", 1)[1].split(
            "static void ala_e2_covariance", 1)[0]
        self.assertEqual(body.count("free(o);"), 1)

    def test_production_path_is_default_and_diagnostic_state_is_restored(self):
        stokes = (ROOT / "lib" / "Stokes_flow_Incomp.c").read_text()
        f1a_source = (ROOT / "lib" / "Strict_ala_stage_f1a.inc").read_text()
        self.assertIn("strict_ala_f1a_patch_stage = -1", stokes)
        self.assertIn("strict_ala_f1a_patch_stage=-1;", f1a_source)
        self.assertIn("E2_REPLAY_LAYOUT_MISMATCH", f1a_source)

    def test_closed_budgets_and_safety(self):
        self.assertEqual(f1a.MAX_PATCHES, 128)
        self.assertEqual(f1a.MAX_TIGHT_SOLVES, 40)
        self.assertFalse(f1a.thresholds().get("automatic_solver_change_authorized", False))

    def test_operator_and_inverse_are_distinct_fields(self):
        self.assertIn("local_solve_relative_residual", f1a.LOCAL_COLUMNS)
        self.assertIn("local_action_error", f1a.LOCAL_COLUMNS)

    def test_complete_synthetic_analysis_is_fail_closed_and_runnable(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); output = root / "analysis"
            def write(name, columns, rows):
                with (root / name).open("w", newline="") as stream:
                    writer = csv.DictWriter(stream, fieldnames=columns)
                    writer.writeheader(); writer.writerows(rows)
            telescoping = rows_from_energies(
                [[1, .9, .9, .9, .9, .9], [1, .9, .9, .9, .9, .9]],
                weights=[.5, .5])
            write("strict_ala_stage_F1a_telescoping.csv", f1a.TELESCOPING_COLUMNS,
                  telescoping)
            write("strict_ala_stage_F1a_tight_solves.csv", f1a.TIGHT_COLUMNS,
                  [{"call_id": "1", "role": "test", "POD_mode": "1",
                    "rhs_norm": "1", "requested_relative_tolerance": "1e-10",
                    "target_residual": "1e-10", "achieved_relative_residual": "1e-11",
                    "cycles": "1", "max_cycles": "2000", "converged": "1"}])
            write("strict_ala_stage_F1a_patch_selection.csv", f1a.PATCH_COLUMNS,
                  [{"patch_ID": "1", "selection_rule": "TOP", "POD_mode": "1",
                    "depth_stratum": "0_200", "contribution_quantile": "top",
                    "deterministic_order": "1", "outlier_guard": "0"}])
            write("strict_ala_stage_F1a_local_operator.csv", f1a.LOCAL_COLUMNS,
                  [{column: ("1" if column in ("patch_or_bin_ID", "POD_mode", "valid") else "0")
                    for column in f1a.LOCAL_COLUMNS}])
            write("strict_ala_stage_F1a_nonlocality.csv", f1a.NONLOCAL_COLUMNS,
                  [{column: ("1" if column in ("source_vector_or_bin", "valid") else "0")
                    for column in f1a.NONLOCAL_COLUMNS}])
            projected = [{"matrix": matrix, "row_mode": str(i), "column_mode": str(j),
                          "value": "0"} for matrix in ("M_B", "M_S", "M_cfg", "H_B", "H_S", "H_cfg")
                         for i in (1, 2) for j in (1, 2)]
            write("strict_ala_stage_F1a_projected_matrices.csv", f1a.PROJECTED_COLUMNS, projected)
            write("strict_ala_stage_F1a_patch_structure.csv", f1a.STRUCTURE_COLUMNS,
                  [{"supported_DOFs": "1", "multiplicity_min": "1",
                    "multiplicity_max": "1", "invalid_multiplicity": "0",
                    "top_W": "1", "mid_W": "1", "transition_W": "1",
                    "partition_formula": "test", "valid": "1"}])
            write("strict_ala_stage_F1a_support.csv", f1a.SUPPORT_COLUMNS,
                  [{"POD_mode": str(mode), "q_energy_in_support_fraction": "1",
                    "eB_energy_in_support_fraction": "1", "valid": "1"}
                   for mode in (1, 2)])
            e2_columns = ("mode", "map", "mode_energy_fraction", "cumulative_energy",
                "E_P", "cosine", "alpha_opt", "Pq_norm", "SPq_norm", "qTPq",
                "qTSPq", "positive_qTPq", "positive_qTSPq", "tight_solve_achieved",
                "map_semantics")
            archived = []
            current = {(int(row["POD_mode"]), row["stage"]): row for row in telescoping}
            for mode in (1, 2):
                for mapping in ("BPI", "CONFIGURED"):
                    row = {column: "0" for column in e2_columns}
                    row.update({"mode": str(mode), "map": mapping,
                        "E_P": current[mode, mapping]["E_P"],
                        "cosine": current[mode, mapping]["cosine"],
                        "qTSPq": current[mode, mapping]["qTSPq"]})
                    archived.append(row)
            write("e2.csv", e2_columns, archived)
            e2_reanalysis = root / "e2.json"
            e2_reanalysis.write_text(json.dumps({"next_authorized_path": "LOCAL_SCHWARZ_PATH",
                "forensic_path_authorized": True,
                "production_schwarz_modification_authorized": False}))
            result = f1a.analyze(root, output, root / "e2.csv", e2_reanalysis)
            self.assertTrue(result["experiment_valid"])
            self.assertFalse(result["automatic_solver_change_authorized"])


if __name__ == "__main__":
    unittest.main()
