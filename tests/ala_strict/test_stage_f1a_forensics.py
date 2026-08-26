import importlib.util
import csv
import json
import shutil
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

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
    def test_normalized_linearity_warns_but_hard_failure_is_rejected(self):
        base = {"role": "CONFIGURED", "difference_norm": "1e-7",
                "combined_norm": "1", "sum_norm": "1",
                "relative_defect": "1e-7", "solve_floor": "1e-10",
                "warn_limit": "1e-8", "hard_limit": "1e-5",
                "status": "WARN", "valid": "1"}
        self.assertEqual(f1a.validate_linearity([base])["status"], "WARN")
        failed = dict(base, relative_defect="1e-4", status="FAIL", valid="0")
        with self.assertRaises(ValueError):
            f1a.validate_linearity([failed])

    def test_configured_replay_uses_validated_linearity_envelope_only(self):
        with tempfile.TemporaryDirectory() as directory:
            archived = Path(directory) / "e2.csv"
            columns = ("mode", "map", "mode_energy_fraction", "cumulative_energy",
                "E_P", "cosine", "alpha_opt", "Pq_norm", "SPq_norm", "qTPq",
                "qTSPq", "positive_qTPq", "positive_qTSPq", "tight_solve_achieved",
                "map_semantics")
            with archived.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=columns); writer.writeheader()
                for mode in (1, 2):
                    for mapping in ("BPI", "CONFIGURED"):
                        row = {column: "0" for column in columns}
                        row.update({"mode": str(mode), "map": mapping,
                                    "E_P": "1", "cosine": "1", "qTSPq": "1"})
                        writer.writerow(row)
            current = rows_from_energies([[1] * 6, [1] * 6], weights=[.5, .5])
            for row in current:
                row.update({"E_P": "1", "cosine": "1", "qTSPq": "1"})
                if row["stage"] == "CONFIGURED": row["cosine"] = "1.00000002"
            linearity = {"relative_defect": "3e-7", "hard_limit": "1e-5"}
            replay = f1a.validate_operator_replay(current, archived, linearity)
            configured = [row for row in replay if row["map"] == "CONFIGURED"]
            self.assertTrue(all(row["effective_relative_tolerance"] == 3e-7
                                for row in configured))
            with self.assertRaises(ValueError):
                f1a.validate_operator_replay(current, archived,
                    {"relative_defect": "1e-9", "hard_limit": "1e-5"})
            for row in current:
                if row["stage"] == "BPI": row["cosine"] = "1.00000002"
            with self.assertRaises(ValueError):
                f1a.validate_operator_replay(current, archived, linearity)

    def test_rank0_log_supplies_missing_tight_solve_cycles(self):
        with tempfile.TemporaryDirectory() as directory:
            log = Path(directory) / "model.log"
            log.write_text("ALA COUPLED INNER VELOCITY summary status=converged "
                "cycles=481 max_cycles=2000 residual=9.91e-10 initial=1.0e1 "
                "target=1.0e-9 relative=9.91e-11 seconds=3.0e2\n")
            row = {"call_id": "1", "role": "BPI", "POD_mode": "1",
                   "rhs_norm": "10", "requested_relative_tolerance": "1e-10",
                   "target_residual": "1e-9", "achieved_relative_residual": "9.9e-11",
                   "cycles": "0", "max_cycles": "2000", "converged": "1"}
            evidence = f1a.validate_tight_solves([row], log)
            self.assertEqual(evidence["source"], "rank0_model_log")
            self.assertEqual(evidence["total_cycles"], 481)
            log.write_text("")
            with self.assertRaises(ValueError):
                f1a.validate_tight_solves([row], log)

    def test_f1a_cfg_contract_accepts_only_configured_schwarz_difference(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            base = ROOT.parents[1] / "runs" / "cmbhf_ALA_strict.cfg"
            c0, c1, diff = root / "c0.cfg", root / "c1.cfg", root / "diff.txt"
            text = base.read_text()
            c0.write_text(text.replace("ala_shallow_patch_preconditioner   = on",
                                       "ala_shallow_patch_preconditioner   = off"))
            shutil.copyfile(base, c1)
            self.assertEqual(f1a.check_cfg_contract(c0, c1, diff), 0)
            self.assertEqual(diff.read_text(),
                             "ala_shallow_patch_preconditioner: off -> on\n")
            c1.write_text(c1.read_text().replace("ala_pcg_restart_interval       = 50",
                                                 "ala_pcg_restart_interval       = 49"))
            with self.assertRaises(ValueError):
                f1a.check_cfg_contract(c0, c1, diff)

    def test_nonlocal_sources_encode_mode_independently_of_row_order(self):
        self.assertEqual(f1a.nonlocal_mode("1000007_mode_12"), 12)
        with self.assertRaises(ValueError):
            f1a.nonlocal_mode("1000007")

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

    def test_hpc_launcher_selects_only_valid_closed_e2_evidence(self):
        launcher = (ROOT.parents[1] / "runs" /
                    "cmbhf_ALA_strict_stage_F1a.lsf").read_text()
        self.assertIn("select_e2_root", launcher)
        self.assertIn("verify_closure_inputs", launcher)
        self.assertIn("current_binary.sha256", launcher)
        self.assertIn("E2_snapshots.post.sha256", launcher)
        self.assertIn("STRICT_ALA_CASE=F1A", launcher)
        self.assertIn("rank_log_contract", launcher)
        self.assertIn("verify_e2_binary_lineage", launcher)
        self.assertIn("F1a_gated_instrumentation_only_rebuild", launcher)
        self.assertIn('--model-log "${MODEL_LOG}"', launcher)
        self.assertIn("lib/Strict_ala_stage_f1a.inc|tools/analyze_strict_ala_stage_F1a.py|tests/ala_strict/test_stage_f1a_forensics.py", launcher)
        self.assertNotIn("STRICT_STAGE_E2_12110470", launcher)

    def test_postprocess_launcher_reuses_only_frozen_analysis_only_evidence(self):
        launcher = (ROOT.parents[1] / "runs" /
                    "cmbhf_ALA_strict_stage_F1a_postprocess.lsf").read_text()
        self.assertIn("#BSUB -q ser", launcher)
        self.assertIn("#BSUB -n 1", launcher)
        self.assertIn("verify_evidence_root", launcher)
        self.assertIn("relation=F1a_analysis_only_replay", launcher)
        self.assertIn("numerical_evidence.pre.sha256", launcher)
        self.assertIn("E2_snapshots.post.sha256", launcher)
        self.assertIn("postprocess-audit", launcher)
        self.assertNotIn("pycitcoms --pyre-start", launcher)

    def test_postprocess_audit_requires_unchanged_frozen_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            def dump(name, value):
                path = root / name; path.write_text(json.dumps(value)); return path
            decision = dump("decision.json", {"experiment_valid": True,
                "automatic_solver_change_authorized": False,
                "production_schwarz_modification_authorized": False,
                "next_authorized_task": "F1B_LOCAL_CANDIDATE"})
            provenance = {"provenance_complete": True, "source": {
                "branch": "cmbhf_ALA_strict", "branch_verified": True,
                "scientific_dirty": []}}
            numerical = dump("numerical.json", provenance)
            analysis = dump("analysis.json", provenance)
            e2 = dump("e2.json", {"experiment_valid": True,
                "numerical_validation_pass": True,
                "observation_only_trajectory_pass": True,
                "next_authorized_path": "LOCAL_SCHWARZ_PATH"})
            log = root / "model.log"; log.write_text("STRICT_ALA_STAGE_F1A_COMPLETE\n")
            lineage = root / "lineage.txt"
            lineage.write_text("relation=F1a_analysis_only_replay\n"
                "tools/analyze_strict_ala_stage_F1a.py\n")
            pairs = {}
            for stem in ("evidence", "analysis", "source", "snapshot"):
                before, after = root / (stem + ".pre"), root / (stem + ".post")
                before.write_text(stem); after.write_text(stem); pairs[stem] = (before, after)
            args = SimpleNamespace(decision=decision, numerical_manifest=numerical,
                analysis_manifest=analysis, evidence_pre=pairs["evidence"][0],
                evidence_post=pairs["evidence"][1], analysis_pre=pairs["analysis"][0],
                analysis_post=pairs["analysis"][1], source_pre=pairs["source"][0],
                source_post=pairs["source"][1], snapshot_pre=pairs["snapshot"][0],
                snapshot_post=pairs["snapshot"][1], lineage=lineage, e2_audit=e2,
                model_log=log, evidence_root=root, output=root / "final.json")
            self.assertEqual(f1a.postprocess_audit(args), 0)
            pairs["evidence"][1].write_text("tampered")
            self.assertEqual(f1a.postprocess_audit(args), 1)

    def test_linearity_guard_is_relative_reported_and_fail_closed(self):
        source = (ROOT / "lib" / "Strict_ala_stage_f1a.inc").read_text()
        guard = source.split("One independently solved dominant-mode linearity guard.", 1)[1]
        self.assertIn("relative_defect=difference_norm/scale", guard)
        self.assertIn("strict_ala_stage_F1a_linearity.csv", source)
        self.assertIn("STRICT_ALA_STAGE_F1A_LINEARITY", guard)
        self.assertIn("relative_defect<=hard_limit", guard)
        self.assertNotIn("sqrt(global_pdot(E,y,y,lev))>1.0e-8", guard)

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
                  [{column: ("1000001_mode_1" if column == "source_vector_or_bin" else
                             "1" if column == "valid" else "0")
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
            write("strict_ala_stage_F1a_linearity.csv", f1a.LINEARITY_COLUMNS,
                  [{"role": "CONFIGURED", "difference_norm": "1e-11",
                    "combined_norm": "1", "sum_norm": "1",
                    "relative_defect": "1e-11", "solve_floor": "1e-10",
                    "warn_limit": "1e-8", "hard_limit": "1e-5",
                    "status": "PASS", "valid": "1"}])
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
