import csv
import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "causal", ROOT / "tools/analyze_strict_ala_stage_F2g_causality.py")
CAUSAL = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CAUSAL)


def write_csv(path, fields, rows):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)


class CausalityTests(unittest.TestCase):
    ITER_FIELDS = (
        "case", "iteration", "krylov_recursive", "krylov_explicit",
        "continuity_relative", "momentum_relative", "cumulative_inner_solves",
        "cumulative_inner_cycles", "cumulative_K_gamma_applications",
        "cumulative_schur_actions", "elapsed_seconds", "final_iterate")
    INNER_FIELDS = (
        "case", "call_id", "outer_iteration", "requested_relative_tolerance",
        "target_absolute", "achieved_absolute", "achieved_relative", "cycles",
        "max_cycles", "seconds", "status")

    def test_cfg_diff_allows_only_inner_cap_and_output_paths(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            specs = []
            for case, cap in CAUSAL.CAPS.items():
                path = root / f"{case}.cfg"
                path.write_text("[CitcomS]\nsteps=1\n"
                                "[CitcomS.solver.mesher]\n"
                                "nprocx=4\nnprocy=4\nnprocz=2\n"
                                "nodex=129\nnodey=129\nnodez=65\nlevels=5\n"
                                "[CitcomS.solver.vsolver]\n"
                                f"ala_inner_accuracy_max={cap}\n"
                                "ala_inner_accuracy_factor=1e-2\n"
                                "piterations=60\ntole_compressibility=1e-02\n"
                                "ala_augmented_lagrangian_gamma=10.0\n"
                                "ala_element_vanka_smoother=on\n"
                                "ala_element_vanka_damping=0.2\n"
                                "ala_pcg_restart_interval=50\n"
                                "ala_outer_solver=fgmres\n"
                                "ala_shallow_patch_preconditioner=on\n"
                                "ala_shallow_patch_horizontal_elements=6\n"
                                "ala_shallow_patch_horizontal_stride=3\n"
                                "ala_shallow_patch_mpi_overlap=3\n"
                                "ala_shallow_patch_velocity_solver=element_vanka\n"
                                "ala_stage_abc_production_logging=on\n"
                                "ala_stage_e_diagnostic=off\n"
                                "[CitcomS.solver]\n"
                                f"datadir={root}/{case}/DATA/%RANK\n"
                                f"datadir_old={root}/{case}/Restart\n")
                specs.append(f"{case}={path}")
            output = root / "diff.json"
            args = SimpleNamespace(cfg=specs, output=output)
            self.assertEqual(CAUSAL.cfg_diff(args), 0)
            value = json.loads(output.read_text())
            self.assertTrue(value["valid"])
            self.assertEqual(value["unauthorized_differences"], [])

    def test_historical_zero_envelope_is_derived_not_invented(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fields = ("case", "iteration", "continuity_relative",
                      "momentum_relative")
            rows = [{"case": "E1", "iteration": i,
                     "continuity_relative": 1 / (i + 1),
                     "momentum_relative": 1e-4} for i in range(1, 61)]
            write_csv(root / "a.csv", fields, rows)
            write_csv(root / "b.csv", fields, rows)
            args = SimpleNamespace(reference_a=root / "a.csv",
                                   reference_b=root / "b.csv",
                                   output=root / "envelope.json")
            self.assertEqual(CAUSAL.derive_envelope(args), 0)
            value = json.loads(args.output.read_text())
            self.assertTrue(value["valid"])
            self.assertEqual(value["frozen_allowed_relative_difference"], 0.0)

    def test_lsf_has_three_cold_processes_and_frozen_thresholds(self):
        runs = ROOT.parents[1] / "runs"
        text = (runs / "cmbhf_ALA_strict_stage_F2g_causality.lsf").read_text()
        thresholds = json.loads((runs /
            "cmbhf_ALA_strict_stage_F2g_causality_thresholds.json").read_text())
        self.assertIn("CASE_NAMES=(INNER_1E2 INNER_3E3 INNER_1E3)", text)
        self.assertIn("for index in 0 1 2", text)
        self.assertIn("STRICT_STAGE_F2G_CAUSALITY_${LSB_JOBID}", text)
        self.assertIn("F2G_POSTPROCESS_V2", text)
        self.assertIn("python/pythia-0.8.1.15-py2.6.egg", text)
        self.assertIn('export PYTHONPATH="${PYRE_START}" &&', text)
        launch = text.rindex('"${CODE_DIR}/bin/pycitcoms"')
        launch_prefix = text[max(0, launch - 250):launch]
        self.assertIn('export PYTHONPATH="${PYRE_START}"', launch_prefix)
        self.assertNotIn("F2H_", text)
        self.assertEqual(thresholds["primary_late_window"], [20, 40])
        self.assertEqual(thresholds["reproduction_gate"]
                         ["allowed_relative_difference"], 0.0)
        self.assertFalse(thresholds["production_default_change_authorized"])

    def make_analysis(self, directory, factors):
        root = Path(directory)
        iterations = []; inner = []; completions = []
        reference_rows = []
        for case in CAUSAL.CASES:
            case_root = root / case; case_root.mkdir()
            iter_path = case_root / "iterations.csv"
            inner_path = case_root / "inner.csv"
            rows = []
            for k in range(1, 61):
                base = 1.0 / (k + 1)
                cont = base * factors[case]
                rows.append({
                    "case": case, "iteration": k,
                    "krylov_recursive": cont * 10,
                    "krylov_explicit": cont * 10,
                    "continuity_relative": cont, "momentum_relative": 5e-4,
                    "cumulative_inner_solves": k,
                    "cumulative_inner_cycles": k * {
                        "INNER_1E2": 10, "INNER_3E3": 15,
                        "INNER_1E3": 20}[case],
                    "cumulative_K_gamma_applications": k,
                    "cumulative_schur_actions": k, "elapsed_seconds": k,
                    "final_iterate": int(k == 60)})
                if case == "INNER_1E2":
                    reference_rows.append({
                        "case": "E1", "iteration": k,
                        "continuity_relative": cont,
                        "momentum_relative": 5e-4})
            write_csv(iter_path, self.ITER_FIELDS, rows)
            cap = CAUSAL.CAPS[case]
            inner_rows = [{
                "case": case, "call_id": k, "outer_iteration": k,
                "requested_relative_tolerance": cap,
                "target_absolute": cap, "achieved_absolute": cap * .9,
                "achieved_relative": cap * .9, "cycles": 10,
                "max_cycles": 2000, "seconds": 1, "status": "CONVERGED"}
                for k in range(1, 61)]
            write_csv(inner_path, self.INNER_FIELDS, inner_rows)
            completion_path = case_root / "completion.json"
            artifacts = {str(p.resolve()): hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in (iter_path, inner_path)}
            completion_path.write_text(json.dumps({
                "complete": True, "valid": True, "exit_status": 0,
                "artifacts": artifacts}))
            iterations.append(f"{case}={iter_path}")
            inner.append(f"{case}={inner_path}")
            completions.append(f"{case}={completion_path}")
        reference = root / "reference.csv"
        write_csv(reference, ("case", "iteration", "continuity_relative",
                              "momentum_relative"), reference_rows)
        thresholds = ROOT.parents[1] / "runs/cmbhf_ALA_strict_stage_F2g_causality_thresholds.json"
        envelope = root / "envelope.json"
        envelope.write_text(json.dumps({"valid": True,
                                        "frozen_allowed_relative_difference": 0.0}))
        cfg_diff = root / "cfg.json"
        cfg_diff.write_text(json.dumps({"valid": True}))
        return SimpleNamespace(
            iterations=iterations, inner=inner, completion=completions,
            thresholds=thresholds, reproduction_envelope=envelope,
            cfg_diff=cfg_diff, reference=reference,
            aligned_csv=root / "aligned.csv",
            matched_work_csv=root / "matched.csv", output=root / "decision.json")

    def test_machine_classification_controls_outer_continuity(self):
        with tempfile.TemporaryDirectory() as directory:
            args = self.make_analysis(directory, {
                "INNER_1E2": 1.0, "INNER_3E3": .8, "INNER_1E3": .64})
            self.assertEqual(CAUSAL.analyze(args), 0)
            value = json.loads(args.output.read_text())
            self.assertEqual(value["decision"],
                             "INNER_ACCURACY_CONTROLS_OUTER_CONTINUITY")
            self.assertTrue(value["observed_systematic_improvement"])
            self.assertFalse(value["independent_repeatability_claimed"])

    def test_machine_classification_does_not_control(self):
        with tempfile.TemporaryDirectory() as directory:
            args = self.make_analysis(directory, {
                "INNER_1E2": 1.0, "INNER_3E3": .96, "INNER_1E3": .95})
            self.assertEqual(CAUSAL.analyze(args), 0)
            value = json.loads(args.output.read_text())
            self.assertEqual(value["decision"],
                             "INNER_ACCURACY_DOES_NOT_CONTROL_OUTER_CONTINUITY")
            self.assertFalse(value["observed_systematic_improvement"])


if __name__ == "__main__":
    unittest.main()
