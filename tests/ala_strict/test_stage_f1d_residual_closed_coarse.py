#!/usr/bin/env python3
import csv
import importlib.util
import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "f1d", ROOT / "tools/analyze_strict_ala_stage_F1d.py")
F1D = importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(F1D)


class Args:
    def __init__(self, **values): self.__dict__.update(values)


def write_csv(path, fields, rows):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)


class F1dContractTests(unittest.TestCase):
    def test_runtime_is_opt_in_and_changes_only_the_coarse_basis(self):
        source = (ROOT / "lib/Strict_ala_stage_f1d.inc").read_text()
        driver = (ROOT / "lib/Stokes_flow_Incomp.c").read_text()
        self.assertIn('getenv("STRICT_ALA_STAGE_F1D_REQUIRED")', source)
        self.assertIn('"operator_consistent"', source)
        self.assertIn("RESIDUAL_CLOSURE_SR", source)
        self.assertIn("RESIDUAL_CLOSED_CANDIDATE", source)
        self.assertIn("for(pass=0;pass<2;pass++)", source)
        self.assertIn("ala_schur_preconditioner_action(", source)
        self.assertIn("fine_rhs,fine_z,work,aux,lev,3,cache", source)
        self.assertIn("F1b_operator_consistent_frozen", source)
        self.assertIn("production_default_change_authorized", source)
        self.assertIn('#include "Strict_ala_stage_f1d.inc"', driver)
        self.assertIn("strict_ala_stage_f1d_run", driver)

    def test_hpc_job_is_single_launch_rank_closed_and_resumable(self):
        lsf = ROOT.parents[1] / "runs/cmbhf_ALA_strict_stage_F1d.lsf"
        text = lsf.read_text()
        self.assertIn("#BSUB -n 400", text)
        self.assertIn("--nodes=384", text)
        self.assertIn("STRICT_STAGE_F1D_${LSB_JOBID}", text)
        self.assertIn("STRICT_ALA_STAGE_F1D_RESUME_ROOT", text)
        self.assertIn("export STRICT_ALA_CASE=F1D", text)
        self.assertIn("numerical_completion.json", text)
        self.assertIn("analysis_completion.json", text)
        self.assertIn("final_audit_completion.json", text)
        self.assertIn("STRICT_ALA_STAGE_F1D_RESIDUAL_CLOSURE_V1", text)
        self.assertEqual(text.count("CitcomS.SimpleApp:SimpleApp"), 1)
        self.assertIn("! -path '*/DATA/0/*'", text)
        checked = subprocess.run(["bash", "-n", str(lsf)],
                                 text=True, capture_output=True)
        self.assertEqual(checked.returncode, 0, checked.stderr)

    def test_analyzer_accepts_only_closed_improving_evidence(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            baseline_fields = ("POD_mode", "POD_energy_weight", "E_P", "operator")
            candidate_fields = (
                "POD_mode", "POD_energy_weight", "E_P", "cosine", "qTPq",
                "qTSPq", "original_projection_max", "enriched_projection_max",
                "tight_solve_achieved", "valid")
            basis_fields = (
                "basis_index", "basis_type", "source_mode", "preorthogonal_norm",
                "postorthogonal_norm", "accepted", "maximum_orthogonality",
                "S_energy", "valid")
            projected_fields = ("matrix", "row_mode", "column_mode", "value")
            tight_fields = (
                "call_id", "role", "POD_mode", "rhs_norm",
                "requested_relative_tolerance", "target_residual",
                "achieved_relative_residual", "cycles", "max_cycles", "converged")
            weights = (0.85, 0.10, 0.05)
            write_csv(td/"baseline.csv", baseline_fields, [
                {"POD_mode": i+1, "POD_energy_weight": weights[i], "E_P": 1.0,
                 "operator": "operator_consistent"} for i in range(3)])
            write_csv(td/"candidate.csv", candidate_fields, [
                {"POD_mode": i+1, "POD_energy_weight": weights[i], "E_P": 0.2,
                 "cosine": 1.0, "qTPq": 1.0, "qTSPq": 1.0,
                 "original_projection_max": 1e-12,
                 "enriched_projection_max": 1e-12,
                 "tight_solve_achieved": 1e-11, "valid": 1}
                for i in range(3)])
            write_csv(td/"basis.csv", basis_fields, [
                {"basis_index": i+4, "basis_type": "coarse_residual",
                 "source_mode": i+1, "preorthogonal_norm": 1.0,
                 "postorthogonal_norm": 1.0, "accepted": 1,
                 "maximum_orthogonality": 1e-12, "S_energy": 1.0, "valid": 1}
                for i in range(3)])
            projected = []
            for matrix, size in (("T_enriched", 6), ("M_balanced", 3),
                                 ("H_balanced", 3)):
                for i in range(size):
                    for j in range(size):
                        projected.append({"matrix": matrix, "row_mode": i+1,
                                          "column_mode": j+1,
                                          "value": 1.0 if i == j else 0.0})
            write_csv(td/"projected.csv", projected_fields, projected)
            write_csv(td/"tight.csv", tight_fields, [
                {"call_id": i+1, "role": "TEST", "POD_mode": i % 3 + 1,
                 "rhs_norm": 1.0, "requested_relative_tolerance": 1e-10,
                 "target_residual": 1e-10, "achieved_relative_residual": 1e-11,
                 "cycles": 10, "max_cycles": 2000, "converged": 1}
                for i in range(9)])
            (td/"runtime.json").write_text(json.dumps({
                "candidate": "residual_closed_balanced_actual_mode_coarse",
                "fine_map": "F1b_operator_consistent_frozen",
                "production_default_change_authorized": False,
                "selected_modes": 3, "selected_energy": 1.0,
                "accepted_residual_modes": 3, "enriched_basis_dimension": 6,
                "enriched_matrix_minimum_pivot_ratio": 0.5,
                "enriched_matrix_symmetry_defect": 1e-12,
                "tight_solve_count": 9}))
            (td/"f1c.json").write_text(json.dumps({
                "experiment_evidence_valid": True,
                "decision": "F1C_BALANCED_COARSE_REJECTED",
                "next_authorized_task": "F1D_COARSE_BASIS_REDESIGN",
                "production_default_change_authorized": False,
                "metrics": {"E_P_RMS_F1c": 3.0}}))
            args = Args(f1b_true_mode=td/"baseline.csv",
                        f1c_decision=td/"f1c.json",
                        candidate_true_mode=td/"candidate.csv",
                        basis=td/"basis.csv", projected=td/"projected.csv",
                        tight=td/"tight.csv", runtime=td/"runtime.json",
                        output=td/"decision.json")
            self.assertEqual(F1D.analyze(args), 0)
            result = json.loads((td/"decision.json").read_text())
            self.assertEqual(result["decision"],
                             "F1D_RESIDUAL_CLOSED_COARSE_QUALIFIED")
            self.assertEqual(result["next_authorized_task"],
                             "G0_GENERALIZABLE_COARSE_BASIS")
            with (td/"candidate.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            rows[0]["E_P"] = "2.0"
            write_csv(td/"candidate.csv", candidate_fields, rows)
            F1D.analyze(args)
            result = json.loads((td/"decision.json").read_text())
            self.assertEqual(result["decision"],
                             "F1D_RESIDUAL_CLOSED_COARSE_REJECTED")
            self.assertEqual(result["next_authorized_task"],
                             "F2_AVV_PRECONDITIONER_REVIEW")


if __name__ == "__main__":
    unittest.main()
