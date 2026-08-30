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
    "f1c", ROOT / "tools/analyze_strict_ala_stage_F1c.py")
F1C = importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(F1C)


class Args:
    def __init__(self, **values): self.__dict__.update(values)


def write_csv(path, fields, rows):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)


class F1cContractTests(unittest.TestCase):
    def test_runtime_is_opt_in_and_preserves_f1b_map(self):
        source = (ROOT / "lib/Strict_ala_stage_f1c.inc").read_text()
        driver = (ROOT / "lib/Stokes_flow_Incomp.c").read_text()
        self.assertIn('getenv("STRICT_ALA_STAGE_F1C_REQUIRED")', source)
        self.assertIn('"operator_consistent"', source)
        self.assertIn("ala_schur_preconditioner_action(", source)
        self.assertIn("fine_rhs,fine_z,work,aux,lev,3,cache", source)
        self.assertIn("value-=Sq[i][m][e]*c[i]", source)
        self.assertIn("value+=q[i][m][e]*(c[i]-d[i])", source)
        self.assertIn('#include "Strict_ala_stage_f1c.inc"', driver)
        self.assertIn("strict_ala_stage_f1c_run", driver)
        self.assertIn("production_default_change_authorized", source)

    def test_hpc_job_is_one_rank_closed_resumable_root(self):
        lsf = ROOT.parents[1] / "runs/cmbhf_ALA_strict_stage_F1c.lsf"
        text = lsf.read_text()
        self.assertIn("#BSUB -n 400", text)
        self.assertIn("--nodes=384", text)
        self.assertIn("STRICT_STAGE_F1C_${LSB_JOBID}", text)
        self.assertIn("STRICT_ALA_STAGE_F1C_RESUME_ROOT", text)
        self.assertIn("completion_valid", text)
        self.assertIn("balanced_coarse_completion.json", text)
        self.assertIn("analysis_completion.json", text)
        self.assertIn("final_audit_completion.json", text)
        self.assertEqual(text.count("CitcomS.SimpleApp:SimpleApp"), 1)
        self.assertIn("! -path '*/DATA/0/*'", text)
        self.assertIn("strict_ala_stage_F1c_final_audit.json", text)
        checked = subprocess.run(["bash", "-n", str(lsf)],
                                 text=True, capture_output=True)
        self.assertEqual(checked.returncode, 0, checked.stderr)

    def test_balanced_formula_annihilates_coarse_residual_projection(self):
        # Diagonal S and diagonal fine P give an exact arithmetic reference
        # for C + (I-CS)P(I-SC) on the first coarse basis vector.
        s = (2.0, 3.0)
        p = (0.4, 0.2)
        q = (1.0, 0.0)
        t = s[0]
        c = (q[0] / t, 0.0)
        fine_rhs = (q[0] - s[0] * c[0], q[1])
        fine_z = (p[0] * fine_rhs[0], p[1] * fine_rhs[1])
        d = s[0] * fine_z[0] / t
        z = (fine_z[0] + c[0] - d, fine_z[1])
        residual = (q[0] - s[0] * z[0], q[1] - s[1] * z[1])
        self.assertEqual(z, (0.5, 0.0))
        self.assertEqual(residual, (0.0, 0.0))

    def test_analyzer_accepts_only_closed_balanced_evidence(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            true_fields = ("POD_mode", "POD_energy_weight", "E_P", "cosine",
                           "qTPq", "qTSPq", "tight_solve_achieved", "operator")
            write_csv(td/"f1b.csv", true_fields, [
                {"POD_mode": str(mode), "POD_energy_weight": "0.5",
                 "E_P": "1.0", "cosine": "0.1", "qTPq": "1",
                 "qTSPq": "0.1", "tight_solve_achieved": "1e-10",
                 "operator": "operator_consistent"} for mode in (1, 2)])
            candidate_fields = ("POD_mode", "POD_energy_weight", "E_P",
                "cosine", "qTPq", "qTSPq",
                "coarse_residual_projection_max", "tight_solve_achieved", "valid")
            write_csv(td/"candidate.csv", candidate_fields, [
                {"POD_mode": str(mode), "POD_energy_weight": "0.5",
                 "E_P": "0.2", "cosine": "0.99", "qTPq": "1",
                 "qTSPq": "0.9", "coarse_residual_projection_max": "1e-10",
                 "tight_solve_achieved": "1e-10", "valid": "1"}
                for mode in (1, 2)])
            (td/"runtime.json").write_text(json.dumps({
                "selected_modes": 2, "selected_energy": 0.97,
                "coarse_matrix_symmetry_defect": 1e-12,
                "coarse_matrix_minimum_pivot_ratio": 0.2,
                "tight_solve_count": 4,
                "candidate": "balanced_actual_mode_coarse",
                "production_default_change_authorized": False}))
            projected = []
            for matrix in ("T", "M_balanced", "H_balanced"):
                for i in (1, 2):
                    for j in (1, 2):
                        projected.append({"matrix": matrix, "row_mode": i,
                                          "column_mode": j,
                                          "value": 1.0 if i == j else 0.0})
            write_csv(td/"projected.csv",
                      ("matrix", "row_mode", "column_mode", "value"), projected)
            tight_fields = ("call_id", "role", "POD_mode", "rhs_norm",
                            "requested_relative_tolerance", "target_residual",
                            "achieved_relative_residual", "cycles",
                            "max_cycles", "converged")
            write_csv(td/"tight.csv", tight_fields, [{
                "call_id": i, "role": "test", "POD_mode": i % 2 + 1,
                "rhs_norm": "1", "requested_relative_tolerance": "1e-10",
                "target_residual": "1e-10",
                "achieved_relative_residual": "9e-11", "cycles": "5",
                "max_cycles": "2000", "converged": "1"} for i in range(4)])
            args = Args(f1b_true_mode=td/"f1b.csv",
                        candidate_true_mode=td/"candidate.csv",
                        runtime=td/"runtime.json", projected=td/"projected.csv",
                        tight=td/"tight.csv", output=td/"decision.json")
            F1C.analyze(args)
            result = json.loads((td/"decision.json").read_text())
            self.assertEqual(result["decision"],
                             "F1C_PRESSURE_PHYSICS_QUALIFIED")
            self.assertEqual(result["next_authorized_task"],
                             "G0_GENERALIZABLE_COARSE_BASIS")
            self.assertFalse(result["production_default_change_authorized"])

            with (td/"candidate.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            rows[0]["coarse_residual_projection_max"] = "1e-4"
            write_csv(td/"candidate.csv", candidate_fields, rows)
            args.output = td/"rejected.json"
            F1C.analyze(args)
            rejected = json.loads((td/"rejected.json").read_text())
            self.assertEqual(rejected["decision"],
                             "F1C_BALANCED_COARSE_REJECTED")
            self.assertFalse(rejected["gates"]["coarse_residual_projection"])


if __name__ == "__main__": unittest.main()
