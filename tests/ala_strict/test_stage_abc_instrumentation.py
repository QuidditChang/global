import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "stage_abc", ROOT / "tools" / "strict_ala_stage_abc.py")
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


class StageABCInstrumentationTest(unittest.TestCase):
    def test_cfg_generator_changes_only_schwarz_switch(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td); base = td / "base.cfg"; c0 = td / "c0.cfg"; c1 = td / "c1.cfg"
            base.write_text("datadir = /same\nala_shallow_patch_preconditioner = on\n")
            MOD.cfg_variant(base, c0, False); MOD.cfg_variant(base, c1, True)
            self.assertEqual([d[0] for d in MOD.normalized_cfg_diff(c0, c1)],
                             ["ala_shallow_patch_preconditioner"])

    def test_cfg_whitelist_exposes_extra_scientific_change(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td); a = td / "a"; b = td / "b"
            a.write_text("ala_shallow_patch_preconditioner = off\ngamma = 10\n")
            b.write_text("ala_shallow_patch_preconditioner = on\ngamma = 20\n")
            self.assertEqual([x[0] for x in MOD.normalized_cfg_diff(a, b)],
                             ["ala_shallow_patch_preconditioner", "gamma"])

    def test_stage_b_exact_near_null_and_broken_transpose(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td); adj = td / "adj.csv"; gauge = td / "gauge.csv"; out = td / "out.json"
            fields = ["layout", "operator", "absolute_defect", "delta_scale",
                      "scale_floor", "denom_lr"]
            with adj.open("w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
                for op in ("D", "C", "B"):
                    w.writerow({"layout": "4x4x2", "operator": op,
                                "absolute_defect": "1e-16", "delta_scale": "1e-16",
                                "scale_floor": "1e-15", "denom_lr": "1"})
            gauge.write_text("probe,BTq_norm,S_gamma_q_norm,scale_floor\nconstant,1e-16,1e-16,1e-15\n")
            args = type("A", (), {"adjoint": str(adj), "gauge": str(gauge),
                                   "output": str(out)})
            self.assertEqual(MOD.stage_b_decision(args), 0)
            self.assertEqual(json.loads(out.read_text())["D"], "PASS")
            with adj.open() as f: rows = list(csv.DictReader(f))
            rows[0]["absolute_defect"] = "1e-3"; rows[0]["delta_scale"] = "1e-2"
            with adj.open("w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)
            self.assertEqual(MOD.stage_b_decision(args), 1)
            self.assertEqual(json.loads(out.read_text())["D"], "TRUE_DEFECT")

    def test_c_schemas_and_default_off_contract_are_present(self):
        stokes = (ROOT / "lib" / "Stokes_flow_Incomp.c").read_text()
        gm = (ROOT / "lib" / "General_matrix_functions.c").read_text()
        ins = (ROOT / "lib" / "Instructions.c").read_text()
        self.assertIn("continuity_numerator,continuity_denominator", stokes)
        self.assertIn("requested_relative_tolerance,target_absolute", gm)
        self.assertIn('input_boolean("ala_stage_abc_production_logging"', ins)
        self.assertIn("ala_stage_abc_production_logging = 0", ins)
        self.assertIn(": solve_del2_u(E,tmpU,tmpF,inner_accuracy,lev)", stokes)

    def test_observer_uses_independent_force_storage(self):
        stokes = (ROOT / "lib" / "Stokes_flow_Incomp.c").read_text()
        element = (ROOT / "lib" / "Element_calculations.c").read_text()
        self.assertIn("assemble_forces_into(E,0,force_work)", stokes)
        helper = element[element.index("void assemble_forces_into"):]
        helper = helper[:helper.index("void assemble_forces_pseudo_surf")]
        self.assertNotIn("E->F[", helper)
        self.assertNotIn("get_buoyancy(", helper)


if __name__ == "__main__":
    unittest.main()
