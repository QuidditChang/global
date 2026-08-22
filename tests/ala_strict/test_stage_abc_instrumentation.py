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
    def write_complete_stage_b(self, adj, gauge):
        adj_fields = ["layout", "operator", "probe", "velocity_probe"] + \
            list(MOD.ADJOINT_NUMERIC_FIELDS)
        with adj.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=adj_fields); w.writeheader()
            for op in MOD.EXPECTED_OPERATORS:
                for probe in MOD.EXPECTED_PRESSURE_PROBES:
                    for velocity in MOD.EXPECTED_VELOCITY_PROBES:
                        row = {k: "1e-16" for k in MOD.ADJOINT_NUMERIC_FIELDS}
                        row.update({"layout": "4x4x2", "operator": op,
                                    "probe": probe, "velocity_probe": velocity,
                                    "denom_lr": "1", "denom_action": "1",
                                    "scale_floor": "1e-15"})
                        w.writerow(row)
        gauge_fields = ["layout", "probe"] + list(MOD.GAUGE_NUMERIC_FIELDS)
        with gauge.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=gauge_fields); w.writeheader()
            for probe in MOD.EXPECTED_GAUGE_PROBES:
                row = {k: "1e-16" for k in MOD.GAUGE_NUMERIC_FIELDS}
                row.update({"layout": "4x4x2", "probe": probe,
                            "q_norm": "1", "DTq_norm": "1", "CTq_norm": "1",
                            "BTq_norm": "1e-16", "S_gamma_q_norm": "1e-16",
                            "ratio_B": "1e-16", "scale_floor": "1e-15"})
                w.writerow(row)

    def test_cfg_input_discovery_ignores_non_file_controls(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            exact = td / "refstate.txt"
            exact.write_text("reference state\n")
            prefix = td / "velocity"
            (td / "velocity.0").write_text("boundary data\n")
            cfg = td / "case.cfg"
            cfg.write_text(
                "profile_optional = temperature,pressure,velocity\n"
                "datafile = global\n"
                "logfile = global_AhatP%%P\n"
                "file_vbcs = 1\n"
                "refstate_file = refstate.txt\n"
                "vel_bound_file = %s\n" % prefix
            )
            self.assertEqual(
                set(MOD.cfg_input_files(cfg)),
                {exact.resolve(), (td / "velocity.0").resolve()},
            )

    def test_cfg_input_discovery_reports_missing_required_prefix(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td); cfg = td / "case.cfg"; missing = []
            cfg.write_text("plate_file = absent_prefix\n")
            self.assertEqual(MOD.cfg_input_files(cfg, missing), [])
            self.assertEqual(missing, [str(td / "absent_prefix")])

    def test_provenance_distinguishes_generated_and_scientific_dirty_files(self):
        status = (" M Makefile.in\n M lib/Stokes_flow_Incomp.c\n"
                  "?? m4/libtool.m4\n?? tools/new_audit.py\n")
        self.assertEqual(
            MOD.scientific_dirty_paths(status),
            ["lib/Stokes_flow_Incomp.c", "tools/new_audit.py"],
        )

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
            self.write_complete_stage_b(adj, gauge)
            args = type("A", (), {"adjoint": str(adj), "gauge": str(gauge),
                                   "expected_layout": ["4x4x2"],
                                   "output": str(out)})
            self.assertEqual(MOD.stage_b_decision(args), 0)
            self.assertEqual(json.loads(out.read_text())["D"], "PASS")
            with adj.open() as f: rows = list(csv.DictReader(f))
            rows[0]["absolute_defect"] = "1e-3"; rows[0]["delta_scale"] = "1e-2"
            with adj.open("w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
            self.assertEqual(MOD.stage_b_decision(args), 1)
            self.assertEqual(json.loads(out.read_text())["D"], "TRUE_DEFECT")

    def test_stage_b_fails_closed_on_missing_duplicate_inf_and_bt_split(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td); adj = td / "adj.csv"; gauge = td / "gauge.csv"; out = td / "out.json"
            args = type("A", (), {"adjoint": str(adj), "gauge": str(gauge),
                                   "expected_layout": ["4x4x2"], "output": str(out)})
            self.write_complete_stage_b(adj, gauge)
            with adj.open() as f: rows = list(csv.DictReader(f))
            with adj.open("w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows[:-1])
            self.assertEqual(MOD.stage_b_decision(args), 1)
            self.assertFalse(json.loads(out.read_text())["adjoint_complete"])
            self.write_complete_stage_b(adj, gauge)
            with adj.open() as f: rows = list(csv.DictReader(f))
            rows[0]["scale_floor"] = "inf"
            with adj.open("w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
            self.assertEqual(MOD.stage_b_decision(args), 1)
            self.write_complete_stage_b(adj, gauge)
            with gauge.open() as f: grows = list(csv.DictReader(f))
            grows[0]["BT_split_absolute"] = "1e-8"
            grows[0]["BT_split_relative"] = "1e-6"
            with gauge.open("w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=grows[0].keys()); w.writeheader(); w.writerows(grows)
            self.assertEqual(MOD.stage_b_decision(args), 1)
            self.assertEqual(json.loads(out.read_text())["B"], "TRUE_DEFECT")

    def test_stage_c_requires_both_terminal_cases_and_compares_cost(self):
        fields = ["case", "iteration", "final_iterate", "krylov_recursive",
                  "krylov_explicit", "krylov_drift", "continuity_numerator",
                  "continuity_denominator", "continuity_relative",
                  "momentum_numerator", "momentum_denominator",
                  "momentum_relative", "cumulative_inner_solves",
                  "cumulative_inner_cycles", "elapsed_seconds"]
        with tempfile.TemporaryDirectory() as td:
            td = Path(td); it = td / "iterations.csv"; inner = td / "inner.csv"; out = td / "out.json"
            rows = []
            for case, cycles, seconds in (("C0", 100, 100), ("C1", 80, 80)):
                row = {k: "1e-4" for k in fields}
                row.update({"case": case, "iteration": "1", "final_iterate": "1",
                            "continuity_relative": "1e-4", "momentum_relative": "1e-4",
                            "cumulative_inner_solves": "1",
                            "cumulative_inner_cycles": str(cycles),
                            "elapsed_seconds": str(seconds)})
                rows.append(row)
            with it.open("w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)
            inner.write_text("case\nC0\nC1\n")
            args = type("A", (), {"iterations": str(it), "inner": str(inner),
                                   "output": str(out)})
            self.assertEqual(MOD.stage_c_decision(args), 0)
            self.assertEqual(json.loads(out.read_text())["decision"], "CONFIGURED_WINS")
            with it.open("w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerow(rows[0])
            inner.write_text("case\nC0\n")
            self.assertEqual(MOD.stage_c_decision(args), 1)
            self.assertEqual(json.loads(out.read_text())["decision"], "INVALID_EXPERIMENT")

    def test_final_audit_rejects_missing_c1_and_nan_inner_evidence(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            iteration = {
                "case": "C0", "iteration": "1", "final_iterate": "1",
                "krylov_recursive": "1e-4", "krylov_explicit": "1e-4",
                "krylov_drift": "1", "continuity_numerator": "1e-4",
                "continuity_denominator": "1", "continuity_relative": "1e-4",
                "momentum_numerator": "1e-4", "momentum_denominator": "1",
                "momentum_relative": "1e-4", "cumulative_inner_solves": "1",
                "cumulative_inner_cycles": "10", "cumulative_K_gamma_applications": "2",
                "cumulative_schur_actions": "1", "cumulative_preconditioner_applications": "1",
                "bpi_construction_seconds": "1", "schwarz_construction_seconds": "0",
                "elapsed_seconds": "2",
            }
            iterations = td / "iterations.csv"
            with iterations.open("w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=iteration.keys()); w.writeheader(); w.writerow(iteration)
            inner_row = {
                "case": "C0", "call_id": "1", "rhs_rms": "1",
                "requested_relative_tolerance": "1e-2", "target_absolute": "1e-2",
                "achieved_absolute": "1e-3", "achieved_relative": "1e-3",
                "cycles": "10", "max_cycles": "100", "seconds": "nan",
                "status": "CONVERGED",
            }
            inner = td / "inner.csv"
            with inner.open("w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=inner_row.keys()); w.writeheader(); w.writerow(inner_row)
            manifest = td / "manifest.json"
            manifest.write_text(json.dumps({"provenance_complete": True,
                                            "files": [], "missing": []}))
            stage_b = td / "stage_b.json"
            stage_b.write_text(json.dumps({"D": "PASS", "C": "PASS", "B": "PASS",
                                           "stage_C": "ALLOWED", "adjoint_complete": True,
                                           "gauge_complete": True}))
            decision = td / "decision.json"
            decision.write_text(json.dumps({"C0": {"viable": True},
                                            "C1": {"viable": False},
                                            "decision": "INVALID_EXPERIMENT",
                                            "structure_valid": False}))
            cfg_diff = td / "cfg.diff"
            cfg_diff.write_text("ala_shallow_patch_preconditioner: off -> on\n")
            pre = td / "pre"; post = td / "post"; pre.write_text("same\n"); post.write_text("same\n")
            logs = []
            for case in ("C0", "C1"):
                log = td / (case + ".log")
                log.write_text("STRICT_ALA_STAGE_C_CASE_COMPLETE case=%s status=CONVERGED\n" % case)
                logs.append(str(log))
            output = td / "final.json"
            args = type("A", (), {
                "iterations": str(iterations), "inner": str(inner),
                "cost": str(td / "cost.csv"), "decision": str(decision),
                "manifest": str(manifest), "stage_b_decision": str(stage_b),
                "cfg_diff": str(cfg_diff), "binary_pre": str(pre),
                "binary_post": str(post), "case_log": logs,
                "output": str(output),
            })
            self.assertEqual(MOD.final_audit(args), 1)
            self.assertEqual(json.loads(output.read_text())["decision"],
                             "INVALID_EXPERIMENT")
            iteration_c1 = dict(iteration); iteration_c1["case"] = "C1"
            with iterations.open("w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=iteration.keys()); w.writeheader()
                w.writerows([iteration, iteration_c1])
            inner_row["seconds"] = "1"
            inner_c1 = dict(inner_row); inner_c1["case"] = "C1"
            with inner.open("w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=inner_row.keys()); w.writeheader()
                w.writerows([inner_row, inner_c1])
            decision.write_text(json.dumps({"C0": {"viable": True},
                                            "C1": {"viable": True},
                                            "decision": "TIE_NEEDS_REPEAT",
                                            "structure_valid": True}))
            self.assertEqual(MOD.final_audit(args), 0)
            self.assertTrue(json.loads(output.read_text())["structural_checks_pass"])

    def test_c_schemas_and_default_off_contract_are_present(self):
        stokes = (ROOT / "lib" / "Stokes_flow_Incomp.c").read_text()
        gm = (ROOT / "lib" / "General_matrix_functions.c").read_text()
        ins = (ROOT / "lib" / "Instructions.c").read_text()
        inventory = (ROOT / "CitcomS" / "Components" / "Stokes_solver" /
                     "Incompressible.py").read_text()
        bridge = (ROOT / "module" / "setProperties.c").read_text()
        self.assertIn("continuity_numerator,continuity_denominator", stokes)
        self.assertIn("requested_relative_tolerance,target_absolute", gm)
        for control in ("ala_stage_abc_adjoint_diagnostic",
                        "ala_stage_abc_production_logging"):
            self.assertIn("%s = prop.bool(" % control, inventory)
            self.assertIn('"%s", default=False' % control, inventory)
            self.assertIn('getIntProperty(properties, "%s"' % control, bridge)
            self.assertIn('input_boolean("%s"' % control, ins)
            self.assertIn("%s = 0" % control, ins)
        self.assertIn(": solve_del2_u(E,tmpU,tmpF,inner_accuracy,lev)", stokes)

    def test_observer_uses_independent_force_storage(self):
        stokes = (ROOT / "lib" / "Stokes_flow_Incomp.c").read_text()
        element = (ROOT / "lib" / "Element_calculations.c").read_text()
        self.assertIn("assemble_forces_into(E,0,P,force_work)", stokes)
        helper = element[element.index("void assemble_forces_into"):]
        helper = helper[:helper.index("void assemble_forces_pseudo_surf")]
        self.assertNotIn("E->F[", helper)
        self.assertNotIn("get_buoyancy(", helper)
        self.assertIn("pressure[m][e]-E->P[m][e]", helper)

    def test_stage_b_probe_is_global_and_combined_coefficients_are_promoted(self):
        stokes = (ROOT / "lib" / "Stokes_flow_Incomp.c").read_text()
        probe = stokes[stokes.index("static void ala_stage_b_build_velocity"):]
        probe = probe[:probe.index("static double ala_stage_b_velocity_difference")]
        self.assertIn("E->X[lev][m][1][node]", probe)
        self.assertNotIn("78.233", probe)
        header = (ROOT / "lib" / "global_defs.h").read_text()
        self.assertIn("((double)(g) + (double)(c))", header)
        for path in (ROOT / "lib" / "Element_calculations.c",
                     ROOT / "lib" / "Construct_arrays.c", ROOT / "lib" / "Stokes_flow_Incomp.c"):
            self.assertIn("ALA_COMBINED_PRESSURE_COEFFICIENT", path.read_text())

    def test_final_audit_is_evidence_driven_and_lsf_aggregates_exactly_one(self):
        tool = (ROOT / "tools" / "strict_ala_stage_abc.py").read_text()
        final = tool[tool.index("def final_audit"):tool.index("def main")]
        self.assertNotIn('"provenance_complete": True', final)
        self.assertNotIn('"hidden_fallback_detected": False', final)
        self.assertIn("case_logs_complete", final)
        lsf = (ROOT.parents[1] / "runs" / "cmbhf_ALA_strict_stage_ABC.lsf").read_text()
        self.assertIn("find_exactly_one", lsf)
        self.assertNotIn("B_ALTERNATE", lsf)
        self.assertIn("--expected-layout 4x4x2", lsf)
        self.assertIn("STRICT_ALA_STAGE_C_CASE_COMPLETE", stokes :=
                      (ROOT / "lib" / "Stokes_flow_Incomp.c").read_text())
        self.assertIn("Strict-ALA Stage-C Schwarz fallback is forbidden", stokes)
        velocity_wrapper = stokes[
            stokes.index("static int strict_ala_stage_c_velocity_solve", 1000):
            stokes.index("static void strict_ala_stage_c_momentum_metrics")]
        self.assertNotIn("ala_stage_abc_schur_action_count++", velocity_wrapper)

    def test_stage_abc_creates_logs_on_rank_zero_only(self):
        instructions = (ROOT / "lib" / "Instructions.c").read_text()
        tracer = (ROOT / "lib" / "Full_tracer_advection.c").read_text()
        guard = 'getenv("STRICT_ALA_CASE") != NULL && E->parallel.me != 0'
        self.assertIn(guard, instructions)
        self.assertIn('E->fp = output_open("/dev/null", "w")', instructions)
        self.assertIn(guard, tracer)
        self.assertIn('sprintf(output_file,"/dev/null")', tracer)


if __name__ == "__main__":
    unittest.main()
