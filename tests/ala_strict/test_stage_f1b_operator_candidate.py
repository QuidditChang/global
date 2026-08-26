#!/usr/bin/env python3
import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT.parents[1] / "runs"
SPEC = importlib.util.spec_from_file_location(
    "f1b", ROOT / "tools/analyze_strict_ala_stage_F1b.py")
F1B = importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(F1B)


class Args:
    def __init__(self, **values): self.__dict__.update(values)


def write_cfg(path, selector):
    values = {
        "ala_shallow_patch_local_operator": selector,
        "ala_shallow_patch_preconditioner": "on",
        "ala_shallow_patch_velocity_solver": "element_vanka",
        "ala_augmented_lagrangian_gamma": "10.0",
        "ala_outer_solver": "fgmres", "ala_pcg_restart_interval": "50",
        "ala_stage_abc_production_logging": "on", "piterations": "30",
        "nprocx": "4", "nprocy": "4", "nprocz": "2",
        "nodex": "129", "nodey": "129", "nodez": "65", "steps": "1"}
    Path(path).write_text("".join(f"{key} = {value}\n" for key, value in values.items()))


def write_stage(path, case, terminal_cont=0.5, terminal_mom=1e-4):
    fields = ("case", "iteration", "true_continuity_relative",
              "true_momentum_relative", "recursive_residual",
              "explicit_residual", "velocity_MG_cycles",
              "K_gamma_operator_applications", "elapsed_seconds")
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader()
        for iteration in range(31):
            fraction = iteration / 30.0
            writer.writerow({"case": case, "iteration": iteration,
                "true_continuity_relative": 1.0 + fraction * (terminal_cont - 1.0),
                "true_momentum_relative": terminal_mom,
                "recursive_residual": 1.0 - 0.5 * fraction,
                "explicit_residual": 1.0 - 0.5 * fraction,
                "velocity_MG_cycles": 100 + 10 * iteration,
                "K_gamma_operator_applications": iteration,
                "elapsed_seconds": 1 + iteration})


def write_memory(path, case, rss):
    Path(path).write_text(
        "case,rank_rss_max_kib,rank_rss_sum_kib,ranks,source\n"
        f"{case},{rss},{rss*384},384,getrusage_RUSAGE_SELF_MPI_reduce\n")


class F1bContractTests(unittest.TestCase):
    @staticmethod
    def audit_args(td, offline, output):
        evidence = {}
        for name in ("binary", "inputs", "source"):
            before, after = td/f"{name}.pre", td/f"{name}.post"
            before.write_text(f"frozen-{name}\n"); after.write_bytes(before.read_bytes())
            evidence[f"{name}_pre"] = before
            evidence[f"{name}_post"] = after
        return Args(offline=offline, short=None, output=output, **evidence)

    def test_source_keeps_legacy_default_and_outer_projection(self):
        instructions = (ROOT / "lib/Instructions.c").read_text()
        source = (ROOT / "lib/Stokes_flow_Incomp.c").read_text()
        candidate = (ROOT / "lib/Strict_ala_stage_f1b.inc").read_text()
        self.assertIn('ala_shallow_patch_local_operator,"legacy"', instructions)
        self.assertIn('"operator_consistent"', instructions)
        self.assertIn("ala_f1b_build_operator_consistent_patch", source)
        self.assertIn("direct-face plus two-hop corner element records", candidate)
        self.assertIn("Kgamma_replay_relative_defect", candidate)
        # Existing pre/post Q_i and multiplicity path remains in the common
        # apply routine; the candidate changes only cached A_i factors.
        self.assertGreaterEqual(source.count("rhs[i] -= sum*constant_mode[i]"), 2)
        self.assertGreaterEqual(source.count("solution[i] -= sum*constant_mode[i]"), 2)

    def test_cfg_diff_is_exactly_one_selector(self):
        with tempfile.TemporaryDirectory() as td:
            base, cand, diff = (Path(td) / name for name in ("base", "cand", "diff"))
            write_cfg(base, "legacy"); write_cfg(cand, "operator_consistent")
            self.assertEqual(F1B.check_cfg(Args(base=base, candidate=cand, diff=diff)), 0)
            self.assertEqual(diff.read_text().strip(),
                             "ala_shallow_patch_local_operator: legacy -> operator_consistent")

    def test_lsf_times_executable_and_preserves_manifest_templates(self):
        lsf = (RUNS / "cmbhf_ALA_strict_stage_F1b.lsf").read_text()
        self.assertNotIn('/usr/bin/time -v -o "${B}/case_time.txt" launch_case', lsf)
        self.assertIn('/usr/bin/time -v -o "${time_file}"', lsf)
        self.assertIn('/python/pythia-0.8.1.15-py2.6.egg:', lsf)
        self.assertIn('if [ "${rc}" -ne 0 ]; then return "${rc}"; fi', lsf)
        self.assertIn('stage_inputs "${O}" "${OFFLINE_CFG}"', lsf)
        self.assertNotIn('stage_inputs "${O}" "${CAND_CFG}"', lsf)
        self.assertIn('completion_valid "${B}/phase3_base_completion.json"', lsf)
        self.assertIn('--binary-pre "${M}/binary.pre.sha256"', lsf)

    def test_hpc_evidence_root_is_valid_lineage(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td); evidence = td/"f1a"
            (evidence/"01_operator_replay").mkdir(parents=True)
            (evidence/"01_operator_replay/strict_ala_stage_F1a_projected_matrices.csv").write_text("x\n")
            final = {"experiment_valid": True,
                     "primary_defect_class": "LOCAL_OPERATOR_MISMATCH",
                     "dominant_damage_stage": "RAW_LOCAL", "F1b_authorized": True,
                     "next_authorized_task": "F1B_LOCAL_CANDIDATE",
                     "numerical_evidence_root": str(evidence),
                     "weighted_telescoping_damage": {"RAW_LOCAL": 1.0}}
            (td/"final.json").write_text(json.dumps(final))
            F1B.verify_lineage(Args(final_audit=td/"final.json",
                evidence_root=evidence, expected_hpc_root=str(evidence),
                output=td/"lineage.json"))
            self.assertTrue(json.loads((td/"lineage.json").read_text())["valid"])

    def test_short_routes_use_discrete_same_state_rows(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td); write_stage(td/"base.csv", "BASE", 0.5)
            write_stage(td/"cand.csv", "CAND", 0.35)
            write_memory(td/"bm.csv", "BASE", 1000)
            write_memory(td/"cm.csv", "CAND", 1100)
            (td/"bw").write_text("40\n"); (td/"cw").write_text("42\n")
            F1B.short(Args(base_iterations=td/"base.csv", cand_iterations=td/"cand.csv",
                base_memory=td/"bm.csv", cand_memory=td/"cm.csv",
                base_wall=td/"bw", cand_wall=td/"cw", output_dir=td/"out"))
            result = json.loads((td/"out/strict_ala_stage_F1b_short_decision.json").read_text())
            self.assertTrue(result["route_A_pass"])
            self.assertTrue(result["memory_gate_pass"])
            self.assertTrue(result["short_AB_pass"])
            with (td/"out/strict_ala_stage_F1b_short_iterations.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(rows[0]["fgmres_candidate_iteration"], "0")
            self.assertEqual(rows[-1]["fgmres_candidate_iteration"], "30")

    def test_scientific_rejection_is_not_invalid_experiment(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            offline = {"experiment_evidence_valid": True, "offline_gate_pass": False,
                       "gates": {"candidate_isolatable": True, "memory": True,
                                 "assembly": True, "local_RMS": True,
                                 "local_heavy_mode": True, "H_RAW": False}}
            (td/"offline.json").write_text(json.dumps(offline))
            F1B.final_audit(self.audit_args(td, td/"offline.json", td/"audit.json"))
            result = json.loads((td/"audit.json").read_text())
            self.assertTrue(result["experiment_valid"])
            self.assertFalse(result["candidate_valid"])
            self.assertEqual(result["decision"], "CANDIDATE_REJECTED_RAW_LOCAL")

    def test_provenance_drift_fails_closed(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            offline = {"experiment_evidence_valid": True, "offline_gate_pass": False,
                       "gates": {"candidate_isolatable": True, "memory": True,
                                 "assembly": True, "local_RMS": True,
                                 "local_heavy_mode": True, "H_RAW": False}}
            (td/"offline.json").write_text(json.dumps(offline))
            args = self.audit_args(td, td/"offline.json", td/"audit.json")
            Path(args.binary_post).write_text("drift\n")
            F1B.final_audit(args)
            result = json.loads((td/"audit.json").read_text())
            self.assertFalse(result["experiment_valid"])
            self.assertEqual(result["decision"], "INVALID_EXPERIMENT")
            self.assertEqual(result["next_authorized_task"],
                             "REPEAT_F1B_VALID_EXPERIMENT")


if __name__ == "__main__": unittest.main()
