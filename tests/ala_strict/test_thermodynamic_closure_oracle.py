#!/usr/bin/env python3
"""Tests for the strict-ALA Stage 0 thermodynamic closure oracle."""

from __future__ import annotations

import math
import subprocess
import sys
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from thermodynamic_closure_oracle import (
    BETA_SOURCE_OCCURRENCE_CONTRACT,
    BOUNDARIES_KM,
    BLOCKED,
    FAIL,
    IMPLEMENTED,
    NOT_IMPLEMENTED,
    PASS,
    TRAJECTORY_CLOSURE_TOLERANCE_K,
    TRAJECTORY_INTEGRATION_TOLERANCE_K,
    GLOBAL_ROOT,
    StrictALAClosureOracle,
    _compiled_library_sources,
    scan_beta_source_occurrences,
)


class ThermodynamicClosureOracleTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.oracle = StrictALAClosureOracle()
        cls.report = cls.oracle.run()

    def test_production_source_contracts_map_external_temperature_mutation(self) -> None:
        self.assertEqual(
            self.report["source_contracts_scope"],
            "PREREQUISITE STATIC ENERGY MAP; no phase parameter loading, "
            "phase-law evaluation, or phase experiment",
        )
        self.assertTrue(all(self.report["source_contracts"].values()))
        self.assertTrue(
            self.report["source_contracts"]["rheo_dat_can_mutate_temperature"]
        )
        self.assertTrue(
            self.report["source_contracts"]["SUPG_active_with_zero_diffusion"]
        )

    def test_manufactured_velocity_conserves_reference_mass(self) -> None:
        result = self.report["manufactured_mass_flux"]
        self.assertLessEqual(
            result["maximum_absolute_residual"], 4.0 * np.finfo(float).eps
        )

    def test_temperature_and_phase_ownership_are_explicit(self) -> None:
        ownership = self.report["temperature_ownership"]
        self.assertEqual(len(ownership["T_total"]), len(self.oracle.radius))
        self.assertLessEqual(max(map(abs, ownership["deltaT"])), np.finfo(float).eps)
        self.assertIsNone(ownership["phases"])
        self.assertEqual(ownership["phase_diagnostics_status"], BLOCKED)
        self.assertEqual(
            set(ownership["energy_contributions"]),
            {
                "R_storage",
                "R_advection",
                "R_adiabatic",
                "R_phase",
                "R_diffusion",
                "R_viscous",
                "R_internal",
                "R_boundary",
                "R_stabilization",
                "R_total",
            },
        )
        self.assertIn("NO FE/SUPG", ownership["scope"])

    def test_reduced_smooth_trajectories_evolve_and_fail_closure(self) -> None:
        self.assertEqual(
            self.report["gates"]["smooth_reduced_trajectory"], FAIL
        )
        for result in self.report["smooth_closure"]:
            self.assertEqual(result["implementation_status"], IMPLEMENTED)
            self.assertEqual(result["result"], FAIL)
            self.assertGreater(result["elapsed_nondimensional_time"], 0.0)
            self.assertGreater(result["state_changed_K"], 1.0)
            self.assertGreater(
                result["maximum_deltaT_K"], TRAJECTORY_CLOSURE_TOLERANCE_K
            )
            self.assertLessEqual(
                result["integration_error_estimate_K"],
                TRAJECTORY_INTEGRATION_TOLERANCE_K,
            )
            self.assertGreater(result["reduced_residual_norms"]["RMS"], 0.0)
            self.assertIn(
                result["largest_reduced_component_by_unweighted_RMS"],
                result["reduced_component_norms"],
            )

    def test_rk4_trajectory_integrator_matches_analytic_state_change(self) -> None:
        decay = 0.7
        initial = 2.5
        start = 0.6
        stop = 0.9
        with mock.patch.object(
            self.oracle,
            "_material_temperature_gradient",
            side_effect=lambda radius, temperature: -decay * temperature,
        ):
            final, path = self.oracle._integrate_segment(
                start, stop, initial, substeps=400
            )
        expected = initial * math.exp(-decay * (stop - start))
        self.assertGreater(abs(final - initial), 0.1)
        self.assertAlmostEqual(final, expected, delta=2.0e-14)
        self.assertEqual(len(path), 400)

    def test_smooth_trajectories_use_complete_same_branch_elements(self) -> None:
        expected_nodes = ([42, 64], [40, 41], [38, 39], [2, 37])
        expected_depths = (
            (398.697180, 14.564106),
            (497.912763, 445.116286),
            (622.618717, 557.080240),
            (2867.695407, 694.534565),
        )
        for result, nodes, depths in zip(
            self.report["smooth_closure"], expected_nodes, expected_depths
        ):
            self.assertEqual(result["mesh_node_indices_one_based"], list(nodes))
            self.assertEqual(result["selected_phase_boundaries_km"], [])
            self.assertEqual(
                result["trajectory_selection"],
                "COMPLETE CONTIGUOUS SAME-BRANCH MESH ELEMENTS",
            )
            self.assertAlmostEqual(result["trajectory_depth_range_km"][0], depths[0], places=5)
            self.assertAlmostEqual(result["trajectory_depth_range_km"][1], depths[1], places=5)
            first, last = (value - 1 for value in nodes)
            selected_depth = (
                self.oracle.depth[first : last + 1]
                * self.oracle.radius_m
                / 1000.0
            )
            for boundary in BOUNDARIES_KM:
                self.assertFalse(selected_depth[-1] <= boundary <= selected_depth[0])
            self.assertIn(
                "sum_e R_e*Delta_r_e",
                result["integrated_reduced_residual_definition"],
            )

    def test_clean_smooth_results_replace_contaminated_values(self) -> None:
        maxima = [item["maximum_deltaT_K"] for item in self.report["smooth_closure"]]
        self.assertTrue(all(np.isfinite(maximum) for maximum in maxima))
        self.assertTrue(all(0.0 < maximum < 0.2 for maximum in maxima))

    def test_reduced_projection_does_not_claim_production_mpi(self) -> None:
        audit = self.report["thermal_buoyancy_equivalence"]
        for result in audit["results"]:
            self.assertLessEqual(result["difference"]["Linf"], 2.0e-15)
        production = audit["production_MPI"]
        self.assertEqual(production["implementation_status"], NOT_IMPLEMENTED)
        self.assertIsNone(production["result"])
        self.assertIsNone(production["measurements"])
        self.assertEqual(production["required_world_ranks"], 384)
        self.assertEqual(production["required_horizontal_comm_size"], 192)

    def test_zero_reference_reports_raw_and_normalized_scales(self) -> None:
        audit = self.report["zero_reference_buoyancy"]
        self.assertEqual(audit["implementation_status"], NOT_IMPLEMENTED)
        self.assertEqual(audit["execution_status"], BLOCKED)
        self.assertIsNone(audit["measurements"])
        self.assertEqual(
            self.report["gates"]["serial_zero_reference_coefficient_space"],
            BLOCKED,
        )

    def test_beta_serialization_and_exact_runtime_authority_close(self) -> None:
        audit = self.report["beta_consistency"]
        self.assertEqual(audit["selected_cfg_source"], "interval")
        self.assertLessEqual(
            audit["serialized_vs_regenerated"]["Linf"], 2.0e-14
        )
        self.assertTrue(all(audit["source_contracts"].values()))
        self.assertEqual(audit["unexpected_runtime_authorities"], [])
        self.assertEqual(audit["missing_expected_occurrences"], [])
        self.assertEqual(
            audit["source_scan"]["expected_occurrence_count"],
            len(BETA_SOURCE_OCCURRENCE_CONTRACT),
        )
        self.assertEqual(
            audit["source_scan"]["actual_occurrence_count"],
            len(BETA_SOURCE_OCCURRENCE_CONTRACT),
        )
        self.assertTrue(audit["source_scan"]["contract_matches"])
        self.assertTrue(audit["source_scan"]["declaration_inventory_matches"])
        self.assertTrue(
            all(
                occurrence["allowlisted_role"] is not None
                for occurrence in audit["source_scan"]["occurrences"]
            )
        )
        fields = {item["field"]: item for item in audit["inventory"]}
        self.assertEqual(
            set(fields),
            {
                "beta_ala",
                "ala_beta_supplied",
                "ala_beta_density",
                "ala_beta_interval",
                "ala_beta",
                "gamma_eff",
            },
        )
        self.assertTrue(fields["ala_beta_interval"]["selected"])
        self.assertTrue(fields["ala_beta"]["selected"])
        self.assertFalse(fields["gamma_eff"]["selected"])

    def test_beta_occurrence_contract_rejects_new_use_in_allowed_file(self) -> None:
        sources = {
            str(path.relative_to(GLOBAL_ROOT)): path.read_text(encoding="utf-8")
            for path in _compiled_library_sources()
        }
        sources["lib/Material_properties.c"] += (
            "\ndouble hidden = E->refstate.ala_beta[1];\n"
        )
        result = scan_beta_source_occurrences(sources)
        self.assertFalse(result["contract_matches"])
        self.assertEqual(
            result["actual_occurrence_count"],
            len(BETA_SOURCE_OCCURRENCE_CONTRACT) + 1,
        )
        self.assertEqual(result["missing_occurrences"], [])
        self.assertEqual(len(result["unexpected_occurrences"]), 1)

    def test_beta_occurrence_contract_is_multiset_and_line_number_independent(self) -> None:
        line = "beta += E->refstate.ala_beta[fine_nz] * dr;"
        baseline = scan_beta_source_occurrences(
            {"allowed.c": "\n\n" + line + "\n"},
            expected_contract=(
                (
                    "ala_beta",
                    "allowed.c",
                    "9cb7075ceb63446dab74fb1b80a3323a0b939991659e908b1426a3c2102e5d76",
                ),
            ),
        )
        self.assertTrue(baseline["contract_matches"])
        duplicate = scan_beta_source_occurrences(
            {"allowed.c": line + "\n" + line + "\n"},
            expected_contract=(
                (
                    "ala_beta",
                    "allowed.c",
                    "9cb7075ceb63446dab74fb1b80a3323a0b939991659e908b1426a3c2102e5d76",
                ),
            ),
        )
        self.assertFalse(duplicate["contract_matches"])
        self.assertEqual(duplicate["unexpected_occurrences"][0]["count"], 1)

    def test_beta_occurrence_contract_reports_missing_and_ignores_comments(self) -> None:
        expected = (
            (
                "ala_beta",
                "allowed.c",
                "9cb7075ceb63446dab74fb1b80a3323a0b939991659e908b1426a3c2102e5d76",
            ),
        )
        missing = scan_beta_source_occurrences(
            {"allowed.c": "/* E->refstate.ala_beta[1] */\n"},
            expected_contract=expected,
        )
        self.assertFalse(missing["contract_matches"])
        self.assertEqual(missing["actual_occurrence_count"], 0)
        self.assertEqual(missing["missing_occurrences"][0]["count"], 1)

    def test_phase_experiments_are_blocked_without_placeholder_measurements(self) -> None:
        self.assertEqual(self.report["gates"]["phase_crossings"], BLOCKED)
        ledger = self.report["phase_energy_ledger"]
        self.assertEqual(ledger["implementation_status"], NOT_IMPLEMENTED)
        self.assertEqual(ledger["execution_status"], BLOCKED)
        self.assertIsNone(ledger["values"])
        for collection in (
            self.report["phase_closure"],
            self.report["count_once"],
            self.report["forward_reverse"],
        ):
            for experiment in collection.values():
                self.assertEqual(experiment["implementation_status"], NOT_IMPLEMENTED)
                self.assertEqual(experiment["execution_status"], BLOCKED)
                self.assertIsNone(experiment["measurements"])
                self.assertIsNone(experiment["result"])

    def test_smooth_failure_short_circuits_phase_ledger(self) -> None:
        smooth = [dict(item, result=FAIL) for item in self.report["smooth_closure"]]
        with mock.patch.object(
            self.oracle, "smooth_closure", return_value=smooth
        ), mock.patch.object(
            self.oracle,
            "_load_phase_parameters",
            side_effect=AssertionError("phase parameters must not be loaded"),
        ) as phase_source, mock.patch.object(
            self.oracle,
            "phase_fraction",
            side_effect=AssertionError("phase laws must not execute"),
        ):
            report = self.oracle.run()
        phase_source.assert_not_called()
        self.assertEqual(
            report["phase_energy_ledger"]["execution_status"],
            BLOCKED,
        )

    def test_smooth_gate_does_not_sample_phase_thermodynamics(self) -> None:
        with mock.patch.object(
            self.oracle,
            "phase_fraction",
            side_effect=AssertionError("smooth gate must not sample phase laws"),
        ):
            results = self.oracle.smooth_closure()
        self.assertEqual(len(results), 4)
        for result in results:
            self.assertGreater(
                result["trajectory_depth_range_km"][0],
                result["trajectory_depth_range_km"][1],
            )

    def test_manifest_freezes_inputs_and_full_sphere_topology(self) -> None:
        manifest = self.report["manifest"]
        self.assertEqual(manifest["topology"]["caps"], 12)
        self.assertEqual(manifest["topology"]["solver"], "full")
        self.assertIn("FullSphere.py", manifest["topology"]["provenance"])
        self.assertEqual(manifest["topology"]["world_ranks"], 384)
        self.assertEqual(manifest["topology"]["horizontal_comm_size"], 192)
        self.assertEqual(manifest["topology"]["horizontal_communicators"], 2)
        self.assertEqual(len(manifest["source"]["commit"]), 40)
        self.assertEqual(len(manifest["runs"]["commit"]), 40)
        self.assertEqual(len(manifest["source"]["tracked_diff_sha256"]), 64)
        self.assertEqual(len(manifest["runs"]["tracked_diff_sha256"]), 64)
        for value in manifest["files"].values():
            self.assertEqual(len(value["sha256"]), 64)
            self.assertTrue(Path(value["path"]).is_file())
        self.assertGreater(len(manifest["inspected_source_sha256"]), 0)
        self.assertGreater(len(manifest["test_artifacts"]), 0)

    def test_final_verdict_does_not_overclaim_missing_evidence(self) -> None:
        self.assertEqual(self.report["verdict"], "UNRESOLVED")
        self.assertEqual(
            self.report["gates"]["smooth_production_FE_trajectory"],
            NOT_IMPLEMENTED,
        )
        self.assertEqual(
            self.report["gates"]["production_MPI_projection"], NOT_IMPLEMENTED
        )
        self.assertIsNone(self.report["internal_gate_state"]["full_fe_smooth_pass"])
        for evidence in self.report["unavailable_production_evidence"].values():
            self.assertEqual(evidence["implementation_status"], NOT_IMPLEMENTED)
            self.assertEqual(evidence["execution_status"], "NOT-RUN")
            self.assertIsNone(evidence["result"])
            self.assertIsNone(evidence["measurements"])
        self.assertFalse(self.report["production_files_modified"])
        self.assertIn("test-only", self.report["recommended_next_action"].lower())

    def test_cli_fails_closed_unless_unresolved_report_is_explicitly_allowed(self) -> None:
        script = Path(__file__).with_name("thermodynamic_closure_oracle.py")
        failed = subprocess.run(
            (sys.executable, "-B", str(script)), capture_output=True, text=True
        )
        self.assertEqual(failed.returncode, 2)
        self.assertIn("verdict: UNRESOLVED", failed.stdout)
        allowed = subprocess.run(
            (sys.executable, "-B", str(script), "--allow-unresolved"),
            capture_output=True,
            text=True,
        )
        self.assertEqual(allowed.returncode, 0)

    def test_cli_classifies_input_failure_as_unresolved_without_traceback(self) -> None:
        script = Path(__file__).with_name("thermodynamic_closure_oracle.py")
        missing = Path(__file__).with_name("missing-stage0-mesh.dat")
        result = subprocess.run(
            (sys.executable, "-B", str(script), "--mesh", str(missing)),
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("verdict: UNRESOLVED", result.stdout)
        self.assertIn("oracle error: FileNotFoundError", result.stdout)
        self.assertNotIn("Traceback", result.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)
