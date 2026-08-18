#!/usr/bin/env python3
"""Static contracts for the frozen-state strict-ALA Schur diagnostic."""

from __future__ import annotations

import re
import unittest
from pathlib import Path


GLOBAL_ROOT = Path(__file__).resolve().parents[2]
PROJECT_ROOT = Path(__file__).resolve().parents[4]
LIB_ROOT = GLOBAL_ROOT / "lib"


class SchurDiagnosticArchitectureTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.stokes = (LIB_ROOT / "Stokes_flow_Incomp.c").read_text()
        cls.element = (LIB_ROOT / "Element_calculations.c").read_text()
        cls.viscosity = (LIB_ROOT / "Viscosity_structures.c").read_text()
        cls.drive = (LIB_ROOT / "Drive_solvers.c").read_text()
        cls.instructions = (LIB_ROOT / "Instructions.c").read_text()
        cls.matrix = (LIB_ROOT / "General_matrix_functions.c").read_text()
        cls.global_defs = (LIB_ROOT / "global_defs.h").read_text()
        cls.component = (
            GLOBAL_ROOT / "CitcomS/Components/Stokes_solver/Incompressible.py"
        ).read_text()
        cls.launcher = (
            PROJECT_ROOT / "runs/cmbhf_ALA_schur_diagnostic.lsf"
        ).read_text()
        cls.config = (PROJECT_ROOT / "runs/cmbhf_ALA_strict.cfg").read_text()

    def test_diagnostic_is_default_off_and_audit_only_is_guarded(self) -> None:
        self.assertIn("E->control.ala_schur_diagnostic = 0", self.instructions)
        self.assertIn("E->control.ala_schur_diagnostic_only = 0", self.instructions)
        self.assertIn(
            "ala_schur_diagnostic_only requires ala_schur_diagnostic=on",
            self.instructions,
        )

    def test_exact_split_transposes_share_exchange_and_bc_stripping(self) -> None:
        start = self.element.index("void assemble_grad_c_p(")
        block = self.element[start:self.element.index(
            "double assemble_dAhatp_entry", start
        )]
        self.assertIn("E->elt_c[lev][m][e].c", block)
        self.assertIn("(E->solver.exchange_id_d)(E,gradP,lev)", block)
        self.assertIn("strip_bcs_from_residual(E,gradP,lev)", block)

    def test_new_old_rebuilds_and_production_mode_restores_new(self) -> None:
        start = self.stokes.index("void strict_ala_schur_diagnostic(")
        block = self.stokes[start:self.stokes.index("static void", start)]
        self.assertIn("for(state=0;state<=1;state++)", block)
        self.assertIn("get_system_viscosity(E,1", block)
        self.assertIn("construct_stiffness_B_matrix(E)", block)
        self.assertIn("production_state_restored=%s", block)
        self.assertIn("not_required_audit_only", block)
        self.assertLess(
            block.index("ala_schur_run_state(E,csv,state ? \"OLD\" : \"NEW\"") ,
            block.index("MPI_Finalize()"),
        )

    def test_p0_is_the_new_operator_momentum_consistent_residual(self) -> None:
        self.assertIn(
            "initial_vel_residual(E,initial_velocity,E->P,initial_force",
            self.stokes,
        )
        self.assertIn(
            "assemble_div_rho_u(E,initial_velocity,frozen_q0,lev)",
            self.stokes,
        )

    def test_all_probe_families_are_internal_to_one_state_loop(self) -> None:
        self.assertIn("#define ALA_SCHUR_PROBE_COUNT 19", self.stokes)
        for name in (
            "P0_initial_continuity", "P1_fixed_random", "P5_horizontal_checkerboard",
            "P6_degree_1", "P6_degree_2", "P6_degree_4", "P6_degree_8",
            "P7_depth_0_200", "P11_depth_1000_cmb", "P12_patch_scale",
            "P13_longer_than_patch", "P14_constant", "P15_density_gauge",
        ):
            self.assertIn(name, self.stokes)
        self.assertIn(
            "for(probe=first_probe;probe<ALA_SCHUR_PROBE_COUNT;probe++)",
            self.stokes,
        )

    def test_two_split_solves_reconstruct_four_schur_terms(self) -> None:
        self.assertIn("assemble_grad_p(E,q,dT,lev)", self.stokes)
        self.assertIn("assemble_grad_c_p(E,q,cT,lev)", self.stokes)
        self.assertIn("duD[m][e]+cuD[m][e]+duC[m][e]+cuC[m][e]", self.stokes)
        self.assertIn("sum_defect", self.stokes)
        self.assertIn("adj_d", self.stokes)
        self.assertIn("adj_c", self.stokes)
        self.assertIn("adj_b", self.stokes)

    def test_preconditioner_variants_and_depth_bands_are_measured(self) -> None:
        for name in (
            "BPI_only", "Schwarz_only", "combined_unscaled", "configured"
        ):
            self.assertIn(name, self.stokes)
        for boundary in ("depth<200.0", "depth<410.0", "depth<660.0", "depth<1000.0"):
            self.assertIn(boundary, self.stokes)
        self.assertIn("ala_schur_write_depth_rows", self.stokes)

    def test_driver_runs_before_production_solve_and_can_exit_cleanly(self) -> None:
        diagnostic = self.drive.index("strict_ala_schur_diagnostic(E)")
        solve = self.drive.index("solve_constrained_flow_iterative(E)")
        self.assertLess(diagnostic, solve)
        self.assertIn("MPI_Barrier(E->parallel.world)", self.stokes)
        self.assertIn("MPI_Finalize()", self.stokes)

    def test_csv_is_incremental_and_probe_progress_is_logged(self) -> None:
        self.assertIn("setvbuf(csv,NULL,_IOLBF,0)", self.stokes)
        header = self.stokes.index('fprintf(csv,"row_type,state,probe')
        self.assertLess(header, self.stokes.index("fflush(csv)", header))
        self.assertIn("STRICT_ALA_SCHUR_DIAGNOSTIC_PROBE_BEGIN", self.stokes)
        self.assertIn("STRICT_ALA_SCHUR_DIAGNOSTIC_PROBE_COMPLETE", self.stokes)

    def test_csv_writes_do_not_hide_mpi_collectives_on_rank_zero(self) -> None:
        self.assertIn("action_norm2=global_pdot(E,y,y,lev)", self.stokes)
        self.assertIn("sensitivity_energy=global_pdot(E,q,y,lev)", self.stokes)
        self.assertNotIn("global_pdot(E,y,y,lev),(qp>0.0", self.stokes)
        self.assertNotIn("global_pdot(E,q,y,lev));", self.stokes)

    def test_nonfinite_fields_are_caught_before_vanka(self) -> None:
        self.assertIn("ala_schur_check_field_finite", self.stokes)
        self.assertIn('"tight_solve_rhs"', self.stokes)
        self.assertIn('variants[variant]', self.stokes)

    def test_zero_rhs_is_a_valid_zero_velocity_action(self) -> None:
        self.assertIn("if(rhs2==0.0)", self.stokes)
        self.assertIn("ala_schur_zero_velocity(E,u,lev)", self.stokes)

    def test_destructive_multigrid_rhs_is_restored_before_audit(self) -> None:
        solve = self.stokes.index("valid=solve_del2_u_bounded")
        restore = self.stokes.index("rhs[m][i]=work[m][i]", solve)
        residual = self.stokes.index("ala_schur_velocity_residual", restore)
        self.assertLess(solve, restore)
        self.assertLess(restore, residual)

    def test_zero_upward_action_cannot_create_nan_alpha(self) -> None:
        self.assertIn("if(AudotAu<=1.0e-300)", self.matrix)
        self.assertIn("alpha=0.0", self.matrix)
        self.assertIn("Invalid multigrid upward action product", self.matrix)

    def test_exact_operator_invariants_fail_before_later_probes(self) -> None:
        self.assertIn("STRICT_ALA_SCHUR_DIAGNOSTIC_INVARIANT_FAILURE", self.stokes)
        self.assertIn("rd>1.0e-8 || rc>1.0e-8", self.stokes)
        self.assertIn("adj_b>1.0e-8", self.stokes)

    def test_csv_rows_match_the_24_column_header(self) -> None:
        formats = re.findall(r'fprintf\(csv,"([^"]+)', self.stokes)
        rows = {
            row.split(",", 1)[0]: row
            for row in formats
            if row.split(",", 1)[0] in {
                "viscosity", "viscosity_depth", "depth", "schur",
                "preconditioner", "inner_sensitivity",
            }
        }
        self.assertEqual(len(rows), 6)
        for row_type, row in rows.items():
            with self.subTest(row_type=row_type):
                self.assertEqual(row.count(",") + 1, 24)

    def test_default_off_returns_before_large_field_allocations(self) -> None:
        start = self.stokes.index("void strict_ala_schur_diagnostic(")
        block = self.stokes[start:self.stokes.index("static void", start)]
        guard = block.index("if(!E->control.ala_schur_diagnostic) return")
        allocation = block.index("frozen_q0=ala_schur_alloc_pressure")
        self.assertLess(guard, allocation)

    def test_resume_control_is_plumbed_and_default_off(self) -> None:
        self.assertIn("int ala_schur_diagnostic_resume", self.global_defs)
        self.assertIn('input_boolean("ala_schur_diagnostic_resume"', self.instructions)
        self.assertIn('"ala_schur_diagnostic_resume", default=False', self.component)
        self.assertRegex(
            self.config,
            r"(?m)^ala_schur_diagnostic_resume\s*=\s*off\s*$",
        )

    def test_checkpoint_advances_only_after_complete_probe(self) -> None:
        complete = self.stokes.index("STRICT_ALA_SCHUR_DIAGNOSTIC_PROBE_COMPLETE")
        checkpoint = self.stokes.index("ala_schur_write_checkpoint", complete)
        self.assertLess(complete, checkpoint)
        self.assertIn("version=1\\nstate=%d\\nprobe=%d\\ncsv_bytes=%lld", self.stokes)
        self.assertIn('ftruncate(fileno(csv),(off_t)checkpoint.csv_bytes)', self.stokes)
        self.assertIn("first_probe=resume_probe+1", self.stokes)

    def test_stage_reports_are_atomic_and_named(self) -> None:
        self.assertIn("strict_ala_schur.stage_%02d.json", self.stokes)
        self.assertIn('{"setup","new","old","overall"}', self.stokes)
        self.assertIn("rename(tmp,path)", self.stokes)
        self.assertIn("global.strict_ala_schur.stage_*.json", self.launcher)
        for call in (
            "ala_schur_write_stage_report(E,csv,0",
            "ala_schur_write_stage_report(E,csv,state ? 2 : 1",
            "ala_schur_write_stage_report(E,csv,3",
        ):
            self.assertIn(call, self.stokes)

    def test_launcher_collects_partial_artifacts_and_preserves_exit_status(self) -> None:
        self.assertIn("ALA_SCHUR_RUN_ROOT", self.launcher)
        self.assertIn("ALA_SCHUR_RESUME", self.launcher)
        self.assertIn("global.strict_ala_schur_diagnostic.checkpoint", self.launcher)
        self.assertIn("global.strict_ala_schur.stage_*.json", self.launcher)
        self.assertIn("APP_STATUS=$?", self.launcher)
        self.assertIn('exit "${APP_STATUS}"', self.launcher)

    def test_old_rheology_switch_is_internal_only(self) -> None:
        self.assertIn("ala_schur_diagnostic_viscosity_mode==1", self.viscosity)
        self.assertNotIn(
            'input_int("ala_schur_diagnostic_viscosity_mode"', self.instructions
        )


if __name__ == "__main__":
    unittest.main()
