#!/usr/bin/env python3
"""Static contracts for the frozen-state strict-ALA Schur diagnostic."""

from __future__ import annotations

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
            "for(probe=0;probe<ALA_SCHUR_PROBE_COUNT;probe++)", self.stokes
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

    def test_old_rheology_switch_is_internal_only(self) -> None:
        self.assertIn("ala_schur_diagnostic_viscosity_mode==1", self.viscosity)
        self.assertNotIn(
            'input_int("ala_schur_diagnostic_viscosity_mode"', self.instructions
        )


if __name__ == "__main__":
    unittest.main()
