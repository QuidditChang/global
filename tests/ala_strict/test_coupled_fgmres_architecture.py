#!/usr/bin/env python3
"""Static contracts for the Stage 9d monolithic strict-ALA prototype."""

from __future__ import annotations

import unittest
from pathlib import Path


GLOBAL_ROOT = Path(__file__).resolve().parents[2]
PROJECT_ROOT = Path(__file__).resolve().parents[4]
LIB_ROOT = GLOBAL_ROOT / "lib"


def _between(text: str, start: str, end: str) -> str:
    return text[text.index(start):text.index(end, text.index(start))]


class CoupledFGMRESArchitectureTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.stokes = (LIB_ROOT / "Stokes_flow_Incomp.c").read_text()
        cls.core = _between(
            cls.stokes,
            "static float solve_ala_coupled_fgmres_core(",
            "static float solve_ala_fgmres_core(",
        )

    def test_selector_is_explicit_and_legacy_fgmres_remains(self) -> None:
        self.assertIn('"coupled_fgmres"', self.stokes)
        self.assertIn("solve_ala_coupled_fgmres_core(", self.stokes)
        self.assertIn("solve_ala_fgmres_core(", self.stokes)
        pyre = (
            GLOBAL_ROOT
            / "CitcomS/Components/Stokes_solver/Incompressible.py"
        ).read_text()
        self.assertIn('"coupled_fgmres"', pyre)
        for path in (LIB_ROOT / "Instructions.c", GLOBAL_ROOT / "module/setProperties.c"):
            self.assertIn('"coupled_fgmres"', path.read_text())

    def test_full_explicit_residual_uses_pressure_independent_rhs(self) -> None:
        self.assertIn("assemble_ala_pressure_independent_force(", self.core)
        self.assertIn("assemble_ala_coupled_residual(", self.core)
        residual = _between(
            self.stokes,
            "static void assemble_ala_coupled_residual(",
            "static void apply_ala_coupled_block_preconditioner(",
        )
        self.assertIn("apply_ala_coupled_operator(", residual)
        self.assertIn("rhs->velocity[m][i]", residual)
        self.assertIn("residual->pressure[m][e]=-action->pressure[m][e]", residual)

    def test_arnoldi_metric_matches_algebraic_schur_inner_product(self) -> None:
        vectors = (LIB_ROOT / "ALA_block_vector.c").read_text()
        self.assertNotIn("/E->ECO[level][m][i].area", vectors)
        self.assertIn("pressure_algebraic", self.core)
        self.assertIn(
            "velocity_weight=1.0/max(force_norm*force_norm", self.core
        )

    def test_right_preconditioner_is_pressure_first_triangular(self) -> None:
        preconditioner = _between(
            self.stokes,
            "static void apply_ala_coupled_triangular_once(",
            "static void apply_ala_coupled_block_preconditioner(",
        )
        pressure = preconditioner.index("apply_ala_pressure_preconditioner(")
        gradient = preconditioner.index("assemble_grad_rho_p(")
        velocity = preconditioner.index("valid=solve_del2_u_bounded(")
        self.assertLess(pressure, gradient)
        self.assertLess(gradient, velocity)
        self.assertIn("apply_ala_pressure_preconditioner(", preconditioner)
        self.assertIn(
            "correction->pressure[m][e]=-correction->pressure[m][e]",
            preconditioner,
        )
        self.assertIn(
            "velocity_work[m][i]=residual->velocity[m][i]",
            preconditioner,
        )

    def test_block_defect_correction_is_multiplicative_and_reversible(self) -> None:
        preconditioner = _between(
            self.stokes,
            "static void apply_ala_coupled_block_preconditioner(",
            "static void strict_ala_coupled_preconditioner_audit(",
        )
        action = preconditioner.index("apply_ala_coupled_operator(")
        defect = preconditioner.index(
            "action_work->pressure[m][e]=residual->pressure[m][e]"
        )
        correction = preconditioner.rindex(
            "ala_block_vector_axpy(E,1.0,correction_work,correction)"
        )
        self.assertLess(action, defect)
        self.assertLess(defect, correction)
        self.assertIn("ala_coupled_defect_corrections", preconditioner)
        self.assertIn("defect_corrections=%d", self.core)
        for path in (
            LIB_ROOT / "Instructions.c",
            GLOBAL_ROOT / "module/setProperties.c",
        ):
            self.assertIn(
                "ala_coupled_defect_corrections requires", path.read_text()
            )

    def test_factor2_coarse_correction_projects_true_shallow_defect(self) -> None:
        preconditioner = _between(
            self.stokes,
            "static void apply_ala_coupled_block_preconditioner(",
            "static void strict_ala_coupled_preconditioner_audit(",
        )
        start = preconditioner.index(
            "if(E->control.ala_coupled_factor2_coarse_correction)"
        )
        block = preconditioner[start:preconditioner.index(
            "/* Multiplicative block defect correction:", start
        )]
        self.assertIn("assemble_ala_coupled_block_defect(", block)
        self.assertIn("projected_layers=2*((min(", block)
        self.assertIn("mean /= 8.0", block)
        self.assertIn("mean += vanka_delta->pressure[m][e]", block)
        self.assertIn("action_work->velocity[m][i]=0.0", block)
        self.assertIn("apply_ala_coupled_triangular_once(", block)
        self.assertIn(
            "ala_block_vector_axpy(E,1.0,correction_work,correction)", block
        )
        self.assertLess(
            block.index("assemble_ala_coupled_block_defect("),
            block.index("apply_ala_coupled_triangular_once("),
        )
        self.assertIn("projection=orthogonal_mean", block)
        self.assertIn("coarse_solver=triangular", block)
        for path in (
            LIB_ROOT / "Instructions.c",
            GLOBAL_ROOT / "module/setProperties.c",
        ):
            text = path.read_text()
            self.assertIn(
                "ala_coupled_factor2_coarse_correction requires", text
            )
            self.assertIn(
                "ala_coupled_shallow_vanka_layers > 0", text
            )

    def test_coupled_velocity_solve_is_bounded_and_observable(self) -> None:
        matrix = (LIB_ROOT / "General_matrix_functions.c").read_text()
        self.assertIn("int solve_del2_u_bounded(", matrix)
        self.assertIn(
            "while (!valid && (max_cycles<=0 || counts<max_cycles))",
            matrix,
        )
        self.assertIn("ALA COUPLED INNER VELOCITY progress", matrix)
        self.assertIn("ALA COUPLED INNER VELOCITY summary", matrix)
        self.assertIn(
            "return solve_del2_u_internal(E,d0,F,acc,high_lev,0,0)",
            matrix,
        )

    def test_arnoldi_reconstruction_updates_both_fields(self) -> None:
        self.assertIn("struct ala_block_vector **vb,**zb", self.core)
        self.assertIn("apply_ala_coupled_operator(", self.core)
        self.assertIn("ala_block_vector_dot(", self.core)
        self.assertIn("V[m][e] += delta*zb[i]->velocity[m][e]", self.core)
        self.assertIn("P[m][e] += delta*zb[i]->pressure[m][e]", self.core)
        self.assertIn("ala_block_vector_copy(E,explicit_r,r)", self.core)
        self.assertIn("drift=%e", self.core)

    def test_acceptance_is_physical_and_joint(self) -> None:
        self.assertIn("strict_ala_continuity_metrics(", self.core)
        self.assertIn("strict_ala_momentum_residual_audit(", self.core)
        self.assertIn(
            "converged=(continuity_converged && momentum_converged)",
            self.core,
        )
        self.assertIn("ALA_COUPLED_FEASIBILITY_SUMMARY", self.core)
        self.assertIn("parallel_process_termination()", self.core)

    def test_restart_and_final_momentum_audits_are_current(self) -> None:
        self.assertIn('"coupled_restart"', self.core)
        self.assertIn('"coupled_final"', self.core)
        self.assertIn("raw_momentum_relative=%e breakdown=%d", self.core)

    def test_preconditioned_block_action_is_audited(self) -> None:
        audit = _between(
            self.stokes,
            "static void strict_ala_coupled_preconditioner_audit(",
            "static float solve_ala_coupled_fgmres_core(",
        )
        self.assertIn("r_components=(velocity:%e,pressure:%e)", audit)
        self.assertIn("z_components=(velocity:%e,pressure:%e)", audit)
        self.assertIn("Az_components=(velocity:%e,pressure:%e)", audit)
        self.assertIn("defect_to_r=(velocity:%e,pressure:%e)", audit)
        self.assertIn("cosine=(velocity:%e,pressure:%e,block:%e)", audit)
        self.assertIn(
            "optimal_scale=(velocity:%e,pressure:%e,block:%e)", audit
        )
        self.assertIn(
            "projected_defect=(velocity:%e,pressure:%e,block:%e)", audit
        )
        self.assertIn("pressure_optimal_scale=rap/max(ap*ap", audit)
        self.assertIn("strict_ala_pressure_mode_audit(", audit)
        self.assertIn(
            "ALA PRESSURE MODE AUDIT iteration=%d offset=%d", self.stokes
        )
        self.assertIn(
            "value_r=residual->pressure[m][e]-mean_r", self.stokes
        )
        self.assertIn("strict_ala_coupled_preconditioner_audit(", self.core)

    def test_active_cfg_is_clean_current_rheology_diagnostic(self) -> None:
        cfg = (PROJECT_ROOT / "runs/cmbhf_ALA_strict.cfg").read_text()
        self.assertIn("ala_outer_solver                = coupled_fgmres", cfg)
        self.assertIn(
            "ala_coupled_initial_velocity_relative_tolerance = 0.0", cfg
        )
        self.assertIn(
            "ala_coupled_inner_relative_tolerance = 1e-2", cfg
        )
        self.assertIn("ala_coupled_inner_max_cycles       = 200", cfg)
        self.assertIn(
            "ala_coupled_inner_progress_interval = 20", cfg
        )
        self.assertIn("ala_coupled_defect_corrections       = 0", cfg)
        self.assertIn("ala_pressure_multigrid                  = off", cfg)
        self.assertIn("ala_pressure_multigrid_galerkin         = off", cfg)
        self.assertIn("ala_augmented_lagrangian_gamma = 10.0", cfg)
        self.assertIn("ala_coupled_debug_stop_iteration            = 50", cfg)
        self.assertIn("ala_element_vanka_damping      = 0.2", cfg)
        self.assertIn(
            "ala_coupled_shallow_vanka_layers             = 0", cfg
        )
        self.assertIn(
            "ala_coupled_shallow_vanka_sweeps             = 0", cfg
        )
        self.assertIn("ala_shallow_patch_preconditioner   = on", cfg)
        self.assertIn("ala_shallow_patch_weight           = 1.0", cfg)
        self.assertIn(
            "ala_shallow_patch_horizontal_elements = 6", cfg
        )
        self.assertIn(
            "ala_shallow_patch_horizontal_stride   = 3", cfg
        )
        self.assertIn("ala_shallow_patch_mpi_overlap      = 3", cfg)
        self.assertIn(
            "ala_shallow_patch_velocity_solver  = element_vanka", cfg
        )
        self.assertIn(
            "ala_coupled_factor2_coarse_correction = off", cfg
        )

    def test_coupled_path_prebalances_every_momentum_residual(self) -> None:
        dispatch = _between(
            self.stokes,
            'if(strcmp(E->control.ala_outer_solver,"fgmres")==0 ||',
            "/* FF contains the current -C^T*P forcing.",
        )
        prebalance = dispatch.index("ALA COUPLED MOMENTUM PREBALANCE")
        velocity_solve = dispatch.index("initial_vel_residual(", prebalance)
        coupled_solve = dispatch.index("solve_ala_coupled_fgmres_core(")
        legacy = dispatch[coupled_solve:]
        self.assertIn(
            "ala_coupled_initial_velocity_relative_tolerance", dispatch
        )
        self.assertIn("solution_cycle=%d relative_tolerance=%e", dispatch)
        self.assertIn("E->monitor.solution_cycles", dispatch)
        self.assertIn("initial_velocity_solve=skipped", dispatch)
        self.assertLess(prebalance, velocity_solve)
        self.assertLess(velocity_solve, coupled_solve)
        self.assertIn("initial_vel_residual(E,V,P,F,imp);", legacy)
        self.assertIn("solve_ala_fgmres_core(", legacy)


if __name__ == "__main__":
    unittest.main(verbosity=2)
