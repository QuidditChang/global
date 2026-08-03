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

    def test_right_preconditioner_has_correct_saddle_point_sign(self) -> None:
        preconditioner = _between(
            self.stokes,
            "static void apply_ala_coupled_block_preconditioner(",
            "static float solve_ala_coupled_fgmres_core(",
        )
        self.assertEqual(preconditioner.count("solve_del2_u("), 2)
        self.assertIn("apply_ala_pressure_preconditioner(", preconditioner)
        self.assertIn(
            "correction->pressure[m][e]=-correction->pressure[m][e]",
            preconditioner,
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

    def test_active_cfg_is_isolated_stage9d_ab(self) -> None:
        cfg = (PROJECT_ROOT / "runs/cmbhf_ALA_strict.cfg").read_text()
        self.assertIn("ala_outer_solver                = coupled_fgmres", cfg)
        self.assertIn("ala_pressure_multigrid                  = off", cfg)


if __name__ == "__main__":
    unittest.main(verbosity=2)
