#!/usr/bin/env python3
"""Static contracts for the strict-ALA monolithic-operator foundation."""

from __future__ import annotations

import unittest
from pathlib import Path


GLOBAL_ROOT = Path(__file__).resolve().parents[2]
LIB_ROOT = GLOBAL_ROOT / "lib"


def _function_body(source: str, signature: str) -> str:
    start = source.index(signature)
    opening = source.index("{", start)
    depth = 0
    for index in range(opening, len(source)):
        if source[index] == "{":
            depth += 1
        elif source[index] == "}":
            depth -= 1
            if depth == 0:
                return source[start : index + 1]
    raise AssertionError(f"Unterminated function {signature}")


class CoupledOperatorArchitectureTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.element = (LIB_ROOT / "Element_calculations.c").read_text()
        cls.header = (LIB_ROOT / "ala_coupled_operator.h").read_text()

    def test_public_block_contract_contains_complete_strict_operator(self) -> None:
        self.assertIn("K_gamma  (D+C)^T", self.header)
        self.assertIn("D+C        0", self.header)
        self.assertIn("apply_ala_coupled_operator", self.header)

    def test_block_action_uses_full_transpose_and_continuity(self) -> None:
        body = _function_body(
            self.element, "void apply_ala_coupled_operator("
        )
        self.assertIn(
            "assemble_del2_u(E,velocity,momentum,level,1);", body
        )
        self.assertIn(
            "assemble_grad_rho_p(E,pressure,velocity_work,level);", body
        )
        self.assertIn(
            "assemble_div_rho_u(E,velocity,continuity,level);", body
        )
        self.assertNotIn("assemble_forces", body)
        self.assertLess(
            body.index("assemble_del2_u"), body.index("assemble_grad_rho_p")
        )
        self.assertLess(
            body.index("assemble_grad_rho_p"),
            body.index("assemble_div_rho_u"),
        )

    def test_block_action_rejects_destructive_aliases(self) -> None:
        body = _function_body(
            self.element, "void apply_ala_coupled_operator("
        )
        for alias in (
            "velocity==momentum",
            "velocity==velocity_work",
            "momentum==velocity_work",
            "pressure==continuity",
        ):
            self.assertIn(alias, body)

    def test_pressure_independent_rhs_removes_legacy_ctp_force(self) -> None:
        body = _function_body(
            self.element,
            "void assemble_ala_pressure_independent_force(",
        )
        full_gradient = "assemble_grad_rho_p(E,E->P,velocity_work,level);"
        ordinary_gradient = "assemble_grad_p(E,E->P,velocity_work,level);"
        self.assertIn("assemble_forces(E,0);", body)
        self.assertIn(full_gradient, body)
        self.assertIn(ordinary_gradient, body)
        self.assertIn("force[m][i]=E->F[m][i]+velocity_work[m][i];", body)
        self.assertIn("force[m][i]-=velocity_work[m][i];", body)
        self.assertLess(body.index(full_gradient), body.index(ordinary_gradient))

    def test_existing_production_solver_does_not_select_new_path(self) -> None:
        stokes = (LIB_ROOT / "Stokes_flow_Incomp.c").read_text()
        drive = (LIB_ROOT / "Drive_solvers.c").read_text()
        self.assertNotIn("apply_ala_coupled_operator", stokes)
        self.assertNotIn("apply_ala_coupled_operator", drive)
        self.assertNotIn("assemble_ala_pressure_independent_force", stokes)
        self.assertNotIn("assemble_ala_pressure_independent_force", drive)


if __name__ == "__main__":
    unittest.main(verbosity=2)
