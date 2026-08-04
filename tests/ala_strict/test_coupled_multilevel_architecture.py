#!/usr/bin/env python3
"""Static contracts for Stage 9f.0 coupled multilevel preparation."""

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


class CoupledMultilevelArchitectureTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.element = (LIB_ROOT / "Element_calculations.c").read_text()
        cls.header = (LIB_ROOT / "ala_coupled_operator.h").read_text()
        cls.stokes = (LIB_ROOT / "Stokes_flow_Incomp.c").read_text()
        cls.defs = (LIB_ROOT / "global_defs.h").read_text()
        cls.instructions = (LIB_ROOT / "Instructions.c").read_text()
        cls.pyre = (
            GLOBAL_ROOT / "CitcomS/Components/Stokes_solver/Incompressible.py"
        ).read_text()
        cls.properties = (GLOBAL_ROOT / "module/setProperties.c").read_text()

    def test_every_coupled_solve_enters_one_time_hierarchy_audit(self) -> None:
        core = _function_body(
            self.stokes, "static float solve_ala_coupled_fgmres_core("
        )
        self.assertIn("audit_ala_coupled_multilevel_contracts(E);", core)
        audit = _function_body(
            self.element, "void audit_ala_coupled_multilevel_contracts("
        )
        self.assertIn("for(level=E->mesh.levmin;level<=E->mesh.levmax;level++)", audit)
        self.assertIn("static int completed=0;", audit)

    def test_audit_only_terminates_after_five_levels_before_fgmres(self) -> None:
        core = _function_body(
            self.stokes, "static float solve_ala_coupled_fgmres_core("
        )
        audit_call = core.index("audit_ala_coupled_multilevel_contracts(E);")
        audit_gate = core.index("ala_coupled_multilevel_audit_only")
        barrier = core.index("MPI_Barrier(E->parallel.world);")
        finalization = core.index("MPI_Finalize();")
        termination = core.index("exit(EXIT_SUCCESS);")
        basis_allocation = core.index("rhs=ala_block_vector_create(E,lev);")
        self.assertLess(audit_call, audit_gate)
        self.assertLess(audit_gate, barrier)
        self.assertLess(barrier, finalization)
        self.assertLess(finalization, termination)
        self.assertLess(termination, basis_allocation)
        self.assertIn("terminating before FGMRES iteration 1", core)

    def test_audit_only_configuration_is_default_off_and_validated(self) -> None:
        self.assertIn("int ala_coupled_multilevel_audit_only;", self.defs)
        self.assertIn(
            'input_boolean("ala_coupled_multilevel_audit_only",',
            self.instructions,
        )
        self.assertIn(
            'E->control.ala_coupled_multilevel_audit_only = 0;',
            self.instructions,
        )
        self.assertIn(
            '"ala_coupled_multilevel_audit_only", default=False', self.pyre
        )
        self.assertIn(
            'getIntProperty(properties, "ala_coupled_multilevel_audit_only",',
            self.properties,
        )
        for source in (self.instructions, self.properties):
            self.assertIn(
                "ala_coupled_multilevel_audit_only requires ", source
            )
            self.assertIn('ala_outer_solver=coupled_fgmres', source)

    def test_level_audit_checks_complete_operator_contracts(self) -> None:
        audit = _function_body(
            self.element, "void audit_ala_coupled_multilevel_contracts("
        )
        for token in (
            "assemble_del2_u(E,x->velocity,Ax->velocity,level,1);",
            "assemble_grad_rho_p(E,x->pressure,Ax->velocity,level);",
            "assemble_div_rho_u(E,y->velocity,Ay->pressure,level);",
            "K_symmetry_defect=%e",
            "G_adjoint_defect=%e",
            "G_element_adjoint_defect=%e",
            "G_multiplicity_adjoint_defect=%e",
            "K_bilinear=(%e,%e)",
            "G_bilinear=(%e,%e)",
            "G_element_bilinear=(%e,%e)",
            "G_multiplicity_bilinear=(%e,%e)",
            "pressure_mass_range=[%e,%e]",
            "duplicate_velocity_dofs=%d",
        ):
            self.assertIn(token, audit)
        scaling = _function_body(
            self.element, "static double ala_norm_scaled_adjoint_defect("
        )
        self.assertIn("left_vector_norm*left_action_norm", scaling)
        self.assertIn("right_vector_norm*right_action_norm", scaling)
        self.assertNotIn("fabs(left)+fabs(right)", scaling)
        self.assertIn("output=(output_index==0) ? E->fp : stderr;", audit)
        self.assertIn("assemble_grad_rho_p_local_terms(", audit)
        self.assertIn("MPI_Allreduce(&local_g_element_right", audit)
        self.assertIn(
            "(E->solver.exchange_id_d)(E,multiplicity->velocity,level);",
            audit,
        )
        self.assertIn("/multiplicity->velocity[m][e]", audit)
        self.assertIn("MPI_Allreduce(&local_g_multiplicity_right", audit)

    def test_production_gradient_wraps_same_local_transpose_terms(self) -> None:
        gradient = _function_body(
            self.element, "void assemble_grad_rho_p("
        )
        self.assertIn(
            "assemble_grad_rho_p_local_terms(E,P,gradP,lev);", gradient
        )
        self.assertIn("(E->solver.exchange_id_d)(E, gradP, lev);", gradient)
        self.assertIn("strip_bcs_from_residual(E,gradP,lev);", gradient)

    def test_probe_is_directly_conforming_not_additively_assembled(self) -> None:
        probe = _function_body(
            self.element, "static void ala_fill_level_probe("
        )
        self.assertIn("radius=E->SX[level][m][3][node];", probe)
        self.assertIn("probe->velocity[m][eq]=value;", probe)
        self.assertNotIn("E->NODE[level][m][node] & SKIP", probe)
        self.assertNotIn("exchange_id_d", probe)
        self.assertIn("velocity_probe=direct_conforming_radial", self.element)

    def test_coarse_beta_is_restricted_from_authoritative_fine_field(self) -> None:
        audit = _function_body(
            self.element, "void audit_ala_coupled_multilevel_contracts("
        )
        beta = _function_body(self.element, "static double ala_restricted_beta(")
        self.assertIn("E->refstate.ala_beta[fine_nz]*dr", beta)
        self.assertIn("return beta/dr_total;", beta)
        self.assertIn("nested_beta += child_beta*child_width;", audit)
        self.assertIn("beta_nested_defect=%e", audit)
        self.assertIn("beta_range=[%e,%e]", audit)

    def test_pressure_transfer_is_constant_p_with_exact_transpose(self) -> None:
        prolong = _function_body(
            self.element, "void ala_coupled_prolong_pressure_p0("
        )
        restrict = _function_body(
            self.element, "void ala_coupled_restrict_pressure_p0_transpose("
        )
        self.assertIn("fine[m][e]=coarse[m][ce];", prolong)
        self.assertIn("coarse[m][ce] += fine[m][e];", restrict)
        self.assertNotIn("ALA_velocity_BI", prolong + restrict)

    def test_velocity_transfer_is_exposed_and_numerically_audited(self) -> None:
        prolong = _function_body(
            self.element, "void ala_coupled_prolong_velocity("
        )
        restrict = _function_body(
            self.element, "void ala_coupled_restrict_velocity("
        )
        self.assertIn("interp_vector(E,coarse_level,coarse,fine);", prolong)
        self.assertIn("project_vector(E,fine_level,fine,coarse,1);", restrict)
        self.assertIn("velocity_transfer_adjoint_defect=%e", self.element)
        self.assertIn("pressure_transfer_adjoint_defect=%e", self.element)

    def test_stage9f1_interface_is_mixed_and_observable(self) -> None:
        for token in (
            "struct ala_coupled_patch_spec",
            "const int *elements",
            "double pressure_regularization",
            "struct ala_coupled_patch_stats",
            "double pivot_ratio",
            "int fallback_count",
            "size_t workspace_bytes",
            "typedef int (*ala_coupled_patch_solver)",
            "const struct ala_block_vector *residual",
            "struct ala_block_vector *correction",
        ):
            self.assertIn(token, self.header)


if __name__ == "__main__":
    unittest.main(verbosity=2)
