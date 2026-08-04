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
        self.assertIn(
            "int ala_coupled_first_preconditioner_audit_only;", self.defs
        )
        self.assertIn(
            'input_boolean("ala_coupled_first_preconditioner_audit_only",',
            self.instructions,
        )
        self.assertIn(
            "E->control.ala_coupled_first_preconditioner_audit_only = 0;",
            self.instructions,
        )
        self.assertIn(
            '"ala_coupled_first_preconditioner_audit_only", default=False',
            self.pyre,
        )
        self.assertIn(
            'getIntProperty(properties, '
            '"ala_coupled_first_preconditioner_audit_only",',
            self.properties,
        )
        self.assertIn("int ala_coupled_debug_stop_iteration;", self.defs)
        self.assertIn(
            'input_int("ala_coupled_debug_stop_iteration",',
            self.instructions,
        )
        self.assertIn(
            '"ala_coupled_debug_stop_iteration", default=0', self.pyre
        )
        self.assertIn(
            'getIntProperty(properties, "ala_coupled_debug_stop_iteration",',
            self.properties,
        )
        self.assertIn("int ala_coupled_element_vanka;", self.defs)
        self.assertIn(
            'input_boolean("ala_coupled_element_vanka",', self.instructions
        )
        self.assertIn(
            '"ala_coupled_element_vanka", default=False', self.pyre
        )
        self.assertIn(
            'getIntProperty(properties, "ala_coupled_element_vanka",',
            self.properties,
        )
        self.assertIn("int ala_coupled_multilevel_vcycle;", self.defs)
        self.assertIn(
            'input_boolean("ala_coupled_multilevel_vcycle",',
            self.instructions,
        )
        self.assertIn(
            '"ala_coupled_multilevel_vcycle", default=False', self.pyre
        )
        self.assertIn(
            'getIntProperty(properties, "ala_coupled_multilevel_vcycle",',
            self.properties,
        )
        self.assertIn(
            "int ala_coupled_multilevel_coarse_sweeps;", self.defs
        )
        self.assertIn(
            'input_int("ala_coupled_multilevel_coarse_sweeps",',
            self.instructions,
        )
        self.assertIn(
            '"ala_coupled_multilevel_coarse_sweeps", default=2', self.pyre
        )
        self.assertIn(
            'getIntProperty(properties, '
            '"ala_coupled_multilevel_coarse_sweeps",',
            self.properties,
        )
        self.assertIn(
            "double ala_coupled_multilevel_coarse_weight;", self.defs
        )
        self.assertIn(
            'input_double("ala_coupled_multilevel_coarse_weight",',
            self.instructions,
        )
        self.assertIn(
            '"ala_coupled_multilevel_coarse_weight", default=1.0', self.pyre
        )
        self.assertIn(
            'getDoubleProperty(properties, '
            '"ala_coupled_multilevel_coarse_weight",',
            self.properties,
        )
        self.assertIn("int ala_coupled_shallow_vanka_layers;", self.defs)
        self.assertIn("int ala_coupled_shallow_vanka_sweeps;", self.defs)
        self.assertIn(
            'input_int("ala_coupled_shallow_vanka_layers",',
            self.instructions,
        )
        self.assertIn(
            'input_int("ala_coupled_shallow_vanka_sweeps",',
            self.instructions,
        )
        self.assertIn(
            '"ala_coupled_shallow_vanka_layers", default=0', self.pyre
        )
        self.assertIn(
            '"ala_coupled_shallow_vanka_sweeps", default=0', self.pyre
        )
        self.assertIn(
            'getIntProperty(properties, "ala_coupled_shallow_vanka_layers",',
            self.properties,
        )
        self.assertIn(
            'getIntProperty(properties, "ala_coupled_shallow_vanka_sweeps",',
            self.properties,
        )

    def test_first_preconditioner_audit_stops_after_action_audit(self) -> None:
        core = _function_body(
            self.stokes, "static float solve_ala_coupled_fgmres_core("
        )
        application = core.index("apply_ala_coupled_block_preconditioner(")
        action = core.index("apply_ala_coupled_operator(", application)
        audit = core.index("strict_ala_coupled_preconditioner_audit(", action)
        gate = core.index(
            "ala_coupled_first_preconditioner_audit_only", audit
        )
        termination = core.index("exit(EXIT_SUCCESS);", gate)
        orthogonalization = core.index("for(i=0;i<=j;i++)", termination)
        self.assertLess(application, action)
        self.assertLess(action, audit)
        self.assertLess(audit, gate)
        self.assertLess(gate, termination)
        self.assertLess(termination, orthogonalization)

    def test_short_run_stop_follows_iteration_and_momentum_audit(self) -> None:
        core = _function_body(
            self.stokes, "static float solve_ala_coupled_fgmres_core("
        )
        iteration_line = core.index("ALA COUPLED FGMRES iteration=%d")
        gate = core.index("ala_coupled_debug_stop_iteration>0", iteration_line)
        momentum = core.index('"coupled_debug_stop"', gate)
        marker = core.index("ALA COUPLED DEBUG ITERATION STOP COMPLETE", momentum)
        termination = core.index("exit(EXIT_SUCCESS);", marker)
        acceptance = core.index("ALA_COUPLED_FEASIBILITY_SUMMARY", termination)
        self.assertLess(iteration_line, gate)
        self.assertLess(gate, momentum)
        self.assertLess(momentum, marker)
        self.assertLess(marker, termination)
        self.assertLess(termination, acceptance)

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
            "G_exchange_adjoint_defect=%e",
            "G_stripped_adjoint_defect=%e",
            "G_multiplicity_adjoint_defect=%e",
            "K_bilinear=(%e,%e)",
            "G_bilinear=(%e,%e)",
            "G_element_bilinear=(%e,%e)",
            "G_exchange_bilinear=(%e,%e)",
            "G_stripped_bilinear=(%e,%e)",
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
        self.assertLess(
            audit.index("(E->solver.exchange_id_d)(E,Ax->velocity,level);"),
            audit.index("strip_bcs_from_residual(E,Ax->velocity,level);"),
        )
        self.assertIn("&g_exchange_right", audit)
        self.assertIn("&g_stripped_right", audit)
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
        self.assertIn("const unsigned int vbc_flag[4]={0,VBX,VBY,VBZ};", probe)
        self.assertIn("E->NODE[level][m][node] & vbc_flag[d]", probe)
        self.assertIn("value=0.0;", probe)
        self.assertNotIn("E->NODE[level][m][node] & SKIP", probe)
        self.assertNotIn("exchange_id_d", probe)
        self.assertIn(
            "velocity_probe=direct_conforming_radial_bc_projected",
            self.element,
        )

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
        self.assertIn(
            "strip_bcs_from_residual(E,fine,coarse_level+1);", prolong
        )
        self.assertNotIn("project_vector", restrict)
        self.assertIn("from_rtf_to_xyz(E,fine_level,fine,E->temp);", restrict)
        self.assertIn("E->NODE[fine_level][m][node] & SKIP", restrict)
        self.assertIn(
            "(E->solver.exchange_id_d)(E,E->temp1,coarse_level);", restrict
        )
        self.assertIn(
            "from_xyz_to_rtf(E,coarse_level,E->temp1,coarse);", restrict
        )
        self.assertIn(
            "strip_bcs_from_residual(E,coarse,coarse_level);", restrict
        )
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

    def test_stage9f1_element_patch_solves_both_fields(self) -> None:
        smoother = _function_body(
            self.stokes,
            "static void apply_ala_coupled_element_vanka_region(",
        )
        for token in (
            "E->ALA_vanka_chol[lev][m]",
            "E->elt_del[lev][m][e].g[i][0]",
            "E->elt_c[lev][m][e].c[i][0]",
            "pressure_solution /= schur;",
            "velocity_base[i]",
            "-velocity_pressure[i]*pressure_solution",
            "correction->pressure[m][e]",
            "ALA COUPLED ELEMENT VANKA APPLICATION",
            "local_schur_range=[%e,%e]",
            "fallback_count=0",
        ):
            self.assertIn(token, smoother)
        block = _function_body(
            self.stokes, "static void apply_ala_coupled_block_preconditioner("
        )
        self.assertIn("if(E->control.ala_coupled_element_vanka)", block)
        self.assertIn("apply_ala_coupled_element_vanka_once(", block)

    def test_vanka_diagnostics_do_not_reduce_on_every_sweep(self) -> None:
        smoother = _function_body(
            self.stokes,
            "static void apply_ala_coupled_element_vanka_region(",
        )
        self.assertIn("reported_cycle[MAX_LEVELS]", smoother)
        self.assertIn("reported_valid[MAX_LEVELS]", smoother)
        self.assertIn("report_diagnostics", smoother)
        self.assertIn("if(report_diagnostics) {", smoother)
        self.assertIn("diagnostic_scope=first_application_per_cycle", smoother)
        self.assertEqual(smoother.count("MPI_Allreduce("), 3)

    def test_coupled_vcycle_skips_unused_legacy_pressure_cache(self) -> None:
        solve = self.stokes[
            self.stokes.rindex("static float solve_Ahat_p_fhat_ALA_PCG("):
        ]
        self.assertIn("coupled_self_contained_preconditioner", solve)
        self.assertIn("ala_coupled_multilevel_vcycle ||", solve)
        self.assertIn("ala_coupled_element_vanka", solve)
        self.assertIn(
            "if(!coupled_self_contained_preconditioner &&\n"
            "       E->control.ala_shallow_patch_preconditioner)",
            solve,
        )
        self.assertIn("self_contained_skip_legacy_pressure_cache", solve)
        self.assertIn('? "coupled_multilevel_vcycle"', solve)

    def test_expensive_action_audit_is_targeted_not_periodic(self) -> None:
        core = _function_body(
            self.stokes, "static float solve_ala_coupled_fgmres_core("
        )
        audit = core[core.index("strict_ala_coupled_preconditioner_audit(") - 300:]
        self.assertIn("count==0 || j==restart-1", audit)
        self.assertIn(
            "count+1==E->control.ala_coupled_debug_stop_iteration", audit
        )
        self.assertNotIn("((count+1)%5)==0", audit)

    def test_coupled_iterations_report_pressure_residual_by_level(self) -> None:
        core = _function_body(
            self.stokes, "static float solve_ala_coupled_fgmres_core("
        )
        self.assertIn(
            "strict_ala_coarse_residual_diagnostics(\n"
            "        E,r->pressure,lev,0);",
            core,
        )
        self.assertIn(
            "strict_ala_coarse_residual_diagnostics(\n"
            "                E,explicit_r->pressure,lev,count);",
            core,
        )
        self.assertIn(
            "strict_ala_depth_diagnostics(\n"
            "        E,w->pressure,operator_work->pressure,lev,0);",
            core,
        )
        self.assertIn(
            "strict_ala_depth_diagnostics(\n"
            "                E,w->pressure,operator_work->pressure,lev,count);",
            core,
        )

    def test_stage9f2_vcycle_uses_both_exact_transfer_pairs(self) -> None:
        vcycle = _function_body(
            self.stokes,
            "static void apply_ala_coupled_multilevel_vcycle(",
        )
        for token in (
            "apply_ala_coupled_element_vanka_once(",
            "assemble_ala_coupled_block_defect(",
            "ala_coupled_restrict_velocity(",
            "ala_coupled_restrict_pressure_p0_transpose(",
            "apply_ala_coupled_multilevel_vcycle(",
            "ala_coupled_prolong_velocity(",
            "ala_coupled_prolong_pressure_p0(",
            "E->control.ala_coupled_multilevel_coarse_weight",
            "E->control.ala_coupled_multilevel_coarse_sweeps",
            "E->control.ala_coupled_shallow_vanka_layers",
            "E->control.ala_coupled_shallow_vanka_sweeps",
            "apply_ala_coupled_element_vanka_region(",
        ):
            self.assertIn(token, vcycle)
        block = _function_body(
            self.stokes, "static void apply_ala_coupled_block_preconditioner("
        )
        self.assertIn("if(E->control.ala_coupled_multilevel_vcycle)", block)
        self.assertIn("ALA COUPLED MULTILEVEL VCYCLE APPLICATION", block)
        self.assertIn("coarse_weight=%e", block)
        self.assertIn("shallow_layers=%d shallow_sweeps=%d", block)

    def test_stage10c_shallow_smoother_uses_selected_patch_partition(self) -> None:
        region = _function_body(
            self.stokes,
            "static void apply_ala_coupled_element_vanka_region(",
        )
        for token in (
            "radial_element=(e-1)%elz",
            "global_radial_element=E->lmesh.EZS[lev]+radial_element",
            "E->mesh.ELZ[lev]-selected_layers",
            "correction->velocity[m][eq] += 1.0",
            "exchange_id_d)(E,correction->velocity,lev)",
            "1.0/correction->velocity[m][eq]",
            'shallow ? "shallow_finest" : "full"',
        ):
            self.assertIn(token, region)


if __name__ == "__main__":
    unittest.main(verbosity=2)
