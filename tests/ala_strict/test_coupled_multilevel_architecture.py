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
        cls.construct = (LIB_ROOT / "Construct_arrays.c").read_text()
        cls.defs = (LIB_ROOT / "global_defs.h").read_text()
        cls.instructions = (LIB_ROOT / "Instructions.c").read_text()
        cls.output = (LIB_ROOT / "Output_h5.c").read_text()
        cls.pyre = (
            GLOBAL_ROOT / "CitcomS/Components/Stokes_solver/Incompressible.py"
        ).read_text()
        cls.properties = (GLOBAL_ROOT / "module/setProperties.c").read_text()
        cls.strict_cfg = (
            GLOBAL_ROOT.parents[1] / "runs/cmbhf_ALA_strict.cfg"
        ).read_text()

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
        for token in (
            "ala_viscosity_spectrum_diagnostics",
            "ala_viscosity_spectrum_interval",
        ):
            self.assertIn(token, self.defs)
            self.assertIn(token, self.instructions)
            self.assertIn(token, self.pyre)
            self.assertIn(token, self.properties)
            self.assertIn(token, self.output)
        self.assertIn("int ala_coupled_shallow_vanka_layers;", self.defs)
        self.assertIn(
            "int ala_coupled_shallow_vanka_core_layers;", self.defs
        )
        self.assertIn(
            "int ala_coupled_shallow_vanka_band_sweeps;", self.defs
        )
        self.assertIn("int ala_coupled_shallow_vanka_sweeps;", self.defs)
        self.assertIn(
            "int ala_coupled_shallow_vanka_warm_sweeps;", self.defs
        )
        self.assertIn(
            'input_int("ala_coupled_shallow_vanka_layers",',
            self.instructions,
        )
        self.assertIn(
            'input_int("ala_coupled_shallow_vanka_core_layers",',
            self.instructions,
        )
        self.assertIn(
            'input_int("ala_coupled_shallow_vanka_band_sweeps",',
            self.instructions,
        )
        self.assertIn(
            'input_int("ala_coupled_shallow_vanka_sweeps",',
            self.instructions,
        )
        self.assertIn(
            'input_int("ala_coupled_shallow_vanka_warm_sweeps",',
            self.instructions,
        )
        self.assertIn(
            '"ala_coupled_shallow_vanka_layers", default=0', self.pyre
        )
        self.assertIn(
            '"ala_coupled_shallow_vanka_core_layers", default=-1',
            self.pyre,
        )
        self.assertIn(
            '"ala_coupled_shallow_vanka_band_sweeps", default=0',
            self.pyre,
        )
        self.assertIn(
            '"ala_coupled_shallow_vanka_sweeps", default=0', self.pyre
        )
        self.assertIn(
            '"ala_coupled_shallow_vanka_warm_sweeps", default=-1',
            self.pyre,
        )
        self.assertIn(
            'getIntProperty(properties, "ala_coupled_shallow_vanka_layers",',
            self.properties,
        )
        self.assertIn(
            'getIntProperty(properties, '
            '"ala_coupled_shallow_vanka_core_layers",',
            self.properties,
        )
        self.assertIn(
            'getIntProperty(properties, '
            '"ala_coupled_shallow_vanka_band_sweeps",',
            self.properties,
        )
        self.assertIn(
            'getIntProperty(properties, "ala_coupled_shallow_vanka_sweeps",',
            self.properties,
        )
        self.assertIn(
            'getIntProperty(properties, '
            '"ala_coupled_shallow_vanka_warm_sweeps",',
            self.properties,
        )
        for source in (self.instructions, self.properties):
            self.assertIn(
                "ala_coupled_shallow_vanka_warm_sweeps must be -1 or in ",
                source,
            )
            self.assertIn(
                "ala_coupled_shallow_vanka_core_layers must be -1 or in ",
                source,
            )
            self.assertIn(
                "ala_coupled_shallow_vanka_band_sweeps must be in [0,8]",
                source,
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

    def test_galerkin_rap_audit_uses_complete_coupled_action(self) -> None:
        element = self.element
        rap = _function_body(
            element, "static void ala_apply_coupled_galerkin_rap("
        )
        for token in (
            "ala_coupled_prolong_velocity(",
            "ala_coupled_prolong_pressure_p0(",
            "apply_ala_coupled_operator(",
            "ala_coupled_restrict_velocity(",
            "ala_coupled_restrict_pressure_p0_transpose(",
        ):
            self.assertIn(token, rap)
        audit = _function_body(
            element, "void audit_ala_coupled_multilevel_contracts("
        )
        self.assertIn("ALA COUPLED GALERKIN RAP AUDIT", audit)
        self.assertIn("rediscretized_action_difference", audit)
        self.assertIn("rap_symmetry_defect", audit)
        self.assertIn("action=observe_only", audit)

    def test_galerkin_pressure_schur_inherits_child_sum(self) -> None:
        builder = _function_body(
            self.construct, "void build_ala_element_vanka_factors("
        )
        for token in (
            "ala_element_vanka_galerkin_schur",
            "for(level=E->mesh.gridmax;level>=E->mesh.gridmin+1;level--)",
            "for(fine_y=2*coarse_y-1",
            "for(fine_x=2*coarse_x-1",
            "for(fine_z=2*coarse_z-1",
            "schur +=",
            "E->ALA_vanka_schur[coarse_level][m][coarse_e]",
            "ALA GALERKIN PRESSURE SCHUR SCALE enabled",
            "operator=PpT_Sfine_Pp",
        ):
            self.assertIn(token, builder)

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

    def test_action_audit_is_targeted_not_periodic(self) -> None:
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
            "coarse_weight=(lev==E->mesh.levmax)",
            ": 1.0",
            "E->control.ala_coupled_multilevel_coarse_sweeps",
            "E->control.ala_coupled_shallow_vanka_layers",
            "E->control.ala_coupled_shallow_vanka_core_layers",
            "E->control.ala_coupled_shallow_vanka_band_sweeps",
            "E->control.ala_coupled_shallow_vanka_sweeps",
            "E->control.ala_coupled_shallow_vanka_warm_sweeps",
            "E->monitor.solution_cycles>0",
            "apply_ala_coupled_element_vanka_region(",
        ):
            self.assertIn(token, vcycle)
        block = _function_body(
            self.stokes, "static void apply_ala_coupled_block_preconditioner("
        )
        self.assertIn("if(E->control.ala_coupled_multilevel_vcycle)", block)
        self.assertIn("ALA COUPLED MULTILEVEL VCYCLE APPLICATION", block)
        self.assertIn("coarse_weight=%e", block)
        self.assertIn("nested_coarse_weight=1.000000e+00", block)
        self.assertIn(
            "shallow_layers=%d shallow_core_layers=%d", block
        )
        self.assertIn(
            "shallow_band_window=sine shallow_band_sweeps=%d", block
        )
        self.assertIn("shallow_cold_sweeps=%d shallow_warm_sweeps=%d", block)
        self.assertIn('solution_start=%s', block)

    def test_vanka_application_rebuilds_overlap_consistent_local_schur(self) -> None:
        builder = _function_body(
            self.construct, "void build_ala_element_vanka_factors("
        )
        self.assertIn("E->ALA_vanka_schur[level][m][e]=schur;", builder)
        self.assertIn("E->EVI[level][m][(e-1)*vpts+v]", builder)
        self.assertIn("diag=1.0/E->BI[level][m][eq];", builder)
        self.assertIn("ALA VISCOSITY SPECTRUM level=%d", builder)
        self.assertIn("eta_mean_times_schur_mean=%e", builder)
        region = _function_body(
            self.stokes,
            "static void apply_ala_coupled_element_vanka_region(",
        )
        self.assertIn("1.0/sqrt(correction->velocity[m][eq])", region)
        self.assertIn(
            "sqrt(E->ALA_vanka_overlap_BI[lev][m][eq])", region
        )
        self.assertIn("velocity_rhs[i]=local_weight[i]", region)
        self.assertIn("gradient[i]=local_weight[i]", region)
        self.assertIn("schur += gradient[i]*velocity_pressure[i];", region)
        self.assertIn("region_weight*local_weight[i]", region)
        self.assertIn("overlap_contract=sqrt_partition", region)

    def test_vanka_velocity_factor_can_remove_external_diagonal(self) -> None:
        builder = _function_body(
            self.construct, "void build_ala_element_vanka_factors("
        )
        self.assertIn(
            "E->control.ala_element_vanka_external_diagonal_weight",
            builder,
        )
        self.assertIn("level==E->mesh.gridmax", builder)
        self.assertIn(": 1.0", builder)
        self.assertIn("*max(external,0.0)", builder)
        self.assertIn("finest_external_diagonal_weight=%e", builder)
        self.assertIn(
            "coarse_external_diagonal_weight=1.000000e+00", builder
        )

    def test_pressure_aggregate_uses_current_rheology_vanka_metric(self) -> None:
        self.assertIn("#define ALA_PATCH_RADIAL_ELEMENTS 2", self.stokes)
        self.assertIn(
            "static double ala_element_vanka_schur_entry(", self.stokes
        )
        self.assertIn("ala_solve_cached_element_k(chol,rhs2,sol2)", self.stokes)
        self.assertIn(
            '"element_vanka")==0', self.stokes
        )
        self.assertIn(
            "=ala_element_vanka_schur_entry(E,e1,e2,lev,m);", self.stokes
        )
        self.assertIn(
            "ala_shallow_patch_horizontal_elements = 6", self.strict_cfg
        )
        self.assertIn(
            "ala_shallow_patch_horizontal_stride   = 3", self.strict_cfg
        )
        self.assertIn(
            "ala_shallow_patch_mpi_overlap      = 3", self.strict_cfg
        )

    def test_pressure_aggregate_interface_carries_vanka_factors(self) -> None:
        for token in (
            "ALA_HALO_ELEMENT_RECORD_BASE_DOUBLES+ALA_VANKA_CHOL_SIZE",
            "ALA_vanka_chol[lev][m]+e*ALA_VANKA_CHOL_SIZE",
            "static double ala_halo_element_vanka_coupling(",
            "patch_records,n,i,j",
            '"halo_element_vanka"',
        ):
            self.assertIn(token, self.stokes)
        interface = _function_body(
            self.stokes, "static double ala_halo_element_vanka_coupling("
        )
        self.assertIn("support[50+i]", interface)
        self.assertIn("rhs_left[i]*solution[i]", interface)

    def test_vanka_pressure_aggregate_supports_full_local_schur(self) -> None:
        signature = "static void apply_ala_pressure_preconditioner("
        definition = self.stokes[self.stokes.rindex(signature) :]
        apply = _function_body(definition, signature)
        for token in (
            "constant_mode[i]=1.0",
            "numerator += constant_mode[i]*rhs[i]",
            "rhs[i] -= sum*constant_mode[i]",
            "numerator += constant_mode[i]*solution[i]",
            "solution[i] -= sum*constant_mode[i]",
            "z[m][e] += weight*work[m][e]",
            '"additive_weighted_highpass"',
            '"convex_blend_full_local_schur"',
            "if(element_vanka_highpass)",
        ):
            self.assertIn(token, apply)

        self.assertIn(
            "ala_shallow_patch_highpass         = off", self.strict_cfg
        )

        interface = apply[apply.index("cache->interface_blocks[m]") :]
        backward_solve = interface.index("for(i=n-1;i>=0;i--)")
        solution_projection = interface.index(
            "numerator += constant_mode[i]*solution[i]"
        )
        self.assertLess(backward_solve, solution_projection)

    def test_triangular_map_can_correct_only_its_shallow_defect(self) -> None:
        block = _function_body(
            self.stokes, "static void apply_ala_coupled_block_preconditioner("
        )
        triangular = block.index("apply_ala_coupled_triangular_once(")
        defect = block.index("assemble_ala_coupled_block_defect(", triangular)
        regional = block.index(
            "apply_ala_coupled_element_vanka_region(", defect
        )
        update = block.index(
            "ala_block_vector_axpy(E,1.0,correction_work,correction);",
            regional,
        )
        self.assertLess(triangular, defect)
        self.assertLess(defect, regional)
        self.assertLess(regional, update)
        self.assertIn("ala_coupled_shallow_vanka_layers,-1,0", block)
        self.assertIn("struct ala_block_vector *operator_work", block)
        self.assertIn(
            "action_work,vanka_delta,operator_work,", block
        )
        self.assertNotIn(
            "action_work,vanka_delta,velocity_work,", block
        )
        self.assertIn(
            "upper_block_triangular_plus_shallow_vanka", self.stokes
        )
        for source in (self.instructions, self.properties):
            self.assertIn(
                '"a coupled V-cycle or triangular coupled_fgmres"', source
            )

    def test_unaugmented_current_rheology_vanka_is_allowed(self) -> None:
        for source in (self.instructions, self.properties):
            self.assertIn(
                "ala_element_vanka_smoother requires multigrid, ", source
            )
            self.assertIn("and compressible_formulation=ala", source)
            self.assertNotIn(
                "ala_element_vanka_smoother requires multigrid, "
                '"\n              "compressible_formulation=ala, and positive gamma',
                source,
            )

    def test_schur_coarsening_audit_is_observational_only(
        self,
    ) -> None:
        audit = _function_body(
            self.stokes,
            "static void ala_audit_local_schur_coarsening(",
        )
        for token in (
            "coarse_schur/fine_sum",
            "relative_mismatch_mean=%e",
            "relative_mismatch_max=%e",
            "ALA SCHUR COARSENING AUDIT fine_level=%d",
            "action=observe_only",
        ):
            self.assertIn(token, audit)
        self.assertNotIn("coarse[m][ce] *=", audit)
        vcycle = _function_body(
            self.stokes,
            "static void apply_ala_coupled_multilevel_vcycle(",
        )
        restriction = vcycle.index(
            "ala_coupled_restrict_pressure_p0_transpose("
        )
        audit_call = vcycle.index(
            "ala_audit_local_schur_coarsening("
        )
        recursion = vcycle.index(
            "apply_ala_coupled_multilevel_vcycle(", audit_call + 1
        )
        self.assertLess(restriction, audit_call)
        self.assertLess(audit_call, recursion)

    def test_shallow_vanka_band_is_separate_and_windowed(self) -> None:
        region = _function_body(
            self.stokes,
            "static void apply_ala_coupled_element_vanka_region(",
        )
        for token in (
            "depth_layer=E->mesh.ELZ[lev]-1-global_radial_element",
            "band_only && depth_layer<selected_core",
            "taper_position=(depth_layer-selected_core+0.5)",
            "region_weight=sin(",
            "delta->pressure[m][e]=region_weight",
            "ala_element_vanka_pressure_damping",
            "delta->velocity[m][eq] += region_weight",
            'band_only ? "sine" : "off"',
        ):
            self.assertIn(token, region)
        vcycle = _function_body(
            self.stokes,
            "static void apply_ala_coupled_multilevel_vcycle(",
        )
        self.assertIn(
            "E->control.ala_coupled_shallow_vanka_band_sweeps", vcycle
        )
        self.assertIn("for(band_sweep=0;band_sweep<band_sweeps", vcycle)
        self.assertIn(
            "E->control.ala_coupled_shallow_vanka_core_layers,1", vcycle
        )

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
            "1.0/sqrt(correction->velocity[m][eq])",
            'shallow ? "shallow_core" : "full"',
        ):
            self.assertIn(token, region)


if __name__ == "__main__":
    unittest.main(verbosity=2)
