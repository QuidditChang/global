from __future__ import division

import math
import os
from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[2]
RUNS_ROOT = Path(os.environ.get(
    "CITCOMS_LZ_RUNS_ROOT",
    ROOT.parent / "runs-cmbhf-ALA-Leng_Zhong_2008"))


def read(relative):
    return (ROOT / relative).read_text()


def function_body(source, name):
    match = re.search(r"\b%s\s*\([^;]*?\)\s*(?:\n\s*[^\n{;]+;)*\s*\{" % name,
                      source, re.S)
    if not match:
        raise AssertionError("function %s not found" % name)
    start = match.end() - 1
    depth = 0
    for index in range(start, len(source)):
        if source[index] == "{":
            depth += 1
        elif source[index] == "}":
            depth -= 1
            if depth == 0:
                return source[start:index + 1]
    raise AssertionError("unterminated function %s" % name)


class Stage2ArchitectureTest(unittest.TestCase):

    def test_selector_is_explicit_and_defaults_off(self):
        instructions = read("lib/Instructions.c")
        properties = read("module/setProperties.c")
        component = read(
            "CitcomS/Components/Stokes_solver/Incompressible.py")
        definitions = read("lib/global_defs.h")

        self.assertIn("int ala_leng_zhong_2008;", definitions)
        self.assertIn('input_boolean("ala_leng_zhong_2008"', instructions)
        self.assertIn("E->control.ala_leng_zhong_2008 = 0;", instructions)
        self.assertIn('getIntProperty(properties, "ala_leng_zhong_2008"',
                      properties)
        self.assertRegex(component,
                         r'ala_leng_zhong_2008\s*=\s*prop\.bool\('
                         r'"ala_leng_zhong_2008",\s*default=False\)')
        self.assertRegex(component,
                         r'ala_leng_zhong_stage3\s*=\s*prop\.bool\('
                         r'"ala_leng_zhong_stage3",\s*default=False\)')
        self.assertRegex(component,
                         r'ala_leng_zhong_stage4\s*=\s*prop\.bool\('
                         r'"ala_leng_zhong_stage4",\s*default=False\)')
        self.assertRegex(component,
                         r'ala_leng_zhong_stage5\s*=\s*prop\.bool\('
                         r'"ala_leng_zhong_stage5",\s*default=False\)')
        self.assertRegex(
            component,
            r'ala_leng_zhong_residual_replacement_interval\s*=\s*prop\.int\('
            r'\s*"ala_leng_zhong_residual_replacement_interval",\s*'
            r'default=0\)')
        self.assertRegex(
            component,
            r'ala_leng_zhong_residual_drift_tolerance\s*=\s*prop\.float\('
            r'\s*"ala_leng_zhong_residual_drift_tolerance",\s*'
            r'default=0\.1\)')
        self.assertRegex(
            component,
            r'ala_leng_zhong_radial_line_preconditioner\s*=\s*prop\.bool\('
            r'\s*"ala_leng_zhong_radial_line_preconditioner",\s*'
            r'default=False\)')
        self.assertIn(
            'input_int("ala_leng_zhong_residual_replacement_interval"',
            instructions)
        self.assertIn(
            'getIntProperty(properties,\n'
            '                   "ala_leng_zhong_residual_replacement_interval"',
            properties)
        self.assertIn(
            'input_boolean("ala_leng_zhong_radial_line_preconditioner"',
            instructions)
        self.assertIn(
            '"ala_leng_zhong_radial_line_preconditioner"', properties)

    def test_stage2_has_hard_exclusion_guards(self):
        for relative in ("lib/Instructions.c", "module/setProperties.c"):
            source = read(relative)
            self.assertIn("Stage 2 requires gruneisen=0", source)
            self.assertIn("E->control.augmented_Lagr", source)
            self.assertIn("E->control.ala_augmented_lagrangian_gamma", source)
            self.assertIn("!E->control.precondition", source)
            self.assertIn("E->control.ala_pressure_multigrid", source)
            self.assertIn("E->control.ala_shallow_patch_preconditioner",
                          source)
            self.assertIn("E->control.ala_coupled_multilevel_vcycle", source)
            self.assertIn("Stage 2 requires uzawa=cg", source)

    def test_d_only_bpi_is_separate_from_full_b_bpi(self):
        element = read("lib/Element_calculations.c")
        definitions = read("lib/global_defs.h")
        builder = function_body(
            element, "build_diagonal_of_leng_zhong_Ahat")
        entry = function_body(
            element, "assemble_leng_zhong_dAhatp_entry")

        self.assertIn("double *LZ_BPI[MAX_LEVELS][NCS];", definitions)
        self.assertIn("LZ_BPI", builder)
        self.assertIn("assemble_leng_zhong_dAhatp_entry", builder)
        self.assertIn("isfinite", builder)
        self.assertIn("elt_del", entry)
        self.assertIn("ALA_velocity_BI", entry)
        self.assertNotIn("elt_c", entry)
        self.assertNotIn("ala_pressure_buoyancy", entry)

    def test_stage2_dispatches_only_to_ordinary_cg(self):
        stokes = read("lib/Stokes_flow_Incomp.c")
        dispatch = function_body(stokes, "solve_Ahat_p_fhat")
        cg = function_body(stokes, "solve_Ahat_p_fhat_CG")
        lz_preconditioner = function_body(
            stokes, "apply_leng_zhong_pressure_preconditioner")
        stage2 = dispatch[dispatch.index("if(E->control.ala_leng_zhong_2008)"):
                          dispatch.index("else if(E->control.inv_gruneisen")]

        self.assertIn("solve_Ahat_p_fhat_CG", stage2)
        self.assertNotIn("ALA_PCG", stage2)
        self.assertNotIn("BiCG", stage2)
        self.assertIn("apply_leng_zhong_pressure_preconditioner", cg)
        self.assertIn("E->LZ_BPI", lz_preconditioner)
        self.assertIn("assemble_grad_p", cg)
        self.assertIn("solve_del2_u", cg)
        self.assertIn("assemble_div_u", cg)
        self.assertNotIn("assemble_grad_rho_p", cg)
        self.assertNotIn("assemble_div_rho_u", cg)

    def test_zero_compressible_diagnostic_does_not_dereference_elt_c(self):
        element = read("lib/Element_calculations.c")
        assembler = function_body(element, "assemble_c_u")
        guard = assembler.index("if(E->control.inv_gruneisen == 0)")
        dereference = assembler.index("E->elt_c")
        self.assertLess(guard, dereference)
        self.assertIn("result[m][e] = 0.0", assembler[:dereference])

    def test_stage3_selector_keeps_c_out_of_inner_action(self):
        definitions = read("lib/global_defs.h")
        instructions = read("lib/Instructions.c")
        properties = read("module/setProperties.c")
        stokes = read("lib/Stokes_flow_Incomp.c")
        self.assertIn("ala_leng_zhong_stage3", definitions)
        self.assertIn('input_boolean("ala_leng_zhong_stage3"', instructions)
        self.assertIn('"ala_leng_zhong_stage3"', properties)
        self.assertIn('"LENG_ZHONG_STAGE3"', stokes)
        self.assertIn('operator=G^T_K^-1_G', stokes)
        inner = function_body(stokes, "solve_Ahat_p_fhat_CG")
        self.assertIn("c_rhs", inner)
        self.assertIn("F[m][j] += c_rhs[m][j]", inner)
        self.assertNotIn("assemble_grad_c_p", inner)

    def test_stage3_inner_cg_replaces_only_drifted_explicit_residual(self):
        stokes = read("lib/Stokes_flow_Incomp.c")
        inner = function_body(stokes, "solve_Ahat_p_fhat_CG")

        self.assertIn("frozen_reduction >= E->control.tole_comp", inner)
        self.assertIn("frozen_norm=max", inner.replace("sqrt(", ""))
        self.assertIn("frozen_norm/initial_frozen_norm", inner)
        self.assertIn("r2[m][j]-F[m][j]", inner)
        self.assertIn("r2[m][j]=F[m][j]", inner)
        self.assertIn("restart_direction = 1", inner)
        self.assertIn("ala_leng_zhong_residual_replacement_interval", inner)
        self.assertIn("ala_leng_zhong_residual_drift_tolerance", inner)
        self.assertRegex(
            inner,
            r"ala_leng_zhong_residual_replacement_interval>0\s*&&")
        self.assertIn("LENG_ZHONG_FROZEN_CONTINUITY", inner)
        stage3_loop = inner[inner.index("while("):
                            inner.index("} /* end loop for conjugate gradient */")]
        self.assertIn("E->control.ala_leng_zhong_stage3", stage3_loop)
        self.assertIn("dpressure >= imp", stage3_loop)
        self.assertRegex(stage3_loop,
                         r"ala_leng_zhong_stage3\s*\n\s*\?\s*"
                         r"\(frozen_reduction")

    def test_leng_radial_line_preconditioner_is_d_only_and_spd(self):
        element = read("lib/Element_calculations.c")
        arrays = read("lib/Construct_arrays.c")
        stokes = read("lib/Stokes_flow_Incomp.c")
        definitions = read("lib/global_defs.h")
        builder = function_body(
            element, "build_leng_zhong_radial_line_preconditioner")
        offdiag = function_body(
            element, "assemble_leng_zhong_Ahatp_jacobi_entry")
        apply = function_body(
            stokes, "apply_leng_zhong_pressure_preconditioner")

        self.assertIn("ala_leng_zhong_radial_line_preconditioner",
                      definitions)
        self.assertIn("LZ_BPI", builder)
        self.assertIn("pivot<=pivot_fraction*line_max", builder)
        self.assertIn("ALA_BPI_line_valid", builder)
        self.assertIn("build_leng_zhong_radial_line_preconditioner(E)",
                      arrays)
        self.assertIn("elt_del", offdiag)
        self.assertIn("ALA_velocity_BI", offdiag)
        self.assertNotIn("elt_c", offdiag)
        self.assertIn("ALA_BPI_line_lower", apply)
        self.assertIn("ALA_BPI_line_diag", apply)
        self.assertIn("LZ_BPI", apply)
        self.assertNotIn("BPI[", apply.replace("LZ_BPI[", ""))

    def test_python_bridge_validates_uzawa_after_loading_it(self):
        properties = read("module/setProperties.c")
        load = properties.index(
            'getStringProperty(properties, "uzawa", E->control.uzawa, fp)')
        validation = properties.index(
            'if(E->control.ala_leng_zhong_radial_line_preconditioner &&\n'
            '       strcmp(E->control.uzawa,"cg") != 0)')

        self.assertLess(load, validation)

    def test_stage4_lags_pressure_buoyancy_outside_inner_action(self):
        definitions = read("lib/global_defs.h")
        instructions = read("lib/Instructions.c")
        properties = read("module/setProperties.c")
        stokes = read("lib/Stokes_flow_Incomp.c")
        outer = function_body(stokes, "solve_Ahat_p_fhat_iterCG")
        inner = function_body(stokes, "solve_Ahat_p_fhat_CG")

        self.assertIn("int ala_leng_zhong_stage4;", definitions)
        self.assertIn('input_boolean("ala_leng_zhong_stage4"', instructions)
        self.assertIn('"ala_leng_zhong_stage4"', properties)
        self.assertIn("LENG_ZHONG_STAGE4", stokes)
        snapshot = outer.index("old_p[m][i] = P[m][i]")
        force = outer.index("assemble_forces(E,0)")
        solve = outer.index("solve_Ahat_p_fhat_CG")
        self.assertLess(snapshot, force)
        self.assertLess(force, solve)
        self.assertNotIn("assemble_grad_rho_p", inner)

    def test_stage5_audits_complete_fixed_point_after_outer_iteration(self):
        definitions = read("lib/global_defs.h")
        stokes = read("lib/Stokes_flow_Incomp.c")
        outer = function_body(stokes, "solve_Ahat_p_fhat_iterCG")
        momentum = function_body(
            stokes, "strict_ala_momentum_residual_audit")

        self.assertIn("int ala_leng_zhong_stage5;", definitions)
        self.assertIn("LENG_ZHONG_STAGE5_RESULT", outer)
        loop_end = outer.index("} /* end of while */")
        continuity = outer.index("assemble_div_rho_u(E,V,diff_p,lev)")
        momentum_audit = outer.index("strict_ala_momentum_residual_audit")
        self.assertLess(loop_end, continuity)
        self.assertLess(continuity, momentum_audit)
        self.assertIn("strict_ala_continuity_metrics", outer)
        self.assertIn("calloc(neq+1,sizeof(double))", outer)
        self.assertIn("assemble_unaugmented_del2_u", momentum)
        self.assertIn("assemble_forces(E,0)", momentum)
        self.assertIn("assemble_grad_p", momentum)
        self.assertNotIn("assemble_grad_rho_p", momentum)
        inner = function_body(stokes, "solve_Ahat_p_fhat_CG")
        self.assertNotIn("assemble_grad_rho_p", inner)
        self.assertIn("LENG_ZHONG_STAGE5_CONTINUITY_INNER", inner)
        self.assertIn("r1[m][j]+c_rhs[m][j]", inner)
        self.assertIn("r1[m][j]+r2[m][j]", inner)
        self.assertIn("r2[m][j]-c_rhs[m][j]", inner)

    def test_stage5_run_files_enable_fixed_point_audit(self):
        cfg = (RUNS_ROOT / "cmbhf_ALA_Leng_Zhong_2008.cfg").read_text()
        lsf = (RUNS_ROOT / "cmbhf_ALA_Leng_Zhong_2008.lsf").read_text()

        required = {
            "steps": "1",
            "tole_compressibility": "1e-06",
            "lith_age": "1",
            "lith_age_time": "1",
            "max_plate_age_Ma": "70",
            "gruneisen": "1",
            "compressible_formulation": "ala",
            "ala_leng_zhong_2008": "on",
            "ala_leng_zhong_stage3": "on",
            "ala_leng_zhong_stage4": "on",
            "ala_leng_zhong_stage5": "on",
            "ala_leng_zhong_residual_replacement_interval": "0",
            "ala_leng_zhong_residual_drift_tolerance": "0.1",
            "ala_leng_zhong_radial_line_preconditioner": "on",
            "ala_leng_zhong_stage5_continuity_tolerance": "1e-3",
            "ala_leng_zhong_stage5_momentum_tolerance": "1e-3",
            "ala_leng_zhong_stage5_initial_guess_scale": "1.0",
            "uzawa": "cg",
            "precond": "on",
            "aug_lagr": "off",
            "aug_number": "0.0",
            "ala_augmented_lagrangian_gamma": "0",
            "ala_outer_solver": "pcg",
            "ala_element_vanka_smoother": "off",
            "ala_radial_line_preconditioner": "off",
            "ala_two_level_preconditioner": "off",
            "ala_pressure_multigrid": "off",
            "ala_global_coarse_preconditioner": "off",
            "ala_shallow_patch_preconditioner": "off",
            "ala_geneo_preconditioner": "off",
            "ala_coupled_element_vanka": "off",
            "ala_coupled_multilevel_vcycle": "off",
            "ala_coupled_factor2_coarse_correction": "off",
            "ala_pressure_defect_corrections": "0",
        }
        for key, value in required.items():
            self.assertRegex(cfg, r"(?m)^%s\s*=\s*%s\s*$" %
                             (re.escape(key), re.escape(value)))

        self.assertIn("#BSUB -J CMBHF_LZ08_S5", lsf)
        self.assertIn("cmbhf_ALA_Leng_Zhong_2008.cfg", lsf)
        self.assertIn("builds/global/cmbhf_ALA_Leng_Zhong_2008", lsf)
        self.assertIn("stage5_lz_radial_line", lsf)
        self.assertIn("cmbhf_ALA_Leng_Zhong_2008/DATA/%RANK", cfg)
        self.assertIn("stage5_lz_radial_line_AhatP%P_%T", cfg)
        self.assertNotIn("runtime.cfg", lsf)
        self.assertNotIn("-i.bak", lsf)

    def test_small_spd_schur_pcg_oracle(self):
        # K=diag(2,3,5), G has two pressure columns, and S=G^T K^-1 G.
        k_inv = (0.5, 1.0 / 3.0, 0.2)
        g = ((1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
        s = [[sum(g[row][i] * k_inv[row] * g[row][j]
                  for row in range(3)) for j in range(2)]
             for i in range(2)]
        rhs = [1.25, -0.75]
        x = [0.0, 0.0]
        r = rhs[:]
        z = [r[i] / s[i][i] for i in range(2)]
        direction = z[:]
        rz = sum(r[i] * z[i] for i in range(2))

        for _ in range(2):
            action = [sum(s[i][j] * direction[j] for j in range(2))
                      for i in range(2)]
            curvature = sum(direction[i] * action[i] for i in range(2))
            self.assertGreater(curvature, 0.0)
            alpha = rz / curvature
            x = [x[i] + alpha * direction[i] for i in range(2)]
            r = [r[i] - alpha * action[i] for i in range(2)]
            z_new = [r[i] / s[i][i] for i in range(2)]
            rz_new = sum(r[i] * z_new[i] for i in range(2))
            beta = rz_new / rz
            direction = [z_new[i] + beta * direction[i] for i in range(2)]
            z, rz = z_new, rz_new

        residual = math.sqrt(sum(
            (sum(s[i][j] * x[j] for j in range(2)) - rhs[i]) ** 2
            for i in range(2)))
        self.assertLessEqual(residual, 1.0e-10)

    def test_tala_does_not_require_unused_supplied_beta(self):
        source = read("lib/Material_properties.c")
        initializer = function_body(source, "initialize_ala_beta")
        self.assertIn("density-derived beta", initializer)
        self.assertIn("E->control.ala_pressure_buoyancy", initializer)
        self.assertIn("ala_beta_element_source", initializer)


if __name__ == "__main__":
    unittest.main(verbosity=2)
