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
        stage2 = dispatch[dispatch.index("if(E->control.ala_leng_zhong_2008)"):
                          dispatch.index("else if(E->control.inv_gruneisen")]

        self.assertIn("solve_Ahat_p_fhat_CG", stage2)
        self.assertNotIn("ALA_PCG", stage2)
        self.assertNotIn("BiCG", stage2)
        self.assertIn("E->LZ_BPI", cg)
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
        self.assertIn("LENG_ZHONG_STAGE3 operator=G^T_K^-1_G", stokes)
        inner = function_body(stokes, "solve_Ahat_p_fhat_CG")
        self.assertIn("c_rhs", inner)
        self.assertIn("F[m][j] += c_rhs[m][j]", inner)
        self.assertNotIn("assemble_grad_c_p", inner)

    def test_stage2_run_files_exclude_later_stage_operators(self):
        cfg = (RUNS_ROOT / "cmbhf_ALA_Leng_Zhong_2008.cfg").read_text()
        lsf = (RUNS_ROOT / "cmbhf_ALA_Leng_Zhong_2008.lsf").read_text()

        required = {
            "gruneisen": "0",
            "compressible_formulation": "tala",
            "ala_leng_zhong_2008": "on",
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

        self.assertIn("#BSUB -J CMBHF_LZ08_S2", lsf)
        self.assertIn("cmbhf_ALA_Leng_Zhong_2008.cfg", lsf)
        self.assertIn("builds/global/cmbhf_ALA_Leng_Zhong_2008", lsf)
        self.assertIn("cmbhf_ALA_Leng_Zhong_2008_stage2", lsf)

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
