#!/usr/bin/env python3
"""Static regression checks for the production strict-ALA configuration."""

from __future__ import annotations

import re
import unittest
from pathlib import Path

import numpy as np


GLOBAL_ROOT = Path(__file__).resolve().parents[2]
PROJECT_ROOT = Path(__file__).resolve().parents[4]
RUNS_ROOT = PROJECT_ROOT / "runs"


def _active_cfg_lines(path: Path) -> list[str]:
    return [
        line
        for line in path.read_text(encoding="utf-8").splitlines()
        if not line.lstrip().startswith("#")
    ]


def _without_cfg_keys(lines: list[str], keys: tuple[str, ...]) -> list[str]:
    pattern = re.compile(
        r"^\s*(?:" + "|".join(map(re.escape, keys)) + r")\s*="
    )
    return [line for line in lines if not pattern.match(line)]


def _cfg_scalar(path: Path, name: str) -> float:
    text = path.read_text(encoding="utf-8")
    match = re.search(rf"^\s*{re.escape(name)}\s*=\s*(.*?)\s*$", text, re.M)
    if match is None:
        raise AssertionError(f"Missing {name} in {path.name}")
    return float(match.group(1))


class StrictProductionArchitectureTest(unittest.TestCase):
    def test_cfg_changes_only_strict_inputs_through_stage9d(
        self,
    ) -> None:
        legacy = _active_cfg_lines(RUNS_ROOT / "cmbhf_ALA.cfg")
        strict = _active_cfg_lines(RUNS_ROOT / "cmbhf_ALA_strict.cfg")
        strict_keys = (
            "steps",
            "walltime",
            "kC_ratio",
            "kC_primordial_flavor",
            "output_optional",
            "profile_optional",
            "tracer_reclassify_flavors",
            "bottom_tbl_thickness",
            "Tbottom",
            "deltaT",
            "surfaceT",
            "k0",
            "kT_exponent",
            "rayleigh",
            "phase_Ra",
            "phase_clapeyron",
            "buoyancy_ratio",
            "ala_augmented_lagrangian_gamma",
            "ala_element_vanka_damping",
            "ala_element_vanka_pressure_damping",
            "ala_element_vanka_external_diagonal_weight",
            "ala_element_vanka_galerkin_schur",
            "ala_element_vanka_rebuild_interval",
            "refstate_file",
            "piterations",
            "ala_beta_element_source",
            "ala_beta_interval_file",
            "ala_shallow_patch_velocity_solver",
            "ala_shallow_patch_preconditioner",
            "ala_shallow_patch_weight",
            "ala_shallow_patch_horizontal_elements",
            "ala_shallow_patch_horizontal_stride",
            "ala_shallow_patch_mpi_overlap",
            "ala_pcg_restart_interval",
            "ala_outer_solver",
            "ala_coupled_initial_velocity_relative_tolerance",
            "ala_coupled_inner_relative_tolerance",
            "ala_coupled_inner_max_cycles",
            "ala_coupled_inner_progress_interval",
            "ala_coupled_defect_corrections",
            "ala_coupled_factor2_coarse_correction",
            "ala_coupled_multilevel_audit_only",
            "ala_coupled_first_preconditioner_audit_only",
            "ala_coupled_debug_stop_iteration",
            "ala_coupled_element_vanka",
            "ala_coupled_multilevel_vcycle",
            "ala_coupled_multilevel_coarse_sweeps",
            "ala_coupled_multilevel_coarse_weight",
            "ala_viscosity_spectrum_diagnostics",
            "ala_viscosity_spectrum_interval",
            "ala_coupled_shallow_vanka_layers",
            "ala_coupled_shallow_vanka_core_layers",
            "ala_coupled_shallow_vanka_band_sweeps",
            "ala_coupled_shallow_vanka_sweeps",
            "ala_coupled_shallow_vanka_warm_sweeps",
            "ala_depth_diagnostics",
            "ala_depth_diagnostic_interval",
            "ala_coarse_residual_diagnostics",
            "ala_coarse_residual_interval",
            "ala_coarse_residual_levels",
            "ala_two_level_preconditioner",
            "ala_two_level_coarse_solver",
            "ala_two_level_coarse_iterations",
            "ala_two_level_coarse_damping",
            "ala_two_level_velocity_solver",
            "file_vbcs",
            "lith_age",
            "lith_age_time",
            "ala_unaugmented_momentum_tolerance",
            "ala_geneo_preconditioner",
            "ala_geneo_basis_type",
            "ala_pressure_multigrid",
            "ala_pressure_multigrid_galerkin",
            "ala_pressure_bpi_weight",
            "ala_pressure_factor2_coarse_action_scale",
            "ala_pressure_multigrid_min_level",
            "ala_pressure_multigrid_pre_smooth",
            "ala_pressure_multigrid_post_smooth",
            "ala_pressure_multigrid_coarse_iterations",
            "ala_pressure_multigrid_damping",
            "ala_pressure_multigrid_weight",
        )
        self.assertEqual(
            _without_cfg_keys(legacy, strict_keys),
            _without_cfg_keys(strict, strict_keys),
        )
        strict_text = "\n".join(strict)
        self.assertRegex(
            strict_text,
            r"(?m)^\s*refstate_file\s*=\s*refstate_ALA_strict\.txt\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_beta_element_source\s*=\s*interval\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_beta_interval_file\s*="
            r"\s*interval_ALA_strict\.txt\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_coupled_multilevel_audit_only\s*=\s*off\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_coupled_first_preconditioner_audit_only"
            r"\s*=\s*off\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_coupled_debug_stop_iteration\s*=\s*20\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_depth_diagnostics\s*=\s*on\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_coarse_residual_diagnostics\s*=\s*on\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*piterations\s*=\s*220\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_coupled_element_vanka\s*=\s*off\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_coupled_multilevel_vcycle\s*=\s*off\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_coupled_multilevel_coarse_weight\s*=\s*1\.0\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_viscosity_spectrum_diagnostics\s*=\s*on\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_shallow_patch_velocity_solver\s*="
            r"\s*element_vanka\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_shallow_patch_preconditioner\s*=\s*on\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_shallow_patch_regularization\s*=\s*1\.0e-3\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_pcg_restart_interval\s*=\s*50\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_coarse_residual_levels\s*=\s*4\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*file_vbcs\s*=\s*1\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*lith_age\s*=\s*1\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*lith_age_time\s*=\s*1\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*steps\s*=\s*1\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*walltime\s*=\s*168\s*\*hour\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*kC_ratio\s*=\s*0\.8\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*kC_primordial_flavor\s*=\s*24\s*$",
        )
        profile = re.search(
            r"(?m)^\s*profile_optional\s*=\s*(.*?)\s*$", strict_text
        )
        self.assertIsNotNone(profile)
        self.assertIn("kC", profile.group(1).split(","))
        output = re.search(
            r"(?m)^\s*output_optional\s*=\s*(.*?)\s*$", strict_text
        )
        self.assertIsNotNone(output)
        self.assertIn("k", output.group(1).split(","))
        interfaces = re.search(
            r"(?m)^\s*z_interface\s*=\s*(.*?)\s*$", strict_text
        )
        self.assertIsNotNone(interfaces)
        self.assertEqual(
            [float(value) for value in interfaces.group(1).split(",")],
            [0.585] * 24,
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*tracer_reclassify_flavors\s*=\s*off\s*$",
        )
        self.assertAlmostEqual(
            _cfg_scalar(RUNS_ROOT / "cmbhf_ALA_strict.cfg",
                        "bottom_tbl_thickness") * 6371.0,
            500.0,
            places=6,
        )
        lsf = (RUNS_ROOT / "cmbhf_ALA_strict.lsf").read_text(
            encoding="utf-8"
        )
        self.assertRegex(lsf, r"(?m)^#BSUB\s+-W\s+168:00\s+")
        self.assertRegex(
            strict_text,
            r"(?m)^\s*monitoringFrequency\s*=\s*50\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*checkpointFrequency\s*=\s*50\s*$",
        )

        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_element_vanka_rebuild_interval\s*=\s*1\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_coupled_multilevel_coarse_sweeps\s*=\s*2\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_coupled_shallow_vanka_layers\s*=\s*0\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_coupled_shallow_vanka_core_layers\s*=\s*-1\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_coupled_shallow_vanka_band_sweeps\s*=\s*0\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_coupled_shallow_vanka_sweeps\s*=\s*0\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_coupled_shallow_vanka_warm_sweeps\s*=\s*-1\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_outer_solver\s*=\s*coupled_fgmres\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_coupled_initial_velocity_relative_tolerance"
            r"\s*=\s*0\.0\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_augmented_lagrangian_gamma\s*=\s*10\.0\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_element_vanka_damping\s*=\s*0\.2\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_element_vanka_pressure_damping\s*=\s*0\.25\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_element_vanka_external_diagonal_weight\s*=\s*1\.0\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_coupled_inner_relative_tolerance\s*=\s*1e-2\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_coupled_inner_max_cycles\s*=\s*200\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_coupled_inner_progress_interval\s*=\s*20\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_coupled_defect_corrections\s*=\s*0\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_coupled_factor2_coarse_correction\s*=\s*off\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_unaugmented_momentum_tolerance\s*=\s*1e-3\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_geneo_preconditioner\s*=\s*off\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_geneo_basis_type\s*="
            r"\s*radial_partition\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_pressure_multigrid\s*=\s*off\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_two_level_preconditioner\s*=\s*off\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_two_level_coarse_solver\s*=\s*scaled_diagonal\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_two_level_coarse_iterations\s*=\s*8\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_two_level_coarse_damping\s*=\s*0\.07\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_two_level_velocity_solver\s*=\s*diagonal\s*$",
        )
        damping_match = re.search(
            r"(?m)^\s*ala_two_level_coarse_damping\s*=\s*([0-9.eE+-]+)",
            strict_text,
        )
        self.assertIsNotNone(damping_match)
        self.assertGreater(float(damping_match.group(1)), 0.0)
        self.assertLess(float(damping_match.group(1)), 2.0 / 27.0)
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_pressure_multigrid_galerkin\s*=\s*off\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_pressure_bpi_weight\s*=\s*1\.0\s*$",
        )
        self.assertRegex(
            strict_text,
            r"(?m)^\s*ala_pressure_multigrid_min_level\s*=\s*0\s*$",
        )

    def test_kappa0_is_authoritative_and_k0_is_closed(self) -> None:
        path = RUNS_ROOT / "cmbhf_ALA_strict.cfg"
        text = path.read_text(encoding="utf-8")
        kappa0 = 1.0e-6
        rho0 = _cfg_scalar(path, "rho0")
        cp0 = _cfg_scalar(path, "Cp0")
        k0 = _cfg_scalar(path, "k0")
        self.assertIn("# @physical kappa0_m2_s = 1.0e-6", text)
        self.assertAlmostEqual(k0, kappa0 * rho0 * cp0, places=12)

    def test_strict_reference_schema_and_tref_endpoints(self) -> None:
        path = RUNS_ROOT / "refstate_ALA_strict.txt"
        comments = [
            line[1:].strip()
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.startswith("#")
        ]
        self.assertEqual(
            tuple(comments[1].split()),
            ("rho", "g", "Tref", "alpha", "Cp", "beta", "Gamma_eff", "Ks"),
        )
        table = np.loadtxt(path, comments="#")
        self.assertEqual(table.shape[1], 8)
        self.assertGreater(float(np.min(table[:, 7])), 0.0)
        self.assertGreater(abs(float(table[0, 2]) - 1.0), 1.0e-6)
        self.assertGreater(abs(float(table[-1, 2])), 1.0e-6)

    def test_gamma_effective_closes_with_cfg_di(self) -> None:
        table = np.loadtxt(
            RUNS_ROOT / "refstate_ALA_strict.txt", comments="#"
        )
        di = _cfg_scalar(
            RUNS_ROOT / "cmbhf_ALA_strict.cfg", "dissipation_number"
        )
        beta_check = (
            di * table[:, 3] * table[:, 1] / (table[:, 4] * table[:, 6])
        )
        relative = np.abs(beta_check - table[:, 5]) / np.abs(table[:, 5])
        self.assertLessEqual(float(np.max(relative)), 8.0 * np.finfo(float).eps)

    def test_interval_beta_is_serialized_density_log_secant(self) -> None:
        interval_path = RUNS_ROOT / "interval_ALA_strict.txt"
        comments = [
            line[1:].strip()
            for line in interval_path.read_text(encoding="utf-8").splitlines()
            if line.startswith("#")
        ]
        self.assertEqual(
            tuple(comments[1].split()),
            ("element_index", "r_inner", "r_outer", "beta_interval"),
        )
        intervals = np.loadtxt(interval_path, comments="#")
        refstate = np.loadtxt(
            RUNS_ROOT / "refstate_ALA_strict.txt", comments="#"
        )
        radii = np.loadtxt(
            RUNS_ROOT / "GLB.coor.global.dat", skiprows=1, usecols=1
        )
        expected = -np.diff(np.log(refstate[:, 0])) / np.diff(radii)
        residual = (
            intervals[:, 3] * np.diff(radii)
            + np.diff(np.log(refstate[:, 0]))
        )
        self.assertEqual(intervals.shape, (len(refstate) - 1, 4))
        self.assertTrue(
            np.array_equal(intervals[:, 0], np.arange(1, len(refstate)))
        )
        self.assertTrue(np.array_equal(intervals[:, 1], radii[:-1]))
        self.assertTrue(np.array_equal(intervals[:, 2], radii[1:]))
        self.assertTrue(np.array_equal(intervals[:, 3], expected))
        self.assertLessEqual(float(np.max(np.abs(residual))), 1.0e-14)

    def test_reader_retains_legacy_seven_column_compatibility(self) -> None:
        legacy = np.loadtxt(RUNS_ROOT / "refstate_ALA.txt", comments="#")
        material = (GLOBAL_ROOT / "lib" / "Material_properties.c").read_text()
        self.assertEqual(legacy.shape[1], 7)
        self.assertIn("columns != 7 && columns != 8", material)
        self.assertIn(
            "E->refstate.Ks[i] = columns == 8 ? values[7] : 0.0;",
            material,
        )

    def test_interval_reader_feeds_single_element_beta_authority(self) -> None:
        material = (GLOBAL_ROOT / "lib" / "Material_properties.c").read_text()
        element = (GLOBAL_ROOT / "lib" / "Element_calculations.c").read_text()
        instructions = (GLOBAL_ROOT / "lib" / "Instructions.c").read_text()
        self.assertIn("read_ala_beta_intervals(E);", material)
        self.assertIn(
            "beta = E->refstate.ala_beta_interval[nz];", material
        )
        self.assertIn("E->refstate.ala_beta[nz] = beta;", material)
        self.assertIn("E->refstate.ala_beta[fine_nz] * dr", element)
        self.assertIn('"interval") != 0', instructions)

    def test_interval_mesh_identity_allows_legacy_float_roundoff(self) -> None:
        material = (GLOBAL_ROOT / "lib" / "Material_properties.c").read_text()
        intervals = np.loadtxt(
            RUNS_ROOT / "interval_ALA_strict.txt", comments="#"
        )
        serialized_radii = np.concatenate(
            (intervals[:1, 1], intervals[:, 2])
        )
        runtime_radii = serialized_radii.astype(np.float32).astype(np.float64)
        maximum_roundoff = float(
            np.max(np.abs(serialized_radii - runtime_radii))
        )
        tolerance = 4.0 * np.finfo(np.float32).eps

        self.assertGreater(maximum_roundoff, 1.0e-10)
        self.assertLessEqual(maximum_roundoff, tolerance)
        self.assertIn("#include <float.h>", material)
        self.assertIn(
            "const double radius_tolerance = 4.0 * FLT_EPSILON;",
            material,
        )

    def test_pyre_bindings_cover_strict_runtime_options(self) -> None:
        param = (
            GLOBAL_ROOT / "CitcomS" / "Components" / "Param.py"
        ).read_text()
        incompressible = (
            GLOBAL_ROOT
            / "CitcomS"
            / "Components"
            / "Stokes_solver"
            / "Incompressible.py"
        ).read_text()
        properties = (GLOBAL_ROOT / "module" / "setProperties.c").read_text()

        self.assertIn('ala_beta_interval_file = pyre.inventory.str(', param)
        self.assertIn('"interval"]))', incompressible)
        for name in (
            "ala_shallow_patch_horizontal_elements",
            "ala_shallow_patch_horizontal_stride",
        ):
            self.assertIn(f'{name} = prop.int(', incompressible)
            self.assertIn(f'getIntProperty(properties, "{name}"', properties)
        self.assertIn(
            'getStringProperty(properties, "ala_beta_interval_file"',
            properties,
        )
        self.assertIn('"density_log_secant") != 0 &&', properties)
        self.assertIn('"interval") != 0)', properties)
        self.assertIn(
            "twice ala_shallow_patch_mpi_overlap must not exceed",
            properties,
        )
        self.assertIn(
            'ala_shallow_patch_velocity_solver = prop.str(',
            incompressible,
        )
        self.assertIn(
            '["diagonal", "node_block", "element_vanka"]',
            incompressible,
        )
        self.assertIn(
            'getStringProperty(properties, '
            '"ala_shallow_patch_velocity_solver"',
            properties,
        )
        self.assertIn(
            'strcmp(E->control.ala_shallow_patch_velocity_solver,'
            '"element_vanka") != 0',
            properties,
        )
        self.assertIn(
            'ala_geneo_basis_type = prop.str(',
            incompressible,
        )
        self.assertIn(
            'getStringProperty(properties, "ala_geneo_basis_type"',
            properties,
        )
        self.assertIn(
            'ala_unaugmented_momentum_tolerance = prop.float(',
            incompressible,
        )
        self.assertIn(
            'getDoubleProperty(properties, '
            '"ala_unaugmented_momentum_tolerance"',
            properties,
        )
        for name in (
            "ala_coupled_inner_max_cycles",
            "ala_coupled_inner_progress_interval",
            "ala_coupled_defect_corrections",
        ):
            self.assertIn(f'{name} = prop.int(', incompressible)
            self.assertIn(f'getIntProperty(properties, "{name}"', properties)
        self.assertIn(
            "ala_coupled_factor2_coarse_correction = prop.bool(",
            incompressible,
        )
        self.assertIn(
            'getIntProperty(properties, '
            '"ala_coupled_factor2_coarse_correction"',
            properties,
        )
        for name in (
            "ala_coupled_initial_velocity_relative_tolerance",
            "ala_coupled_inner_relative_tolerance",
        ):
            self.assertIn(f"{name} = prop.float(", incompressible)
            self.assertIn("getDoubleProperty(", properties)
            self.assertIn(f'properties, "{name}"', properties)
        for name in (
            "ala_pressure_multigrid",
            "ala_pressure_multigrid_galerkin",
            "ala_pressure_multigrid_min_level",
            "ala_pressure_multigrid_pre_smooth",
            "ala_pressure_multigrid_post_smooth",
            "ala_pressure_multigrid_coarse_iterations",
        ):
            self.assertIn(f'{name} = prop.', incompressible)
            self.assertIn(f'getIntProperty(properties, "{name}"', properties)
        for name in (
            "ala_pressure_bpi_weight",
            "ala_pressure_factor2_coarse_action_scale",
            "ala_pressure_multigrid_damping",
            "ala_pressure_multigrid_weight",
        ):
            self.assertIn(f'{name} = prop.float(', incompressible)
            self.assertIn(
                f'getDoubleProperty(properties, "{name}"', properties
            )

    def test_pressure_bpi_weight_scales_only_the_legacy_schur_map(
        self,
    ) -> None:
        stokes = (
            GLOBAL_ROOT / "lib" / "Stokes_flow_Incomp.c"
        ).read_text(encoding="utf-8")
        self.assertIn(
            "z[m][j] *= E->control.ala_pressure_bpi_weight;",
            stokes,
        )
        self.assertIn("BPI_weight=%e", stokes)

    def test_scaled_diagonal_coarse_map_is_recurrence_free(self) -> None:
        instructions = (
            GLOBAL_ROOT / "lib" / "Instructions.c"
        ).read_text(encoding="utf-8")
        pyre = (
            GLOBAL_ROOT
            / "CitcomS/Components/Stokes_solver/Incompressible.py"
        ).read_text(encoding="utf-8")
        properties = (
            GLOBAL_ROOT / "module" / "setProperties.c"
        ).read_text(encoding="utf-8")
        stokes = (
            GLOBAL_ROOT / "lib" / "Stokes_flow_Incomp.c"
        ).read_text(encoding="utf-8")
        self.assertIn('"scaled_diagonal") != 0', instructions)
        self.assertIn('["jacobi", "chebyshev", "scaled_diagonal"]', pyre)
        self.assertIn('"scaled_diagonal") != 0', properties)
        self.assertIn(
            "scaled_diagonal coarse solver requires diagonal pressure BPI",
            instructions,
        )
        branch = stokes[stokes.index('"scaled_diagonal")==0') :]
        branch = branch[: branch.index("else if")]
        self.assertIn("ala_two_level_coarse_iterations*damping", branch)
        self.assertIn("*cache->coarse_bpi[m][ce]*coarse_rhs[m][ce]", branch)
        self.assertIn("n=factor*factor*factor", branch)
        self.assertIn("-coarse_residual[m][ce]/(double)n", branch)
        self.assertIn("+coarse_Ax[m][ce]/(double)n", branch)
        self.assertNotIn("apply_ala_galerkin_fixed_schur", branch)

    def test_stage8_pressure_multigrid_is_recursive_strict_ala(
        self,
    ) -> None:
        stokes = (
            GLOBAL_ROOT / "lib" / "Stokes_flow_Incomp.c"
        ).read_text(encoding="utf-8")
        instructions = (
            GLOBAL_ROOT / "lib" / "Instructions.c"
        ).read_text(encoding="utf-8")

        operator = stokes[
            stokes.index("apply_ala_pressure_multigrid_operator"):
            stokes.index("smooth_ala_pressure_multigrid_level")
        ]
        self.assertIn("assemble_grad_rho_p", operator)
        self.assertIn("E->ALA_velocity_BI[level]", operator)
        self.assertIn("assemble_div_rho_u", operator)
        self.assertIn("pressure_mg_galerkin_p", operator)
        self.assertIn(
            "apply_ala_pressure_multigrid_operator(\n"
            "            E,cache->pressure_mg_galerkin_p,",
            operator,
        )
        self.assertIn(
            "Ap[m][ce] += cache->pressure_mg_galerkin_Ap[m][e]\n"
            "                            /sqrt((double)(factor*factor*factor));",
            operator,
        )
        vcycle = stokes[
            stokes.index("static void ala_pressure_multigrid_vcycle"):
            stokes.index("static void apply_ala_pressure_multigrid_correction")
        ]
        self.assertIn("ala_pressure_multigrid_vcycle(E,level-1,cache);", vcycle)
        self.assertIn(
            "coarse_rhs[m][ce] += residual[m][e]/sqrt(8.0);", vcycle
        )
        self.assertIn(
            "x[m][e] += coarse_x[m][ce]/sqrt(8.0);", vcycle
        )
        self.assertIn("ala_pressure_multigrid_post_smooth", vcycle)
        self.assertIn(
            "mpi_overlap_schwarz_plus_pressure_vcycle", stokes
        )
        self.assertIn(
            "ala_pressure_multigrid requires precond=on, strict ALA",
            instructions,
        )
        audit_start = stokes.index(
            "static void audit_ala_shallow_patch_preconditioner"
        )
        audit = stokes[
            audit_start:
            stokes.index("/* Master loop for pressure", audit_start)
        ]
        disable = audit.index("E->control.ala_pressure_multigrid=0;")
        first_apply = audit.index("apply_ala_pressure_preconditioner(", disable)
        second_apply = audit.index(
            "apply_ala_pressure_preconditioner(", first_apply + 1
        )
        restore = audit.index(
            "E->control.ala_pressure_multigrid=pressure_multigrid_enabled;",
            second_apply,
        )
        self.assertLess(disable, first_apply)
        self.assertLess(first_apply, second_apply)
        self.assertLess(second_apply, restore)
        self.assertIn("scope=base_without_pressure_mg", audit)

    def test_stage6c_shallow_schur_is_spd_gram_with_rollback(self) -> None:
        stokes = (
            GLOBAL_ROOT / "lib" / "Stokes_flow_Incomp.c"
        ).read_text(encoding="utf-8")
        instructions = (
            GLOBAL_ROOT / "lib" / "Instructions.c"
        ).read_text(encoding="utf-8")
        definitions = (
            GLOBAL_ROOT / "lib" / "global_defs.h"
        ).read_text(encoding="utf-8")

        self.assertIn("ala_build_velocity_node_factor", stokes)
        self.assertIn("ala_velocity_node_feature", stokes)
        self.assertIn(
            "value += left[74+p1+d]*right[74+p2+d];",
            stokes,
        )
        self.assertIn(
            "=assemble_Ahatp_jacobi_entry(E,e1,e2,lev,m);",
            stokes,
        )
        self.assertIn("matrix[i][j]=patch_records[i][1];", stokes)
        self.assertIn(
            "bi=0.5*(left[50+p1+d]+right[50+p2+d]);",
            stokes,
        )
        self.assertIn(
            "operator=principal((D+C)Mv^-1(D+C)^T)",
            stokes,
        )
        self.assertIn("velocity_block_fallbacks=%d", stokes)
        self.assertIn(
            'strcpy(E->control.ala_shallow_patch_velocity_solver,'
            '"diagonal");',
            instructions,
        )
        self.assertIn(
            "char ala_shallow_patch_velocity_solver[20];",
            definitions,
        )

    def test_stage6d_cross_rank_radial_partition_is_fixed_galerkin(
        self,
    ) -> None:
        stokes = (
            GLOBAL_ROOT / "lib" / "Stokes_flow_Incomp.c"
        ).read_text(encoding="utf-8")
        instructions = (
            GLOBAL_ROOT / "lib" / "Instructions.c"
        ).read_text(encoding="utf-8")
        definitions = (
            GLOBAL_ROOT / "lib" / "global_defs.h"
        ).read_text(encoding="utf-8")

        self.assertIn("ala_build_radial_partition_shapes", stokes)
        self.assertIn("(j%rbins==mode) ? 1.0 : 0.0", stokes)
        radial_builder = stokes[
            stokes.index("ala_build_radial_partition_shapes"):
            stokes.index("ala_geneo_jacobi_eigensolve")
        ]
        self.assertNotIn("active_map", radial_builder)
        self.assertNotIn(
            "if(!isfinite(denominator) || denominator<=1.0e-30)",
            radial_builder,
        )
        self.assertIn(
            "? ala_build_radial_partition_shapes(",
            stokes,
        )
        self.assertIn(
            "desired_modes=(shallow_layers>0) ? rbins : 0;",
            stokes,
        )
        self.assertIn(
            "accepted=(desired_modes>0)",
            stokes,
        )
        collective = stokes.index(
            "MPI_Bcast(&group_modes[m],1,MPI_INT,0,group_comm[m]);"
        )
        synchronized_error = stokes.index(
            "ALA Stage 6d radial partition basis is incomplete",
            collective,
        )
        self.assertLess(collective, synchronized_error)
        self.assertIn("geometric_radial_partition", stokes)
        self.assertIn(
            "mpi_overlap_schwarz_plus_cross_rank_radial_aggregate",
            stokes,
        )
        self.assertIn("apply_ala_galerkin_fixed_schur(", stokes)
        self.assertIn(
            "radial_partition requires a cross-rank GenEO group",
            instructions,
        )
        self.assertIn("char ala_geneo_basis_type[24];", definitions)

    def test_stage7_fgmres_reconstruction_and_residual_audits(
        self,
    ) -> None:
        stokes = (
            GLOBAL_ROOT / "lib" / "Stokes_flow_Incomp.c"
        ).read_text(encoding="utf-8")
        dispatch = stokes[
            stokes.index(
                'if(strcmp(E->control.ala_outer_solver,"fgmres")==0 ||'
            ):
            stokes.index(
                "/* FF contains the current -C^T*P forcing.",
            )
        ]
        self.assertLess(
            dispatch.index("initial_vel_residual(E,V,P,F,imp);"),
            dispatch.index("solve_ala_fgmres_core("),
        )
        fgmres = stokes[
            stokes.index("static float solve_ala_fgmres_core"):
            stokes.index("static float solve_Ahat_p_fhat_iterCG")
        ]
        self.assertIn(
            "algebraic_relative=explicit_norm/initial_rnorm;",
            fgmres,
        )
        self.assertIn(
            "residual_drift=residual_est\n"
            "                /max(algebraic_relative,1.0e-300);",
            fgmres,
        )
        self.assertNotIn(
            "residual_est/max(cancellation_l2,1.0e-300)",
            fgmres,
        )
        self.assertIn("mass_norm=%e mass_relative=%e", fgmres)
        self.assertIn("Q=%e Q_relative=%e", fgmres)
        self.assertIn(
            "strict_ala_momentum_residual_audit(",
            fgmres,
        )
        self.assertIn(
            "tmpF[m]=(double *)calloc(neq+1,sizeof(double));",
            fgmres,
        )
        self.assertIn(
            "tmpU[m]=(double *)calloc(neq+1,sizeof(double));",
            fgmres,
        )
        self.assertIn(
            "ub[j][m]=(double *)calloc(neq+1,sizeof(double));",
            fgmres,
        )
        candidate = fgmres.index(
            "if(continuity_converged && momentum_gate)"
        )
        candidate_audit = fgmres.index(
            "strict_ala_momentum_residual_audit(", candidate
        )
        joint_acceptance = fgmres.index(
            "converged=(continuity_converged && momentum_converged);",
            candidate_audit,
        )
        self.assertLess(candidate_audit, joint_acceptance)
        self.assertIn("momentum_target_not_reached", fgmres)
        self.assertIn("joint_target_reached", fgmres)
        self.assertIn(
            'E,V,P,tmpF,tmpU,lev,"startup",0,',
            fgmres,
        )
        self.assertIn(
            'E,V,P,tmpF,tmpU,lev,"restart",count,',
            fgmres,
        )
        self.assertIn("if(converged || arnoldi_breakdown)", fgmres)
        self.assertLess(
            fgmres.rindex("ALA FGMRES momentum audit"),
            fgmres.index("Strict ALA FGMRES failed joint acceptance"),
        )

    def test_stage7c_momentum_decomposition_is_low_memory_and_exact(
        self,
    ) -> None:
        stokes = (
            GLOBAL_ROOT / "lib" / "Stokes_flow_Incomp.c"
        ).read_text(encoding="utf-8")
        decomposition_start = stokes.index(
            "static void strict_ala_momentum_decomposition_audit(",
            stokes.index("Audit the original, unaugmented momentum"),
        )
        decomposition = stokes[
            decomposition_start:
            stokes.index(
                "static void strict_ala_depth_diagnostics",
                decomposition_start,
            )
        ]
        self.assertIn("assemble_del2_u(E,V,work_a,lev,1);", decomposition)
        self.assertIn(
            "assemble_unaugmented_del2_u(E,V,work_b,lev,1);",
            decomposition,
        )
        self.assertIn(
            "assemble_ala_augmented_u(E,V,work_b,lev,1);",
            decomposition,
        )
        self.assertIn("raw_penalty_cosine", decomposition)
        self.assertIn("split_defect", decomposition)
        self.assertIn("residual_norm_defect", decomposition)
        self.assertNotIn("calloc", decomposition)
        self.assertNotIn("malloc", decomposition)

    def test_stage7d_refreshes_strict_force_after_stiffness_update(
        self,
    ) -> None:
        drive = (GLOBAL_ROOT / "lib" / "Drive_solvers.c").read_text(
            encoding="utf-8"
        )
        solver_start = drive.index(
            "void general_stokes_solver(struct All_variables *E)"
        )
        solver = drive[
            solver_start:
            drive.index("void general_stokes_solver_pseudo_surf", solver_start)
        ]
        stiffness = solver.index("construct_stiffness_B_matrix(E);")
        strict_guard = solver.index(
            "if(E->control.ala_pressure_buoyancy &&", stiffness
        )
        gamma_guard = solver.index(
            "E->control.ala_augmented_lagrangian_gamma > 0.0",
            strict_guard,
        )
        refreshed_force = solver.index("assemble_forces(E,0);", gamma_guard)
        solve = solver.index("solve_constrained_flow_iterative(E);")
        self.assertLess(stiffness, strict_guard)
        self.assertLess(strict_guard, gamma_guard)
        self.assertLess(gamma_guard, refreshed_force)
        self.assertLess(refreshed_force, solve)
        self.assertEqual(solver.count("assemble_forces(E,0);"), 2)
        self.assertIn(
            "ALA strict force reassembled after stiffness update",
            solver,
        )

    def test_stage9g_lags_only_vanka_preconditioner_factors(self) -> None:
        construct = (GLOBAL_ROOT / "lib" / "Construct_arrays.c").read_text(
            encoding="utf-8"
        )
        instructions = (GLOBAL_ROOT / "lib" / "Instructions.c").read_text(
            encoding="utf-8"
        )
        properties = (
            GLOBAL_ROOT / "module" / "setProperties.c"
        ).read_text(encoding="utf-8")
        pyre = (
            GLOBAL_ROOT
            / "CitcomS/Components/Stokes_solver/Incompressible.py"
        ).read_text(encoding="utf-8")
        output = (GLOBAL_ROOT / "lib" / "Output_h5.c").read_text(
            encoding="utf-8"
        )

        self.assertIn(
            'input_int("ala_element_vanka_rebuild_interval",', instructions
        )
        self.assertIn(
            'input_double("ala_element_vanka_pressure_damping",', instructions
        )
        self.assertIn(
            'input_double("ala_element_vanka_external_diagonal_weight",',
            instructions,
        )
        self.assertIn(
            'input_boolean("ala_element_vanka_galerkin_schur",', instructions
        )
        self.assertIn(
            '"ala_element_vanka_rebuild_interval", default=1', pyre
        )
        self.assertIn(
            '"ala_element_vanka_pressure_damping", default=1.0', pyre
        )
        self.assertIn(
            '"ala_element_vanka_external_diagonal_weight", default=1.0', pyre
        )
        self.assertIn(
            '"ala_element_vanka_galerkin_schur", default=False', pyre
        )
        self.assertIn(
            'getIntProperty(properties, "ala_element_vanka_rebuild_interval",',
            properties,
        )
        self.assertIn(
            'getDoubleProperty(properties, "ala_element_vanka_pressure_damping",',
            properties,
        )
        self.assertIn(
            'getDoubleProperty(properties, '
            '"ala_element_vanka_external_diagonal_weight",', properties
        )
        self.assertIn(
            'getIntProperty(properties, "ala_element_vanka_galerkin_schur",',
            properties,
        )
        self.assertIn("ala_element_vanka_rebuild_interval", output)
        self.assertIn("ala_element_vanka_pressure_damping", output)
        self.assertIn("ala_element_vanka_external_diagonal_weight", output)
        self.assertIn("ala_element_vanka_last_build_cycle", construct)
        self.assertIn("factor_age>=0", construct)
        self.assertIn("build_ala_element_vanka_factors(E);", construct)
        self.assertIn("action=reuse", construct)
        self.assertIn("scope=preconditioner_only", construct)
        self.assertIn("current_operator=reassembled", construct)
        self.assertIn("build_seconds_max=%e", construct)
        for source in (instructions, properties):
            self.assertIn(
                "ala_element_vanka_rebuild_interval must be in [1,100]",
                source,
            )
            self.assertIn(
                "effective element Vanka pressure damping must be in (0,1]",
                source,
            )
            self.assertIn(
                "ala_element_vanka_damping *",
                source,
            )
            self.assertIn(
                "ala_element_vanka_external_diagonal_weight must be in [0,1]",
                source,
            )

    def test_stage10_continuity_audit_localizes_linf_and_rejects_fake_div_v(
        self,
    ) -> None:
        drive = (GLOBAL_ROOT / "lib" / "Drive_solvers.c").read_text(
            encoding="utf-8"
        )
        start = drive.index("static void write_ala_residual(")
        residual = drive[start:drive.index(
            "static void write_ala_unaugmented_momentum_residual(", start
        )]
        self.assertIn("MPI_DOUBLE_INT", residual)
        self.assertIn("MPI_MAXLOC", residual)
        self.assertIn("MPI_Bcast(global_linf_location", residual)
        for metric in (
            "Linf_signed",
            "Linf_owner_rank",
            "Linf_owner_cap",
            "Linf_owner_element",
            "Linf_global_radial_element",
            "Linf_depth_km",
        ):
            self.assertIn(metric, residual)
        self.assertIn(
            'strcmp(E->control.ala_outer_solver,"coupled_fgmres") == 0',
            residual,
        )
        self.assertIn('"solver_div_v", "N/A"', residual)
        self.assertIn('"solver_div_v_status"', residual)
        self.assertIn('"not_applicable"', residual)
        power_row = drive[
            drive.index("static void print_power_row("):
            drive.index("static void write_mechanical_power(")
        ]
        self.assertIn("isfinite(value.maximum)", power_row)
        self.assertIn('"N/A", "N/A"', power_row)

    def test_stage6c_failure_path_audits_physical_scales_and_momentum(
        self,
    ) -> None:
        stokes = (
            GLOBAL_ROOT / "lib" / "Stokes_flow_Incomp.c"
        ).read_text(encoding="utf-8")
        failure = stokes.index("Strict ALA PCG failed physical continuity")
        final_audit = stokes.index("ALA PCG momentum audit final")
        self.assertLess(final_audit, failure)
        self.assertIn("mass_norm=%e", stokes)
        self.assertIn("Q=%e Q_relative=%e", stokes)
        self.assertIn(
            "strict_ala_momentum_residual_audit(\n"
            "        E,V,P,F,E->u1",
            stokes,
        )

    def test_total_temperature_is_the_only_temperature_state(self) -> None:
        definitions = (GLOBAL_ROOT / "lib" / "global_defs.h").read_text()
        instructions = (GLOBAL_ROOT / "lib" / "Instructions.c").read_text()
        initial = (GLOBAL_ROOT / "lib" / "Initial_temperature.c").read_text()
        advection = (GLOBAL_ROOT / "lib" / "Advection_diffusion.c").read_text()
        self.assertIn("double *T[NCS],*Tdot[NCS],*buoyancy[NCS];", definitions)
        self.assertNotIn("DataT", definitions)
        self.assertNotIn("DataT", instructions)
        self.assertNotIn("DataT", initial)
        self.assertIn("E->T", advection)
        self.assertIn("E->assim_delta_T", advection)

    def test_strict_rheology_uses_refstate_temperature_profile(self) -> None:
        viscosity = (
            GLOBAL_ROOT / "lib" / "Viscosity_structures.c"
        ).read_text(encoding="utf-8")
        helper_start = viscosity.index(
            "static double strict_rheology_reference_temperature(",
            viscosity.index("void visc_from_mat("),
        )
        helper_end = viscosity.index("void visc_from_T(", helper_start)
        helper = viscosity[helper_start:helper_end]
        case3_start = viscosity.index("case 3:")
        case3_end = viscosity.index("case 4:", case3_start)
        case3 = viscosity[case3_start:case3_end]

        self.assertIn("E->refstate.Tref[radial_node]", helper)
        self.assertIn("E->N.vpt[GNVINDEX(a,gp)]", helper)
        self.assertIn("E->control.ala_pressure_buoyancy", case3)
        self.assertIn(
            "strict_rheology_reference_temperature(E,m,i,jj)", case3
        )

    def test_kc_tracks_the_explicit_primordial_flavor(self) -> None:
        advection = (GLOBAL_ROOT / "lib" / "Advection_diffusion.c").read_text()
        composition = (GLOBAL_ROOT / "lib" / "Composition_related.c").read_text()
        material = (GLOBAL_ROOT / "lib" / "Material_properties.c").read_text()
        properties = (GLOBAL_ROOT / "module" / "setProperties.c").read_text()
        inventory = (
            GLOBAL_ROOT
            / "CitcomS/Components/Advection_diffusion/Advection_diffusion.py"
        ).read_text()
        profile = (GLOBAL_ROOT / "lib" / "Profile_output.c").read_text()

        self.assertIn('input_int("kC_primordial_flavor",', advection)
        self.assertIn(
            'getIntProperty(properties, "kC_primordial_flavor",', properties
        )
        self.assertIn(
            'prop.int("kC_primordial_flavor", default=1)', inventory
        )
        self.assertIn("E->control.kC_primordial_flavor - 1", material)
        self.assertNotIn("E->composition.ncomp-1", material)
        self.assertIn("Active kC requires tracers, composition", composition)
        self.assertIn("kC_primordial_flavor=%d is invalid", composition)
        self.assertIn('{"kC", "1", ELEMENT_PROFILE', profile)

    def test_layered_composition_disables_legacy_flavor_reclassification(
        self,
    ) -> None:
        tracer = (GLOBAL_ROOT / "lib" / "Tracer_setup.c").read_text()
        definitions = (GLOBAL_ROOT / "lib" / "tracer_defs.h").read_text()
        properties = (GLOBAL_ROOT / "module" / "setProperties.c").read_text()
        inventory = (GLOBAL_ROOT / "CitcomS/Components/Tracer.py").read_text()

        self.assertIn("int reclassify_flavors;", definitions)
        self.assertIn('input_boolean("tracer_reclassify_flavors",', tracer)
        self.assertIn("if(!E->trace.reclassify_flavors)", tracer)
        self.assertIn(
            'getIntProperty(properties, "tracer_reclassify_flavors",',
            properties,
        )
        self.assertIn(
            '"tracer_reclassify_flavors", default=True', inventory
        )

    def test_runtime_audits_are_absent(self) -> None:
        material = (GLOBAL_ROOT / "lib" / "Material_properties.c").read_text()
        initial = (GLOBAL_ROOT / "lib" / "Initial_temperature.c").read_text()
        instructions = (GLOBAL_ROOT / "lib" / "Instructions.c").read_text()
        combined = material + initial + instructions
        self.assertNotIn("validate_strict_reference_state", combined)
        self.assertNotIn("audit_strict_temperature_reference", combined)
        self.assertNotIn("strict_reference_audit", combined)

    def test_phase_temperature_authorities_are_separate(self) -> None:
        phase = (GLOBAL_ROOT / "lib" / "Phase_change.c").read_text()
        self.assertRegex(
            phase,
            re.compile(
                r"phase_change_reference_fraction.*?"
                r"E->refstate\.Tref\[nz\]",
                re.S,
            ),
        )
        self.assertIn(
            "phase->clapeyron * (E->T[m][i] - phase->transT)", phase
        )
        self.assertIn("phase->Ra * (B[m][i] - Xref)", phase)


if __name__ == "__main__":
    unittest.main(verbosity=2)
