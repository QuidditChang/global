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
    def test_cfg_changes_only_strict_inputs_and_stage6c_metric(self) -> None:
        legacy = _active_cfg_lines(RUNS_ROOT / "cmbhf_ALA.cfg")
        strict = _active_cfg_lines(RUNS_ROOT / "cmbhf_ALA_strict.cfg")
        strict_keys = (
            "refstate_file",
            "ala_beta_element_source",
            "ala_beta_interval_file",
            "ala_shallow_patch_velocity_solver",
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
            r"(?m)^\s*ala_shallow_patch_velocity_solver\s*="
            r"\s*node_block\s*$",
        )

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
            'getStringProperty(properties, '
            '"ala_shallow_patch_velocity_solver"',
            properties,
        )

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
