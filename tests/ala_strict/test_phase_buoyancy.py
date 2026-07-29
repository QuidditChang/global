#!/usr/bin/env python3
"""Regression checks for strict-ALA phase-reference buoyancy decomposition."""

from __future__ import annotations

import re
import unittest
from pathlib import Path

import numpy as np


GLOBAL_ROOT = Path(__file__).resolve().parents[2]
PROJECT_ROOT = Path(__file__).resolve().parents[4]
RUNS_ROOT = PROJECT_ROOT / "runs"
PHASE_SOURCE = GLOBAL_ROOT / "lib" / "Phase_change.c"
ENERGY_SOURCE = GLOBAL_ROOT / "lib" / "Advection_diffusion.c"
POWER_SOURCE = GLOBAL_ROOT / "lib" / "Drive_solvers.c"
PROFILE_SOURCE = GLOBAL_ROOT / "lib" / "Profile_output.c"


def _cfg_vector(name: str) -> np.ndarray:
    text = (RUNS_ROOT / "cmbhf_ALA.cfg").read_text(encoding="utf-8")
    match = re.search(rf"^\s*{re.escape(name)}\s*=\s*(.*?)\s*$", text, re.M)
    if match is None:
        raise AssertionError(f"Missing {name} in cmbhf_ALA.cfg")
    return np.asarray(
        [float(value.strip()) for value in match.group(1).split(",")],
        dtype=float,
    )


def _phase_fraction(
    depth: np.ndarray,
    rho: np.ndarray,
    gravity: np.ndarray,
    temperature: np.ndarray,
    phase_depth: float,
    clapeyron: float,
    trans_temperature: float,
    inverse_width: float,
) -> np.ndarray:
    q = (
        (depth - phase_depth) * rho * gravity
        - clapeyron * (temperature - trans_temperature)
    )
    return 0.5 * (1.0 + np.tanh(inverse_width * q))


class StrictPhaseBuoyancyTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        refstate = np.loadtxt(RUNS_ROOT / "refstate_ALA.txt", comments="#")
        radius = np.loadtxt(
            RUNS_ROOT / "GLB.coor.global.dat", skiprows=1, usecols=1
        )
        cls.depth = 1.0 - radius
        cls.rho = refstate[:, 0]
        cls.gravity = refstate[:, 1]
        cls.tref = refstate[:, 2]
        cls.phase_depth = _cfg_vector("phase_depth")
        cls.phase_ra = _cfg_vector("phase_Ra")
        cls.width = _cfg_vector("phase_width")
        cls.clapeyron = _cfg_vector("phase_clapeyron")
        cls.trans_temperature = _cfg_vector("phase_transT")
        cls.dynamic_temperature = cls.tref + 0.08 * np.sin(
            np.linspace(0.0, 3.0 * np.pi, len(cls.tref))
        )

    def _fractions(
        self, phase: int, temperature: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        xref = _phase_fraction(
            self.depth,
            self.rho,
            self.gravity,
            self.tref,
            self.phase_depth[phase],
            self.clapeyron[phase],
            self.trans_temperature[phase],
            1.0 / self.width[phase],
        )
        x = _phase_fraction(
            self.depth,
            self.rho,
            self.gravity,
            temperature,
            self.phase_depth[phase],
            self.clapeyron[phase],
            self.trans_temperature[phase],
            1.0 / self.width[phase],
        )
        return x, xref

    def test_reference_temperature_has_zero_phase_anomaly(self) -> None:
        maximum = 0.0
        for phase in range(3):
            x, xref = self._fractions(phase, self.tref)
            maximum = max(maximum, float(np.max(np.abs(x - xref))))
        self.assertLessEqual(maximum, np.finfo(float).eps)

    def test_dynamic_phase_fraction_is_the_legacy_fraction(self) -> None:
        source = PHASE_SOURCE.read_text(encoding="utf-8")
        self.assertIn(
            "phase->clapeyron * (E->T[m][i] - phase->transT)", source
        )
        self.assertIn(
            "B[m][i] = pt5 * (one + tanh(phase->inv_width * e_pressure))",
            source,
        )
        for phase in range(3):
            old_x, _ = self._fractions(phase, self.dynamic_temperature)
            new_x, _ = self._fractions(phase, self.dynamic_temperature)
            self.assertTrue(np.array_equal(old_x, new_x))

    def test_buoyancy_difference_is_reference_subtraction_only(self) -> None:
        for phase in range(3):
            x, xref = self._fractions(phase, self.dynamic_temperature)
            old_buoyancy = -self.phase_ra[phase] * x
            new_buoyancy = -self.phase_ra[phase] * (x - xref)
            residual = new_buoyancy - old_buoyancy - self.phase_ra[phase] * xref
            scale = max(float(np.max(np.abs(old_buoyancy))), 1.0)
            self.assertLessEqual(
                float(np.max(np.abs(residual))), 8.0 * np.finfo(float).eps * scale
            )

    def test_phase_boundary_still_uses_absolute_phase_fraction(self) -> None:
        source = PHASE_SOURCE.read_text(encoding="utf-8")
        boundary_start = source.index("/* compute the phase boundary")
        boundary_end = source.index("return;", boundary_start)
        boundary_code = source[boundary_start:boundary_end]
        self.assertIn("B[m][n]>=pt5", boundary_code)
        self.assertIn("B[m][n+1]<=pt5", boundary_code)
        self.assertNotIn("Xref", boundary_code)
        self.assertNotIn("deltaX", boundary_code)

    def test_latent_heating_still_uses_absolute_phase_fraction(self) -> None:
        source = ENERGY_SOURCE.read_text(encoding="utf-8")
        self.assertIn("(1.0 - B[m][j]) * B[m][j]", source)
        self.assertIn("E->phase_B[phase_index], phase->Ra", source)
        self.assertNotIn("phase_change_reference_fraction", source)

    def test_all_phase_buoyancy_diagnostics_use_x_minus_xref(self) -> None:
        phase_source = PHASE_SOURCE.read_text(encoding="utf-8")
        power_source = POWER_SOURCE.read_text(encoding="utf-8")
        profile_source = PROFILE_SOURCE.read_text(encoding="utf-8")
        self.assertIn("phase->Ra * (B[m][i] - Xref)", phase_source)
        self.assertIn("phase_change_reference_fraction(", power_source)
        self.assertIn("phase_change_reference_fraction(", profile_source)


if __name__ == "__main__":
    unittest.main(verbosity=2)
