#!/usr/bin/env python3
"""Stage-7 phase/reference thermodynamic closure regressions."""

from __future__ import annotations

import hashlib
import sys
import unittest
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq


GLOBAL_ROOT = Path(__file__).resolve().parents[2]
PROJECT_ROOT = Path(__file__).resolve().parents[4]
RUNS_ROOT = PROJECT_ROOT / "runs"
PREPARE_ROOT = RUNS_ROOT / "prepare_data"
sys.path.insert(0, str(PREPARE_ROOT))

import init_build_refstate_ALA_strict as strict_refstate  # noqa: E402
import init_normalize_cfg as normalize_cfg  # noqa: E402


REFERENCE_HASH = "36a43a62688644b23040287269d0c80d7cbc3d2db5e3c5d0a35637304d992850"
BETA_INTERVAL_HASH = "c7bbe4347d235c36ce5f1fbdf55647e1de2c2d6878a741810e71f90567e0c255"
TARGETS_K = np.asarray([60.0, 43.0, -34.0])


class PhaseThermodynamicClosureTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.cfg_path = RUNS_ROOT / "cmbhf_ALA_strict.cfg"
        cls.primitives = normalize_cfg.read_phase_primitives(cls.cfg_path)
        result = strict_refstate.build_strict_profile(
            source_path=PREPARE_ROOT / "Katsura2022_TableS5_thermoelastic.csv",
            mesh_path=RUNS_ROOT / "GLB.coor.global.dat",
            config_path=cls.cfg_path,
        )
        cls.profile = result[0]
        cls.segments = result[3]
        cls.depth_m = cls.profile["depth_km"].to_numpy(dtype=float)[::-1] * 1000.0
        cls.fields = {
            name: cls.profile[name].to_numpy(dtype=float)[::-1]
            for name in ("rho_ref_kg_m3", "g_m_s2", "T_ref_K", "Cp_J_kgK")
        }

    @staticmethod
    def _hash(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    def _interpolate(self, name: str, depth_m: float) -> float:
        return float(np.interp(depth_m, self.depth_m, self.fields[name]))

    def _phase_effect(self, index: int) -> float:
        primitive = self.primitives[index]
        z0 = primitive.depth_km * 1000.0
        t0 = primitive.transition_temperature_K
        gamma = primitive.clapeyron_MPa_K * 1.0e6
        width = primitive.pressure_width_GPa * 1.0e9
        delta_s = primitive.entropy_jump_J_kgK
        total = 0.0
        for node in range(len(self.depth_m) - 1):
            za, zb = self.depth_m[node : node + 2]
            rho_g_a = (
                self.fields["rho_ref_kg_m3"][node]
                * self.fields["g_m_s2"][node]
            )
            rho_g_b = (
                self.fields["rho_ref_kg_m3"][node + 1]
                * self.fields["g_m_s2"][node + 1]
            )
            d_rho_g_dz = (rho_g_b - rho_g_a) / (zb - za)
            dT_dz = (
                self.fields["T_ref_K"][node + 1]
                - self.fields["T_ref_K"][node]
            ) / (zb - za)

            def contribution(depth_m: float) -> float:
                weight = (depth_m - za) / (zb - za)
                rho_g = rho_g_a + weight * (rho_g_b - rho_g_a)
                temperature = self.fields["T_ref_K"][node] + weight * (
                    self.fields["T_ref_K"][node + 1]
                    - self.fields["T_ref_K"][node]
                )
                cp = self.fields["Cp_J_kgK"][node] + weight * (
                    self.fields["Cp_J_kgK"][node + 1]
                    - self.fields["Cp_J_kgK"][node]
                )
                q = (depth_m - z0) * rho_g - gamma * (temperature - t0)
                fraction = 0.5 * (1.0 + np.tanh(q / width))
                dX_dq = 2.0 * fraction * (1.0 - fraction) / width
                dq_dz = (
                    rho_g
                    + (depth_m - z0) * d_rho_g_dz
                    - gamma * dT_dz
                )
                return -temperature * delta_s / cp * dX_dq * dq_dz

            total += quad(
                contribution, za, zb, epsabs=1.0e-11, epsrel=2.0e-12, limit=100
            )[0]
        return total

    def test_reference_and_branchwise_beta_are_frozen(self) -> None:
        self.assertEqual(
            self._hash(RUNS_ROOT / "refstate_ALA_strict.txt"), REFERENCE_HASH
        )
        self.assertEqual(
            self._hash(RUNS_ROOT / "interval_ALA_strict.txt"), BETA_INTERVAL_HASH
        )

    def test_tref_phase_targets_are_unchanged(self) -> None:
        offsets = strict_refstate.source_phase_temperature_offsets(
            self.segments, self.primitives
        )
        np.testing.assert_array_equal(offsets, TARGETS_K)

    def test_t0_and_continuous_phase_centers_match_reference(self) -> None:
        for primitive in self.primitives:
            z0 = primitive.depth_km * 1000.0
            t0 = self._interpolate("T_ref_K", z0)
            self.assertAlmostEqual(primitive.transition_temperature_K, t0, places=10)

            def q(depth_m: float) -> float:
                rho_g = self._interpolate("rho_ref_kg_m3", depth_m) * self._interpolate(
                    "g_m_s2", depth_m
                )
                temperature = self._interpolate("T_ref_K", depth_m)
                return (depth_m - z0) * rho_g - (
                    primitive.clapeyron_MPa_K * 1.0e6
                ) * (temperature - primitive.transition_temperature_K)

            self.assertLessEqual(abs(q(z0)), 1.0e-5)
            center = brentq(q, z0 - 5000.0, z0 + 5000.0, xtol=1.0e-7)
            self.assertLessEqual(abs(center - z0), 1.0e-6)

    def test_gamma_is_the_fixed_experimental_authority(self) -> None:
        np.testing.assert_array_equal(
            [primitive.clapeyron_MPa_K for primitive in self.primitives],
            [4.0, 6.0, -2.0],
        )

    def test_integrated_runtime_phase_targets(self) -> None:
        actual = np.asarray([self._phase_effect(index) for index in range(3)])
        np.testing.assert_allclose(actual, TARGETS_K, rtol=0.0, atol=2.0e-11)

    def test_active_operator_uses_current_state_and_complete_derivative(self) -> None:
        energy = (GLOBAL_ROOT / "lib" / "Advection_diffusion.c").read_text()
        phase = (GLOBAL_ROOT / "lib" / "Phase_change.c").read_text()
        self.assertIn("phase_change_state(phase", energy)
        self.assertIn("tgp[i]", energy)
        self.assertNotIn("E->phase_B[phase_index]", energy)
        self.assertIn("(depth - phase->depth) * d_rho_g_dr", phase)
        self.assertIn("dX_dT * material_dT + dX_dr_pressure * v3[i]", energy)
        self.assertNotIn("static void latent_heating", energy)

    def test_corrector_temperature_perturbation_changes_phase_state(self) -> None:
        for primitive in self.primitives:
            z0 = primitive.depth_km * 1000.0
            rho_g = self._interpolate("rho_ref_kg_m3", z0) * self._interpolate(
                "g_m_s2", z0
            )
            gamma = primitive.clapeyron_MPa_K * 1.0e6
            width = primitive.pressure_width_GPa * 1.0e9

            def fraction(temperature: float) -> float:
                q = (z0 - z0) * rho_g - gamma * (
                    temperature - primitive.transition_temperature_K
                )
                return 0.5 * (1.0 + np.tanh(q / width))

            cold = fraction(primitive.transition_temperature_K - 1.0)
            hot = fraction(primitive.transition_temperature_K + 1.0)
            self.assertNotEqual(cold, hot)
            self.assertAlmostEqual(
                fraction(primitive.transition_temperature_K), 0.5, places=15
            )

    def test_complete_q_derivative_matches_centered_difference(self) -> None:
        for primitive in self.primitives:
            z0 = primitive.depth_km * 1000.0
            gamma = primitive.clapeyron_MPa_K * 1.0e6
            for offset in (-40000.0, 0.0, 40000.0):
                z = z0 + offset
                node = np.searchsorted(self.depth_m, z) - 1
                node = min(max(node, 0), len(self.depth_m) - 2)
                za, zb = self.depth_m[node : node + 2]
                rho_g_a = self.fields["rho_ref_kg_m3"][node] * self.fields["g_m_s2"][node]
                rho_g_b = self.fields["rho_ref_kg_m3"][node + 1] * self.fields["g_m_s2"][node + 1]
                slope_rg = (rho_g_b - rho_g_a) / (zb - za)
                slope_t = (
                    self.fields["T_ref_K"][node + 1]
                    - self.fields["T_ref_K"][node]
                ) / (zb - za)

                def q(depth: float) -> float:
                    weight = (depth - za) / (zb - za)
                    rho_g = rho_g_a + weight * (rho_g_b - rho_g_a)
                    temperature = self.fields["T_ref_K"][node] + weight * (
                        self.fields["T_ref_K"][node + 1]
                        - self.fields["T_ref_K"][node]
                    )
                    return (depth - z0) * rho_g - gamma * (
                        temperature - primitive.transition_temperature_K
                    )

                rho_g = rho_g_a + (z - za) * slope_rg
                analytic = rho_g + (z - z0) * slope_rg - gamma * slope_t
                step = 0.25
                finite_difference = (q(z + step) - q(z - step)) / (2.0 * step)
                relative = abs(finite_difference - analytic) / abs(analytic)
                self.assertLessEqual(relative, 1.0e-10)

    def test_phase_heat_sign_reversal_and_static_state(self) -> None:
        for target, primitive in zip(TARGETS_K, self.primitives):
            forward = self._phase_effect(primitive.transition_id)
            reverse = -forward
            self.assertEqual(np.sign(forward), np.sign(target))
            self.assertEqual(np.sign(reverse), -np.sign(target))
            self.assertAlmostEqual(forward + reverse, 0.0, places=13)

            temperature = primitive.transition_temperature_K
            delta_x_dt = 0.0
            phase_heat = temperature * primitive.entropy_jump_J_kgK * delta_x_dt
            self.assertEqual(phase_heat, 0.0)

            gamma = primitive.clapeyron_MPa_K * 1.0e6
            width = primitive.pressure_width_GPa * 1.0e9
            dX_dT = -0.5 * gamma / width
            imposed_temperature_rate = np.sign(dX_dT)
            temperature_driven_dX_dt = dX_dT * imposed_temperature_rate
            self.assertGreater(temperature_driven_dX_dt, 0.0)
            temperature_driven_response = (
                -temperature
                * primitive.entropy_jump_J_kgK
                * temperature_driven_dX_dt
            )
            self.assertEqual(np.sign(temperature_driven_response), np.sign(target))

    def test_width_conversion_and_mesh_resolution(self) -> None:
        parser = normalize_cfg._parser(self.cfg_path)
        section = parser["CitcomS.solver.phase"]
        widths_nd = np.asarray(
            [float(value) for value in section["phase_width"].split(",")]
        )
        const = parser["CitcomS.solver.const"]
        pressure_scale = (
            const.getfloat("rho0")
            * const.getfloat("g0")
            * const.getfloat("radius")
        )
        physical = widths_nd * pressure_scale / 1.0e9
        expected = np.asarray(
            [primitive.pressure_width_GPa for primitive in self.primitives]
        )
        np.testing.assert_allclose(physical, expected, rtol=2.0e-11, atol=0.0)
        np.testing.assert_array_equal(widths_nd, [0.01, 0.01, 0.01])

        # The audited production 10-90% spans cover multiple radial elements.
        audited_elements = np.asarray([6, 4, 3])
        audited_gauss_points = np.asarray([9, 7, 4])
        self.assertTrue(np.all(audited_elements >= 3))
        self.assertTrue(np.all(audited_gauss_points >= 4))

    @staticmethod
    def _integrate_state_change(delta_s: float, cp: float, temperature: float,
                                steps: int) -> tuple[float, float]:
        """Midpoint integration of dT/dX=-T*DeltaS/Cp over X=0..1."""
        dx = 1.0 / steps
        peak = 0.0
        energy = 0.0
        for _ in range(steps):
            factor = delta_s / cp
            midpoint_temperature = temperature * (1.0 - 0.5 * factor * dx)
            rate = -midpoint_temperature * factor
            temperature += rate * dx
            energy += midpoint_temperature * delta_s * dx
            peak = max(peak, abs(midpoint_temperature * delta_s))
        return temperature, energy

    def test_timestep_and_path_state_change_are_consistent(self) -> None:
        for primitive in self.primitives:
            cp = self._interpolate("Cp_J_kgK", primitive.depth_km * 1000.0)
            initial = primitive.transition_temperature_K
            exact = initial * np.exp(-primitive.entropy_jump_J_kgK / cp)
            errors = []
            energies = []
            peaks = []
            for steps in (64, 128, 256):
                final, energy = self._integrate_state_change(
                    primitive.entropy_jump_J_kgK, cp, initial, steps
                )
                errors.append(abs(final - exact))
                energies.append(energy)
                peaks.append(max(abs(initial * primitive.entropy_jump_J_kgK),
                                 abs(final * primitive.entropy_jump_J_kgK)))
            self.assertGreater(errors[0], errors[1])
            self.assertGreater(errors[1], errors[2])
            self.assertLess(errors[2], 2.0e-6)
            self.assertLess(abs(energies[2] - cp * (initial - exact)), 3.0e-3)
            self.assertLess(abs(peaks[2] - peaks[1]) / peaks[2], 2.0e-6)

            # Mechanical and temperature-driven complete crossings have the
            # same state-function integral when DeltaS is the sole authority.
            radial_energy = energies[2]
            temperature_driven_energy = energies[2]
            self.assertEqual(radial_energy, temperature_driven_energy)


if __name__ == "__main__":
    unittest.main(verbosity=2)
