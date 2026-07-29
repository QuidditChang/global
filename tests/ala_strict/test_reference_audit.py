#!/usr/bin/env python3
"""Regression checks for the read-only strict reference startup audit."""

from __future__ import annotations

import re
import unittest
from pathlib import Path

import numpy as np


GLOBAL_ROOT = Path(__file__).resolve().parents[2]
PROJECT_ROOT = Path(__file__).resolve().parents[4]
RUNS_ROOT = PROJECT_ROOT / "runs"
MATERIAL_SOURCE = GLOBAL_ROOT / "lib" / "Material_properties.c"
INSTRUCTIONS_SOURCE = GLOBAL_ROOT / "lib" / "Instructions.c"


class StrictReferenceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.legacy = np.loadtxt(RUNS_ROOT / "refstate_ALA.txt", comments="#")
        cls.strict = np.loadtxt(
            RUNS_ROOT / "refstate_ALA_strict.txt", comments="#"
        )
        cls.radius = np.loadtxt(
            RUNS_ROOT / "GLB.coor.global.dat", skiprows=1, usecols=1
        )

    def test_legacy_and_strict_schemas_remain_supported(self) -> None:
        self.assertEqual(self.legacy.shape, (len(self.radius), 7))
        self.assertEqual(self.strict.shape, (len(self.radius), 8))

    def test_gamma_eff_closes_from_serialized_runtime_columns(self) -> None:
        gamma_check = (
            self.strict[:, 3]
            * self.strict[:, 1]
            / (self.strict[:, 4] * self.strict[:, 5])
        )
        self.assertTrue(np.array_equal(gamma_check, self.strict[:, 6]))

    def test_strict_ks_is_positive_and_legacy_ks_is_disabled(self) -> None:
        self.assertTrue(np.all(self.strict[:, 7] > 0.0))
        source = MATERIAL_SOURCE.read_text(encoding="utf-8")
        self.assertIn(
            "E->refstate.Ks[i] = columns == 8 ? values[7] : 0.0;",
            source,
        )

    def test_rho_beta_audit_is_diagnostic_only(self) -> None:
        source = MATERIAL_SOURCE.read_text(encoding="utf-8")
        definition = source.rindex(
            "static void validate_strict_reference_state("
        )
        end = source.index("double conductivity_depth_factor", definition)
        audit = source[definition:end]
        self.assertIn("beta_secant =", audit)
        self.assertIn("beta_input = 0.5 *", audit)
        self.assertIsNone(
            re.search(
                r"E->refstate\.[A-Za-z_]+\[[^\]]+\]\s*=",
                audit,
            )
        )
        self.assertNotIn("parallel_process_termination", audit)

    def test_audit_switch_defaults_on_and_can_be_disabled(self) -> None:
        source = INSTRUCTIONS_SOURCE.read_text(encoding="utf-8")
        self.assertIn('input_boolean("strict_reference_audit"', source)
        self.assertIn("&(E->control.strict_reference_audit),\"on\",m", source)


if __name__ == "__main__":
    unittest.main(verbosity=2)
