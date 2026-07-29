#!/usr/bin/env python3
"""Static regression checks for the strict temperature semantic audit."""

from __future__ import annotations

import re
import unittest
from pathlib import Path


GLOBAL_ROOT = Path(__file__).resolve().parents[2]
SOURCE = GLOBAL_ROOT / "lib" / "Initial_temperature.c"


class StrictTemperatureReferenceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        source = SOURCE.read_text(encoding="utf-8")
        start = source.index(
            "static void audit_strict_temperature_reference("
        )
        start = source.index(
            "static void audit_strict_temperature_reference(", start + 1
        )
        end = source.index("\n\nvoid debug_tic", start)
        cls.audit = source[start:end]
        cls.source = source

    def test_audit_runs_after_initial_temperature_construction(self) -> None:
        function = self.source[
            self.source.index("void convection_initial_temperature(") :
            self.source.index("\n\n/* Read-only semantic audit",)
        ]
        self.assertGreater(
            function.index("audit_strict_temperature_reference(E);"),
            function.index("(E->solver.construct_tic_from_input)(E);"),
        )

    def test_delta_temperature_is_read_only_difference(self) -> None:
        self.assertIn(
            "E->T[cap][node]-E->refstate.temperature[k]",
            self.audit,
        )
        self.assertIsNone(
            re.search(r"E->T\[[^\n;]+\]\s*=", self.audit)
        )
        self.assertIsNone(
            re.search(
                r"E->refstate\.[A-Za-z_]+\[[^\]]+\]\s*=",
                self.audit,
            )
        )

    def test_adiabatic_check_uses_dimensional_temperature(self) -> None:
        self.assertIn("log(temperature1)-log(temperature0)", self.audit)
        self.assertIn("-E->control.disptn_number", self.audit)
        self.assertIn("relative RMS:", self.audit)
        self.assertIn("relative MAX:", self.audit)

    def test_semantics_and_phase_paths_are_reported(self) -> None:
        self.assertIn("Katsura 2022 Table S5 temperature", self.audit)
        self.assertIn(
            "four phase branches; endpoint-constrained piecewise cubic fit",
            self.audit,
        )
        self.assertIn(
            "Xref uses E->refstate.temperature",
            self.audit,
        )
        self.assertIn("dynamic X uses E->T", self.audit)

    def test_restart_is_not_mislabeled_as_initial_anomaly(self) -> None:
        self.assertIn("E->convection.tic_method == -1", self.audit)
        self.assertIn("E->T is a restart state, skip", self.audit)


if __name__ == "__main__":
    unittest.main(verbosity=2)
