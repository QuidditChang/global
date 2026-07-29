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


class StrictProductionArchitectureTest(unittest.TestCase):
    def test_cfg_changes_only_reference_state_input(self) -> None:
        legacy = _active_cfg_lines(RUNS_ROOT / "cmbhf_ALA.cfg")
        strict = _active_cfg_lines(RUNS_ROOT / "cmbhf_ALA_strict.cfg")
        self.assertEqual(len(legacy), len(strict))
        differences = [
            (old, new)
            for old, new in zip(legacy, strict)
            if old != new
        ]
        self.assertEqual(
            differences,
            [
                (
                    "refstate_file                = refstate_ALA.txt",
                    "refstate_file                = refstate_ALA_strict.txt",
                )
            ],
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

    def test_reader_retains_legacy_seven_column_compatibility(self) -> None:
        legacy = np.loadtxt(RUNS_ROOT / "refstate_ALA.txt", comments="#")
        material = (GLOBAL_ROOT / "lib" / "Material_properties.c").read_text()
        self.assertEqual(legacy.shape[1], 7)
        self.assertIn("columns != 7 && columns != 8", material)
        self.assertIn(
            "E->refstate.Ks[i] = columns == 8 ? values[7] : 0.0;",
            material,
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
