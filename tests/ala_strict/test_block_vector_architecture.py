#!/usr/bin/env python3
"""Static contracts for strict-ALA mixed velocity-pressure vectors."""

from __future__ import annotations

import unittest
from pathlib import Path


GLOBAL_ROOT = Path(__file__).resolve().parents[2]
LIB_ROOT = GLOBAL_ROOT / "lib"


class BlockVectorArchitectureTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.source = (LIB_ROOT / "ALA_block_vector.c").read_text()
        cls.header = (LIB_ROOT / "ala_block_vector.h").read_text()

    def test_vector_owns_velocity_pressure_and_level(self) -> None:
        self.assertIn("double **velocity;", self.header)
        self.assertIn("double **pressure;", self.header)
        self.assertIn("int level;", self.header)
        self.assertIn("ala_block_vector_create", self.header)
        self.assertIn("ala_block_vector_destroy", self.header)

    def test_storage_preserves_legacy_sentinels(self) -> None:
        self.assertIn("neq+1", self.source)
        self.assertIn("npno+1", self.source)
        self.assertIn("destination->velocity[m][neq]=0.0;", self.source)
        self.assertIn("destination->pressure[m][0]=0.0;", self.source)
        self.assertIn("vector->velocity[m][neq]=0.0;", self.source)
        self.assertIn("vector->pressure[m][0]=0.0;", self.source)

    def test_algebra_rejects_cross_level_vectors(self) -> None:
        self.assertIn("vector->level!=level", self.source)
        self.assertIn(
            "Mismatched level in strict-ALA block vector operation",
            self.source,
        )

    def test_metric_uses_one_mpi_reduction_and_pressure_dual_mass(self) -> None:
        self.assertIn("E->parallel.Skip_neq[level][m]", self.source)
        self.assertIn("E->parallel.Skip_id[level][m][i]", self.source)
        self.assertIn("/E->ECO[level][m][i].area", self.source)
        self.assertIn("Nonpositive pressure volume", self.source)
        self.assertIn(
            "MPI_Allreduce(local,global,2,MPI_DOUBLE,MPI_SUM", self.source
        )
        self.assertIn("velocity_weight*component_dot[0]", self.source)
        self.assertIn("pressure_weight*component_dot[1]", self.source)
        self.assertIn("metric weights must be positive", self.source)
        self.assertIn("ala_block_vector_component_norms", self.source)

    def test_legacy_makefile_template_links_new_module(self) -> None:
        makefile = (LIB_ROOT / "Makefile.in").read_text()
        self.assertIn("libCitcomS_a-ALA_block_vector.$(OBJEXT)", makefile)
        self.assertIn("ALA_block_vector.lo", makefile)
        self.assertIn("libCitcomS_a-ALA_block_vector.o:", makefile)

    def test_strict_vector_does_not_apply_legacy_pressure_projection(self) -> None:
        combined = self.header + self.source
        self.assertIn("C^T generally", self.header)
        self.assertNotIn("remove_pressure_mean", combined)
        self.assertNotIn("project_pressure_nullspace", combined)

    def test_module_is_selected_only_by_explicit_coupled_solver(self) -> None:
        stokes = (LIB_ROOT / "Stokes_flow_Incomp.c").read_text()
        drive = (LIB_ROOT / "Drive_solvers.c").read_text()
        self.assertIn("solve_ala_coupled_fgmres_core", stokes)
        self.assertIn("ala_block_vector_create", stokes)
        self.assertIn('"coupled_fgmres"', stokes)
        self.assertNotIn("ala_block_vector_create", drive)


if __name__ == "__main__":
    unittest.main(verbosity=2)
