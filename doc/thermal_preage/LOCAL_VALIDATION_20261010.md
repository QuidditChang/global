# Local validation, 2026-10-10

Executed again on the user's macOS x86_64 host, against the thermal-preage
patch based on `bbb22a924618f358a7c7b7b1b9b65a5bceafedd0`. This is separate
from the supplied bundle's cloud validation report.

- Compiler: Apple clang 21.0.0; MPI: Open MPI 5.0.9.
- Python 3.10, NumPy 2.2.6, SciPy 1.15.3.
- All real C solver sources compiled and linked into the initialization-only
  driver. Apple ld lacks GNU `--wrap`; the test builder now uses per-file
  symbol definitions on Darwin, without changing production sources.
- `MPIEXEC_FLAGS=--oversubscribe` regression: **40/40 PASS** (36 MPI cases,
  four convergence summaries). See `LOCAL_VALIDATION_20261010.json`.
- TA regression: 3/3 PASS; phase-energy geometry: 3/3 PASS; Qvis: 6/6 PASS.
- New module compiled with `-Wall -Wextra -Werror` plus
  `-Wno-deprecated-non-prototype` for the existing `get_global_shape_fn`
  declaration used throughout the solver.
- `git diff --check` passed.

Commands executed from the solver worktree:

```sh
python3 tests/thermal_preage/build_validation.py \
  --build-dir /tmp/citcoms-preage-local-20261010 --jobs 4
MPIEXEC_FLAGS=--oversubscribe python3 tests/thermal_preage/run_validation.py \
  --build-dir /tmp/citcoms-preage-local-20261010 \
  --output /tmp/citcoms-preage-validation-20261010-mpi
python3 tests/thermal_assimilation/test_ta.py
python3 tests/phase_energy/test_phase_energy_geometry.py
python3 tests/qvis/test_qvis.py
```

Open MPI required local communication sockets unavailable inside the sandbox;
the successful run used the approved local execution outside that sandbox.
The earlier sandbox attempt failed before launching solver ranks.

The production Python2/Pyre application, Intel MPI build, geological forcing,
and 384-rank target mesh were not executed here. Python 3's parser also rejects
pre-existing mixed indentation in `Param.py`; this is not a Python2/Pyre runtime
test. No HPC job was submitted. Small-mesh convergence does not establish
target-mesh accuracy, initialization time, or memory cost.
