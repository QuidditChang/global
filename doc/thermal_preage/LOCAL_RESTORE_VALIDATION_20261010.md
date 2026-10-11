# Restoration on the ASCII baseline: local validation

Restored the opt-in thermal-preage implementation on
`22399749866f452b3c06960e3d00eea5333a37f0`, preserving its ASCII output changes.
No production thermal age was selected. Default mode remains off; the runs
production cfg and 400-core/PostProc launcher are unchanged.

Actual local execution on macOS x86_64, Apple clang 21.0.0, Open MPI 5.0.9,
Python 3.10, NumPy 2.2.6 and SciPy 1.15.3:

- All real solver C sources compiled and linked into the initialization-only
  driver (Apple ld uses per-file symbol definitions instead of GNU `--wrap`).
- **41/41 PASS** in the 12/24-rank MPI suite: 37 MPI cases and four convergence
  summaries. Full result: `LOCAL_RESTORE_VALIDATION_20261010.json`.
- New ASCII test calls the real field writer after initialization, before
  Stokes/time marching. Its final temperatures exactly match the gzip-mode
  initialization and agree with the saved dimensional ASCII temperatures to
  the writer's precision.
- TA regression 3/3, phase geometry 3/3 and Qvis 6/6 PASS.
- The build helper avoids `Path.is_relative_to` so it also supports Python 3.8
  in the cluster's Anaconda environment. Python 3.8 itself was not run here.

Executed commands from the solver worktree:

```sh
python3 tests/thermal_preage/build_validation.py \
  --build-dir /tmp/citcoms-preage-restored-20261010 --jobs 4
MPIEXEC_FLAGS=--oversubscribe python3 tests/thermal_preage/run_validation.py \
  --build-dir /tmp/citcoms-preage-restored-20261010 \
  --output /tmp/citcoms-preage-restored-validation-20261010
python3 tests/thermal_assimilation/test_ta.py
python3 tests/phase_energy/test_phase_energy_geometry.py
python3 tests/qvis/test_qvis.py
```

The MPI suite needed approved local communication sockets outside the sandbox.
This is new execution evidence, not a reused cloud test report.

The separate runs tools have local tests for initialization-only control flow,
clock/restart guards, independent age/step-cap configuration, preserved
production files, input checksums, directory reuse rejection, and comparison
rejection for changed composition or physics. Those use mocks for Pyre; they do
not certify the cluster Python2/Pyre application or Intel MPI launch.

No HPC job was submitted. The full mesh, geological forcing, target-mesh cost,
memory use, physical age selection and 5/2.5/1.25 Ma target-mesh convergence
remain to be evaluated on the cluster.
