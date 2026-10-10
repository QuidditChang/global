# Conductive pre-age regression tests

These tests run the real full-sphere MPI solver initialization with small,
self-contained synthetic inputs. They do not need the geological reconstruction
or change production run configurations. The generated main exits immediately
after `initial_conditions`, before Stokes or any production timestep.

## Run

```sh
export MPICC=mpicc MPIEXEC=mpiexec
# OpenMPI on a machine with fewer than 24 slots also needs:
# export MPIEXEC_FLAGS=--oversubscribe
python tests/thermal_preage/build_validation.py --build-dir /tmp/preage-build
python tests/thermal_preage/run_validation.py --build-dir /tmp/preage-build \
  --output /tmp/preage-validation
```

Use a fresh output directory each time. `--prepare-only` creates sample 12-rank
and 24-rank fixtures without launching MPI. Python requires NumPy and SciPy.
Python 3.8 or newer is supported by the build and test scripts.
`MPICC` and `MPIEXEC` must refer to the same MPI implementation. The build links
all real C solver sources and the real parameter parser; it has no parser patch
or mock MPI layer. Build provenance is saved in `VALIDATION_BUILD.txt` and
`SOURCE_SHA256.json`.

On macOS, where Apple ld lacks GNU `--wrap`, the test build renames the
implementation and wrapper symbols with per-file compiler definitions. Solver
source files are not rewritten; production builds are unchanged.

## Coverage

- Default-off and explicit-off fields agree exactly; enabling a zero duration
  gives the same final field when the legacy bottom layer is disabled. Enabling
  pre-age also bypasses a deliberately nonzero legacy bottom-layer thickness.
- A linker wrapper records the real call's input and output. Eight independent
  byte hashes check that clocks, controls, tracer state, composition, PICES state,
  reference profiles, physical boundary metadata, and other production fields
  remain unchanged. This includes flavor counts, phase caches and heating arrays.
  Only temperature may change.
- Pre-age uses CMB=1 and surface=Tref(surface), then the real shallow overlay
  restores the production surface=0. Deep pre-aged nodes survive the overlay.
- `DataT` is final T minus Tref, and particle temperatures interpolate final T.
- Independent global SciPy sparse assembly and direct solves reproduce constant,
  nonlinear k(T), and laterally varying composition cases. Quadrature geometry
  and the depth/composition conductivity prefactor come from production
  helpers; no initializer K, mass, linear
  solve, MPI exchange, or Picard code is reused in the reference solve.
- The lateral-composition fixture also computes an operator without angular
  derivatives. Its distinct answer checks that a radial-column implementation
  cannot accidentally pass the full 3-D reference test.
- Three successively halved timesteps check first-order convergence in constant,
  nonlinear-conductivity and lateral-composition cases.
- A phase-active case checks the frozen latent capacity against an independent
  formula and refines the timestep. Sharp-transition cases require rejection
  of overly large steps, including a case with constant k and capacity that
  isolates the phase-fraction safeguard.
- Setting production Q0=30 leaves pre-age T exactly unchanged and Q0 itself
  untouched. A 25-flavor fixture checks the actual primordial flavor24 mapping.
- Identical 12-rank and 24-rank globally sampled fields check radial interfaces
  and global radial boundary indexing.
- Negative/nonfinite durations, zero/nonfinite timesteps and an incompatible
  solver configuration, the P5 benchmark override, and nonfinite/nonpositive
  initial absolute temperatures must fail before a successful snapshot.
- An ASCII initialization case calls the real production field writer before
  any Stokes solve or production timestep, checks saved dimensional temperature,
  and verifies identical initialized fields to the gzip-mode baseline.
- Frozen-capacity boundary/storage balance is checked separately from nonlinear
  thermodynamic enthalpy; the test does not claim exact enthalpy conservation.

The optional lateral-composition fixture assigns tracer flavors to deterministic
elementwise hemispheres immediately before the pre-age call, then uses the real
recount and composition reconstruction to make all fields and counts consistent.
This is test-only setup, is identical across radial decompositions, and is never
used in a production configuration. Tracer immutability is measured
across pre-age itself; final Tp sampling happens later in normal initialization.

`summary.json`, stderr, per-rank nodes, particle samples, and before/after state
hashes retain the evidence for each run. The independent oracle also retains
quadrature files. Small fixtures verify numerical plumbing and convergence;
physical accuracy still requires refinement on the intended production mesh.

## Restricted-container MPI transport

If the MPICH runtime cannot use UCX sockets or cross-process memory reads, the
following process-local options select portable shared-memory transport without
changing system settings:

```sh
export MPIR_CVAR_CH4_NETMOD=ofi FI_PROVIDER=shm FI_SHM_DISABLE_CMA=1
export MPIR_CVAR_CH4_CMA_ENABLE=0
```

These are runtime-specific options, not required by CitcomS or this test suite.
A successful `mpiexec hostname` is not an MPI transport check: run the test
executable, which calls `MPI_Init` and collective communication.
