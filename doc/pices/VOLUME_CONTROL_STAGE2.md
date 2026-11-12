# EBA composition volume controls, stage two

This opt-in experiment diagnoses particle-count-ratio reconstruction and
transport errors. It does not apply a conservation correction. It initializes
the preserved RH20 configuration from zero once, freezes particle labels and
thermal evolution, and runs five independent controls from the same state.

| Case | Particle path | Steps |
| --- | --- | --- |
| A_static | No motion; repeat composition reconstruction | N |
| B_exact | Exact volume-preserving map; production host search/MPI exchange | N |
| B_dt | Production grid interpolation and predictor/corrector | N |
| C_half | Same as B_dt | 2N |
| C_quarter | Same as B_dt | 4N |

The requested N is 100 and the requested end time is 4e-6 nondimensional.
The core increases N when needed to keep the analytic maximum speed times dt
below 0.25 times the shortest mesh edge. The requested and realized parameters
are recorded. The realized end time is N times the float-representable base
dt; half and quarter dt are exact binary subdivisions. The trajectory clock
is double and the legacy audit clock is float. The report checks the latter
against the actual float32 rounding of the former.

## Analytic transport

Let c=0.64(1,1,0.2)/sqrt(2.04), q=x-c, R=0.08 and
omega(q)=(pi/T)(1-|q|^2/R^2)^2 inside the support ball, zero outside.
The velocity is omega(q) ez cross q. Since q dot velocity=0, |q| and
omega remain constant on a trajectory. The exact map is a rotation of q
around ez through angle omega(q)t, added back to c. Its divergence is zero;
the support lies strictly inside the shell. The eccentric centre permits
radial movement relative to the Earth's centre, including through the top
of the initial C24 layer. Exact coordinates are evaluated from initial
coordinates at total elapsed time rather than accumulated rotations.

The normal cases fill the actual spherical velocity array used by production
particle interpolation. They do not substitute analytic particle velocities
for the grid interpolation. Guards abort if the existing speed limiter or
shell boundary correction would modify the prescribed experiment. An empty
active C24 population also aborts, avoiding an uninformative control.

## State, isolation and acceptance

`PICES_VOLUME_CONTROL=1` reserves four migrating double particle attributes:
initial Cartesian coordinates and immutable initial flavor. Registration
rejects restart, non-EBA and conflicting motion-audit modes. At the end of
PICES initialization, before either driver advances heat or Stokes, the
controlled driver performs the five cases and exits collectively.

`Initial/baseline.RANK.bin` stores all populated particle quantities, labels,
hosts, element/node composition, flavor counts and nodal temperature. These
native files are private snapshots for this executable and MPI topology,
not production restart checkpoints. Every case restores and checks the same
SHA256-verified snapshot, including the fallback for empty cells. Counts of
all 25 flavors and two permutation-independent property fingerprints check
MPI migration at every step. Nodal temperature is checked unchanged at each
case end. The extra four attributes require about 10.9 GB across 384 ranks
for 339 million particles; the retained snapshot occupies tens of GB.

A fails if any species' volume changes by more than 1e-10 of shell volume,
or if any particle moves. B_exact fails if the path differs by more than
1e-11 nondimensional from its analytic reference. All cases require unchanged
labels, counts and temperature properties. Finite volume drift in a moving
case is an output, not a failed control. Existing pointwise nodal fraction
warnings remain in the solver logs; a closed volume ledger does not certify
pointwise bounds.

Each `Cases/NAME/Diagnostic/VolumeAudit` directory holds the unchanged stage-one
ledger. Its accepted marker means a completed controlled transport step here;
`control_case.json` records this semantic distinction from a production
Stokes acceptance. `trajectory.csv` includes errors for all particles, active
particles, and active C24 particles. The active metrics avoid dilution by
the stationary majority. `cases.csv` and `control_complete.json` certify all
five completed cases. The separate report validates their raw histories and
produces one comparison figure, CSV and JSON.

B_exact drift demonstrates that exact incompressible trajectories do not make
count-ratio reconstruction conservative. B_dt minus B_exact is a controlled
comparison, not an additive percentage attribution to the full RH20 run.
Reducing dt can leave a spatial-interpolation error floor, and reconstructed
volume need not converge monotonically. This compact synthetic flow cannot
predict the amount of drift in the evolving geological velocity field.

## Local validation, 2026-11-12

The full standalone solver was built with the system MPI compiler. The
production mathematical map passed a 1,600-point finite-difference test:
maximum |det(J)-1|=1.97e-9, normalized divergence=2.2e-10 and inverse mapping
error=2.22e-16. These are finite-difference checks, not new runtime tolerances.

Real 12-rank and 24-rank small-mesh controls completed all five cases with
393,216 particles, unchanged all-flavor counts and identical initial volumes.
The static drift was zero. A nonuniform radial 24-rank test placed the MPI
interface within the vortex and verified 202 outward and 31 inward departure
events using production migration logs. Its largest normalized ledger
residual was 1.283e-16; the comparison report passed.

With controls disabled, 96 physical checkpoint payloads at steps 0 and 2
matched the committed pre-control solver (7a2d427e4bb281a99b147624d14e48727dae018b),
both with and without the stage-one audit. Unused index-zero storage was masked
by the existing checkpoint comparison helper. Launcher/report tests cover
incomplete cases, failed static gates, altered baselines, clocks, labels,
frozen runtime, affinity and process cleanup.

Reproduce with a fresh temporary build/output directory:

```bash
MPICC=/usr/local/bin/mpicc python3 tests/cbf/build_validation.py \
  --build-dir /tmp/pices-controls-build --optimization O1
python3 tests/pices/check_volume_control_math.py
MPIEXEC=/usr/local/bin/mpiexec python3 tests/pices/check_volume_control_integration.py \
  /tmp/pices-controls-build ../runs-cmbhf_EBA /tmp/pices-controls-12 \
  --report ../../scripts/PostProc.volume_control.py
MPIEXEC=/usr/local/bin/mpiexec python3 tests/pices/check_volume_control_integration.py \
  /tmp/pices-controls-build ../runs-cmbhf_EBA /tmp/pices-controls-24 \
  --radial-partitions 2 --report ../../scripts/PostProc.volume_control.py
MPIEXEC=/usr/local/bin/mpiexec MPICC=/usr/local/bin/mpicc \
  python3 tests/pices/check_volume_control_production_regression.py \
  /tmp/pices-controls-build ../runs-cmbhf_EBA /tmp/pices-controls-regression \
  --baseline-ref 7a2d427e4bb281a99b147624d14e48727dae018b
```

The HPC launcher retains the production 384-rank layout in a 400-slot medium
allocation, without an online postprocessor. Its 48-hour limit is a requested
ceiling, not a measured production runtime. The local synchronized RH20 folder
lacks the original grid file, so a real-input `--check-only` must succeed on
HPC before submission. The separate scripts repository contains the launch
recipe in `PICES.RH20_volume_control.md`.
