# P6: opt-in bounded consistent projection

Status (2026-10-06): HPC job 12259598 passes all 46 matrix cases and coupled
restart checks. The original large stationary projection defect is substantially
reduced. Controlled production-configuration pilot may proceed, but formal
production remains on hold: assimilation vrms differs from PG by 15.96% at base
dt (4.37% at dt/4), the launcher/library MPI providers are mixed, and long-time
real-configuration validation remains outstanding. HPC tests did not activate
projection bounds; that branch currently has local sharp-profile coverage only.
See runs/PICES_P6_HPC_AUDIT_12259598.md and its JSON for full evidence.

## Problem and method

P5 frozen-state isolation showed that lumped P(Q(T)) changes even an unmoved,
finite-element grid temperature. The assimilation baseline loses about 1.29% of
its initial FE heat proxy before any physical operator runs. More particles or
smaller timesteps do not remove this consistency defect.

`pices_projection=bounded_consistent` selects a matrix-free least-squares
projection with the existing interpolation weights W. Absolute temperature solves
min ||W g - Tp||² / 2, with existing fixed nodal temperatures and global bounds
spanning current particle and previous nodal temperatures. It first solves the
unconstrained free-node system with Jacobi-preconditioned CG, then uses projected
Jacobi iteration to satisfy the bounded quadratic problem. The row sums of WᵀW
provide the diagonal scaling. The reported residual is the maximum scaled
projected-gradient step (KKT residual), not an energy error.

Signed subgrid/TA increments and diagnostic projections use the unbounded linear
operator. Applying absolute-temperature bounds to those increments would corrupt
cooling and violate P(Q(delta)) = delta. Duplicate MPI interface nodes are counted
once in global dot products. Every local node must have positive particle support;
loss of support or convergence aborts rather than silently substituting a method.

The default remains `lumped`. The new method contributes an explicit version tag
to checkpoint physics fingerprints; switching methods across restart is rejected.
Old lumped fingerprint bytes and checkpoint layout are unchanged.

## Scope and limitations

The particle-jump probe intentionally replaces particle temperatures while retaining
the original grid as a boundary reference. Its energy difference from that grid
is not a remap-conservation error; use its bounds and independent KKT check.

This is a consistent reconstruction with global bounds. It does not prove local
monotonicity, physical-particle-mass conservation, or arbitrary-decomposition
restart. No global energy rescaling is applied. When bounds activate, the
reconstruction can differ from an unconstrained solution and its FE heat proxy
can change. Moving-particle remap energy remains a quantity to audit.

The matrix-free implementation caches 8 node indices and 8 weights per local
particle for each projection. It performs MPI exchanges in iterative solves.
Memory and wall cost must be measured on production-sized meshes. Rank-deficient
particle arrangements are not regularized; a converged least-squares residual
alone cannot establish uniqueness for an arbitrary future particle distribution.

Standalone C compilation and MPI execution are tested locally. The Pyre property
is wired but its Python 2 runtime is only available in the HPC build environment.

## Reproducible local checks

- `tests/pices/run_consistent_projection.py`: 18 frozen cases; independent SciPy
  assembly verifies distributed unconstrained residual, bounded KKT residual,
  fixed nodes, global bounds, FE reconstruction, and signed increments.
- `tests/pices/run_p2_smoke.py --stage P4`: continuous/split/restart, decoded
  fields, CBF, and binary live state, for both defaults and new projection.
- `tests/pices/run_p4_guards.py`: rejects a changed projection, TA law/targets,
  forcing, or corrupted CBF state before advancement.
- `tests/pices/run_p5_local.py`: selected 37/46 local cases with PIC opt-in; PG
  remains the unchanged comparator. The full 46 cases are mandatory on HPC. `verify_p6.py` adds projection and restart checks.

## HPC handoff

`prepare_p6.py RUNS` derives the new `pices_p6` inputs from reviewed P5/P4 inputs.
`cmbhf_EBA_PICES_P6.lsf` submits directly: ser, 12 cores on one node, 4h limit.
The job executes 46 matrix cases and then continuous/split/restart (4/2/2 steps).
No cluster Python or submission shell wrapper is required. Matrix stdin is
isolated from MPI. Logs, archive, executable dependency listing and launcher path
are captured under the dated model job directory. Download the archive and the
scheduler `.out/.err`; audit with `verify_p6.py` plus scientific comparison to P5.

## Gate to production

1. P6 local and HPC results must remove the static defect and pass temperature,
   heat/TA/CBF, convergence, restart, and resource checks. Shape error against a
   rotated initial field with diffusion active is not an exact-solution error.
2. Run a pilot using the intended real production initial state, viscosity,
   reference state, forcing, grid, partition and duration; compare PG/PICES,
   exercise save/restart, and inspect long-time budgets and particle coverage.
3. Limited production is allowed only for configurations supported by that
   evidence. Formal long runs follow acceptance of the pilot. Passing P6 alone
   does not authorize all PICES configurations or establish a calendar deadline.

## Local evidence already completed

- 18 frozen cases: 17 FE-representable cases recover the previous grid exactly
  from the initial guess; independently solved signed and unconstrained fields
  verify the operator. The sharp particle profile reaches bounds [0,1] with
  independent scaled KKT residual 8.70e-13; it required 42 CG and 139 projected
  iterations. Signed increments retain negative values.
- P4 with new projection: 216 decoded output comparisons, 120 CBF comparisons,
  and binary live-state equality pass. CBF/heat residual is 3.17e-16.
- All 6 rejection guards pass, including changing projection on restart.
- Existing default lumped P4 regression also passes.
- Multiple heat substeps with TA: grid/particle increment errors below 6e-17,
  CBF unchanged across TA, heat/CBF residual 1.07e-16. The signed mapping defect
  is below 3e-13 on both tested grids. A near-roundoff defect need not decrease
  monotonically with refinement.
- The direct matrix LSF loop passes stdin-consumption, truncated-manifest,
  early-exit and MPI-failure regression checks.

## Local time-integration scope and preliminary comparison

37/46 matrix cases completed locally: all 28 cold/hot transport cases, four
assimilation PG comparators (base, half/quarter dt, coarse grid), assimilation PIC
base, and all four real coupled cold/hot cases. The unchanged assimilation PG
fine-grid comparator was interrupted for local runtime cost; it had no numerical
failure at interruption. The remaining assimilation PIC sensitivity cases and
that PG fine-grid comparator must run on HPC. No local 46/46 claim is made.

Cold/hot base final PG/PIC RMS differences are 0.6410/0.4847 K. Assimilation base
is 9.2055 K. Earlier HPC P5 values were approximately 12.8/12.7/71.6 K; these are
method differences, not errors against a known exact solution, and the machines
are different. The local first-step
assimilation remap changes from the prior lumped control -1.293626% to the new
-0.000144509% of initial FE energy proxy. Residual advective and constrained
boundary reconstruction effects remain.

The new matrix-free iteration costs more than lumped averaging. Local timing
includes concurrent validation workloads and is not a production-speed estimate.
HPC launcher now records MPI executable path and executable shared libraries;
module names alone did not establish the actual MPI provider in P5.
