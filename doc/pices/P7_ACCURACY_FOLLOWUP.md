# P7 accuracy follow-up — 2026-10-10

Status: diagnosis and reproducible HPC task prepared; numerical resolution is
not yet established. No C solver or physical production parameters changed.

## Evidence and interpretation

P7 job 12275166 at dt/8 has PG/PICES nodal temperature RMS difference
11.5173 K. The radial RMS profile, bottom to top, is
`0, 19.5312, 9.6929, 6.3592, 6.6122, 6.5640, 11.7638, 21.2787, 0` K.
The strongest difference lies next to fixed temperature boundaries. This is
localization evidence, not proof of a boundary implementation defect.

The old velocity percentage measures the ratio of scalar volume RMS speeds,
not the difference of velocity fields. PG is not an exact solution. P6 already
shows sensitivity to particle count; shrinking dt alone cannot establish
spatial convergence. The final PICES `mismatch` also compares an unconstrained
projection with a field evolved from a Dirichlet-constrained reconstruction;
it is not, by itself, an error against an exact temperature solution.

## Minimal discriminating experiment

`prepare_p7_accuracy.py RUNS` generates ten cases using the P7 dt/8 physical
configuration and current runtime. All reach nondimensional time 4e-7:

- PG/PICES pairs on nested 5, 9 and 17 node grids per cap direction; PICES uses
  64 particles per element on average; 32 steps at dt=1.25e-8 (six cases).
- On the 9-node grid, PICES additionally uses 32 and 128 particles per element
  (two cases). This checks successive particle sensitivity, not random-seed
  uncertainty; existing tracer placement is retained.
- Both methods on the 17-node grid also take 64 steps at dt=6.25e-9 (two cases).
  This measures temporal contamination of the fine-grid comparison.

Reference profiles and age forcing use the same analytic definitions at each
resolution. The inner Stokes budget is 2000 for all ten cases; tolerance and
algorithm are unchanged. Files are prepared locally; the HPC LSF directly
launches each case, without a Python dependency or new shell wrapper.

`verify_p7_accuracy.py JOB --summary report.json` reuses the existing field,
TA, heat/CBF and provenance checks, and adds projection/coverage checks and:
paired radial RMS; successive own-method spatial differences at common nodes;
32→64→128 particle sensitivity; finest-grid time sensitivity. Spatial common
nodes are sampled, not volume-weighted; boundary copies are retained.

Interpret results together. If time sensitivity remains material, the grid
experiment is temporally unresolved. If particle sensitivity remains material,
fix/investigate sampling before attributing the gap to the grid. If both methods
converge spatially and their gap decreases, a coarse-grid method difference is
supported. If a gap persists, isolate transport versus heat/TA and boundary
reconstruction next. Do not tune either method to force agreement or relax the
production gate to make this test pass.

## Validation and execution

Generator determinism, matrix/config invariants, Python compilation and Bash
syntax are checked locally. MPI numerical execution has not been completed for
this follow-up; the previous approval service usage-limit rejection prevents
local escalated MPI checks. HPC evidence is required before claiming resolution.

Direct LSF: `cmbhf_EBA_PICES_P7_accuracy.lsf`, ser, 12 ranks on one node,
8-hour limit, current dated scratch directory. Download the resulting tar.gz
and scheduler .out/.err. Existing P7 and P8c outputs remain the comparison record.
