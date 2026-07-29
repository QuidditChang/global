# ALA strict regression baseline

## Baseline identity

- Legacy branch: `cmbhf_ALA`
- Legacy commit: `deeccebf22e55a6d608163aa5bfc8be8129e4d0b`
- Strict branch start: the same commit
- Configuration: `runs/cmbhf_ALA.cfg`
- Reference state: `runs/refstate_ALA.txt`
- Configuration SHA-256:
  `c4a59b6e07bcfd22771a7eecaa14698d3ee68d459e48423fc4296b105d710b21`
- Reference-state SHA-256:
  `af8967ab45dd84ccc3deb4042c8b2926984b46226c4e21c8fab619fe21945874`
- Reference-state shape: 65 numeric rows by 7 columns

Commit 0 does not copy or rewrite either runtime input.

## Frozen solver options

The active configuration uses:

- `compressible_formulation=ala`;
- `uzawa=ala_cg`;
- PCG as the active strict-ALA pressure iteration;
- `ala_inner_accuracy_max=1e-4`;
- physical cancellation tolerance `1e-2`;
- augmented-Lagrangian and multigrid infrastructure already present on the
  legacy branch.

The code also retains:

- FGMRES through the complete ALA operator path;
- BiCGStab through `uzawa=bicg`.

Commit 0 does not select, tune, or modify any solver.

## Archived numerical evidence

The repository contains documentation summaries rather than the complete raw
MPI logs. The following table records only evidence already present in the
baseline documentation.

| Quantity | Archived value | Interpretation |
|---|---:|---|
| Iteration count | 55 | Last reported feasibility iteration, not a convergence count |
| Physical cancellation | `0.7777` | `global_cancellation_L2` at iteration 55 in the documented `1e-6` diagonal-BPI feasibility run |
| Schur symmetry defect | `1.820353e-5` maximum | 37 samples, below the documented `1e-3` feasibility threshold |
| Nonpositive curvature count | 0 | Reported for the symmetry feasibility evidence |
| Residual depth attribution | `92.09%` above 410 km | Iteration-55 residual-energy fraction |
| Algebraic residual | Not archived | No complete raw log is committed |
| Velocity norm | Not archived | No trustworthy value is available in the committed evidence |
| Thermal budget | Not closed or archived | Existing diagnostics do not constitute a complete thermal balance |

Missing values are intentionally marked as unavailable. They must not be
reconstructed or invented from partial summaries.

## Solver comparison record

| Solver | Frozen status | Commit 0 evidence |
|---|---|---|
| PCG | Active experimental path | Partial iteration-55 and symmetry evidence above |
| FGMRES | Implemented baseline option | No complete current MPI log archived |
| BiCGStab | Preserved rollback baseline | No complete current MPI log archived |

## Commit 0 verification

The Commit 0 verification was run on 2026-07-29.

- Build: **PASS** in an isolated temporary source copy.
  - The source tree had no configured `Makefile`, so the temporary copy was
    configured before running the requested `make clean` and `make`.
  - The local host lacks the legacy Python 2/Automake environment. The pure
    C/MPI build therefore disabled Pyre, HDF5, and Exchanger and used the
    legacy-compatible `gnu89` language mode.
  - Both `CitcomSFull` and `CitcomSRegional` were produced.
- `make check`: **PASS**, with no registered C test targets in `lib` or `bin`.
- Static reference-state structural check: **PASS**, reporting `(65, 7)`.
- Configuration and reference-state hashes after verification: **UNCHANGED**.
- Source equivalence: **PASS**. There is no diff from `cmbhf_ALA` in C headers,
  C sources, build files, solver files, configuration, or reference data.
- Fresh `cmbhf_ALA.cfg` MPI short run: **NOT RUNNABLE LOCALLY**. The production
  case requires 384 solver ranks and external run-site inputs.
- Temporary MPI smoke tests:
  - 12-rank `CitcomSFull`: local executable terminated with `SIGILL`;
  - one-rank `CitcomSRegional`: local executable terminated with `SIGILL`.

Both smoke failures occurred in the unchanged legacy executable before usable
velocity, temperature, pressure, or residual output was produced. Consequently
Commit 0 does not claim a fresh dynamic-field regression pass. It records that
no numerical-code or input difference exists from the frozen branch; a fresh
field comparison remains an HPC validation gate.

No initializer is run during Commit 0 because it rewrites configuration,
reference, mesh, and diagnostic products.
