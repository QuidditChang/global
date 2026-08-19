# Strict ALA Stage A+B+C Implementation Report

## Scope and baseline

- Branch: `cmbhf_ALA_strict`
- Source commit before work: `d15f5ae8fec15e3d9f38fefc8908571c10497081`
- Isolated source tree: `/work/ess-zhangkd/CitcomS/src/global` on HPC; local development used `/Volumes/Pōwehi/CitcomS/worktrees/strict_ala_stage_abc`.
- This implementation adds provenance and observation infrastructure only. It does not change `D`, `C`, `B=D+C`, `K`, `K_gamma`, gamma, BPI, Schwarz, overlap, depth scales, FGMRES restart, solver selection, gauge treatment, viscosity, reference state, phase treatment, or tolerances.

## Modified implementation

- `lib/Stokes_flow_Incomp.c`: Stage B D/C/B adjoint, gauge, split, and repeatability diagnostics; centralized Stage C iteration and inner-solve observation.
- `lib/General_matrix_functions.c`: production-gated bounded `K_gamma^-1` contract and one-row-per-call CSV.
- `lib/Element_calculations.c`: caller-owned physical-force assembly and Stage C K application counting.
- `lib/Drive_solvers.c`: Stage B dispatch and startup control/operator log.
- `lib/Instructions.c`, `lib/global_defs.h`, and `CitcomS/Components/Stokes_solver/Incompressible.py`: controls, defaults, validation, and counters.
- `tools/strict_ala_stage_abc.py`: cfg generation/whitelist, manifest, Stage B decision, Stage C decision, and final integrity audit.
- `cmbhf_ALA_strict_stage_ABC.lsf`: one clean build and one fail-fast A+B+C allocation.
- `tests/ala_strict/test_stage_abc_instrumentation.py`: local infrastructure tests.

## New controls

- `ala_stage_abc_adjoint_diagnostic`, default `off`.
- `ala_stage_abc_production_logging`, default `off`.

The flags are mutually exclusive, require strict ALA and positive gamma, and Stage C additionally requires pressure FGMRES and the physical momentum gate. Startup logs identify both flags, operator, gamma, iteration budgets, and physical tolerances.

## Observation safety

Default-off calls retain the historical unbounded `solve_del2_u()` route. Stage C uses bounded solves only when its explicit flag is on. Physical momentum measurement assembles force into independent storage, never `E->F`, and reuses the already assembled buoyancy instead of mutating model state. Iteration measurements use dedicated work arrays and never feed values into Krylov recurrence or convergence state. Clean outer-budget exhaustion returns only under Stage C logging; inner failure remains a collective nonzero failure.

## Stage A

The job rejects a dirty or wrong source branch, clones the exact source into scratch, performs one clean `config_script` build, and uses only that executable/library set. The manifest records source status/diff, build environment, scheduler context, executable, every `CitcomSLib*.so`, canonical cfg, reference state, beta interval, and file-valued cfg dependencies discoverable on the compute node. Binary/library hashes are checked before and after all runs.

The cfg generator creates C0 and C1 from one canonical cfg and then checks the final patched job cfgs. The normalized scientific difference must be exactly `ala_shallow_patch_preconditioner: off -> on`.

## Stage B

Stage B measures `<u,A^Tq>` and `<q,Au>` independently for `A=D,C,B`, seven normalized pressure probes, and two deterministic velocity probes. CSV rows include signed/absolute defects, action norms, two denominators, repeated actions, repeatability-derived scale floor, and both normalized defects. Gauge rows include low-degree, constant, and density probes, `D^T/C^T/B^T`, `S_gamma`, energy, exact `B^T-(D^T+C^T)` numerator, and repeated tight actions.

The same binary runs baseline and compatible alternate rank decompositions. The machine gate independently classifies D/C/B and blocks Stage C for true defects, decomposition-specific ownership mismatches, or unresolved evidence. It does not repair operators in the job.

## Stage C

Each accepted pressure iteration records recursive and explicit Krylov residuals and drift; independent continuity numerator/denominator/relative metric; independent unaugmented physical momentum numerator/denominator/relative/RMS; cumulative K applications, Schur actions, preconditioner applications, inner calls/cycles, and elapsed time.

Every production `K_gamma^-1` call records a unique id, iteration, role, RHS norm, requested relative tolerance, absolute target, achieved absolute/relative residual, cycles, maximum cycles, time, and status. A call is valid only if finite, at/below target, and below the cycle limit. There is no fallback.

C0 and C1 are separate executable launches in the same allocation. Each launch reads a separate cfg and separate copies of the immutable cold-start inputs; C1 consumes no C0 output. Numerical non-convergence is retained as a result so C1 still runs, while crashes, NaN/Inf, invalid inner solves, provenance mismatches, and structural failures abort.

## HPC control flow

`stage_A -> stage_B baseline -> stage_B alternate -> Stage B gate -> C0 fresh launch -> C1 fresh launch -> final audit`.

Runtime root:

`/scratch/<date>/ess-zhangkd/STRICT_STAGE_ABC_<JOBID>/`

Generated cfgs:

- `02_C0_BPI/cmbhf_ALA_strict_stage_C0_BPI.cfg`
- `03_C1_CONFIGURED/cmbhf_ALA_strict_stage_C1_configured.cfg`

## Local validation

- New Stage A/B/C tests: 5 passed.
- Full strict-ALA discovery: 133 tests, 131 passed, one pre-existing failure and one pre-existing setup error.
- Pre-existing failure: canonical `runs/cmbhf_ALA_strict.cfg` lacks `ala_pressure_defect_corrections = 1`, required by `test_active_cfg_is_clean_current_rheology_diagnostic` before this task.
- Pre-existing setup error: `test_thermodynamic_closure_oracle` requires uncommitted `tests/ala_strict/STAGE0_THERMODYNAMIC_CLOSURE_REPORT.md`, which is absent from the clean baseline.
- `mpicc -fsyntax-only` passed for all five affected C compilation units. Repository legacy prototype warnings remain; no syntax errors were reported.
- `bash -n cmbhf_ALA_strict_stage_ABC.lsf`, Python byte compilation, and `git diff --check` passed.
- A complete local autotools build was unavailable because the workstation `autoreconf` cannot execute `aclocal`. The LSF performs the required clean build with the documented HPC modules before any diagnostic launch.

No production numerical diagnostic was run locally.

## Submission and return files

Submit from the committed source checkout:

```bash
bsub < cmbhf_ALA_strict_stage_ABC.lsf
```

Return the entire `STRICT_STAGE_ABC_<JOBID>` directory or its failure archive. At minimum retain the Stage A manifest/cfg diff/build log/pre/post hashes, both Stage B CSVs and decision, both raw Stage C logs and cfgs, combined iteration/inner/cost CSVs, Stage C decision, final audit JSON, scheduler accounting, and `STRICT_ALA_STAGE_ABC_COMPLETE` when present.
