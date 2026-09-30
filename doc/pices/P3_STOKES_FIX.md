# P3 Stokes tolerance repair after HPC job 12249186

2026-09-30. The HPC job failed its continuous step 3 at 500 pressure iterations.
Its inner velocity tolerance stayed at 4.085e-3, set by accuracy times the
initial force norm. At outer iteration 300 the current velocity RHS norm was
6.722e-3, and a single inner iteration returned residual 1.491e-3. At iteration
500 the corresponding numbers were 8.689e-3 and 2.257e-3. These are poor
relative solves of the increasingly small Schur-action RHS. The outer residual
stagnated and then increased. This is a concrete precision defect; the archive
alone cannot establish why one MPI/compiler trajectory stalls and another does
not.

## Change

Only when `pices.enabled && pices.eba`, use
`max(1e-14, 1e-6 * RMS(current_velocity_rhs))` for the initial momentum correction
and each CG Schur-action velocity solve. The outer accuracy remains 1e-6 and
the pressure iteration limit remains 500; no acceptance threshold is relaxed.
There is no new cfg option. PG and P1/P2 keep their original tolerance arithmetic.
An exploratory 1e-8 inner relative target failed a local inner solve; the tested
1e-6 policy is used, rather than accepting an unconverged inner result.

P3 immediately calls MPI_Abort through pices_fail if an inner solve fails, or
if the outer solve returns without satisfying the existing OR stopping rule.
The outer guard also rejects nonfinite diagnostics and negative increment
ratios as evidence of convergence. Rank 0 writes PICES_STOKES with step, count,
inner policy, outer tolerance and PASS/FAIL. This prevents the old behavior of
advancing after the pressure iteration cap and publishing later checkpoints.

The LSF additionally requires these in-process PASS records at every requested
step, alongside its existing independent AhatP and heat-ledger checks. An old
P3 binary without those records is rejected.

## Validation

- Local 12-rank continuous/split/restart test: 55/48/48/48 pressure iterations,
  maximum CFL 0.03883063749; 216 corresponding decoded outputs and effective
  binary checkpoint arrays match exactly. Step 1 meets dp/p<1e-6; subsequent
  steps meet both increment bounds. These are the unchanged OR criteria.
- Requested a different Open MPI allreduce algorithm (recursive doubling via
  coll_tuned settings) and repeated all three cases including the revised LSF
  gates: PASS, same output comparison and iteration counts. This is local
  robustness evidence, not an Intel MPI test or proof that every collective
  used a different implementation.
- `piterations=52`: initialization converges, step 1 fails its outer limit,
  exit 72, no step 2 advancement and no failed-step checkpoint/manifest.
- `vlowstep=1`: initial momentum solve fails, exit 72, no checkpoint published.
- P2 three-case regression passes. Its continuous 180 decoded velocity/
  temperature, particle and viscosity files match the pre-fix run bit-for-bit;
  constant-physics PICES takes the unchanged branch.
- The new LSF record gate accepts the repaired result and rejects the old
  locally successful P3 output that lacks in-process guards. Bash syntax and
  git diff whitespace checks pass.

See P3_STOKES_FIX_RESULT.json for the machine-readable local results. Rebuild
from the new commit on HPC and run all three cases from fresh initial state.
Do not resume the failed job's step 4. P3 remains pending HPC acceptance, and P4
has not started. Thermal cfg, timestep, mesh, viscosity and phase inputs are
unchanged.

```bash
python3 tests/pices/run_p2_smoke.py BUILD ../runs-cmbhf_EBA --stage P3 --output NEW_OUTPUT
python3 tests/pices/run_p3_stokes_guards.py BUILD ../runs-cmbhf_EBA --output NEW_GUARD_OUTPUT
```
