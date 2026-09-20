> Naming update: native files now use expanded rank datadir/q.botm.<rank>.<step>
> and q.surf.<rank>.<step>; see README.md. Earlier paths below are historical.

> Current update: runtime GRD output has been replaced by native Q1 rank files.
> See README.md for the format and one-step restart benchmark. GRD results below
> describe the earlier output implementation, not the current runtime format.
> Native builds require no NetCDF; 12/24-rank one-step local runs and native
> node-equation/face-integral validation passed. HPC checkpoint benchmark has
> been configured but not run locally (remote checkpoint data unavailable).

# CBF implementation status (2026-09-19)

Authorized implementation and local validation are complete on `cmbhf_EBA`.
Read [the validation report](VALIDATION_REPORT.md) for evidence and limitations,
and [runtime/build instructions](README.md) for usage.

- Appendix C method audited; Q1 GLL boundary mass and physical Galerkin residual implemented.
- Both boundaries output runtime NetCDF GRD in PostProc/HF_CBF at cfg monitoring intervals.
- Numerical kernels, phase terms, MPI radial exchange, spherical mapping, native integrals and non-mutation validated.
- User explicitly authorized the shared parser fix: two invalid pointer writes corrected, comma fields parsed correctly; 8 parser cases pass.
- Final build /tmp/cbf-final-build compiles repository sources without a temporary parser patch.
- Final 12/24-rank smoke runs both reach cycles=1 and write step 0/1 top/bottom grids; all eight grids pass independent size/finite/seam/pole checks.
- Earlier manufactured real-mesh 12/24-rank outputs are bitwise identical; 5/9/17-node refinement decreases local and integral errors; 96 regular gzip products are unchanged with CBF enabled.
- Full production HPC runs and complete literature benchmark suites are outside completed local validation. Output-state/solver-Tdot and point-sampling limitations remain documented.

No pending approval or parser startup blocker remains. Pause hourly automation `cbf`; do not restart completed work.
