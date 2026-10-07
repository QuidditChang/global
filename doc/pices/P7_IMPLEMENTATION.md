# P7 targeted preproduction gates

Target: `cmbhf_EBA_Q0_30_rheol7_scold1.0_LLSVPsV0.04_B24_0.3_TA`.
P7 is not yet this production model. Its input package pins a reference-only copy
and SHA256 of the exact target. Matrix execution never submits that snapshot.

## Changes

- Direct LSF selects `mpiexec.hydra` from the installation of the executable's
  resolved Intel `libmpi.so.12`, rejects an unexpected provider, and archives
  the path, version and shared-library listing. No PATH fallback to ParaView.
- Eight P6 assimilation PG/PIC controls: dt, dt/2, dt/4, dt/8 at equal final
  time 4e-7. This isolates the remaining velocity sensitivity with matched MPI.
- Two rheol7 probes: target cold_scale=1, temperature dependence and viscosity
  limits, fixed zero upper velocity, 8 steps of 1e-9, synthetic P6 reference
  state and TA. Single-level CG, vlowstep=2000, unchanged accuracy 1e-6.
  They are not the production multigrid mesh or its composition/forcing.
- One explicit `p5_case=sharp` stress: zero prescribed flow, one 1e-10 step,
  radial particle step at 0.77 (away from radial nodes), bounded consistent
  projection. Requires PIC, bounded projection, fresh run and no checkpoint.
  Other benchmark modes and production initialization are unchanged.
- Same-partition endurance: 64 continuous, 32 split, restart 32→64. dt=2.5e-8
  keeps the total physical time within the existing age forcing coverage.
  Every step's fields/CBF and every two-step checkpoint are audited.
- Existing P2/P4 audits accept optional endpoints, keeping defaults 2/4. Float
  elapsed-time comparisons use an accumulation roundoff bound for long runs.
- Sharp's nominal zero boundary conduction uses a 1e-12 nondimensional absolute
  flux floor; all ordinary cases keep the previous relative CBF criterion.

The launch is ser, 12 ranks on one node, 4h limit. It runs 11 matrix cases plus
three endurance segments. Output/err stay under the dated model directory.
The launcher requires actual final-step logs, bound activation in sharp, and
restart field equality; legacy MPI exit 8 alone is not success.

## Local validation

- New rheol7 PG/PIC 8-step controls and sharp execute successfully.
- Sharp: 42 CG + 134 projected iterations, 4621 active free nodes, scaled
  KKT residual 9.768e-13 <= 1e-12. Global reconstruction bounds are [0,1].
- 64-step continuous/split/restart: 2376 decoded outputs, 1560 CBF outputs,
  binary live-state equality PASS.
- The generalized verifier still passes the full archived P6 HPC result.
- Matrix-loop regression passes normal, truncated input, early break, MPI
  failure and stdin-isolation checks. Shell syntax and Python compilation pass.
- Eight temporal controls were subsequently verified on HPC (see below).
- Intel MPI provider selection cannot be executed on this macOS/OpenMPI host;
  the HPC script checks the actual library/launcher installation before runs.

Exploratory rheology setups using free-slip CG and coarse two-level multigrid
failed the existing Stokes guards; no guard or tolerance was weakened. The final
fixed-top single-level control passes. The production multigrid setup remains
unvalidated. A first long-run attempt exceeded the age-file window and exited
with legacy status 8; final step validation catches this, and the final dt keeps
all 64 steps within coverage.

## Production compatibility: still blocked

The named target has 384 ranks, 25 tracer flavors, chemical buoyancy, flavor
reclassification, kC_ratio=0.8, qvis_mode=2, file-driven plate velocities and
rheol7. Current PICES explicitly rejects the first composition/transport/source/
boundary combinations, and checkpoint supports only uniform constant viscosity.
These are implementation gaps, not issues solved by running P7 longer. P7 does
not remove those guards. After auditing P7, development must enable and test
these exact target features, including their checkpoint/fingerprint state,
before a faithful production pilot or formal production can be approved.

Audit downloaded results with `tests/pices/verify_p7.py JOBDIR --summary RESULT`.

## HPC audit: job 12275166 (2026-10-06)

PASS for the 11 targeted cases and same-partition constant-viscosity 64-step
restart: 2376 decoded outputs, 1560 CBF outputs and binary live state match.
Intel launcher/library installation matches. Sharp activates 4622 bounds
with 42 CG and 133 projected iterations, residual 8.52e-13.

Assimilation velocity difference drops from 15.9633% to 1.5233% at dt/8.
However, PG/PIC temperature RMS difference grows from 9.2303 K to 11.5173 K.
This is not evidence of full method equivalence or production readiness.
The rheol7 probe does not validate rheol7 restart or the actual target setup.

The initial launcher wrote a P7 completion value under p6_complete.txt.
Future launchers use p7_complete.txt; audit accepts the legacy name only
for P7 and verifies the exact P7 content. No HPC rerun is needed.
The matrix report also now uses len(rows) for cases_expected, avoiding
reuse of the particle-count variable. Numerical checks are unchanged.
See runs PICES_P7_HPC_AUDIT_12275166.md/.json for full evidence and
shared tracer/physics/checkpoint interface constraints on production work.
