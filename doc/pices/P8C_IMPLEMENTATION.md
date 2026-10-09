# P8c: coupled restart and original-model pilot

Branch: `cmbhf_EBA_PICES`; based on P8b `ad70c94948bfd9511c421921cbbc152e54847fb7`.
HPC acceptance is pending. This phase does not change the Stokes algorithm.

## Accepted-state ownership

The existing legacy checkpoint still owns T/Tdot, U/P, particles/Tp/flavors,
composition and clock. The existing `.pices.state` owns accepted nodal velocity,
P4 heat residual/CBF validity and geological reclassification history. Schema 4
appends finest-level quadrature viscosity (8*nel float32) then nodal viscosity
(nno float32), both without unused index zero. It is selected for rheol7,
plate files, qvis, or automatic timesteps. P4/EBA is required. Older supported
constant-viscosity cases retain schema 1/2/3 and their fingerprint layout.

On restart, restore this state before rebuilding thermal coefficients. Do not
call initial_viscosity over schema-4 accepted viscosity: the next thermal step
must use exactly the viscosity belonging to its accepted Stokes state. The
next normal velocity solve rebuilds rheology, coarse viscosities, matrices and
preconditioners through the existing VISC_UPDATE path. No second restart stack,
particle system or post-restart Stokes solve is introduced.

Supported checkpoint rheology remains system-generated Newtonian rheol1
(uniform, temperature independent) or rheol7, with VISC_UPDATE enabled. Nonlinear,
compositional, frozen and prescribed weak-zone/channel rheologies remain rejected.

## Integrity

Existing payload/state SHA-256, exact sizes, collective manifest, strict compiled
commit, same partition and canonical metadata checks remain mandatory. Schema 4
extends the physics hash with cold_scale, qvis strength and reference pressure,
automatic timestep ceiling, MG controls, viscosity averaging/smoothing/layers,
and element material identity. Restored quadrature/nodal viscosities must be
finite and positive. `accepted_velocity_sha256` keeps its legacy metadata name;
it hashes the *entire* state companion, including new viscosity arrays.

Raw current-bracket plate/age/trench files and current new-trench/transform files
are hashed by content, not absolute path. These checks identify inputs at the
checkpoint age, not all future geological inputs. The two-step production pilot
separately records and checks all its required 249/250 Ma files before and after
MPI. Long production continuations still need an immutable input manifest for
their complete future forcing interval; schema 4 alone does not certify that.

## Gates

`prepare_p8c.py`, `run_p8c_local.py`, `verify_p8c.py` implement a 12-rank 9^3-cap
continuous 0->4 / split 0->2 / restart 2->4 gate, with normal initialization
(p5_case=off), MG, rheol7, 25 flavors, TA, plate files, qvis and automatic dt.
Start 2.13 Ma; dt sequence approximately .05,.05,.03,.05 Ma crosses 2 Ma.
The audit checks exact decoded outputs and binary live state, manifest integrity,
TA once per step, native CBF integration/balance, MG residuals, projection,
heat ledger and particle/primordial conservation. Negative tests change qvis,
cold_scale, maximum dt, plate contents and a state byte; all must stop before
any PICES_STEP. The prior P8a gate covers unchanged older checkpoint schemas.

Local evidence is summarized in runs `PICES_P8C_LOCAL_VALIDATION.json` and
`PICES_P8C_RUNBOOK.md`. Local MPI is Open MPI on macOS; Intel MPI/Pyre on HPC
and actual forcing data are not emulated by these tests.

## Original-model pilot

`prepare_p8c_pilot.py` derives a Pyre cfg directly from
`cmbhf_EBA_Q0_30_rheol7_scold1.0_LLSVPsV0.04_B24_0.3_TA.cfg`.
`PICES_P8C_PILOT_CONFIG_DIFF.json` records every section/key difference.
Grid 129x129x65 per cap, 384 ranks (4x4x2x12), 27 particles/element,
249.9 Ma fresh initialization, Q0=30, rheol7/cold_scale=1, B24=.3,
kC=.8, qvis and all geological paths remain the target's values.
Only PICES controls, heat-stage CBF semantics, job name/queue and two-step
output/checkpoint frequency change. LSF copies reference/coordinate/polygon
inputs and makes run-local path substitutions. No flat-config translation.

The medium pilot requests 384 slots, max 40/node (normally ten nodes), 24h.
It runs the existing Pyre application directly without concurrent visualization.
339,738,624 initial particles make full checkpoints large; require 100 GiB free,
keep checkpoints on HPC and do not automatically tar the whole run.
The 12-rank ser restart gate must first pass HPC audit.

Production release remains separate: measure actual-data memory/runtime and
resolve the outstanding P7 temperature-accuracy trend with targeted spatial/
particle-density comparisons. Passing restart integrity is not a scientific
accuracy certification.
