# PICES route B — P2

2026-09-30: implementation and local validation complete; HPC validation pending. P1 HPC job 12244722 passed. Work remains on solver branch `cmbhf_EBA_PICES` and runs branch `cmbhf_EBA`.

## Checkpoint protocol

`pices_checkpoint=on` explicitly enables P2 writes and same-partition restart. The default is off, retaining P1 behavior. The old PG general header and payload layout are unchanged. Tp is serialized by the existing tracer extraq segment, after flavor when present.

Each rank publishes these files together for an accepted state:

- `PICES_P2.chkpt.RANK.STEP`: legacy binary payload.
- `.pices.json`: canonical schema-1 metadata, magic `CITCOMS_EBA_PICES`, compiled solver commit, accepted phase, step/time/dt/counter, MPI decomposition, local mesh/cap, particle count and attribute registry, temperature normalization, geometry versions, physics/mesh fingerprint, payload and accepted-velocity SHA256.
- `.pices.state`: three float arrays of accepted nodal spherical velocity, omitting unused index zero. This is necessary because the original payload saves equation velocity U, while particle transport uses nodal V after rigid-rotation removal. Reconstructing V from U alone changes the continuation. The filename and ordering are fixed by schema 1; its SHA256 is in metadata.
- `.pices.manifest`: complete-set marker replicated on every rank. All markers contain the same SHA256 of the ordered per-rank metadata hashes (65-byte hex plus NUL records).

Payload/velocity/metadata are written under temporary names and flushed. After every rank completes writing, ranks rename the data files, synchronize, construct the metadata-set digest, and publish completion markers. Missing markers or a mixed rank set fail closed. This is an application-level completion protocol, not a guarantee against storage hardware failure after fsync. Do not manually create completion markers.

Preflight verifies the complete set, file hashes, legacy header, counts, exact lengths and property schema before allocating tracer data. Full validation compares canonical metadata against the restored state and current physics. Edited formatting, unknown schema, another build commit, changed partition, normalization, geometry or supported physical parameters are rejected. Cross-endian/ABI and repartitioned restarts are not supported. Old PG checkpoints without PICES metadata are rejected; there is no implicit grid-to-particle conversion.

Restart restores Tp rather than initializing it from grid temperature, restores the accepted nodal velocity, reconstructs host elements and diffusion caches, and restores total_timesteps=step+1. Zero-flavor restart skips legacy flavor counting, which otherwise incorrectly treats Tp as a flavor. P2 checkpoints currently require uniform, temperature-independent Newtonian viscosity with viscosity rebuilds enabled, in addition to P1 constant thermal physics. The metadata fingerprint includes geometry, local capacity/stiffness, transport scaling, boundary data, viscosity and selected solver controls. No reseeding/merging or particle mass field is introduced.

## Required preconditioner correction

The existing element-based CG assembly added the new diagonal to the previous inverse diagonal without clearing it. A fresh process and a continuing process therefore built different BI arrays, although U, P, F, BPI and viscosity were identical. `Construct_arrays.c/construct_elt_ks` now clears BI before assembly only when PICES is enabled. This makes the preconditioner independent of solve history and allows exact continuation. The PG branch retains its historical arithmetic and matches the P0 output regression. P1 PICES velocities may change with this correction; the relevant comparison is continuous versus restarted runs using the same new executable.

## Local validation

See `P2_LOCAL_RESULT.json`.

- Exact HPC configurations: one four-step continuous run, an independent two-step run, and continuation from its step-2 checkpoint to step 4. All 216 overlapping decoded outputs match. All live binary arrays (excluding unused index zero), including double Tp, T, coordinates, U/P and accepted velocity, match exactly.
- Stokes steps converge under the existing CG OR stopping criterion in 133, 247, 86 and 121 iterations. Maximum CFL is 0.0359454; every output retains 98,304 particles.
- Three prescribed rotation axes, 32 steps, 24,576 labeled particles each: 13,688 / 12,028 / 11,998 particles change owner rank; all 12 origin ranks participate. Persistent Tp labels remain bitwise unchanged through migration. Maximum position errors against exact rigid rotation are 0.01070 / 0.01226 / 0.01234 shell-radius units on this coarse mesh; these are diagnostics, not a spatial-order validation.
- Rotation checkpoint/restart yields bitwise identical final particle snapshots. A separate coupled flavor=1 case preserves slot 0 flavor and slot 1 Tp, with 108 matching decoded outputs and identical live checkpoint state.
- MPI abort-72 rejection tests: changed dt/physics, mismatched flavor layout, missing metadata, missing completion manifest, corrupt main payload, corrupt accepted velocity, and wrong Tp slot after re-signing the metadata manifest. No thermal step is accepted in these failure cases.
- All nine P1 manufactured/guard tests pass. PG's 36 decoded velocity/temperature outputs match P0 exactly, and the P0 native heat-flux/budget audit passes. CBF kernel and P0 auditor unit tests pass. Embedded SHA256 matches Python hashlib on empty, standard, padding-boundary and large binary inputs.

Local build uses the real standalone C parser/operators and an MPI test driver for prescribed fields. Autoconf generation and LSF shell syntax were checked. HPC Autotools/Python 2 installation, actual Intel MPI restart equivalence, and Pyre runtime remain unverified until separately exercised. The compile-time commit records HEAD; local working-tree tests are not represented as tests of a committed HPC binary.

## Reproduction

Use fresh temporary output directories:

```bash
MPICC=/usr/local/bin/mpicc python3 tests/cbf/build_validation.py --build-dir /tmp/pices-p2-build-new
python3 tests/pices/build_manufactured.py /tmp/pices-p2-build-new
python3 tests/pices/run_p2_local.py /tmp/pices-p2-build-new ../runs-cmbhf_EBA/cmbhf_EBA_PICES_P1.cfg --output /tmp/pices-p2-seams-new
python3 tests/pices/run_p2_smoke.py /tmp/pices-p2-build-new ../runs-cmbhf_EBA --output /tmp/pices-p2-smoke-new
```

HPC uses the direct P2 LSF in runs. Download the resulting archive and LSF logs; run `python3 tests/pices/verify_p2.py JOB_DIRECTORY --summary result.json` without --local. Local success does not waive HPC audit. P3 physical extensions, old-PG conversion tooling, long-time/resolution convergence, arbitrary mesh decompositions and reseeding remain outside this test acceptance.
