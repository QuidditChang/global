# PICES startup rheology flags

The downloaded run
`cmbhf_EBA_Q0_30_Vis7_scold1.0_LLSVPsV0.04_B24_0.5_PICES_1.8Ga`
records solver commit `c1af0a8699b71ab15103004c1feb25a5168f4a23` and aborts
before its first accepted step:

```text
PICES_ERROR rank=16 step=0: checkpoint requires supported Newtonian rheology and viscosity rebuilds
application called MPI_Abort(MPI_COMM_WORLD, 72) - process 16
```

This is at lines 3151-3152 of the downloaded root
`stderr_EBA_Q0_30_rheol7_scold1.0_TA.txt`. The configuration has rheol=7,
VISC_UPDATE=on, and SDEPV/PDEPV/CDEPV/low_visc_channel off. The Pyre defaults
select system viscosity and low_visc_wedge off. The postprocessing state has
empty done/failed collections and pPv_enabled=false; no completed output or
checkpoint is present in the downloaded subset.

`citcom_init` allocates `All_variables` with malloc. Three reserved viscosity
fields, FREEZE, weak_blobs and weak_zones, are read by the PICES checkpoint
guard but were never initialized by either input frontend or common defaults.
Nonzero allocator contents therefore reject an otherwise supported model.
The guard already accepts the 7/71/72 family; pPv_enabled is not a condition
in this guard. The old log does not identify which field was nonzero on rank
16, so the exact historical field value cannot be recovered from this log.

The common constructor now initializes these three fields to zero. All
unsupported-rheology checks remain active, and their failure message now
prints every checked value. No physical parameter, checkpoint schema or
fingerprint is changed.

## Verification

`tests/pices/check_checkpoint_rheology_init.py` compiles the production
constructor, coupled-checkpoint selector and validation function. An allocator
stub supplies nonzero bytes for the three reserved flags. Before the fix it
reproduces the same step-zero checkpoint error. After the fix it accepts
rheol=1/7/71/72 with pPv off/on and still rejects each of the eleven unsupported
conditions, including explicitly re-enabled reserved flags.

```bash
python3 -B tests/pices/check_checkpoint_rheology_init.py
python3 -B tests/test_rheol7_checkpoint.py
```

Both pass. The existing fingerprint test's extraction boundary was updated
to use the next code condition instead of an obsolete Stokes-policy comment.
The complete standalone C build also passes. A 12-rank rheol=7 small-grid
MPI run completes step 1 with pPv off, and the matching pPv-on control also
passes. Each writes 24 checkpoint metadata files (steps 0 and 1 on 12 ranks),
with no PICES_ERROR/fatal message. The local summary is
`/private/tmp/citcoms-checkpoint-init-mpi-20261118-approved/summary.json`.
These are standalone/OpenMPI checks, not an Intel MPI/Pyre or full-resolution
production run. HPC must rebuild and install the fix before a fresh rerun.
