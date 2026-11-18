# Rheology family validation — 2026-11-18

Implemented on `cmbhf_EBA_PICES`; companion inputs are in the runs worktree on
`cmbhf_EBA`. Mode 7 is the common Fig.1/EBA baseline, 71 applies positive layer
constants, and 72 applies positive layer quadratics. Mode 72 defaults preserve
the historical M2-A reference profile. See [model and input definitions](rheol7_nuref.md).

## Completed local checks

| Check | Result |
| --- | --- |
| Production C kernel, 65 archived M2-A nodes | PASS; mode 72 maximum relative error 3.775e-13 |
| Baseline against archived M2-A divided by its correction | PASS; maximum relative error 6.521e-13 |
| Identity/custom coefficients, interface ownership, thermal response, interval positivity | PASS |
| Production standalone parser and Pyre property bridge harness | PASS; exact four-value vectors, double precision, inactive fields ignored, malformed input rejected |
| Production dispatch and temperature interpolation | PASS for 7/71/72, including TDEPV=off |
| Two-rank MPI failure handling | PASS; invalid rheology inputs terminate the communicator |
| Production checkpoint fingerprint harness | PASS; family version, active coefficients and cold_scale affect identity; inactive coefficients do not |
| Full standalone solver build | PASS with local OpenMPI and O1 |
| Actual Python 2.6 header compile of module/setProperties.c | PASS; legacy unrelated Python/header format warnings remain |
| Real coupled PICES/pPv MPI smoke test | PASS for all three modes: 12 ranks, one accepted step per rank, checkpoint sidecars written |
| New Vis7 config normalization and LSF syntax | PASS |
| New LSF producer/inventory preflight harness | PASS; rejects stale producer or missing input inventory |
| Repository whitespace checks | PASS |

Commands from the solver worktree:

```bash
python3 tests/test_rheol7_nuref.py /Volumes/Pōwehi/CitcomS/output
python3 tests/test_rheol7_runtime.py
python3 tests/test_rheol7_temperature_clip.py
python3 tests/test_rheol7_checkpoint.py
python3 tests/cbf/build_validation.py \
  --build-dir /private/tmp/rheol-family-build-20261118 --optimization O1
python3 tests/test_rheol7_mpi.py \
  /private/tmp/rheol-family-build-20261118 \
  ../runs-cmbhf_EBA/cmbhf_EBA_PICES_P3.cfg \
  --output /private/tmp/rheol-family-mpi-20261118-r2
```

Choose a new output directory when repeating the MPI test. It uses the matching
local `/usr/local/bin/mpiexec`; set `MPIEXEC` to another compatible launcher if
needed. MPI tests need local socket access. Summary from this execution:
`/private/tmp/rheol-family-mpi-20261118-r2/summary.json`.
The smoke test enables pPv with eta_ratio=1, as in the new experiment; it checks
integration and checkpoint writing, not a non-unit pPv viscosity contrast or a
full restart cycle. Checkpoint identity is tested separately with the production
fingerprint function.

## New experiment and thermal-column diagnostic

The new input is
`cmbhf_EBA_Q0_30_Vis7_scold1.0_LLSVPsV0.04_B24_0.5_PICES_1.8Ga_pPv.cfg`.
Compared with the committed source template, only three identity/path fields
and four explicit coefficient vectors differ. All four vectors are inactive
for mode 7. The rheology change is the new definition of mode 7 itself.
It remains a fresh 1799.9 Ma initialization with cold_scale=1, Q0=30,
LLSVP volume 0.04, B24=0.5 and the original PICES/pPv settings.

```bash
python3 tests/validate_rheol7_initial_profile.py \
  /Volumes/Pōwehi/CitcomS \
  ../runs-cmbhf_EBA/cmbhf_EBA_Q0_30_Vis7_scold1.0_LLSVPsV0.04_B24_0.5_PICES_1.8Ga_pPv.cfg
```

The representative 70 Ma column passed the Kelvin formula, reference identity,
finite-positive and thermal-response checks at 65 nodes and 128 Gauss points.
It uses the actual full assimilation depth, **254.84 km**, and an approximately
500 km basal TBL at the 1% anomaly height. Raw surface and CMB viscosities are
2.51266e41 and 2.29369e19 Pa s; illustrative nodal limited values are 1e23 and
1e20 Pa s. Six Gauss points hit the minimum and seven hit the maximum.
Artifacts are under the parent workspace's
`output/rheol7_initial_validation/rheol_7/cold_scale_1/` followed by the config
stem. This diagnostic excludes pPv, composition, plasticity and filtering and
does not reconstruct the full 3D initial state.

## Migration and deployment boundary

The 49 tracked top-level EBA templates formerly selecting 7 now select 72,
preserving their old profiles and all other existing key values. Thirteen
historical P7/P8 fixtures retain their original text. Old checkpoints are not
automatically migrated: the family version and selected coefficients are part
of the new physics fingerprint.

No commit, push, HPC build, or production submission was performed for this
implementation. The local machine cannot execute the cluster's Linux Python
2.6 binary, so full Pyre installation/execution remains an HPC check. Rebuild
and install the solver and Python components together because the C viscosity
structure gained coefficient arrays. The companion runs
`RHEOL_FAMILY_RUNBOOK.md` gives the exact local staging and HPC commands.
