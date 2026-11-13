# EBA species-volume diagnosis, stage one

This opt-in observer targets the RH20 EBA/PICES ratio-composition experiment.
It changes neither particle motion, geological relabeling, composition,
temperature, heating coefficients, nor checkpoint representation. It does not
repair volume drift. All flavors 0 through 24 are accounted for; spatial
diagnostics focus on C24.

Enable `PICES_VOLUME_AUDIT=1`; use `PICES_VOLUME_AUDIT_INTERVAL=50` for spatial
summaries every 50 steps. Default is disabled and creates no diagnostic files.
The supported configuration is PICES EBA, ratio composition, 25 flavors,
24 stored components, and the ordinary Stokes solver. Unsupported configurations
are rejected when the observer is explicitly enabled.

## Ledger and lifecycle

Owned-element finite-element quadrature supplies geometric volumes and cached
local nodal integration weights. Shared nodes are integrated through their
owned-element contributions, never summed as duplicate whole nodal volumes.
Flavor zero is the complement of stored components 1–24.

- A: before movement, integrate the actual element composition and count fields.
- B: after the final particle corrector and MPI migration, count new hosts with
  old labels and evaluate the ratio without writing the composition. Empty
  cells retain the old composition, matching the real reconstruction.
- S: within the existing unique relabeling pass, record entry-to-exit label
  transitions. The volume-equivalent contribution is `V_element/N_element`.
  This is an accounting weight, not a fixed physical particle volume.
- C: after the real composition reconstruction, integrate actual `comp_el`.
- D: integrate actual `comp_node` independently after its mapping.

For each species, `C - A = (B - A) + S + residual`, where
`residual = C - B - S`. C24 has no classification source; any transition
involving its label is flagged. Movement/reconstruction includes sampling,
trajectory integration and velocity interpolation effects and cannot by itself
separate those causes.

Pending rows are written after composition reconstruction. Accepted rows are
published only after the normal Stokes solve and rotation removal return.
The existing PICES thermal `accepted_step` log alone is not this acceptance
marker. Auxiliary CSVs are flushed before the accepted global table. Readers
must require all 25 unique flavor rows before accepting a step.

Initialization records an explicit baseline at the actual start step. A restored
checkpoint uses its own diagnostic-start baseline; it cannot recover unrecorded
earlier history. Existing diagnostic files are refused rather than silently
overwritten or appended. The RH20 from-zero launcher therefore uses a fresh run
directory, preserving the original experiment.

## Outputs and numerical precision

Only rank zero writes `Diagnostic/VolumeAudit/`:

- `metadata.json`: schema, baseline, MPI size, quadrature/precision conventions.
- `globals.csv`: per-step all-flavor counts, A/B/C/D integrals, conversion sources
  and closure residuals, with pending/accepted status.
- `transitions.csv`: nonzero entry-to-exit classification transfers.
- `radial.csv`: C24 volume, particle sampling and empty volume by radial layer.
- `bins.csv`: occupancy transitions, cell-size, particle-count and density bins.
- `samples.csv`: up to ten globally selected cells with largest absolute C24
  volume changes. Local changes can represent genuine transport; these samples
  are localization evidence, not individually proven conservation errors.

Volume accumulation is double precision and species particle counts are 64-bit
integers. Element budget closure uses a diagnostic threshold of `1e-10` times
shell volume. Existing node-mapping weights (`TWW`, `MASS`) are single precision;
node-integral comparisons use a separate `1e-7` shell-volume threshold and node
fraction bounds use `1e-7`. Raw offsets and bound violations are retained. These
thresholds do not define an acceptable physical C24 drift and do not relax
particle-host validity, Stokes convergence, or existing solver checks.

The EBA runs repository supplies `PostProc.volume_audit.py RUNROOT` to generate
accepted-step CSV/JSON summaries and C24/radial figures. Diagnostic execution
time is reported separately; production overhead must be measured on HPC.
The local checkout is `../runs-cmbhf_EBA` relative to this solver worktree;
the HPC checkout is `/work/ess-zhangkd/CitcomS/runs` on branch `cmbhf_EBA`.
Its `PICES.RH20_volume_audit.md` contains the preparation and submission recipe.

## Validation commands

```bash
python3 tests/pices/check_volume_audit.py
MPICC=/usr/local/bin/mpicc python3 tests/cbf/build_validation.py \
  --build-dir /tmp/rh20-volume-audit-build --optimization O1
MPIEXEC=/usr/local/bin/mpiexec python3 tests/pices/check_volume_audit_integration.py \
  /tmp/rh20-volume-audit-build ../runs-cmbhf_EBA /tmp/rh20-volume-audit-12
MPIEXEC=/usr/local/bin/mpiexec python3 tests/pices/check_volume_audit_integration.py \
  /tmp/rh20-volume-audit-build ../runs-cmbhf_EBA /tmp/rh20-volume-audit-24 \
  --radial-partitions 2
```

The core harness exercises unequal-volume cells, relabeling, empty/recovered
cells, shared-node integration, pending-only rows and bitwise read-only state.
The real MPI regression runs the production EBA/RH20 operators with the audit
off/on, compares physical checkpoint payloads (excluding unused index-zero
storage), checks all-species volume/count budgets and radial partition sums.
These small meshes do not replace the original 384-rank HPC/Pyre experiment.

Local validation on 2026-11-12 passed with the full production C source build:

| Test | Result |
| --- | --- |
| Actual audit module, 1/2 MPI ranks | Unequal-cell drift, conversions, empties, precision and read-only checks passed |
| Actual EBA/RH20 operators, 12 ranks, steps 0–2 | 48 physical checkpoint/state payloads identical off/on; maximum closure/shell `8.90e-17` |
| Actual EBA/RH20 operators, 24 ranks, two radial partitions, steps 0–2 | 96 payloads identical off/on; maximum closure/shell `7.05e-18`; all radial layers integrate to global totals |
| Existing final-host, shared-coordinate and restart metadata regressions | Passed |
| Scripts repository fresh launch and report tests | 8 + 10 passed; LSF shell syntax passed |

The 24-rank fixture retains the original five-node radial mesh/reference state
and sets `mgunitz=2,nprocz=2`. No interpolation of experimental reference data
is performed. HPC Intel/Pyre installation, original-resolution initialization
and production diagnostic overhead remain to be checked by the submitted run.
