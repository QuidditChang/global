# EBA fixed-volume reconstruction experiment (stage three)

Status: opt-in research prototype. **The local C24 shape gate failed.** The
ordinary solver defaults to `off`; the submitted experiment defaults to shadow
diagnosis only, with no active composition feedback.
Conservation, boundedness and successful restart are necessary checks; they do
not prove that the reconstructed field follows the correct local transport
equation. A full-resolution RH20 result must pass its own scientific gates.

## Particle budget and cell capacity

At fresh initialization, every particle in element e receives the immutable
double-precision weight `w_p = V_e / N_e`, using eight-point FE element volume.
Initialization rejects empty elements. The weight travels through the existing
extra-quantity migration and checkpoint machinery. It is never reassigned when
a particle enters a different element.

For every flavor, including the background, the target is its original weighted
volume plus accumulated classification sources. A relabeling contributes
`-w_p` to the old flavor and `+w_p` to the new flavor. C24 relabeling is refused.
Thus reacting flavors need a source budget, not individually constant volumes.
The initial measured C24 volume is retained; the algorithm does not force it to
the nominal 4 percent.

Raw cell inventories are `M[e,f] = sum(w_p)` for particles of flavor f in cell e.
Their total `q[e]` can exceed or fall below the actual cell capacity `V[e]`.
The reconstruction retains the raw mixture up to `min(q,V)`, then transfers
only surplus into deficits. The nonnegative transport matrix has the form
`W[d,s] = a[d] K[d,s] b[s]`, with prescribed deficit and surplus marginals.
The same matrix transports all 25 flavors, without separate global rescaling.

`K = diag(V) H^hops`, where H is positive cell-to-shared-node-to-cell averaging
with double-precision FE integrals. It is applied without a stored global graph;
existing double MPI node exchange handles cap and radial partition boundaries.
Alternating row/column scaling solves the capacity constraints. Only after
convergence is floating-point row closure permitted; both its size and the
resulting species-budget discrepancy are recorded.

This maximum-retention, surplus-only construction can have insufficient support:
neutral cells cannot relay capacity. Increasing iterations cannot guarantee a
feasible transport matrix. The code never silently expands the support or
normalizes global species volumes. Shadow mode records failure; active mode
stops before accepting an infeasible field.

## Field projection, switches and restart

The active cell field is the corrected inventory divided by cell volume.
Nodal composition uses a positive double-precision lumped-volume projection.
Thermal projection weights, thermal operators and the ordinary off-mode
composition mapper are unchanged. Active nodal values need not be bitwise equal
to the legacy single-precision mapping, even at initialization. The present
gates check nodal bounds and integrated volume; they do not establish equality
of every local buoyancy or heat-source contribution.

The native and Pyre options are:

| Option | Solver default | RH20 stage-three plan |
| --- | --- | --- |
| `pices_volume_mode` | `off` | `shadow` for controls/E; `active` for F/G |
| `pices_volume_hops` | 2 | 4 |
| `pices_volume_iterations` | 1000 | 1000 |
| `pices_volume_tolerance` | 1e-12 | 1e-12 |

The capacity tolerance is relative to each element volume. The species budget
tolerance is `1e-10` times total shell volume. Fraction/partition checks use
`1e-12`. These bookkeeping tolerances are not physical accuracy targets.

Modes shadow and active use authenticated checkpoint schema 7, with a named
`volume_weight` extra quantity, the 25 original budgets and the 25 cumulative
sources. Weightless checkpoints and mismatched reconstruction settings are
rejected. The off mode preserves its previous schema and physical payload.
Restarts recompute the field from original weights and retain the original
budget history; they do not establish a new zero-drift baseline.

## Controlled and coupled experiment

All experiment launchers and reports live in the `runs` repository on branch
`cmbhf_EBA`, locally `worktrees/runs-cmbhf_EBA`. Shared cluster helpers remain in
`scripts`. Use `PICES.RH20_volume_conserve.py`, its `.cfg` and `.lsf`, and
`PostProc.volume_conserve.py` there. The source RH20 archive is read-only.

One 400-slot medium allocation runs 384 solver ranks, with 16 unused solver
slots; no online CBF/VTU pipeline is started. The implemented full sequence is:

1. A static; B exact and numerical transport; C half and quarter timesteps;
   D exact forward-return. All share one initial snapshot. The expanded flow
   uses 152 nonoverlapping compact vortices at radius 0.64, support radius 0.08.
   Their minimum center spacing is about 0.160453, greater than 0.16. The report
   requires at least 10 percent of weighted initial C24 inside the active flow.
2. E: original RH20 physics, with conservative reconstruction observed in shadow,
   from step 0 through 100.
3. F: active reconstruction, from step 0 through 100, compared with E at common
   physical times without extrapolation.
4. G: active 0–50 and restart 50–100, compared with the continuous F checkpoint.

The current explicit `diagnostic_only = true` plan runs only items 1 and 2.
After a successful controls MPI execution, a recorded scientific-gate failure
does not prevent the independent shadow observation of original RH20 physics.
The failure is retained, and the suite returns nonzero if either scientific
report fails. MPI failures, missing reports or corrupt certificates still stop
execution. F/G are not scheduled in this mode. They remain implemented for
future algorithm validation; switching off this restriction is not recommended
on the present local evidence.

A–C share the same final physical time; D lasts twice as long. C24 shape uses
eight-point quadrature of the inverse analytic image of the initial layer.
This is an approximate geometric reference, not an exact cut-cell integral.
The declared gates allow 10 percent L1 degradation relative to the larger of
the old-method error and initial quadrature floor, plus `1e-8 * shell_volume`;
the centroid allowance is a quarter mean cell width. D additionally checks
recovery of the actual initial cell field, avoiding geometric-reference error.
The redistribution RMS distance is reported, not hidden behind budget closure.

All accepted steps and 25 flavor rows must be present and feasible. Pending-only
rows cannot certify completion. The original ratio-specific audit is used only
in controls and E. F/G use the weighted budget. In active diagnostics,
`original_volume_nd` denotes the pre-update conservative field; it is not a
fresh evaluation of the legacy ratio reconstruction.

## Local verification and practical limits

Reproducible tests are under `tests/pices/`:

- `check_composition_volume.py`: real C kernel on small synthetic MPI geometries,
  including positivity, weighted reactions, empty destinations, unreachable
  capacity, double projection and shadow isolation.
- `check_volume_checkpoint_config.py`: actual configuration/checkpoint C code,
  schema compatibility and authenticated budget persistence.
- `check_volume_conserve_integration.py`: full standalone EBA solver, controls,
  coupled short segments and physical restart comparison. Native fixtures are
  manufactured inputs, not a substitute for the archived RH20 grid.
- `check_volume_conserve_regression.py`: complete independent build of the
  committed baseline and cross-version physical checkpoint comparison.
- `check_volume_shadow_regression.py`: same-build off/shadow comparison, with
  field-by-field native and persistent-state parsing. Only the new named weight,
  its count and authenticated ledger are excluded from the comparison.

The native controls include an initially motionless radial element next to the
CMB. Without this buffer, a very coarse FE velocity interpolation reaches the
boundary cushion despite compact analytic support, invalidating the no-boundary-
projection test. The solver refuses such a test instead of silently clipping.
Inadequate initial tracer sampling likewise fails rather than fabricating
particle weights in empty initial cells.

The initial two-hop coupled test did not converge within 1000 iterations;
10000 iterations still failed its tolerance. A four-hop, 12-rank coupled test
passed two steps (at most 492 iterations). This is the explicit reason for the
four-hop experiment setting, not proof of full-resolution feasibility. A larger
support also permits longer redistribution; quality must be assessed alongside
conservation. Iterative node exchanges and the extra cell/node flavor arrays
add measurable runtime and memory. Each particle adds one double weight.

The disabled-method regression against commit
`8234a71c753e9483602739be0a3a87a186689467` passed: all 96 cross-version physical
checkpoint payloads were identical, ignoring only unused native index-zero
storage; auxiliary state bytes were identical. Real 12-rank active 0–1–2
restart also matched continuous 0–2 in all 24 physical payloads, and passed the
complete restart report with original budgets and source history intact.
The complete off/shadow test also passed: 72 payloads at steps 0, 1 and 2 across
12 ranks, including 1728 main physical fields and 468 original persistent-state
fields, were bitwise equal. All 53 input files and checkpoint publications were
verified. This manufactured RH20 fixture uses rheol=1; it does not certify the
full-resolution rheol=7 HPC model.

### Recorded research outcomes, 2026-11-13

The archived production radii are `ri=3480/6371` and
`r_C24=0.5813265437084149`. Exact sphere-intersection geometry gives only
9.99158 percent C24 coverage for 150 vortices. The final design therefore uses
152, giving about 10.1248 percent analytic coverage without overlapping supports;
the actual particle-weighted coverage still has to pass the unchanged 10 percent
gate. The earlier 150-vortex local results below remain explicitly identified.

The complete 9-node, 12-rank, 150-vortex fixture covered 10.5486 percent of
initial weighted C24. All six controls completed with feasible reconstruction
at every step; maximum iterations were 816, maximum species error
`3.849e-16 * shell_volume`, and maximum node-integral error
`7.698e-16 * shell_volume`. D returned to the actual initial C24 field with
absolute L1 error `5.546e-18` nondimensional volume.

| Control | New C24 shape L1 relative to original ratio method | 10% gate |
| --- | ---: | --- |
| B exact | +7.6973% | pass |
| B numerical | +25.5282% | fail |
| C half timestep | +25.5286% | fail |
| C quarter timestep | +25.5286% | fail |

The lack of improvement under timestep refinement is evidence of a spatial
transport/reconstruction limitation, rather than a time-step-only error. The
algorithm conserves global inventory but does not minimize interface error:
it distributes donor mixtures to multiple capacity deficits. This result must
not be described as successful physical validation or repaired production
volume drift. The 10% gate was not relaxed after observing this result.

A 17-node fixture with 64 nominal tracers per element initialized successfully
and covered 10.7862 percent of initial weighted C24. Static reconstruction passed
116 steps; the first moving step with 4 hops and 1000 iterations had row error
`6.006e-9`, above the unchanged `1e-12` tolerance. This establishes nonconvergence
within the iteration budget, not mathematical infeasibility. An 8-hop probe
converged for one step in 692 iterations, but never completed the six controls
or shape gates; it is diagnostic evidence only and was not adopted as a fix.
A separate 4-hop/10000 probe failed while writing initial snapshots because of
local disk capacity; it supplies no convergence evidence.

Full 12- and 24-rank two-step coupled active fixtures passed their volume and
checkpoint reports. Their success does not override the failed shape and finer
grid controls. Subsequent algorithm work should address local/interface
transport fidelity, rather than tuning support or tolerances until a budget
test alone passes.

The final **152-vortex** design was also run through all six controls on the
9-node, 12-rank fixture. Weighted C24 participation was 10.6656 percent; every
reconstruction was feasible, at most 932 iterations, with maximum species error
`3.849e-16 * shell_volume`. D again restored the initial field to `5.546e-18`
absolute L1. B exact shape degradation was 7.9519 percent, B numerical and C half
were 25.3043 percent, and C quarter was 25.3030 percent. The final report therefore
correctly returned scientific FAIL for the three numerical-transport shape
gates. This confirms the diagnostic-only choice for the final submitted design.

Final runs-side preparation, migration, reporting and launcher checks comprised
84 passing Python tests, plus successful syntax and full native C build checks.
They are software checks, distinct from the failed scientific shape validation.
