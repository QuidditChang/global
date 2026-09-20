# EBA thermal log and CBF audit

Scope: `cmbhf_EBA`, `THERMAL_BUDGET` records in rank 0's log (`global.log`
for the production output layout; `DATA/0/log` in the gzip smoke fixture).
These are heat **rates**, not stored energies. Other log sections, including
`MECHANICAL_POWER`, retain their existing normalization and state labels.

## Output contract (schema 2)

`output_time()` evaluates one record for the initial/restart state and every
completed output timestep, after temperature, tracer and Stokes updates. Both
standalone and Pyre use this C entry point. Native CBF files retain their existing
frequency and boundary output switches. No new cfg setting is required.

Each record includes the timestep, nondimensional time, state, conversion
scales, source rows, and `q_surf` / `q_botm`. All ranks participate in reductions;
only rank 0 prints. CBF computes both Dirichlet boundaries every step even when
native files are disabled. Unsupported non-Dirichlet boundaries are `nan`, not
zero. The diagnostic always includes physical advection; native CBF retains
its existing check requiring `CBF_use_advection=on`.

| Rows | TOTAL_W | MAX / MIN |
|---|---|---|
| Volume sources (`Q...`) | Volume integral, W | Element volume-average source density, W/m³ |
| `q_surf`, `q_botm` | Boundary surface integral, W | Native boundary nodal flux, W/m² |

This deliberately changes the old unlabeled nondimensional thermal table to
SI units. Parsers should use `schema=2` and the new column labels. Numeric rows
retain 17 significant digits. Divide TOTAL_W by 1e12 for TW.

Let L be reference radius in metres, and ΔT the reference temperature contrast:
volume-integral scale = k0 ΔT L; volume-density scale = k0 ΔT/L²;
boundary-flux scale = k0 ΔT/L. CBF multiplication promotes to double before
multiplying the float-valued model constants, avoiding an extra float rounding.

`q_surf > 0` leaves the mantle; `q_botm > 0` enters it from the core. Each physical
boundary face is integrated once with its local GLL weights; shared-node copies
are not added as separate physical area. Extremes are global nodal extrema.
The logged total and extrema match native q-file data at saved steps.

## Source audit

| Term | Meaning in the discrete thermal equation | Audit / correction |
|---|---|---|
| Qinternal | Element mean reference density × Q0, with tracer-enrichment mixing when enabled | Same density and composition formula as the residual; correct volume integral |
| Qvisc | Actual viscous heating used by the residual | Recomputed from current output velocity/viscosity; includes applied cap in mode 2 |
| Qvisc_raw | Viscous source before cap | Diagnostic only; must not be added again |
| Qvisc_capped | Source that would remain after cap | Diagnostic in mode 1, applied source in mode 2; must not be added again |
| Qvisc_removed | Raw minus actually used | Zero in diagnostic-only mode; reported separately from sources |
| Qadi | Adiabatic term stored on the left-hand side | Positive is cooling; enters Qtotal with minus sign |
| Qadi_base | Base adiabatic term | Same as Qadi in this branch; diagnostic duplicate, not another source |
| Qphase | Latent/phase term stored on the left-hand side | Positive is a sink; recomputed from output T, solver Tdot and velocity using the production kernel |
| Qassim | Imposed lithosphere temperature increment divided by accepted dt, weighted by rho Cp | Signed input/removal; rejected attempts reset increments; capacity now uses separately interpolated rho and Cp as in the residual |
| Qvisc-Qadi_base | Difference of two existing terms | Diagnostic only |
| Qtotal | Qinternal + Qvisc − Qadi − Qphase + Qassim | No boundary term or duplicate diagnostic added |

Previously the log was printed inside the temperature step, before the updated
Stokes solution, and phase diagnostics could come from the previous corrector
iterate. It now evaluates sources at the same output state used by CBF. Phase
values are averaged with `dOmega × Gauss weight`, replacing the arithmetic
Gauss-point average, so multiplication by the element volume recovers the
quadrature integral. The phase contribution in the solver's residual itself is
unchanged. Temporary CBF source pointers and geometry caches are restored.
Di=0 explicitly zeroes viscous and adiabatic source diagnostics.

## Limits of interpretation

- The production viscous/adiabatic source formulas are element approximations;
  this audit makes the log match those formulas, not a new higher-order physical
  source discretization. Phase extrema are element averages, not GP extrema.
- `Qassim` is an accepted-step average; other sources and CBF are output-state
  diagnostics. Stored solver Tdot is not recomputed after filtering/assimilation
  or the later Stokes update. Initial/restart Tdot may be initialized to zero.
- `NET_SOURCE_PLUS_BOUNDARY_W = Qtotal + q_botm − q_surf` is explicitly **not** a
  measured storage derivative or a numerical energy-closure residual. No U(t)
  differencing, filter/rheo temperature adjustment ledger, or complete advective
  boundary energy accounting is introduced by this change.
- Each step now requires CBF assembly/communication and a volume diagnostic
  residual evaluation. Native output steps also perform their existing CBF pass.
- The legacy `.heating` file is unchanged; its third value `heating_latent=1`
  is not this log's Qphase. Do not interpret that legacy multiplier as power.

## Validation

- Six kernel tests, including an independent nonuniform Gauss-weight phase
  integral and physical-RHS/state-restoration tests.
- Standalone build of all production library sources with the existing local
  validation driver (not the HPC Pyre/Intel production installation).
- 12-rank, 5³-node, two-step EBA fixture: CBF files every two steps. Logs at
  steps 0,1,2; step 1 has no native q files. Source sum and constant-Q0 SI scaling
  checked; saved-step CBF totals/extrema checked against native faces/nodes.
- 24-rank fixture with two radial partitions, nonzero phase entropy/Clapeyron
  slope and mode-2 viscous cap: the same checks pass; phase rates are nonzero
  and raw-minus-used viscous heat matches the removed-heat diagnostic.
- 12-rank Di=0 fixture with both native output switches off and CBF_frequency=0:
  step-0/1 boundary rows still present, no q files created, and viscous/adiabatic
  source rows exactly zero.
- Verification scripts: `tests/cbf/verify_thermal_budget.py` and
  `tests/cbf/verify_native_outputs.py`. The legacy standalone termination routine
  exits with code 8 even after normal completion; validation requires the expected
  completed steps and all numerical checks, not exit-code acceptance alone.

Example validation commands after building and running the fixture:

```sh
python3 tests/cbf/test_cbf_kernel.py
python3 tests/cbf/verify_thermal_budget.py /tmp/CBF-budget-smoke --steps 0 1 2 --constant-q0 .3
python3 tests/cbf/verify_native_outputs.py /tmp/CBF-budget-smoke --step 2
```
