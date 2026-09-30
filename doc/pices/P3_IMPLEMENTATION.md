# P3: variable thermal properties and EBA heat sources

2026-09-30. P2 HPC job 12247638 passed before this implementation. P3 is
implemented and locally validated; HPC acceptance remains pending. Enable with
`pices_eba=on`; the default remains the P1/P2 constant-physics path.

## Discrete heat equation

Each heat substep freezes the Gauss-point coefficients at its starting grid
field. The existing `thermal_transport_at_gp` evaluates k(T,r), reference rho
and Cp. `phase_change_state` supplies the same reduced-pressure phase model as
PG. The lumped capacity is the integral of
`Ceff = rho * (Cp + Tabs * sum(delta_s * dX/dT))` against each shape function.
The pressure source is `-rho*Tabs*sum(delta_s*dX/dr_pressure)*ur`, using the
accepted physical radial velocity. No grid advection or SUPG term is added.

Internal heating uses the existing PG element-average rho times Q0; adiabatic
and viscous heating reuse `CBF_heat_sources`, which wraps the production EBA
source functions. This does not call a CBF residual or enable legacy CBF output.
The heat-stage temperature view is restored immediately after that call.
Velocity and viscosity remain frozen over the outer step. The phase temperature
term enters capacity exactly once; the legacy latent multiplier is not used.

Initialization and restore build capacity without reading velocity. The source
calculation reads velocity only during heat advancement. A 9^3-node test caught
and now covers the otherwise uninitialized-velocity path.

Particle movement and P projection occur once per outer step. Heat substeps
retain fixed particle positions and the existing Tp/subgrid correction. For P3,
each element uses max_GP(k/Ceff) in tau=Le^2/kappa. No particle mass, reseeding,
repartitioned restart, assimilation, composition coupling or strict ALA is added.
The checkpoint viscosity remains uniform Newtonian. Legacy CBF and heat-flux
output remain explicitly rejected; P4 must implement its own derivative semantics.
`qvis_mode` must be 0; enriched and legacy time-dependent heat sources are rejected.

## Step controls and diagnostics

The diffusion bound remains 0.8 divided by the largest absolute element-row
sum / assembled lumped mass. A source/diffusion rate bound limits max nodal
|deltaT| to 0.01 nondimensional per heat substep. Trial steps are halved until
all GP relative k and Ceff changes are <=10%, and phase-fraction changes <=0.05.
The trial never changes Tp. Nonpositive/nonfinite capacity or conductivity aborts
immediately, including if encountered at a trial state. Invalid absolute
trial temperature retries; invalid Tp aborts. Limits are fixed for this stage,
with max 50 retries and the configured pices_max_substeps. Particle CFL<=0.25
remains mandatory. Failure does not publish a completed outer-step checkpoint.

`PICES_EBA` reports frozen-capacity storage, internal/adiabatic/viscous/phase
pressure contributions and an independently calculated Dirichlet boundary
reaction. Shared boundary-node reaction is apportioned by local/global lumped
mass. Their sum must close. This is a discrete stage ledger, not an exact
nonlinear enthalpy integral or a CBF heat flux. `PICES_STEP heat_energy` is the
sum of these storage increments; `remap_energy` uses pre-advection frozen
capacity. Neither demonstrates physical particle energy conservation.

## Restart

P3 fingerprints static physical inputs, Cp, phase parameters and mesh rather
than temperature-dependent cached K/M. k, capacity and relaxation caches are
rebuilt after restore and each heat stage. P2 retains its former fingerprint.
The existing canonical metadata, compiled solver commit, payload/state SHA256
and collective completion manifest remain mandatory. Only the same executable
commit, native ABI and MPI partition can resume a checkpoint.

## Local evidence

- Exact P3 HPC configurations, 12 ranks, continuous 4 / independent 2 / restart
  2-to-4: 216 decoded outputs and all effective binary state arrays identical.
  Stokes iterations 148/182/138/169 satisfy the existing stopping criterion;
  maximum CFL=0.0388426140. All four heat-source contributions are nonzero and
  independent boundary-reaction ledgers close. The LSF gates were executed.
- Uniform Q/Cp=1 without diffusion: multi-substep solution at dt=0.025 has
  maximum free-node temperature error 2.220446049250313e-16.
- Exact spherical steady nonlinear k(T) solution with radially varying rho/Cp:
  one full remap + heat step at 5^3 and 9^3 nodes/cap gives RMS deviations
  0.01848780724 and 0.005676359998, ratio 3.25698. This is a two-grid consistency
  check of the combined method, not an established convergence order or a
  long-time phase-boundary benchmark.
- Three phase-kernel tests cover capacity and signed radial source against
  independent finite-difference derivatives, zero entropy, and initialization
  without allocated velocity. Three integration guards reject invalid capacity,
  exhausted heat substeps and composition-dependent conductivity.
- Four intact-checkpoint tests reject changed conductivity, phase entropy,
  internal source and reference Cp before any accepted thermal step.
- P1 nine-test integration suite passes (independent FE/P/Q particle error
  4.44e-16); P2 exact configurations still give 216 identical output comparisons;
  six existing PG/CBF kernel tests pass. PG formula files are unchanged.

Machine-readable results are in P3_LOCAL_RESULT.json. Local Open MPI is not an
Intel MPI/HPC acceptance result. Production material values, long integrations,
phase-width/mesh convergence and calibrated timestep tolerances remain separate
validation work; the P3 cfg is a small nonzero-physics acceptance fixture.

## Reproduce

```bash
MPICC=/usr/local/bin/mpicc python3 tests/cbf/build_validation.py --build-dir /tmp/p3-build
python3 tests/pices/build_manufactured.py /tmp/p3-build
python3 -m unittest discover -s tests/pices -p test_p3_phase.py
python3 tests/pices/run_p3_local.py /tmp/p3-build ../runs-cmbhf_EBA/cmbhf_EBA_PICES_P2.cfg --output /tmp/p3-physics
python3 tests/pices/run_p2_smoke.py /tmp/p3-build ../runs-cmbhf_EBA --stage P3 --output /tmp/p3-smoke
python3 tests/pices/run_p3_restart_guards.py /tmp/p3-build /tmp/p3-smoke --output /tmp/p3-guards
```
