# P8 production adaptation plan

Status: P8a HPC audit PASS, job 12280286, 2026-10-07.
P8b local implementation and validation PASS, 2026-10-09; HPC audit pending.
P8c remains planned. See P8A_IMPLEMENTATION.md and P8B_IMPLEMENTATION.md.
Branch: cmbhf_EBA_PICES. Target remains the exact P7 production snapshot.

## Scope and architecture

Preserve target physics: 25 flavors, immutable primordial flavor24, existing
background reclassification, B24=0.3, kC_ratio=0.8, qvis_mode=2, rheol7,
time-dependent plate velocities and TA. Temperature is an independent extraq
slot on the existing tracer population, not a second particle system.

Pices.c owns the accepted-step particle lifecycle through the shared
tracer_move_particles and tracer_update_composition functions. It reuses
thermal_transport_at_gp and thermal_heat_sources; CBF wraps the same source
function. The outer tracer_advection returns early for PICES, avoiding a
second move without a transient moving flag. Keep one implementation of
transport, reclassification, composition and physical heating.
Make ownership/order explicit without broadly refactoring the PG driver.

A previously omitted production gap is automatic timestep selection:
PICES validation requires fixed_timestep>0; the target omits this setting,
and the parser default is zero. Retain fixed steps for comparisons, and
add a global particle CFL constraint to automatic outer step selection.
Thermal subcycling remains separate. Forcing timestamps and accepted clocks
must be consistent; automatic steps must not skip forcing coverage/end time.

## Ordered delivery gates

1. P8a: tracer/composition and composition-dependent conductivity.
   Reuse migration and extraq serialization; ensure exactly one move and one
   reclassification per accepted step. Preserve flavor24, exclude reserved
   flavors18/19, and do not change Tp merely because a flavor changes.
   Reconstruct composition before the PICES heat-stage coefficient assembly;
   audit old/new time-level use explicitly rather than reproducing an
   accidental ordering. Keep PG semantics unchanged.
   Validate MPI migration of flavor+Tp, total particle count, primordial
   count, composition closure, buoyancy and actual kC=0.8 response.
   Include restart of composition+Tp in the supported constant-viscosity
   subset; expand fingerprints for all newly supported settings/state.

2. P8b: production flow/heat coupling and automatic timesteps.
   Share capped viscous heating and existing boundary-file interpolation.
   Define which accepted velocity/viscosity supplies each thermal stage;
   avoid mixing current temperature with stale or reconstructed viscosity
   silently. Add variable-step particle CFL and forcing-time checks.
   Validate rheol7 plus qvis_mode=2 and a changing plate/age forcing interval,
   then a three-level small-mesh multigrid gate. Fixed-step controls remain.
   The actual 384-rank multigrid mesh is validated in the integrated P8c pilot.

3. P8c: restart and integrated pilot.
   Inventory persistent versus derived fields, especially accepted viscosity,
   velocity, flavor, Tp, TA/CBF history and external forcing identity.
   Extend existing checkpoint schema/manifest, not a second restart stack.
   Rebuild viscosity only if continuation equivalence is demonstrated;
   otherwise persist the minimal accepted state. Do not relax the existing
   constant-viscosity guard until this gate passes.
   Test continuous/split runs across forcing changes using rheol7 and all
   target features. Then run a short pilot on the actual 384-rank mesh,
   measure memory/time/output costs and agree the production acceptance.

Each gate gets its own reviewable changes and direct LSF/cfg in runs,
commit/push commands, then user HPC execution and local audit before advancing.
No new shell wrapper. Reuse the ser 12-rank small-test setup where suitable;
actual 384-rank pilot needs a matching queue (medium under supplied limits).
Production paths use the agreed dated scratch convention; out/err remain
in each model/job directory.

## Numerical issue carried from P7

Velocity differences shrink with dt, but PG/PIC temperature RMS difference
grows to 11.52 K. PG is not an exact reference. Use one targeted spatial /
particle-density refinement and localized error/energy diagnostics to distinguish
projection, boundary and splitting errors; do not repeat the entire P0-P7
matrix. Production approval needs explicit error criteria supported by these
results, not an invented universal temperature tolerance. Bounded projection
alone does not establish physical-energy conservation or local monotonicity.

## Confirmed initialization

User confirmed fresh initialization at 249.9 Ma using the target cfg.
PG-to-PICES checkpoint conversion is outside this adaptation scope.
Normal PICES continuation/restart remains required for production.

Internal implementation uncertainties (resolve in code audit, not delegated
to user): time-level ordering of reclassification/TA/buoyancy, qvis inputs,
viscosity rebuild equivalence and complete checkpoint fingerprint coverage.
