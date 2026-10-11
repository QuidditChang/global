# Composition-aware conductive pre-initialization

This opt-in mode replaces the prescribed basal `erfc` anomaly with a finite-time,
three-dimensional conductive evolution of the **actual initial composition**.
It initializes a state at a chosen thermal age; it is not a steady-state solver,
and is not 1000 Myr of mantle convection or plate reconstruction.

## Configuration

In Pyre `[CitcomS.solver.param]` (or the same bare keys in standalone input):

```ini
thermal_preage = on
thermal_preage_Ma = 1000.0
thermal_preage_max_dt_Ma = 5.0
bottom_tbl_thickness = 0.0
```

Defaults are `off`, `1000.0`, and `5.0`. Age must be finite and nonnegative;
maximum step must be finite and positive. `thermal_preage_Ma=0` constructs the
reference background and then normal initial shallow structure without thermal
pre-evolution. `off` retains the historical initialization path.

The values above illustrate the parameter interface, not a selected production
thermal age. The current ASCII production cfg remains unchanged with this mode
off. The separate runs preparation tool requires an explicit total age and
step cap; physical age scans and timestep refinement must be kept separate.

The initial implementation requires fresh full-sphere PICES EBA, one cap per
MPI rank, `lith_age=on`, a thermodynamic reference-temperature file, and physical
radial Dirichlet values 0/1. Existing PICES restrictions still apply. File-based
`tic_method=-1` and P5 benchmark overrides cannot be combined with this mode.
Full checkpoint restart bypasses fresh initialization and does not repeat it.

The old `bottom_tbl_thickness`/`bottom_tbl_diffusivity_ratio` do not contribute
in this mode. No TBL upper cutoff is imposed. Setting them to 0/1 in the example
cfg makes this explicit. Duration is independent of `start_age` and `Q0`.

## Physics and order

1. Initialize material, tracer positions/flavors, and elemental/nodal composition
   by the existing production path.
2. Initialize the total nondimensional nodal temperature from fixed `Tref(r)`.
3. Conduct for the specified **separate thermal age**, holding composition fixed.
   Conductivity uses the same production elemental `kC(Cprim)`, depth dependence,
   temperature dependence, spherical geometry, and quadrature. Heat capacity
   includes the existing temperature-dependent latent-phase contribution.
4. Impose CMB temperature 1 and temporary surface temperature
   `refstate.temperature_surface` in an isolated boundary view. The physical
   `Ttop`, `Tbottom`, normalization, boundary arrays and node flags are unchanged.
5. No velocity, advection, radiogenic `Q0`, viscous/adiabatic/phase-pressure source,
   reconstruction, tracer transport/reclassification, or plate-age assimilation
   is advanced in this stage. Production clocks/counters are not modified.
6. Apply the existing initial-age shallow HSC anomaly and initial-slab overlay
   once. Restore the real radial boundary values through the usual BC routine.
7. Compute `DataT` and interpolate the final nodal temperature to PICES particles.
   Normal production Q0, phase physics, velocities and assimilation then resume.

The temporary hot surface avoids introducing a surface temperature jump to
physical `Ttop` throughout the thermal pre-age. It does **not** freeze the rest
of the mantle: `Tref` is not generally a conductive equilibrium and can diffuse.
The final shallow overlay may have its own adjustment when production starts.
There is no claim of a globally steady initial state or an exact basal-only heat
signal. A threshold on `T-Tref` includes background diffusion, not only CMB heat.

## Numerical method

`lib/Thermal_preage.c` owns temporary scalar arrays. Each step freezes the
production lumped effective heat capacity at the accepted old temperature and
solves backward Euler with new-temperature conductivity. Picard iterations
rebuild conductivity; preconditioned CG solves the SPD mass-plus-diffusion
system using homogeneous Dirichlet corrections. Shared-node residuals and mass
are exchanged; global inner products divide by shared-node multiplicity.
Physical radial boundaries use the **global** radial index, including radial
MPI decomposition.

The method is first-order in time, including the lagged latent capacity. It is
not an exact nonlinear enthalpy discretization. A trial is retried at half the
step if the nonlinear solve fails, conductivity changes by over 10%, effective
capacity by over 5%, or any phase fraction by over 0.05. The last step ends at
the requested thermal age. Temperature and true nonlinear residual checks are
independent of the production time-step state.

The configured step cap is an accuracy control, not just a stability limit.
Before interpreting production science, compare e.g. 5, 2.5, and 1.25 Ma at fixed
thermal age, mesh, composition realization and physics. Compare temperature,
basal flux and any diagnosed TBL thickness, then choose a tolerance appropriate
to the scientific question. Also check spatial and particle/composition
resolution. Cost on the 384-rank production mesh is not established by local
small-mesh tests.

## Diagnostics

Each rank's usual log receives `THERMAL_PREAGE_BEGIN`, periodic
`THERMAL_PREAGE_PROGRESS`, and `THERMAL_PREAGE_END`. They record thermal age,
accepted/rejected steps, CG/Picard iterations, residual, and unchanged production
time. `frozen_storage` and `boundary_input` form the numerical frozen-capacity
backward-Euler ledger. `balance` tests that discrete ledger; it is not exact
nonlinear thermodynamic enthalpy conservation.

No TBL-thickness threshold is hardcoded. Temperature profiles / field output can
be used to diagnose an explicitly defined anomaly threshold after initialization.

## Build and tests

Rebuild both C and Pyre components after adding parameters; source ABI changed.
The supported production `config_script` regenerates Autotools inputs from
`lib/Makefile.am`. The existing LSF script's clean committed build checks remain
unchanged; a local patch is not a submitted or production-ready binary.

Local standalone build:

```sh
MPICC=mpicc python tests/cbf/build_validation.py --build-dir /tmp/preage-build --optimization O2
```

See `tests/thermal_preage/` for real MPI initialization/state-isolation and
convergence tests. `tests/thermal_assimilation/test_ta.py` covers preservation of
the lifted shallow initializer and existing assimilation behavior.
