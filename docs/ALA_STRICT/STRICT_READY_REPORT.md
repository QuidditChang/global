# Strict ALA Production-Readiness Report

## Status

`runs/cmbhf_ALA_strict.cfg` is the production strict-ALA configuration. Its
physical inputs remain aligned with `runs/cmbhf_ALA.cfg`; the active
strict-only differences are:

```text
refstate_file = refstate_ALA_strict.txt
ala_beta_element_source = interval
ala_beta_interval_file = interval_ALA_strict.txt
ala_shallow_patch_velocity_solver = diagonal
ala_geneo_preconditioner = on
ala_geneo_basis_type = radial_partition
```

The first three select the phase-inclusive strict reference state.  The final
three keys select the Stage 6d Schur-preconditioner experiment described
below.  No
energy-equation, phase-equation, boundary-condition, assimilation, timestep,
or viscosity parameter is changed.

## A. Reference-state source

`refstate_ALA_strict.txt` is generated from Katsura (2022) Table S5 by the
strict-only Python pipeline. Density is phase-inclusive. Density and beta are
evaluated directly at the 65 CitcomS radial nodes from the same piecewise
analytical `ln(rho)` model.

The serialized schema is:

```text
rho g Tref alpha Cp beta Gamma_eff Ks
```

The first seven fields retain the legacy nondimensionalization and ordering.
`Ks` is appended in physical GPa. It is stored by the reader but is not used by
the current Stokes or energy formulation.

## B. Tref source and usage

Column 3 is the fixed thermodynamic reference temperature. It comes from the
Katsura Table S5 `T_K` values through the same four-branch endpoint-constrained
cubic fitting and direct radial-node evaluation used for the other fitted
thermal properties.

The 300 K surface and 3800 K CMB total-temperature boundary values are not
written into Tref. With `Ttop=300 K` and `Tbottom=3800 K`, the generated endpoint
values are:

| Endpoint | Tref (K) | serialized Tref* |
|---|---:|---:|
| CMB | 2607.653828 | 0.6593296652 |
| surface | 1618.046307 | 0.3765846591 |

The reader stores column 3 in `E->refstate.temperature`; `E->refstate.Tref` is
a non-owning semantic alias to the same storage. The phase reference uses
`E->refstate.Tref`.

## C. Initial and dynamic temperature

`E->T` remains the sole total-temperature state and the energy-equation
unknown. The initialization path is inherited from `cmbhf_ALA`:

```text
Katsura Tref
  + bottom thermal boundary layer
  + surface half-space cooling / plate-age structure
  + slab and configured perturbations
  -> E->T(x, 0)
  -> total-temperature boundary conditions
```

No `DataT` field is allocated or cached. A temperature anomaly remains
well-defined when needed as the derived quantity `E->T-Tref`; it is not an
independent state variable.

Assimilation continues to update `E->T` through the existing
`assim_delta_T` path. It has no write path to `E->refstate.temperature` or
`E->refstate.Tref`.

## D. Phase-density closure

The strict reference density already contains the equilibrium 410, 520, and
660 km density transitions. The absolute dynamic phase fraction remains
`phase_B = X(E->T)` and is still used for phase-boundary output and latent
heating.

The momentum and matching mechanical-power/profile diagnostic contribution is:

```text
-Ra_phase * (X(E->T) - Xref(Tref))
```

Therefore `Ttotal=Tref` gives zero dynamic phase buoyancy, and the equilibrium
phase density is not counted a second time.

## E. Nondimensional consistency

The production strict cfg preserves `rho0`, `g0`, `alpha0`, `Cp0`, `Ttop`,
`Tbottom`, radius, diffusivity, viscosity, Rayleigh number, and dissipation
number from `cmbhf_ALA.cfg`.

The generator-reader contract is:

| Field | Serialized form |
|---|---|
| rho | `rho/rho0` |
| g | `g/g0` |
| Tref | `(Tref_K-Ttop)/(Tbottom-Ttop)` |
| alpha | `alpha/alpha0` |
| Cp | `Cp/Cp0` |
| beta | `R0*beta_SI` |
| Gamma_eff | `Di*alpha*g/(Cp*beta)` using serialized fields |
| Ks | physical GPa, diagnostic storage only |

The strict generator and validator report exact serialized-profile closure for
rho, beta, Gamma_eff, and Ks.

## F. Stage 6c shallow Schur experiment

Stage 6b.3 used `G diag(K)^-1 G^T` inside every cached shallow pressure
patch.  Stage 6c retains the same `6x6x2` patch geometry, overlap-three MPI
halo, global partition of unity, pressure shift, and PCG outer iteration, but
replaces the scalar velocity diagonal by assembled symmetric `3x3` same-node
velocity blocks.

Each node block is Cholesky factored as `K_n=L_n L_n^T`.  Pressure-gradient
features are formed as `phi=L_n^-1 g`; patch entries are dot products
`phi_i^T phi_j`.  Local and exchanged halo patches are therefore Gram
matrices and remain positive semidefinite before the existing pressure shift.
Nodes whose coupled factor cannot be formed fall back to the exact Stage 6b.3
diagonal metric and are counted in `velocity_block_fallbacks`.  Setting
`ala_shallow_patch_velocity_solver=diagonal` is the complete rollback.

Every PCG iteration now reports the absolute pressure-mass continuity norm
`mass_norm=N_M`, `Q`, and their relative scales.  The final original,
unaugmented momentum residual is audited before a failed continuity exit as
well as after a successful solve.

The 400-rank Stage 6c A/B run was numerically healthy but negative:
`node_block` changed the best cancellation only from `0.01339024` to
`0.01338443` at iteration 42 and ended at `0.01558697`.  It had zero node
fallbacks, positive curvature, unit residual drift, and a `3.0e-14`
preconditioner symmetry defect.  The strict production experiment therefore
returns to the exact `diagonal` rollback.

## G. Stage 6d cross-rank radial aggregate

Stage 6d retains the Stage 6b.3 `6x6x2`, overlap-three shallow Schwarz map and
adds a deterministic coarse space on each `2x2` horizontal MPI-rank group.
With four bins per member rank and two radial bins, each processor aggregate
has `8x8x2` selection bins and crosses four ranks.

Unlike the earlier GenEO experiment, the two coarse shapes are not selected
by an eigenvalue threshold.  They are the disjoint indicators of the two
radial portions of the 0--410 km shallow layer.  Their span contains the
complete aggregate constant while retaining an independent shallow/deep
radial amplitude.  The supports form an exact disjoint partition across
processor aggregates.  Normalization is deliberately deferred until after
transfer to the Galerkin pressure level.  Construction is purely geometric
on all `8x8x2` bins: neither the legacy GenEO `active_map` nor a selection
Rayleigh quotient may suppress a radial indicator.  The authoritative
validity checks are complete cross-rank support after transfer and the final
Galerkin Cholesky factorization.

Processor groups are separated by radial rank coordinate.  Only groups whose
local radial interval intersects the configured 0--410 km shallow layer own
these two modes.  Deeper groups have zero shallow coarse modes; requiring a
nonzero radial indicator on them would manufacture an identically zero basis
vector.

The coarse matrix remains

```text
S_c = P^T (D+C) Kfixed^-1 (D+C)^T P
```

and is assembled by the existing fixed strict-ALA Galerkin action, globally
symmetrized, shifted, and Cholesky factored.  Consequently the additive
coarse correction is fixed and SPD.  `ala_geneo_basis_type=spectral` is the
legacy rollback; `radial_partition` requires a cross-rank group and exactly
one configured mode per radial bin.

## Validation record

- Strict generator: PASS.
- Strict validator: PASS.
- Schema: 65 rows by 8 columns.
- rho analytical-profile relative RMS: 0.
- beta analytical-profile relative RMS: 0.
- interval beta-integral identity maximum: 0.
- Gamma_eff closure relative RMS and maximum: 0.
- Ks fitted-profile relative RMS: 0.
- strict static regressions: 20/20 PASS.
- changed C translation units: `mpicc -fsyntax-only` PASS.
- isolated full source compilation: all C objects and `libCitcomS.a` built;
  `CitcomSFull` linked manually because the imported Automake 1.9 templates do
  not expand `AR`/`LIBTOOL` with the locally available toolchain.

The production Pyre/LSF launcher and its Python 2.6/Linux runtime are not
executable on the local macOS host. The raw legacy C parser is not a substitute
for that launcher and is incompatible with the sectioned Pyre cfg. Consequently
the final 384-rank scheduler startup remains an HPC deployment check; no
reference-column mismatch was observed by generation, schema validation, or C
compilation.
