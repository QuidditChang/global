# Strict ALA Production-Readiness Report

## Status

`runs/cmbhf_ALA_strict.cfg` is the production strict-ALA configuration. It is
an exact physical-configuration copy of `runs/cmbhf_ALA.cfg` except for:

```text
refstate_file = refstate_ALA_strict.txt
```

No Stokes, continuity, multigrid, energy-equation, phase-equation, boundary-
condition, assimilation, timestep, or viscosity parameter is changed.

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
| Gamma_eff | `alpha*g/(Cp*beta)` using serialized fields |
| Ks | physical GPa, diagnostic storage only |

The strict generator and validator report exact serialized-profile closure for
rho, beta, Gamma_eff, and Ks.

## Validation record

- Strict generator: PASS.
- Strict validator: PASS.
- Schema: 65 rows by 8 columns.
- rho analytical-profile relative RMS: 0.
- beta analytical-profile relative RMS: 0.
- interval beta-integral identity maximum: 0.
- Gamma_eff closure relative RMS and maximum: 0.
- Ks fitted-profile relative RMS: 0.
- strict static regressions: 12/12 PASS.
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
