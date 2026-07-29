# Strict ALA thermal, density, and assimilation architecture audit

## Executive judgment

**Result: B — the current operational configuration has an actual reference
state inconsistency and cannot yet be called fully closed.**

The implementation has several correct separations:

- reference arrays are fixed after loading;
- the energy equation solves total temperature `E->T`;
- reference phase uses `Tref`, while dynamic phase uses total temperature;
- phase buoyancy uses `X-Xref`, avoiding double counting when the loaded
  density is phase-inclusive;
- assimilation writes total temperature and never writes the reference state.

The blocking operational inconsistency is configuration selection:

```text
runs/cmbhf_ALA.cfg
    refstate_file = refstate_ALA.txt
```

The selected seven-column file is the phase-smoothed product, while the current
phase code subtracts `Xref`. The phase-inclusive eight-column
`refstate_ALA_strict.txt` exists but is not selected. Consequently the
configured reference density does not contain the equilibrium phase density
that the updated phase force assumes is already in `rho_ref`.
Because the selected file has no eighth-column `Ks`, the current strict
temperature startup audit also classifies it as legacy and skips its
four-depth runtime report.

Two additional limitations are semantic or diagnostic rather than active
solver errors:

- serialized Tref has total-temperature endpoint rows;
- `DataT` is initialized once and becomes stale after thermal evolution or
  assimilation, but no physical operator consumes it.

## Reference and dynamic state ownership

| State | Runtime storage | Mutability | Consumer |
|---|---|---|---|
| reference density | `E->refstate.rho` | fixed | continuity coefficients, thermal inertia, buoyancy scaling, phase pressure approximation |
| reference temperature | `E->refstate.Tref` | fixed alias | reference phase and diagnostics |
| total temperature | `E->T` | evolving | energy equation, dynamic phase, thermal buoyancy, output |
| stored temperature anomaly | `E->DataT` | initialized once | diagnostics only |
| assimilation increment | `E->assim_delta_T` | reset and accumulated per attempted step | assimilation heat accounting |

Reference and dynamic state are separated in the physical paths. The existing
`DataT` cache is not an evolution variable and must not be treated as current
after initialization.

## Initial temperature construction

For the active `lith_age` initialization path:

```text
fitted Katsura background
    with reconstructed background endpoints
        |
        +--> add CMB erfc thermal-boundary-layer anomaly
        |
        +--> subtract shallow half-space-cooling anomaly
        |
        +--> overwrite selected shallow slab nodes with prescribed Ttotal
        |
        +--> update and enforce total-temperature boundary conditions
        |
        v
E->T(x,0) = T_initial
        |
        v
E->DataT(x,0) = E->T(x,0) - E->refstate.Tref(r)
```

The surface cooling and CMB TBL operations are explicitly additive. The slab
path is an absolute assignment, so its provenance is not retained as a
separate named anomaly. Algebraically the final field still satisfies

\[
T_{\mathrm{initial}}=T_{\mathrm{ref}}+\delta T_{\mathrm{initial}}
\]

when the anomaly is defined after construction.

### Current four-depth snapshot

The table uses the nodes nearest the requested depths in the current 65-node
mesh and the available strict eight-column product. The mesh does not contain
nodes exactly at 410 or 660 km. Surface/CMB values are enforced boundary
values. At the selected 398.697 and 694.535 km nodes, the shallow
lithosphere/slab operations do not apply and the configured bottom-TBL
contribution underflows to zero at double precision.

| Requested depth | Selected node depth | Tref* | Tref K | Ttotal* | Ttotal K | deltaT* |
|---:|---:|---:|---:|---:|---:|---:|
| 0 km | 0.000000 km | 0.000000000 | 300.000 | 0.000000000 | 300.000 | 0 |
| 410 km | 398.697180 km | 0.427145927 | 1795.011 | 0.427145927 | 1795.011 | 0 |
| 660 km | 694.534565 km | 0.478247081 | 1973.865 | 0.478247081 | 1973.865 | 0 |
| CMB | 2891.000525 km | 1.000000000 | 3800.000 | 1.000000000 | 3800.000 | 0 |

For these analytically determined samples,

\[
\max|T_{\mathrm{total}}-(T_{\mathrm{ref}}+\delta T)|=0
\]

in the reproduced double-precision calculation. A full horizontal runtime
table was not fabricated: the production initialization requires external
plate-age/reconstruction data and the current cfg does not select the strict
file. The existing startup audit prints the corresponding MPI-reduced values
when that configuration is run.

## Density and buoyancy decomposition

### What the code actually stores

CitcomS does not reconstruct or evolve a nodal `rho_total` field in this path.
`E->rho` is allocated but has no active assignments. The implementation
constructs buoyancy-equivalent contributions:

\[
B_T
=Ra\,\rho_{\mathrm{ref}}\alpha T_{\mathrm{total}},
\]

\[
B_{\mathrm{phase},j}
=-Ra_{\mathrm{phase},j}(X_j-X_{\mathrm{ref},j}),
\]

plus chemical buoyancy when enabled. Each contribution is horizontally
mean-removed, multiplied by reference gravity, and assembled into momentum.

Because `rho_ref`, `alpha`, and `Tref` are radial,

\[
\mathcal{R}_h[\rho_{\mathrm{ref}}\alpha T_{\mathrm{total}}]
=
\mathcal{R}_h[
\rho_{\mathrm{ref}}\alpha(T_{\mathrm{total}}-T_{\mathrm{ref}})
],
\]

where \(\mathcal{R}_h\) denotes horizontal-mean removal. Thus the thermal
forcing is equivalent to using the lateral part of the temperature anomaly
even though it is implemented with total temperature.

A physical-density interpretation may be written diagnostically as

\[
\delta\rho_T
=-\rho_{\mathrm{ref}}\alpha
(T_{\mathrm{total}}-T_{\mathrm{ref}}),
\]

\[
\delta\rho_{\mathrm{phase}}
=\sum_j\Delta\rho_j(X_j-X_{\mathrm{ref},j}),
\]

but those arrays are not explicitly constructed. With the active cfg,
chemical buoyancy is also enabled, so the complete diagnostic decomposition
would require a composition term:

\[
\rho_{\mathrm{diagnostic}}
=\rho_{\mathrm{ref}}
+\delta\rho_T
+\delta\rho_{\mathrm{phase}}
+\delta\rho_{\mathrm{chemical}}.
\]

The requested three-term expression omits an active model contribution.

### Four-depth diagnostic decomposition

At the four initial samples above, `Ttotal=Tref`, hence
\(\delta\rho_T=0\). Dynamic and reference phase fractions are also identical,
hence \(\delta\rho_{\mathrm{phase}}=0\). Ignoring the separately configured
chemical field, the diagnostic total equals strict `rho_ref`:

| Requested depth | Selected node depth | rho_ref* | rho_ref kg m^-3 | delta rho thermal | delta rho phase | diagnostic rho total kg m^-3 |
|---:|---:|---:|---:|---:|---:|---:|
| 0 km | 0.000000 km | 0.692366203 | 3076.875 | 0 | 0 | 3076.875 |
| 410 km | 398.697180 km | 0.768353114 | 3414.561 | 0 | 0 | 3414.561 |
| 660 km | 694.534565 km | 0.967915728 | 4301.417 | 0 | 0 | 4301.417 |
| CMB | 2891.000525 km | 1.213894013 | 5394.545 | 0 | 0 | 5394.545 |

These are diagnostic identities, not an output from a stored total-density
field.

### Phase double-counting result

When `refstate_ALA_strict.txt` is loaded:

```text
raw Katsura phase-inclusive rho
        + Delta rho (X-Xref)
```

is free of static phase double counting. At `Ttotal=Tref`, the dynamic phase
contribution is pointwise zero.

When the currently configured `refstate_ALA.txt` is loaded:

```text
phase-smoothed rho
        + Delta rho (X-Xref)
```

does not double count phase density, but it omits the equilibrium reference
phase density. For scale, the current files differ by about 697 kg m^-3 at the
surface and at the node nearest 410 km; they converge below the 660-km
transition. This is the principal physical closure failure found by the audit.

## Temperature assimilation

The accepted-step assimilation flow is:

```text
energy solver advances E->Ttotal
        |
        +--> rejected attempt:
        |       restore saved E->T and E->Tdot
        |       clear assim_delta_T
        |
        +--> accepted attempt:
                update plate-age thermal boundary data
                        |
                        v
                assimilate_lith_conform_bcs()
                        |
                        v
                Tnew = daf*Told + (1-daf)*T_assimilated
                        |
                        +--> write E->T = Tnew
                        |
                        +--> assim_delta_T += Tnew-Told
                        |
                        v
                measure rho_ref*Cp*assim_delta_T/dt heating
```

Assimilation modifies total temperature, not the reference state and not an
independent anomaly equation. Direct assignment to `E->T` does not invalidate
`Tref`; after the update one may always define

\[
T_{\mathrm{anomaly,new}}=E\mathbin{\to}T_{\mathrm{new}}-T_{\mathrm{ref}}.
\]

The actual issue is narrower: `E->DataT` is not refreshed, so the stored cache
no longer satisfies `E->T=Tref+DataT` after the first thermal or assimilation
update. Since no physics reads `DataT`, this causes no current solver error,
but it makes `DataT` unsuitable as a dynamic diagnostic after initialization.

## Is DataT required?

No independent `DataT` field is required by the present formulation:

- the energy equation solves total temperature;
- restart files store and restore total temperature;
- assimilation supplies and modifies total temperature;
- phase needs total temperature for `X` and fixed `Tref` for `Xref`;
- thermal anomaly can be evaluated on demand as `E->T-Tref`;
- horizontal-mean removal already eliminates any purely radial Tref
  contribution from thermal buoyancy.

The existing `DataT` allocation is therefore redundant for current physics.
It is only an initialization diagnostic cache and should not be promoted to an
authority. A future implementation should either compute the anomaly on
demand or explicitly synchronize a diagnostic cache after every accepted
temperature modification. An independently evolved anomaly is justified only
if a future energy formulation, restart schema, or assimilation algorithm
requires it.

## Closure matrix

| Requirement | Code path with strict file | Current cfg |
|---|---|---|
| reference/dynamic storage separated | pass | pass |
| fixed thermodynamic Tref | partial: internal fit yes, serialized endpoints mixed | partial |
| total temperature solved by energy equation | pass | pass |
| initial temperature = reference + anomaly | pass algebraically | pass algebraically |
| assimilation leaves reference unchanged | pass | pass |
| `Xref=f(Tref)`, `X=f(Ttotal)` | pass | pass |
| phase-inclusive rho plus `Delta rho (X-Xref)` | pass | **fail: cfg selects phase-smoothed rho** |
| stored DataT remains current dynamically | fail, diagnostic cache only | fail, diagnostic cache only |
| explicit total-density reconstruction | not implemented or required by current ALA force path | same |

## Final decision

The temperature evolution and assimilation architecture do not corrupt the
reference state. Solving `E->T=Ttotal` and calculating an anomaly on demand is
sufficient; assimilation does not require an independent anomaly variable.

The current branch is nevertheless category **B**, because the checked runtime
configuration selects a phase-smoothed reference density while phase buoyancy
assumes the equilibrium phase density is already in the reference state.
Selecting and validating the phase-inclusive strict product is required before
claiming operational strict-ALA density closure. Purifying the serialized Tref
endpoints is a separate reference-data issue.
