# Strict ALA Tref usage report

## Audit snapshot

This report audits the working `cmbhf_ALA_strict` implementation without
changing reference data or runtime physics.

| Input | SHA-256 |
|---|---|
| `runs/refstate_ALA.txt` | `af8967ab45dd84ccc3deb4042c8b2926984b46226c4e21c8fab619fe21945874` |
| `runs/refstate_ALA_strict.txt` | `5a118047ae1e03c1449522702c44c030859a95807f55ad29f8df6c7b828cdfe6` |
| `runs/cmbhf_ALA.cfg` | `b4e8f939e107eb9bb2f7f924ddb080cdf5d0365157d92621c36c60654bd07714` |

The active configuration names `refstate_ALA.txt`, not
`refstate_ALA_strict.txt`. Conclusions about the strict eight-column product
therefore describe the available strict path; they are not evidence that the
current cfg selects it.

## Column 3 provenance

The strict generator obtains column 3 through this chain:

```text
Katsura2022_TableS5_thermoelastic.csv / T_K
        |
        v
canonical field Tref_K
        |
        v
four phase branches:
0-410, 410-520, 520-660, 660-CMB
        |
        v
endpoint-constrained cubic fit in each branch
        |
        v
direct evaluation at GLB.coor.global.dat radial nodes
        |
        +--> interior values: fitted Katsura temperature
        |
        +--> CMB serialized row: Tbottom = 3800 K
        |
        +--> surface serialized row: Ttop = 300 K
        |
        v
(T_K-Ttop)/(Tbottom-Ttop)
        |
        v
refstate_ALA_strict.txt column 3
```

For normalized branch coordinate

\[
x=\frac{z-z_0}{z_1-z_0},
\]

the fitted field has the existing form

\[
T_{\mathrm{fit}}(x)
=T_0+(T_1-T_0)x
+a\,x(x-1)
+c\,x(x-1)(2x-1).
\]

The fit coefficients come from least squares while both source endpoints are
retained. `alpha`, `g`, and `Cp` are fitted independently; they are not used to
integrate column 3.

## Serialization and reader semantics

The serialized array is not a pure thermodynamic adiabat:

- internal rows are the phase-branch Katsura fit;
- the surface row is overwritten with 300 K;
- the CMB row is overwritten with 3800 K.

The C reader stores column 3 in `E->refstate.temperature`.
`E->refstate.Tref` is a non-owning pointer alias of the same allocation, so the
alias changes semantics but no value.

For initialization only, the reader separately reconstructs
`temperature_surface` and `temperature_cmb` by quadratic continuation of the
three neighboring internal rows. `Lith_age.c` uses those reconstructed values
as the background endpoints before adding thermal boundary layers. The
serialized 300/3800 K rows remain present in `Tref` itself.

The precise distinction is therefore:

| Quantity | Meaning |
|---|---|
| `Tref_internal` | fitted Katsura branch temperature at internal nodes |
| serialized surface/CMB `Tref` rows | total-temperature Dirichlet values |
| `temperature_surface`, `temperature_cmb` | reconstructed background endpoints used only by initialization |

The whole serialized array should not be described as a semantically pure
thermodynamic reference temperature.

## Complete call classification

No `temperature_ref` or `reference_temperature` runtime variable exists.
Relevant uses are:

| File | Function / location | Class | Purpose | Should use Tref? |
|---|---|---|---|---|
| `Material_properties.c` | `mat_prop_allocate` | storage | allocate legacy column-3 owner and assign `Tref` alias | yes |
| `Material_properties.c` | `read_refstate` | input | parse column 3 into `temperature`; reconstruct separate initialization endpoints | yes for parsed reference; endpoint mixing remains |
| `Material_properties.c` | `reference_state` | diagnostic | print radial reference/background profile | yes, diagnostic |
| `Material_properties.c` | `validate_strict_reference_state` | reference thermodynamics | verify `X(Tref)-Xref` | yes |
| `Material_properties.c` | `adams_williamson_eos` | legacy initialization | set unavailable reference temperature to zero | legacy only |
| `Phase_change.c` | `phase_change_reference_fraction` | reference thermodynamics | calculate `Xref=f(Tref)` | yes |
| `Phase_change.c` | `calc_phase_change` | dynamic state | calculate `X=f(E->T)` | no; total temperature is correct |
| `Lith_age.c` | `lith_age_construct_tic` | initial construction | seed background from column 3, using reconstructed endpoints | yes as initialization background |
| `Initial_temperature.c` | `initialize_temperature_anomaly` | diagnostic storage | initialize `DataT=E->T-Tref` | yes |
| `Initial_temperature.c` | semantic audit | diagnostic | report `Ttotal`, `Tref`, `DataT`, adiabatic representation and reconstruction | yes |
| `Advection_diffusion.c` | temperature solve | dynamic state | evolve `E->T` and use total absolute temperature in adiabatic/latent terms | no; total temperature is correct |
| `Pan_problem_misc_functions.c` | `get_buoyancy` | dynamic forcing | thermal buoyancy from `E->T`, followed by horizontal-mean removal | explicit Tref not required |
| `Output.c` / profile diagnostics | temperature output | output | output dimensional total `E->T` or derived buoyancy fields | no |

### Terms that do not use Tref

- `Gamma_eff` closure uses `alpha`, `g`, `Cp`, and `beta`; it does not use
  temperature.
- Adiabatic heating uses total `E->T`, as required by the present total-
  temperature energy formulation.
- Dynamic pressure buoyancy uses `beta` and dynamic pressure.
- Dynamic pressure work is not implemented.
- Boundary-condition routines modify or constrain `E->T`; they do not write
  the reference-state array.

## Phase confirmation

The phase temperature paths are correctly separated:

\[
X_{\mathrm{ref}}=f(E\mathbin{\to}\mathrm{refstate.Tref}),
\]

\[
X=f(E\mathbin{\to}T).
\]

Strict phase buoyancy uses `X-Xref`. Neither initialization nor assimilation
writes `E->refstate.Tref`.

## Tref judgment

The internal Katsura fit is suitable in role as a thermodynamic-reference
candidate. The serialized array as a whole is not yet a pure thermodynamic
reference because total-temperature boundary values occupy its endpoint rows
and its previously measured discrete adiabatic representation mismatch is
large.

This is partly a semantic/data-product problem rather than evidence that all
internal Katsura temperatures are physically wrong. Nevertheless, strict
reference-state closure cannot be claimed for the complete serialized
temperature array in its current form.
