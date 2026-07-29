# Strict temperature reference semantic audit

## Scope

This audit adds diagnostics only. It does not replace the serialized
reference temperature, change temperature initialization, alter phase
fractions, or modify the energy equation.

The runtime diagnostic is executed immediately after the initial temperature
constructor returns. It reads the completed `E->T` field and the fixed
reference-state arrays without writing to either one.

## Reference-temperature source

Reference-state column 3 originates from the `T_K` values in Katsura 2022
Table S5. The strict generator:

1. separates the source samples into four phase branches;
2. applies the existing endpoint-constrained piecewise cubic fit independently
   to each branch;
3. evaluates that fit directly at CitcomS radial nodes;
4. serializes the surface and CMB rows using the physical total-temperature
   boundary values.

The reader retains the serialized nodal values in
`E->refstate.temperature` and exposes the same allocation through the semantic
alias `E->refstate.Tref`. It separately reconstructs smooth surface and CMB
background endpoints from neighboring interior samples for use by the
initial-temperature builder. Consequently the serialized array has mixed
endpoint semantics even though the initializer already treats boundary-layer
temperature as a superposed anomaly.

## Runtime temperature quantities

The audit names the existing quantities without introducing a new state:

\[
T_{\mathrm{ref}}^* = E\mathbin{\to}\mathrm{refstate.Tref},
\]

\[
E\mathbin{\to}DataT(t=0)
=E\mathbin{\to}T(t=0)-T_{\mathrm{ref}}^*,
\]

and `E->T` remains the dynamic total-temperature field evolved by the existing
energy equation.

The horizontal means of `Ttotal`, `Tref`, and `DataT` are reported at the model nodes nearest
the surface, 410 km, 660 km, and CMB. Both nondimensional values and kelvin
values are printed. Restart input is explicitly skipped because a restart
field is not the model's constructed initial temperature.

## Adiabatic representation diagnostic

For every radial interval the audit reconstructs dimensional temperature,

\[
T_K=T_{\mathrm{top}}+T_{\mathrm{ref}}^*\Delta T,
\]

and compares

\[
\frac{\Delta\ln T_K}{\Delta r^*}
\quad\hbox{with}\quad
-Di\,\overline{\frac{\alpha^*g^*}{C_p^*}}.
\]

It reports the relative RMS, maximum relative difference, and maximum absolute
residual. Serialized endpoints participate intentionally: the result therefore
diagnoses the actual runtime array and exposes the distinction between its
interior Katsura fit and its total-temperature endpoint rows. The diagnostic
does not modify or reject the profile.

For the current 65-node `runs/refstate_ALA_strict.txt` artifact, an offline
reproduction of the runtime calculation gives:

| Metric | Value |
|---|---:|
| relative RMS, all intervals | `5.0206176469e+01` |
| maximum relative difference | `3.5092304583e+02` |
| maximum absolute residual | `7.3733610580e+02` |
| depth of maximum relative difference | `7.282053 km` |
| relative RMS excluding the two endpoint intervals | `5.9044027661e-01` |

These are representation diagnostics for the serialized profile, not failure
criteria. In particular, the all-interval result is dominated by the imposed
total-temperature endpoint treatment. Interior phase-branch joins and nodal
sampling also remain represented in the second value; this audit does not
reinterpret them as a new temperature profile.

## Phase semantics

The phase paths remain unchanged:

- reference phase fraction `Xref` uses `E->refstate.Tref`;
- dynamic absolute phase fraction `X` uses the current total temperature
  `E->T`;
- strict phase buoyancy continues to use `X-Xref`.

## Assessment

The interior Katsura temperature fit is appropriate in purpose as a candidate
thermodynamic reference adiabat. Suitability for strict discrete closure is
determined by the new adiabatic diagnostic, not by assuming that the
serialized surface and CMB total-temperature values are adiabat endpoints.

The current serialized array as a whole is **not semantically suitable as a
pure thermodynamic reference temperature**, because its two endpoint rows are
total-temperature boundary values and its measured representation difference
is large. This does not establish that the interior Katsura adiabat is
physically wrong. It establishes that the interior profile and boundary-layer
closure must not be judged as one undifferentiated `Tref` quantity.

The architecture now names the fixed reference explicitly through `Tref` and
stores an initialization anomaly in `DataT`. `E->T` remains the authoritative
total-temperature field constructed from a Katsura background plus thermal
anomalies. `DataT` is not yet synchronized after subsequent accepted thermal
or assimilation updates, so it is not an independently evolved state.
