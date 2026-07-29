# Strict ALA temperature architecture

## Runtime identities

This migration establishes three explicit temperature meanings without
changing the temperature equation:

| Quantity | Storage | Meaning |
|---|---|---|
| `Tref` | `E->refstate.Tref` | fixed thermodynamic reference temperature |
| `Ttotal` | `E->T` | authoritative total temperature, evolved by the existing energy equation |
| `DataT` | `E->DataT` | derived anomaly cache, `Ttotal-Tref` |

The legacy `E->refstate.temperature` pointer remains the storage owner for
reference-state column 3. `E->refstate.Tref` is a non-owning semantic alias of
the same allocation:

```text
E->refstate.Tref == E->refstate.temperature
```

There is no new reference-state column, interpolation, normalization, or
temperature value.

## Initialization

The existing initializer first constructs `E->T` exactly as before, including
the Katsura background, thermal boundary layers, lithosphere/slab structure,
restart input, and boundary-condition conformance. Only after that process
returns does the new synchronization pass compute

\[
E\mathbin{\to}DataT(\mathbf{x},0)
=E\mathbin{\to}T(\mathbf{x},0)-T_{\mathrm{ref}}(r).
\]

The audit reports horizontal means of `Ttotal`, `Tref`, and `DataT` at the
nodes nearest the surface, 410 km, 660 km, and CMB. It also checks

\[
\max_{\mathbf{x}}
\left|T_{\mathrm{total}}-
\left(T_{\mathrm{ref}}+DataT\right)\right|.
\]

Because `DataT` is constructed by the same floating-point subtraction, the
residual is expected at machine precision. It is diagnostic only and never
changes `E->T`.

## Phase temperature paths

The phase physics and formulas are unchanged. Their temperature inputs are now
explicit:

```text
reference phase: Xref = f(E->refstate.Tref)
dynamic phase:   X    = f(E->T)
phase anomaly:          X - Xref
```

Changing `temperature` to its pointer alias `Tref` in the reference path has
no numerical effect.

## Assimilation interface

The current assimilation path remains unchanged:

1. it computes an assimilated total temperature;
2. it updates `E->T`;
3. it accumulates that change in `E->assim_delta_T` for the existing heating
   accounting.

`assim_delta_T` is an incremental assimilation update, whereas `DataT` is the
full anomaly relative to the fixed thermodynamic reference. They must not be
identified:

\[
\delta T_{\mathrm{assim}}
\ne
T_{\mathrm{total}}-T_{\mathrm{ref}}.
\]

A future migration may synchronize the anomaly cache after an accepted
assimilation/thermal update:

```text
assimilation updates E->T
        |
        +--> existing assim_delta_T accounting
        |
        +--> DataT = E->T - E->refstate.Tref
```

That future synchronization must occur only after an accepted timestep so a
rejected thermal attempt cannot leak state. This commit does not modify
assimilation, advection, diffusion, heating, or boundary-condition physics.

## Current lifecycle limitation

`E->T` remains the sole evolving temperature authority. `DataT` is synchronized
after initial/restart temperature construction and is not yet consumed by any
physical operator. Until a later non-physics synchronization commit is
designed for accepted timesteps, it is an initialization snapshot rather than
an independently evolved field.
