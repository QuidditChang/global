# ALA strict legacy baseline

## Current baseline

- Repository component: `src/global`
- Branch: `cmbhf_ALA`
- Commit: `deeccebf22e55a6d608163aa5bfc8be8129e4d0b`
- Frozen on: 2026-07-29

The `cmbhf_ALA_strict` branch starts from this exact commit. Commit 0 adds
documentation and placeholder directories only. It does not change the
formulation, runtime configuration, reference data, or numerical results.

## Current formulation

The continuity constraint is

\[
G u = 0,
\]

where

\[
G = D + C,
\qquad
D u = \nabla\cdot u,
\qquad
C u = -\beta u_r,
\]

and the intended beta definition is

\[
\beta = -\frac{d\ln(\rho_{ref})}{dr}.
\]

The pressure operator is

\[
G^T p.
\]

The current implementation constructs the continuity and pressure operators
from the same discrete beta field.

## Current reference state

- Input source: Katsura 2022 Table S5
- Runtime file: `runs/refstate_ALA.txt`
- Runtime schema: seven numeric columns
- Fields:
  1. `rho`
  2. `g`
  3. `Tref`
  4. `alpha`
  5. `Cp`
  6. `beta`
  7. `Gamma_eff`

At the frozen baseline, `runs/refstate_ALA.txt` has 65 numeric rows and seven
fields per row.

## Known limitations

1. `rho_ref` was phase-smoothed during preprocessing.
2. `Tref` and `rho_ref` do not represent one thermodynamically closed state.
3. The reference state has no `Xref`.
4. Phase buoyancy uses the absolute phase function rather than `X-Xref`.
5. Latent heat is not represented in a conservative thermodynamic form.
6. Dynamic pressure work is not implemented in the energy equation.
7. The active beta authority is node-based and is mapped to elements by the
   current material-property path.

These limitations are baseline facts. Commit 0 does not correct or otherwise
change them.

## Purpose of the strict branch

The `cmbhf_ALA_strict` branch will introduce, in later reviewed commits:

- ReferenceState v2;
- reference phase fractions;
- explicit density decomposition;
- a classic strict energy formulation;
- an independent validation framework.

The legacy `cmbhf_ALA` branch remains the numerical and scientific comparison
baseline throughout the migration.
