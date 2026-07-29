# ALA strict target design

## Status

This document describes a future architecture. Commit 0 implements none of the
features below and changes no physical or numerical formulation.

## 1. Strict ALA definition

The strict branch retains

\[
G u = 0,
\]

\[
G = D + C,
\]

\[
C u = -\beta u_r,
\]

with

\[
\beta = -\frac{d\ln(\rho_{ref})}{dr}.
\]

The existing Stokes block structure is retained. Future work must preserve the
discrete adjoint relationship between `G` and `G^T`.

## 2. Reference-state goal

ReferenceState v2 will provide a single, versioned source for:

- `rho_ref`;
- `T_ref`;
- `P_ref`;
- `g`;
- `alpha`;
- `Cp`;
- `Ks`;
- `gamma`;
- `beta`;
- reference phase fractions `Xref`.

The state will carry units, provenance, source hashes, generation metadata, and
independent closure diagnostics.

## 3. Architecture principles

### Invariant 1: one reference state

All reference thermodynamic variables must originate from one ReferenceState.
Runtime modules must not independently reconstruct competing reference
profiles.

### Invariant 2: reference and dynamic phase separation

Reference phase enters `rho_ref` and `Xref`. Dynamic phase density uses only

\[
\delta X = X-Xref.
\]

An absolute phase function must not be added to a reference density that
already includes the same transition.

### Invariant 3: one discrete G

Continuity, pressure, augmented Lagrangian, multigrid, and diagnostics must use
an identical discrete `G`. A solver-specific or diagnostic-only beta
reconstruction is not permitted.

## 4. Migration philosophy

The current `cmbhf_ALA` branch remains the legacy numerical baseline.
`cmbhf_ALA_strict` evolves independently through feature-gated, reversible
milestones.

New infrastructure is introduced before it is allowed to affect production
physics. The complete strict path is enabled only after reference-state, phase,
density, energy, operator-adjoint, and multigrid validation gates pass.

Commit 0 establishes documentation and ownership boundaries only.
