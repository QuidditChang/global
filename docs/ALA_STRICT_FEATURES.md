# ALA strict feature-gate plan

## Status

Commit 0 documents future configuration contracts only. No feature gate is
implemented, parsed, or enabled.

## Reference-state schema

Planned values:

```text
reference_state_schema = legacy_v1 | strict_v2
```

- `legacy_v1`: current seven-column reference-state behavior.
- `strict_v2`: future versioned ReferenceState.

Commit 0 effective state: `legacy_v1` by unchanged legacy behavior.

## Phase model

Planned values:

```text
phase_model = legacy | strict
```

- `legacy`: current absolute phase-function behavior.
- `strict`: future `Xref`, `X`, and `X-Xref` behavior.

Commit 0 effective state: `legacy`.

## Density model

Planned values:

```text
density_model = legacy | decomposed
```

- `legacy`: current buoyancy assembly.
- `decomposed`: future thermal, composition, phase, and pressure-diagnostic
  density components.

Commit 0 effective state: `legacy`.

## Energy formulation

Planned values:

```text
energy_formulation = classic_strict | extended_anelastic
```

- `classic_strict`: future conservative classic strict-ALA energy path without
  dynamic pressure work.
- `extended_anelastic`: future extended energy path; it must remain unavailable
  until separately implemented and validated.

Commit 0 does not add either key to a runtime configuration. Existing energy
behavior remains unchanged.

## Gate policy

Future gates must:

- fail fast on unsupported values;
- preserve a complete legacy rollback path;
- never silently select an incomplete strict feature;
- be enabled only with corresponding validation coverage.
