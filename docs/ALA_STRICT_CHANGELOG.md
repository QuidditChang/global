# ALA strict changelog

## Commit 0

Date: 2026-07-29

Changes:

- created the `cmbhf_ALA_strict` branch from the frozen `cmbhf_ALA` baseline;
- recorded the legacy baseline identity and regression evidence;
- added strict-ALA architecture documentation;
- reserved future subsystem and test directories with README files.

Physics changes: **NONE**

Numerical changes: **NONE**

Runtime configuration changes: **NONE**

Reference-data changes: **NONE**

## Temperature reference semantic audit

Changes:

- added a read-only post-initialization audit of
  `Tref`, `deltaT_initial`, and the dynamic total-temperature field;
- added a runtime adiabatic representation diagnostic;
- documented the unchanged reference and dynamic phase-temperature paths.

Physics changes: **NONE**

Numerical changes: **NONE**

Runtime configuration changes: **NONE**

Reference-data changes: **NONE**

## Temperature architecture separation

Changes:

- added `E->refstate.Tref` as a semantic alias of the legacy reference
  temperature allocation;
- added derived `E->DataT` storage initialized as `E->T-Tref`;
- added a machine-precision `Ttotal=Tref+DataT` audit;
- documented the unchanged assimilation interface and future synchronization
  boundary.

Physics changes: **NONE**

Numerical changes: **NONE**

Runtime configuration changes: **NONE**

Reference-data changes: **NONE**

## Thermal, density, and assimilation architecture audit

Changes:

- documented complete Tref provenance and runtime use classification;
- audited initial-temperature, phase-density, thermal-buoyancy, restart, and
  assimilation data flows;
- recorded the active cfg/reference-file mismatch and the diagnostic-only
  lifecycle of `DataT`.

Physics changes: **NONE**

Numerical changes: **NONE**

Runtime configuration changes: **NONE**

Reference-data changes: **NONE**

## Production configuration freeze

Changes:

- added `runs/cmbhf_ALA_strict.cfg` with the strict reference-state input;
- made serialized Tref a pure Katsura fitted thermodynamic profile at all
  nodes, without 300 K/3800 K endpoint overwrites;
- made the reader use the serialized Tref endpoints directly;
- retained `E->T` as the total-temperature state and `assim_delta_T` as the
  existing assimilation increment;
- removed the runtime strict-reference and temperature-semantic audit hooks;
- removed the redundant one-time `DataT` allocation and cache;
- added production static regressions and the readiness report.

Physics changes: strict reference-temperature semantics finalized; no solver,
energy, continuity, or phase-equation change.

Numerical changes: initial total temperature now superposes the unchanged
`cmbhf_ALA` anomaly construction on pure fitted Tref endpoints.

Runtime configuration changes: strict production cfg selects
`refstate_ALA_strict.txt`.

Reference-data changes: column 3 endpoint overwrites removed; schema and all
other columns unchanged.
