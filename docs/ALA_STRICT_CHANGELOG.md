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
