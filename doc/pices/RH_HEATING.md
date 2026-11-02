# Fixed primordial heating contrast

`RH` defaults to 1. Q0 remains the reference uniform source coefficient.
For RH > 1 this release supports PICES with ratio composition, flavor 24,
and PICES checkpointing. Legacy tracer_enriched is incompatible. The common
source evaluator is also used by PG residual and profile paths; enabling an
RH PG experiment requires its own restart-state support first.

After fresh composition initialization, sum the same reference-density-weighted
Gauss volumes used by the EBA source. f0=sum(rho*V*C24)/sum(rho*V).
Set Qbg=Q0/(1+(RH-1)*f0), Qprim=RH*Qbg. These are immutable thereafter.
The elemental source is rho*((1-C24)*Qbg+C24*Qprim). No new markers, changing
particle weights, decay or crust heat budget is introduced.

Schema 6 appends four doubles to the existing authenticated PICES state:
f0, Qbg, Qprim, initial nondimensional integrated power. They are restored
before PICES reconstruction, validated, and included in the physics hash.
RH=1 retains the legacy checkpoint schema and source calculation.

RH_HEATING on rank zero at initialization and storage_spacing reports current
whole-shell component powers in W, initial baseline, drift, density-weighted
fraction, primordial volume (nondimensional), and primordial power assigned to
empty cells. This is instantaneous composition, not an additional heat source.
PICES global/profile heat-stage diagnostics retain their existing stage and
units and use the actual modified source. Empty cells retain composition by
the existing policy; normalization does not hide subsequent inventory drift.

Validation: tests/pices/check_rh_heating.py BUILD RUNS OUTPUT runs real 12-rank
solves using the existing P8a miniature forcing. It compares RH absent/1,
checks initial power and contrast, source-ledger agreement, fixed coefficients,
continuous versus split restart, mismatched RH and corrupted state rejection.
Unused array index zero is excluded from binary comparisons (legacy checkpoint
writer explicitly stores uninitialized index zero); physical entries and state
sidecars are compared exactly.
