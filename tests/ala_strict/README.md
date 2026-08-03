# ALA strict tests

Run the production strict-ALA static regressions with:

```text
python3 -m unittest discover -s tests/ala_strict -p 'test_*.py' -v
```

The phase test verifies:

- `T=Tref` gives `X-Xref=0`;
- the dynamic absolute phase fraction remains the legacy `X`;
- strict buoyancy differs only by the radial reference subtraction;
- phase-boundary output still uses absolute `X`;
- latent heating still receives absolute `phase_B`;
- mechanical-power and buoyancy-profile diagnostics use `X-Xref`.

The production-architecture test also verifies:

- the strict cfg changes only the reference-state inputs and the explicit
  Stage 6d cross-rank aggregate experiment, Stage 8 recursive pressure
  multigrid transfer/operator wiring;
- the eight-column strict reference-state order and pure fitted Tref endpoints;
- `E->T` remains the sole temperature state and assimilation target;
- no `DataT` cache or runtime strict-audit hook remains;
- `Xref` uses `E->refstate.Tref` while dynamic `X` uses `E->T`;
- the Stage 6c node-block shallow Schur map is constructed as an SPD Gram
  operator and retains the exact Stage 6b.3 diagonal rollback;
- the Stage 6d basis consists of deterministic radial partitions on a
  cross-rank aggregate and is solved with the fixed strict Galerkin operator;
- PCG failure output contains absolute `N_M`, `Q`, and the original
  unaugmented momentum audit.

The coupled-operator architecture test verifies the first monolithic-solver
foundation without selecting it in production:

- one allocation-free action returns `K_gamma*u+(D+C)^T*p` and `(D+C)*u`;
- the block action never assembles or consumes the pressure-dependent legacy
  force;
- a finest-level helper reconstructs the pressure-independent momentum RHS as
  `E->F + (D+C)^T*p - D^T*p`;
- destructive field aliases are rejected; and
- Stage 7d/Stage 8 remain the only runtime paths until a coupled Krylov method
  is added and validated.

Future test groups remain reserved for operator adjoint, energy-budget, and
multigrid/Galerkin validation.
