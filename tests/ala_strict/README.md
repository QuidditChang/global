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
  Stage 6c shallow velocity metric;
- the eight-column strict reference-state order and pure fitted Tref endpoints;
- `E->T` remains the sole temperature state and assimilation target;
- no `DataT` cache or runtime strict-audit hook remains;
- `Xref` uses `E->refstate.Tref` while dynamic `X` uses `E->T`;
- the Stage 6c node-block shallow Schur map is constructed as an SPD Gram
  operator and retains the exact Stage 6b.3 diagonal rollback;
- PCG failure output contains absolute `N_M`, `Q`, and the original
  unaugmented momentum audit.

Future test groups remain reserved for operator adjoint, energy-budget, and
multigrid/Galerkin validation.
