# ALA strict tests

Run the phase-reference buoyancy regression with:

```text
python3 tests/ala_strict/test_phase_buoyancy.py
```

The phase test verifies:

- `T=Tref` gives `X-Xref=0`;
- the dynamic absolute phase fraction remains the legacy `X`;
- strict buoyancy differs only by the radial reference subtraction;
- phase-boundary output still uses absolute `X`;
- latent heating still receives absolute `phase_B`;
- mechanical-power and buoyancy-profile diagnostics use `X-Xref`.

Future test groups remain reserved for operator adjoint, energy-budget, and
multigrid/Galerkin validation.
