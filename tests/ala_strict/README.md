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

- the strict cfg changes only the reference-state inputs and explicit solver
  experiments, including the Stage 9d coupled-FGMRES selector and the Stage 8
  pressure-MG rollback;
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
foundation and its explicit Stage 9d selector:

- one allocation-free action returns `K_gamma*u+(D+C)^T*p` and `(D+C)*u`;
- the block action never assembles or consumes the pressure-dependent legacy
  force;
- a finest-level helper reconstructs the pressure-independent momentum RHS as
  `E->F + (D+C)^T*p - D^T*p`;
- destructive field aliases are rejected; and
- the legacy Schur FGMRES remains available as a separate rollback path.

The block-vector architecture test verifies the Stage 9c mixed-vector
foundation:

- velocity and element-pressure fields share one level-tagged lifetime;
- the nodal `neq` and pressure-index-zero sentinels are preserved;
- copy, scale and AXPY reject cross-level use;
- the algebraic Krylov metric combines a duplicate-aware velocity dot product
  with the Euclidean P0 pressure coefficient dot product used by the symmetric
  Schur map; the physical dual mass remains a separate stopping audit;
- no incompressible pressure-mean projection is silently applied to strict
  ALA, where `C^T` generally makes a constant dynamic pressure observable.

The coupled-FGMRES architecture test verifies the Stage 9d feasibility path:

- an explicit `coupled_fgmres` selector leaves legacy `fgmres` intact;
- the full residual is formed from the pressure-independent force and the
  monolithic strict-ALA block action;
- the right block preconditioner uses one velocity solve and the signed
  pressure Schur inverse approximation;
- every Krylov correction updates velocity and pressure with the same basis
  coefficient, followed by explicit block-residual replacement; and
- acceptance requires both physical continuity cancellation and the original
  unaugmented momentum residual when its audit gate is enabled.

Stage 9d.1 additionally verifies that the coupled preconditioner is an upper
block-triangular pressure-first map, and that only its velocity MG call uses a
bounded, periodically reported cycle budget.  The historical velocity-solver
entry point retains its original unbounded behavior for rollback equivalence.

Future test groups remain reserved for operator adjoint, energy-budget, and
multigrid/Galerkin validation.
