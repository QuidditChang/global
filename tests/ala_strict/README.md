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
- strict rheology shape-interpolates `E->refstate.Tref` to each viscosity
  Gauss point instead of constructing a lith-age reference profile;
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

Stage 9d.4 refreshes the original momentum decomposition at every coupled
restart and before any final summary.  Sparse block-action audits report the
velocity and pressure components of `r`, `z`, and `Az`, their componentwise
defects, and alignment, so a restarted-FGMRES tail can be assigned to a block
without changing the Krylov recurrence.

Stage 9d.5 optionally refines the complete triangular block inverse with
`z <- z + P_tri^-1(r-Az)`.  The defect uses the exact coupled operator, so it
retains all `D+C` cross terms; setting the correction count to zero reproduces
the Stage 9d.4 preconditioner exactly.

Stage 9e removes the separate initial `K_gamma` momentum solve only from the
monolithic path.  Coupled FGMRES starts from the supplied `(u,p)` and reduces
the force-scaled momentum and continuity components together.  The Stage 7
pressure-only FGMRES retains its momentum-consistent initialization.

Stage 9d.1 additionally verifies that the coupled preconditioner is an upper
block-triangular pressure-first map, and that only its velocity MG call uses a
bounded, periodically reported cycle budget.  The historical velocity-solver
entry point retains its original unbounded behavior for rollback equivalence.

Stage 9f.0 adds a one-time five-level runtime contract audit for `K_gamma`
symmetry, `G/G^T` adjointness, velocity and P0-pressure transfer adjointness,
coarse beta, pressure mass, and duplicate ownership.  The P0 pair is constant
prolongation with its exact Euclidean transpose.  A mixed patch callback
contract reserves Stage 9f.1 for a real velocity-pressure Vanka kernel; it does
not relabel the existing velocity-only smoother as coupled.

Future test groups remain reserved for energy-budget and coupled-V-cycle
Galerkin validation.
