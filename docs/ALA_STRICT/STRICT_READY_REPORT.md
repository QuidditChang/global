# Strict ALA Production-Readiness Report

## Status

`runs/cmbhf_ALA_strict.cfg` is the production strict-ALA configuration. Its
physical inputs remain aligned with `runs/cmbhf_ALA.cfg`; the active
strict-only differences are:

```text
refstate_file = refstate_ALA_strict.txt
ala_beta_element_source = interval
ala_beta_interval_file = interval_ALA_strict.txt
ala_shallow_patch_velocity_solver = diagonal
ala_pcg_restart_interval = 20
ala_outer_solver = fgmres
ala_unaugmented_momentum_tolerance = 1e-3
ala_geneo_preconditioner = off
ala_geneo_basis_type = radial_partition
```

The first three select the phase-inclusive strict reference state.  The final
five retain the validated Stage 6b.3 Schwarz metric, select the Stage 7
FGMRES path, require Stage 7b joint momentum/continuity acceptance, and
disable the negative Stage 6d coarse-space A/B.
`radial_partition` remains configured only as an inert rollback.  No
energy-equation, phase-equation, boundary-condition, assimilation, timestep,
or viscosity parameter is changed.

## A. Reference-state source

`refstate_ALA_strict.txt` is generated from Katsura (2022) Table S5 by the
strict-only Python pipeline. Density is phase-inclusive. Density and beta are
evaluated directly at the 65 CitcomS radial nodes from the same piecewise
analytical `ln(rho)` model.

The serialized schema is:

```text
rho g Tref alpha Cp beta Gamma_eff Ks
```

The first seven fields retain the legacy nondimensionalization and ordering.
`Ks` is appended in physical GPa. It is stored by the reader but is not used by
the current Stokes or energy formulation.

## B. Tref source and usage

Column 3 is the fixed thermodynamic reference temperature. It comes from the
Katsura Table S5 `T_K` values through the same four-branch endpoint-constrained
cubic fitting and direct radial-node evaluation used for the other fitted
thermal properties.

The 300 K surface and 3800 K CMB total-temperature boundary values are not
written into Tref. With `Ttop=300 K` and `Tbottom=3800 K`, the generated endpoint
values are:

| Endpoint | Tref (K) | serialized Tref* |
|---|---:|---:|
| CMB | 2607.653828 | 0.6593296652 |
| surface | 1618.046307 | 0.3765846591 |

The reader stores column 3 in `E->refstate.temperature`; `E->refstate.Tref` is
a non-owning semantic alias to the same storage. The phase reference uses
`E->refstate.Tref`.

## C. Initial and dynamic temperature

`E->T` remains the sole total-temperature state and the energy-equation
unknown. The initialization path is inherited from `cmbhf_ALA`:

```text
Katsura Tref
  + bottom thermal boundary layer
  + surface half-space cooling / plate-age structure
  + slab and configured perturbations
  -> E->T(x, 0)
  -> total-temperature boundary conditions
```

No `DataT` field is allocated or cached. A temperature anomaly remains
well-defined when needed as the derived quantity `E->T-Tref`; it is not an
independent state variable.

Assimilation continues to update `E->T` through the existing
`assim_delta_T` path. It has no write path to `E->refstate.temperature` or
`E->refstate.Tref`.

## D. Phase-density closure

The strict reference density already contains the equilibrium 410, 520, and
660 km density transitions. The absolute dynamic phase fraction remains
`phase_B = X(E->T)` and is still used for phase-boundary output and latent
heating.

The momentum and matching mechanical-power/profile diagnostic contribution is:

```text
-Ra_phase * (X(E->T) - Xref(Tref))
```

Therefore `Ttotal=Tref` gives zero dynamic phase buoyancy, and the equilibrium
phase density is not counted a second time.

## E. Nondimensional consistency

The production strict cfg preserves `rho0`, `g0`, `alpha0`, `Cp0`, `Ttop`,
`Tbottom`, radius, diffusivity, viscosity, Rayleigh number, and dissipation
number from `cmbhf_ALA.cfg`.

The generator-reader contract is:

| Field | Serialized form |
|---|---|
| rho | `rho/rho0` |
| g | `g/g0` |
| Tref | `(Tref_K-Ttop)/(Tbottom-Ttop)` |
| alpha | `alpha/alpha0` |
| Cp | `Cp/Cp0` |
| beta | `R0*beta_SI` |
| Gamma_eff | `Di*alpha*g/(Cp*beta)` using serialized fields |
| Ks | physical GPa, diagnostic storage only |

The strict generator and validator report exact serialized-profile closure for
rho, beta, Gamma_eff, and Ks.

## F. Stage 6c shallow Schur experiment

Stage 6b.3 used `G diag(K)^-1 G^T` inside every cached shallow pressure
patch.  Stage 6c retains the same `6x6x2` patch geometry, overlap-three MPI
halo, global partition of unity, pressure shift, and PCG outer iteration, but
replaces the scalar velocity diagonal by assembled symmetric `3x3` same-node
velocity blocks.

Each node block is Cholesky factored as `K_n=L_n L_n^T`.  Pressure-gradient
features are formed as `phi=L_n^-1 g`; patch entries are dot products
`phi_i^T phi_j`.  Local and exchanged halo patches are therefore Gram
matrices and remain positive semidefinite before the existing pressure shift.
Nodes whose coupled factor cannot be formed fall back to the exact Stage 6b.3
diagonal metric and are counted in `velocity_block_fallbacks`.  Setting
`ala_shallow_patch_velocity_solver=diagonal` is the complete rollback.

Every PCG iteration now reports the absolute pressure-mass continuity norm
`mass_norm=N_M`, `Q`, and their relative scales.  The final original,
unaugmented momentum residual is audited before a failed continuity exit as
well as after a successful solve.

The 400-rank Stage 6c A/B run was numerically healthy but negative:
`node_block` changed the best cancellation only from `0.01339024` to
`0.01338443` at iteration 42 and ended at `0.01558697`.  It had zero node
fallbacks, positive curvature, unit residual drift, and a `3.0e-14`
preconditioner symmetry defect.  The strict production experiment therefore
returns to the exact `diagonal` rollback.

## G. Stage 6d cross-rank radial aggregate

Stage 6d retains the Stage 6b.3 `6x6x2`, overlap-three shallow Schwarz map and
adds a deterministic coarse space on each `2x2` horizontal MPI-rank group.
With four bins per member rank and two radial bins, each processor aggregate
has `8x8x2` selection bins and crosses four ranks.

Unlike the earlier GenEO experiment, the two coarse shapes are not selected
by an eigenvalue threshold.  They are the disjoint indicators of the two
radial portions of the 0--410 km shallow layer.  Their span contains the
complete aggregate constant while retaining an independent shallow/deep
radial amplitude.  The supports form an exact disjoint partition across
processor aggregates.  Normalization is deliberately deferred until after
transfer to the Galerkin pressure level.  Construction is purely geometric
on all `8x8x2` bins: neither the legacy GenEO `active_map` nor a selection
Rayleigh quotient may suppress a radial indicator.  The authoritative
validity checks are complete cross-rank support after transfer and the final
Galerkin Cholesky factorization.

Processor groups are separated by radial rank coordinate.  Only groups whose
local radial interval intersects the configured 0--410 km shallow layer own
these two modes.  Deeper groups have zero shallow coarse modes; requiring a
nonzero radial indicator on them would manufacture an identically zero basis
vector.

The coarse matrix remains

```text
S_c = P^T (D+C) Kfixed^-1 (D+C)^T P
```

and is assembled by the existing fixed strict-ALA Galerkin action, globally
symmetrized, shifted, and Cholesky factored.  Consequently the additive
coarse correction is fixed and SPD.  `ala_geneo_basis_type=spectral` is the
legacy rollback; `radial_partition` requires a cross-rank group and exactly
one configured mode per radial bin.

The completed 400-rank Stage 6d run constructed all 96 modes in 18.85 seconds
and passed its MPI support, symmetry, Galerkin, and Cholesky audits.  It did
not improve convergence: cancellation was `0.029719` at iteration 30,
`0.013451` at the best iteration 42, and `0.015685` at iteration 50, compared
with `0.029556`, `0.013390`, and `0.015584` for the Stage 6b.3 baseline.
Consequently the production experiment disables this correction.

## H. Stage 6e PCG residual-replacement interval

In both Stage 6b.3 and Stage 6d, explicit PCG replacement produced the largest
single-step improvements: `0.038063 -> 0.028689` after iteration 20 and
`0.019184 -> 0.013780` after iteration 40 in Stage 6d.  Stage 6e changes only
`ala_pcg_restart_interval`, from 20 to 10, with GenEO disabled.  The acceptance
target remains cancellation at most `0.01` within 30 iterations, with positive
curvature, recursive/explicit drift near one, stable velocity norm, and no
degradation of the unaugmented momentum residual.

The completed run rejected this hypothesis.  Although the first iteration
after each replacement improved, cancellation at iteration 30 increased from
`0.029556` to `0.030423`; the best value degraded from `0.013390` to
`0.015231`, and iteration 50 degraded from `0.015584` to `0.020312`.
Consequently fixed ten-step replacement is not retained.

## I. Stage 7 momentum-consistent FGMRES

Stage 7 restores a restart length of 20, keeps GenEO disabled, and changes
only the outer Krylov method to FGMRES.  Before Arnoldi starts, the solver
applies the same momentum-consistent initial velocity correction as PCG.
Every flexible pressure basis vector retains its matching velocity correction,
so the reconstructed `P` and `V` remain coupled.  Each completed restart
assembles the explicit strict continuity residual and audits the original
unaugmented momentum equation using the existing FGMRES work vectors.

The former `drift` diagnostic divided Arnoldi reduction by physical
cancellation even though the two quantities have different denominators.
Stage 7 now defines drift as Arnoldi recursive algebraic residual divided by
the explicitly assembled algebraic residual.  Iteration records also include
`mass_norm`, `mass_relative`, `Q`, `Q_relative`, term strength, and velocity
norm.  Reaching the physical cancellation target now exits the active Arnoldi
cycle immediately after its diagnostics instead of allowing later vectors in
the same cycle to move an already accepted state away from tolerance.  The
run must compare these quantities and the restart momentum audit; cancellation
remains the only physical stopping criterion.

The first 400-rank Stage 7 run reached `cancellation=0.00992642` at iteration
46.  Its recursive and explicit algebraic reductions were both
`2.992918e-03`, `drift=1`, the velocity norm remained at `0.9965` of its
initial value, and the final unaugmented momentum residual was
`3.932077e-06` relative to its audited reference.  The process then aborted
while releasing memory.  This was not an out-of-memory or convergence failure:
nodal velocity assembly writes a zero sentinel at index `neq`, while the new
FGMRES velocity work and correction-basis arrays had been allocated with only
`neq` entries.  They now use the established `neq+1` velocity-vector layout,
and the static regression checks that allocation contract.

## J. Stage 7b joint saddle-point acceptance

The first post-fix evolving run completed the initialization solve and at
least two subsequent Stokes solves without the former allocator abort.  The
initial solve reached cancellation `0.00992642` in 46 iterations with
unaugmented momentum relative residual `3.932077e-06`.  Warm-started solves
then reached cancellation `0.00991258` and `0.00972404` in only 10 and 13
iterations, but their unaugmented momentum residuals remained `2.768549e-02`
and `2.182796e-02`.  Recursive and explicit Schur residuals still agreed and
velocity norms remained stable; the failure was cancellation-only acceptance,
not Arnoldi drift or velocity collapse.

The cause is the state-dependent denominator
`Q=|| |D|+|C| ||`: its relative value was `0.230881` at the end of the
initial solve but `0.932805` and `0.851919` in the warm-started solves.  A
fixed cancellation ratio therefore admitted very different absolute `Gu`
errors and, with `gamma=10`, left a visible augmented-to-unaugmented momentum
defect.

Stage 7b retains `cancellation < 1e-2` and requires the existing audited
unaugmented momentum relative residual to be at most `1e-3`.  Its numerator is
the original-equation defect `||Ku+G^Tp-f||`; its established normalization is
the assembled strict-force norm `||f-C^Tp||`.  FGMRES performs the expensive
original-momentum audit only after cancellation first qualifies.  A failed
momentum candidate continues the current Arnoldi cycle; it is no longer
reported as a converged Stokes state.  The new parameter defaults to zero, so
legacy and `gamma=0` configurations retain their previous stopping behavior
unless they explicitly enable the gate.

## K. Current thresholds, evidence, and production target

The active Stage 7b numerical controls are:

| Control | Value | Role |
|---|---:|---|
| `tole_compressibility` | `1e-2` | physical `D+C` cancellation target |
| `ala_unaugmented_momentum_tolerance` | `1e-3` | original momentum acceptance gate |
| `ala_inner_accuracy_max` | `1e-4` | fixed upper bound for velocity solves |
| `ala_inner_accuracy_factor` | `1e-2` | adaptive inner/outer forcing factor |
| `ala_pcg_restart_interval` | `20` | FGMRES Arnoldi basis length |
| `piterations` | `50` | current diagnostic outer budget |
| `ala_schur_symmetry_tolerance` | `1e-3` | fatal SPD symmetry guard |
| `ala_feasibility_window` | `20` | progress-audit window, not acceptance |
| `ala_feasibility_min_reduction` | `0.02` | progress-audit floor, not acceptance |

The observed Schwarz symmetry defect remains about `1e-14`, recursive and
explicit FGMRES residuals agree to reported precision, drift is one, and the
velocity norm is stable.  The remaining open result is performance under the
new joint gate: the initialization solve already satisfies both equations but
requires 46 iterations, while the former 10--13 iteration warm-start results
must now continue because their momentum defects exceed `1e-3`.

Production acceptance requires both equations to meet their thresholds in
20--30 pressure iterations, positive curvature/SPD diagnostics, recursive and
explicit residual agreement, drift near one, stable velocity scale, auditable
unaugmented momentum balance, and no behavior change when the new gate is zero
or the augmented-Lagrangian path is disabled.

## L. Stage 7c augmented/unaugmented momentum decomposition

The first Stage 7b warm-start run confirmed that the gate works but did not
explain the rejected momentum state.  From iteration 10 through iteration 46,
cancellation decreased from `9.912579e-03` to `4.542411e-03`, mass norm from
`1.912743e+02` to `8.764589e+01`, and algebraic relative residual from
`3.407300e-02` to `1.599477e-02`.  Over the same interval, unaugmented
momentum relative residual remained `2.76855e-02` to reported precision.
Drift remained one and velocity relative norm remained about `0.99987`.

Stage 7c therefore decomposes, at FGMRES startup and restart boundaries,

```text
r_raw = f - G^T p - K u
r_aug = f - G^T p - K_gamma u
r_pen = gamma G^T M_p^-1 G u
r_aug = r_raw - r_pen
```

The log now reports raw, augmented, and penalty RMS/relative norms;
`raw_penalty_cosine`; the independently assembled operator split defect
`||K_gamma u-Ku-r_pen||`; and the residual-norm identity defect.  If
`r_aug` is small, `r_pen` matches `r_raw`, cosine is near one, and both
defects are near roundoff, the warm-start result is an expected finite-`Gu`
augmented state and the Schur solve must reduce the penalty action further.
If `r_aug` is large, pressure/velocity reconstruction is not preserving the
augmented momentum equation.  A large split defect instead identifies an
operator assembly mismatch.

The audit reuses the existing FGMRES `tmpF/tmpU` vectors and allocates no new
global field.  Candidate acceptance retains the cheaper raw audit; the full
decomposition runs only at startup and restart boundaries.  It changes no
equation, Krylov update, preconditioner, threshold, or cfg value.

## M. Stage 7d force/stiffness consistency repair

The Stage 7c production run made the warm-start failure unambiguous.  The
initial solve started with raw and penalty relative residuals of
`2.739295e-3`, an augmented relative residual of `9.865294e-7`, and a split
defect of `1.736952e-14`; it reached joint acceptance at iteration 46.  The
next Stokes solve instead started with raw and augmented relative residuals
both near `2.76855e-2`, while the penalty relative residual was only
`3.727995e-4`.  By restart 40 the penalty had fallen to `1.727965e-6`, but
the raw and augmented defects remained `2.76855e-2`.  The split defect stayed
below `5.5e-11`.  The Schur iteration was therefore reducing continuity while
preserving a different, internally consistent momentum right-hand side.

The cause was the call order in `general_stokes_solver()`.  `assemble_forces()`
ran before a permitted viscosity update and before
`construct_stiffness_B_matrix()`.  The assembled force contains nonzero
velocity-boundary contributions formed from the element stiffness.  Updating
viscosity and rebuilding `K_gamma` consequently paired a stale force with the
new operator during `initial_vel_residual()`.  The Stage 7c audit reassembled
the current force and exposed the resulting fixed augmented residual.

Stage 7d refreshes the force immediately after a stiffness rebuild only when
strict pressure buoyancy and positive augmented-Lagrangian gamma are active.
The refresh precedes `solve_constrained_flow_iterative()` and emits:

```text
ALA strict force reassembled after stiffness update cycle=... gamma=...
```

The guard deliberately preserves the historical call sequence for legacy
solvers, TALA, `gamma=0`, and cycles with no stiffness rebuild.  It adds no
configuration parameter and changes no stopping threshold.  On the next HPC
run, a corrected warm start should have augmented relative residual near the
inner velocity tolerance; raw residual should be explained predominantly by
the AL penalty.  Joint convergence in 20--30 iterations remains a runtime
acceptance requirement, not a result inferred from the static repair.

The subsequent 400-rank Stage 7d run verified the repair.  The initialization
solve retained its 46-iteration trajectory, while four completed evolving
Stokes solves reached joint acceptance in 11, 13, 12, and 11 iterations.
Their startup augmented momentum relative residuals were `9.7e-7` to
`9.9e-7`; raw residuals of `2.0e-4` to `3.6e-4` aligned with the AL penalty.
Accepted unaugmented momentum relative residuals were `2.4e-6` to `5.7e-6`.
The production convergence objective is therefore achieved for these warm
steps.

## N. Stage 8 recursive pressure multigrid A/B

Stage 8 is a performance experiment motivated by the one-to-two pressure
iterations of the former incompressible EBA Stokes solve.  It does not replace
the successful Stage 7d FGMRES baseline.  The new preconditioner traverses the
existing five distributed pressure grids recursively.  At level `l` it
applies the inexpensive complete strict-ALA Gram approximation

```text
A_l = (D_l+C_l) diag(K_gamma,l)^-1 (D_l+C_l)^T.
```

Two damped Jacobi pre- and post-sweeps use the existing level BPI.  Residual
restriction is the sum `P^T` over each `2x2x2` parent, prolongation is the
matching constant `P`, and the lowest `2x2x2` local pressure grid receives 20
Jacobi steps.  Unlike the former two-level Galerkin polynomial, a coarse
operator application stays on its own mesh level and never invokes a fine
velocity MG solve.  Persistent pressure and velocity workspaces avoid
per-iteration global allocation.

The V-cycle is blended at weight `0.5` with the Stage-6b.3 MPI-overlap Schwarz
baseline and is available only to strict ALA with positive gamma and FGMRES.
The latter restriction avoids asserting PCG symmetry for the rediscretized
hierarchy.  Every application audits `r^T z_mg`; a non-finite or nonpositive
value emits `ALA_PRESSURE_MG_ROLLBACK` and preserves the Schwarz result.  At
diagnostic iterations, `ALA_PRESSURE_MG` reports V-cycle/base energy, the
approximate residual reduction, weight, and cumulative level-operator count.

The first HPC A/B must compare both iteration count and wall time against the
verified 11--13 iteration baseline.  Fewer pressure iterations alone is not a
win if the added level operators increase total Stokes time.  The cfg rollback
is a single change: `ala_pressure_multigrid=off`.

## Validation record

- Strict generator: PASS.
- Strict validator: PASS.
- Schema: 65 rows by 8 columns.
- rho analytical-profile relative RMS: 0.
- beta analytical-profile relative RMS: 0.
- interval beta-integral identity maximum: 0.
- Gamma_eff closure relative RMS and maximum: 0.
- Ks fitted-profile relative RMS: 0.
- strict static regressions: 24/24 PASS.
- changed C translation units: `mpicc -fsyntax-only` PASS.
- isolated full source compilation: all C objects and `libCitcomS.a` built;
  `CitcomSFull` linked manually because the imported Automake 1.9 templates do
  not expand `AR`/`LIBTOOL` with the locally available toolchain.

The production Pyre/LSF launcher and its Python 2.6/Linux runtime are not
executable on the local macOS host. The raw legacy C parser is not a substitute
for that launcher and is incompatible with the sectioned Pyre cfg. Consequently
the final 384-rank scheduler startup remains an HPC deployment check; no
reference-column mismatch was observed by generation, schema validation, or C
compilation.
