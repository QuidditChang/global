# EBA mechanical power diagnostics

Enabled automatically after each completed EBA Stokes solve
(`eba_formulation=1`, `ala_pressure_buoyancy=0`). The ALA diagnostic path is
unchanged. This is read-only instrumentation: viscosity, heating and the
momentum/energy equations are not modified.

## Output and timing

`global.log` contains `MECHANICAL_POWER ... formulation=EBA`, followed by
full-precision `TERM TOTAL` rows. These are global integrals, not local maxima.
The snapshot is taken before optional rigid-rotation removal. It uses the
completed solve's U, P, viscosity and body force. The heat-source diagnostic
`THERMAL_BUDGET` is produced elsewhere in the timestep and must not be assumed
to describe this identical state.

When profiles are enabled, `profiles_hist_<step>.npz` (schema version 3) retains
all existing fields and includes:

- `mechanical_valid`: 1 if a Stokes snapshot is available; otherwise 0, with
  no power values written. Unavailable diagnostics are not recorded as zero.
- `mechanical_step`, `mechanical_elapsed_time`: the solve's own timestamp.
- `mechanical_pre_rigid_rotation=1`.
- `mechanical_scale_Di_over_Atemp`.
- `mechanical_<TERM>`: exactly the same total as the corresponding log row.
- `mechanical_<TERM>_shell_integral`: surface-to-CMB shell integrals for Wbody,
  Qvisc, thermal/chemical/phase work and pressure work, aligned with
  `element_depth_km` and `element_depth_edges_km`.

Shell powers are already integrated over their volumes. Sum them directly;
do not multiply by `element_volume_km3` again. Pplate is a surface integral and
is deliberately not assigned an artificial radial volume distribution.
No histograms of mechanical work are introduced. NPZ snapshot metadata remains
that of the solve even if output occurs after the thermal state has advanced.
Old NPZ files cannot acquire these missing quantities by conversion.

## Definitions

All totals are dimensionless powers scaled by Di/Atemp to the thermal source
integral convention. They are not watts or TW. Positive boundary/body powers
supply mechanical energy to the mantle.

Let K be the actual assembled Stokes velocity operator, G the unstripped
pressure-gradient operator, f the assembled body plus applied traction loads,
and s=Di/Atemp. Prescribed-velocity reactions are `K U + G P - f`.

| TERM | Definition |
| --- | --- |
| Pplate | s times U dotted with reaction at prescribed surface velocity DOFs |
| Pother | Same at other prescribed velocity DOFs |
| Wtraction | s times U dotted with applied traction loads |
| Wbody | s times U dotted with the actual element body-force loads |
| Dvisc_operator | s U^T K U; includes any terms present in the chosen operator |
| Qvisc | Recomputed using the EBA thermal-source strain/viscosity formula |
| Wthermal, Wchemical | Thermal and chemical body-force work |
| Wphase_410/520/660 | Individual phase buoyancy work, not latent heat |
| Wphase_total | Sum of the three phase body-force work terms |
| Wpressure | -s U^T G P; no ALA beta term is included in EBA |

The force decomposition follows `get_buoyancy`, including reference phase
subtraction, gravity and horizontal mean removal. Work uses the same spherical
basis rotation and quadrature as `get_elt_f`, not a product of separate radial
averages. Wbody is computed independently from the actual force assembly.

Define `input = Pplate + Pother + Wtraction + Wbody + Wpressure`:

- `Roperator = input - Dvisc_operator`: discrete Stokes power residual.
- `Rmechanical = input - Qvisc`: residual relative to the thermal heating formula.
- `Rbody_split = Wbody - Wthermal - Wchemical - Wphase_total`.
- `Rheating_operator = Dvisc_operator - Qvisc`.

Thus `Rmechanical = Roperator + Rheating_operator`. A nonzero Rmechanical alone
must not be attributed to solver convergence; the heating approximation and
operator may differ. In particular, operator penalty terms are not physical
heating. Relative residuals should use a non-cancelling scale, e.g. the sum of
absolute input terms plus abs(Dvisc_operator), with an appropriate small floor.
Neither Qvisc-Qadi_base nor Qtotal is used as a substitute for Pplate.

## Validation

Run from the EBA submodule:

```sh
python3 -m unittest discover -s tests/mechanical_power -v
python3 -m unittest discover -s tests/phase_energy -v
mpicc -std=gnu99 -Wno-deprecated-non-prototype -Wno-unused-command-line-argument \
  -fsyntax-only -Ilib lib/Drive_solvers.c lib/Profile_output.c
```

The mechanical fixture executes the production force/gradient/divergence
assembly, mechanical diagnostic and NPZ writer with a manufactured element,
a deterministic stiffness action and rotated basis. It tests closed power
balance, deliberate free-DOF imbalance, a pressure offset, stationary flow,
Di=0, absent snapshots, shell-to-global sums and log/NPZ equality. It is a
single-process regression, not a global MPI convergence benchmark.

For a production acceptance check, rebuild the EBA binary and run a short
restart with the same physics and profile output enabled. Check that:

1. Log and NPZ totals agree at `mechanical_step` and `mechanical_elapsed_time`.
2. Each shell-power sum equals its global total.
3. Rbody_split is small relative to the force-work terms.
4. Roperator decreases appropriately with tighter Stokes tolerances.
5. Rheating_operator is inspected separately rather than hidden in Roperator.
6. Compare identical snapshots across MPI decompositions before interpreting
   small residuals as physical effects.

The new code performs extra diagnostic operator applications/quadrature and
MPI reductions once per completed solve. Its production runtime overhead and
global MPI closure have not been measured by the local fixture.
