# EBA Drucker–Prager shear-heating cap

`qvis_mode=0` preserves the original heat source; `1` computes a candidate cap
but uses the original heat; `2` applies the candidate. Defaults are mode 0,
`qvis_cohesion_pa=1e7` and `qvis_friction_angle_rad=0.085`. Both standalone and
Pyre property paths validate these inputs. Modes 1/2 require EBA and the new
six-column reference state. They do not change EVi or the Stokes operator.

The EBA schema is `rho g Tref alpha Cp P_lith_Pa`, in global CMB-to-surface
order. Columns 1–5 retain their normalization; column 6 is SI Pa. The generator
integrates the same fitted phase-inclusive density and gravity from surface
pressure zero, splitting quadrature at 410/520/660 km. Every MPI rank validates
the full pressure profile before reading its own radial slice. Five-column
files remain accepted in mode 0; pressure is never reconstructed at runtime.
ALA and legacy schemas are unchanged.

At the pressure integration point, the existing strain routine returns
`S = 2 eps:eps = 4 eps_II^2`. With the arithmetic mean of the eight viscosities,
`tau_II = eta_mean sqrt(S)`. The yield criterion is interpreted in the sqrt(J2)
stress convention, and the entire stress tensor is radially scaled by
`min(1, sigma_y/tau_II)`. The SI stress scale is `eta0*kappa0/R0^2`.
No componentwise tensor clipping, extra factor of two, or viscosity floor is
applied to this thermal-only effective viscosity. Pressure is averaged over
the two radial nodes, consistent with the element-centered strain measure.
The original heat quadrature is retained so changing the cap does not also
change the baseline dissipation approximation.

## Diagnostics and energy accounting

- Mechanical `Qvisc` remains the original raw dissipation. `Qvisc_raw` is an
  explicit total alias. Existing Roperator/Rmechanical/Rheating_operator keep
  their definitions; they are not forced to close against reduced heat.
- `Qvisc_capped`, `Qvisc_used`, `Qvisc_removed` and
  `Qvisc_potential_removed` have global totals and surface-to-CMB shell
  integrals. Removed is raw minus used; potential removed is raw minus capped.
- `Qvisc_limited_volume` is a dimensionless volume (NOT power), also with
  shell integrals. Its nonzero cells mark where the candidate cap is active.
- Mechanical snapshots store their own mode, step, time and pre-rigid-rotation
  status. NPZ schema version is 4. Mechanical powers remain dimensionless
  Di/Ra-scaled integrals, not W or TW.
- Thermal-budget Qvisc and the qvisc profile are the source actually used.
  Thermal logs additionally retain raw/capped and removal totals. CBF uses
  the same capped source law for its independently evaluated state without
  overwriting the thermal-step diagnostic arrays.

The removed heat is an explicit parameterized energy deficit, not numerical
Stokes error or modeled stored energy. The cap uses total strain and cannot
isolate a purely assimilation-driven component. It does not guarantee
long-term cooling or realistic plate-boundary tractions.

## Local validation

```bash
python3 -m unittest discover -s tests/qvis -v
python3 -m unittest discover -s tests/mechanical_power -v
python3 -m unittest discover -s tests/phase_energy -v
python3 -m unittest discover -s tests/cbf -v
python3 tests/cbf/build_validation.py --build-dir /tmp/citcoms-qvis-build
```

The Qvis tests execute the production heat-source routine and refstate reader:
simple shear and invariant factors, off/diagnose/apply, zero strain/strength,
pressure dependence, invalid pressure/schema and two radial slices. Mechanical
fixtures verify unchanged raw power/operator residuals and matching NPZ/shell
budgets. The standalone build compiles and links the active solver sources;
it is not a production Pyre or 384-rank benchmark run.

Rebuild all C objects and the Python extension because All_variables changed.
The associated runs-branch Qvis LSF script performs two independent one-step
restarts at 13600, verifies stage 1 before launching stage 2, and saves native
CBF output and comparison.json. HPC execution is left to the user.
