# rheol=7/71/72: Steinberger / EBA reference viscosity family

## Selection and scope

All three modes compute the reference even with `TDEPV=off`:

| rheol | Reference correction | Active configuration |
| --- | --- | --- |
| 7 | `f(z)=1`, physical baseline | none |
| 71 | constant positive multiplier in each depth segment | `rheol7_layer_factors` |
| 72 | `f(z)=a_i+b_i*u+c_i*u^2` | `rheol7_poly_a`, `rheol7_poly_b`, `rheol7_poly_c` |

The four fixed depth segments are 0-410, 410-520, 520-660 and 660-2891 km.
`u=(z-segment_top)/(segment_bottom-segment_top)` and the deeper segment owns
each interface. The vectors are independent of `num_mat` and material IDs.
They are stored and transferred in double precision; selected vectors must
contain exactly four finite numbers. Mode 71 requires positive multipliers;
mode 72 additionally checks both endpoints and any interior quadratic minimum
to require `f(z)>0` on the whole segment. Inactive vectors are ignored.

With `TDEPV=on`, each mode
multiplies nuref by exp(s*Az*(1/T_K-1/Tref_K)), where s=cold_scale (default0.5) on the cold
side and 1 on the hot side. Actual nodal temperature is interpolated to the
same quadrature point as Tref; conversion to Kelvin uses Ttop and ref_temperature. Existing rheol options retain their behavior.
The paired runs-cmbhf_EBA/cmbhf_EBA.cfg selects 72 to preserve its historical
M2-A physics. The new `cmbhf_EBA_Q0_30_Vis7_scold1.0_LLSVPsV0.04_B24_0.5_PICES_1.8Ga_pPv`
experiment selects baseline 7, TDEPV=on and cold_scale=1.0.

The base viscosity is evaluated at viscosity quadrature points using interpolated
radius and the existing reference-state temperature:

    Tref_K = Ttop + ref_temperature * Tref_nd
    nuref = (1e21 Pa s / refvisc) * f(z) * exp(Araw(z)/Tref_K + B(z) - C)

`nuref` denotes dimensionless dynamic viscosity eta/refvisc, not kinematic
viscosity eta/rho. No instantaneous temperature, visc0 multiplier, viscE/viscT,
activation-volume or strain-rate law is introduced in this base calculation.
Az equals Araw from the accepted Figure1 fits, excluding the B/C normalization
offsets. Cold scaling affects only the anomaly exponent and leaves nuref unchanged.
The solver clips rheology-input nodal temperatures to the configured bounds
before interpolation (see below); the Kelvin kernel requires positive absolute
temperature and uses no empirical viscT offset. Nonfinite or nonrepresentable float viscosities fail
explicitly; the usual physical min/max caps are applied downstream.

## Frozen stages in lib/Steinberger_nuref.h

1. Figure 1 upper mantle: fifth-degree H polynomial in x=z/660, in kJ/mol.
2. Figure 1 lower mantle: quadratic melting-temperature polynomial in
   q=(z-660)/2231, in K. Upper H is truncated at660; lower fit extends to660.
3. Raw exponent: H*1000/(3.5*8.3144*Tref_K) above660, 12*Tm/Tref_K below660.
4. Lower-mantle log multiplier B=-7.466790490109812, upper B=0;
   log normalization C=13.163586986394256. These are the frozen accepted
   Figure4 reconstruction convention, not a new additive correction to H.
5. Mode 7 omits the empirical factor. Mode 71 defaults to identity factors.
   Mode 72 defaults to the historical M2-A factor: quadratic on 0-410 and
   660-2891 km, constant on 410-520 and 520-660 km. Its original calibration
   excluded lithosphere and basal thermal-boundary-layer anchors. Polynomial
   continuation reaches segment endpoints.

Figure 1 coefficients and B/C normalization remain frozen in the header.
The depth correction coefficients are configurable, with frozen M2-A defaults.
Sources in the parent
workspace: output/steinberger_adiabatic/fit_and_validation.json and
output/steinberger_M2_polynomial_factors/validation.json (results.A).
Original paper: Steinberger & Calderwood (2006), Figures1,4,9.
M2-A is the previously named A vector curve; its individual lmax is unassigned.
Changing Tref changes the computed profile; correction coefficients are not
automatically refitted or renormalized for a different background. The Fig.1
polynomials and chosen normalization define the common EBA baseline, rather
than an exact reproduction of every original-paper temperature model.

## Configuration examples

All entries belong in `[CitcomS.solver.visc]`. Choose one mode; a different
coefficient set for the same correction type does not need another rheol number.

```
# Physical baseline; the example cold-side response matches the new Vis7 run.
rheol = 7
TDEPV = on
cold_scale = 1.0

# For layer constants, replace the rheol line above with 71:
# rheol = 71
rheol7_layer_factors = 1.0,1.0,1.0,1.0

# For historical M2-A, replace the rheol line with 72:
# rheol = 72
rheol7_poly_a = 5.067917152408204,0.22571234789343927,0.3127706714585134,2.073829265582211
rheol7_poly_b = -8.15474461356726,0.0,0.0,-0.7521860487837174
rheol7_poly_c = 4.799095524666932,0.0,0.0,0.09721318357053432
```

For a mode-72 identity correction use a=1,1,1,1 and b=c=0,0,0,0.
Mode 71 with all factors equal to one and mode 72 with identity coefficients
reproduce mode 7, including its shared temperature response.

## Migration and restart identity (2026-11-18)

Before this change, `rheol=7` always included M2-A. Existing experiments must
select `rheol=72` to preserve that profile; merely rebuilding an old input with
`rheol=7` changes its physics. The 49 tracked top-level EBA run templates were
explicitly migrated to 72; their names and other settings are retained. The
13 historical P7/P8 gate fixtures were not rewritten: copy and migrate them to
72 when reproducing their former M2-A physics. The runs worktree records the
changed templates in `RHEOL_FAMILY_MIGRATED_CONFIGS.txt`.

PICES recognizes all three modes as coupled, rebuildable Newtonian rheologies.
Its checkpoint fingerprint now includes `rheol_family_v1` and the coefficients
active in 71 or 72. Changes to the model or coefficients are rejected on restart.
Historical checkpoints are not automatically migrated by changing their cfg to
72; use their original build or a separately audited checkpoint migration.
The new Vis7 LSF is a fresh experiment and checks the installed binary marker
and Python inventory before creating a run directory.

## Numerical and integration limits

The kernel requires depth0–2891 km, positive finite Tref_K and refvisc. It
allows 1e-5 km endpoint roundoff and rejects invalid results. Coefficients are
evaluated in double precision; the existing viscosity array stores floats.
Depth determines the phase segment directly, independent of material IDs.
Existing nodal Tref interpolation is retained: a cell crossing a phase interface
can therefore interpolate its temperature jump. The hardcoded viscosity-factor
jump is not smoothed by this kernel. Exact one-sided thermal evaluation at
quadrature points would be a separate reference-state interpolation change.

The normal get_system_viscosity pipeline still applies any enabled composition,
strain-rate, plasticity, channel corrections and min/max limits afterward.
visc_from_PB currently has no case7, hence its legacy case3 weakening does not
execute for this option. For a pure reference-profile comparison keep CDEPV,
SDEPV,PDEPV and low-viscosity channel/wedge off; inspect limits/smoothing too.

## Validation

    python3 tests/test_rheol7_nuref.py /path/to/CitcomS/output

Compiles the actual production header with warnings-as-errors, compares mode 72
with the 65-node historical M2-A profile, and exercises baseline/identity/custom
corrections, interfaces, temperature response and invalid coefficients. The
runtime tests cover the standalone parser, Pyre property bridge and production
dispatch including TDEPV=off; see [the validation report](rheol7_family_validation.md) for the
latest commands and integration results. A local standalone MPI check is not
a production HPC/Pyre installation or a full convection experiment.

## Initial-temperature validation (70 Ma)

Run `python3 tests/validate_rheol7_initial_profile.py WORKSPACE_ROOT CFG`. It uses the
current cfg and refstate_EBA.txt, the same erfc anomaly construction as
Lith_age.c/init_validate.py, and the compiled production C kernel. Surface
age is70 Ma with the cfg's full lith_age_depth (254.84 km in the new Vis7 experiment); the trench-specific half-depth
branch and initial slab anomalies are outside this representative radial column.
Bottom TBL thickness is approximately500 km at the1% amplitude height, with
Tbottom=3700 K. Geological start_age (1799.9 Ma in the new Vis7 experiment) is distinct from plate thermal age.

The diagnostic evaluates65 nodes and128 radial Gauss points, checks reference
identity, the selected cold-side exponent, hot weakening and the Kelvin formula, and applies
the current depth-dependent viscosity limits. Actual production uses Gauss-point
viscosity; nodal limited values in the plot are illustrative. Outputs are in
output/rheol7_initial_validation/rheol_MODE/cold_scale_VALUE/CFG_STEM.
This is not a full3D initialization or MPI solve.

`cold_scale` is declared in the Pyre Visc inventory, transferred to the C state
by `pyCitcom_Visc_set_properties`, and read by the standalone C parser. It must be
finite and nonnegative; 0 disables cold strengthening, 1 restores full Az.
The hot-side exponent is independent of cold_scale. Rebuild the EBA solver including
the Python component before using the new cfg/lsf; jobs are not auto-submitted.

## Runtime binding and failure diagnostics

Both input paths validate cold_scale and print the selected rheol, TDEPV,
cold_scale, hot_scale, `rheol_family_v1` and the active coefficients at startup.
The property bridge also records cold_scale in its normal
parameter output. This fixes an omitted Pyre-to-C assignment that left the
malloc-allocated cold_scale field uninitialized.

The existing limits remain after temperature/composition/strain-rate/plasticity/
plate-boundary/channel corrections, before the two GP-to-node-to-GP filtering
cycles. The maximum is visc_max where element-center radius > 0.89641 and
5*visc_max below; the minimum is visc_min. These physical limits are not relocated.
For rheol7, finite nodal temperatures used for viscosity are clipped inline to
the nondimensional interval [0,1] before Gauss-point interpolation, matching
rheol3 and the last cluster-tested implementation. E->T and Tref are unchanged.
The dimensional mapping remains cfg-driven:
T_K = Ttop + (Tbottom-Ttop)*T_nd. The bounds therefore map to 300 and 3700 K
for the current experiments, but those Kelvin values are not hardcoded.
Nonfinite temperatures still fail validation. This protects rheology from finite
thermal overshoots; it does not fix overshoots in the evolved temperature field.
Invalid or unrepresentable raw viscosity
still fails before these limits; the kernel's fit and temperature formula are
unchanged.

Rheol7 failures report rank, step, local cap/element/Gauss point, depth, rheology-input
and reference temperature, Az, cold_scale, nuref, raw viscosity and its natural
logarithm, plus element nodal temperatures. Runtime diagnostics go to stderr and
the rank log; MPI_Abort terminates the solver communicator so other ranks do not
remain waiting for the failing rank. Other rheologies retain their error paths.

Failures additionally carry the marker `rheol7_clip_diag_v1`, temperature scales,
boundary types/values and the effective [0,1] clip bounds. Each node reports its
original and clipped temperature plus
the Gauss interpolation weight. A nonfinite interpolated rheology temperature
is identified before calling the viscosity kernel. Finite out-of-range inputs
are clipped; nonfinite inputs remain nonfinite and abort.

Run `python3 tests/test_rheol7_runtime.py` for the production bridge and diagnostic
helpers with mocked Python property APIs and real two-process MPI abort tests.
This does not replace a build and end-to-end run with the cluster's Python2/Pyre.

Temperature scaling is already cfg-driven: DeltaT=Tbottom-Ttop and the offset
is Ttop/DeltaT. The solver uses T_K=Ttop+DeltaT*T_nd for both actual and reference
temperatures; 300, 3400 and 0.0882353 are not production rheol7 constants.
refvisc, cold_scale and radius likewise come from solver inputs. The 1e21 Pa s
factor in nuref is the dimensional calibration of the fitted profile, divided by
the configurable refvisc; replacing that calibration with refvisc would change
the physical fit. The fitted phase depths, Az, B/C and M2 coefficients remain
frozen intentionally, as do the pre-existing downstream viscosity-limit rules.

Run `python3 tests/test_rheol7_temperature_clip.py` to check the production nodal
interpolation block with finite overshoots, unchanged in-range inputs, and
nonfinite temperatures, including the nodal temperatures from the step3 failure.
