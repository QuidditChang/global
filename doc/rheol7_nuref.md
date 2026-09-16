# rheol=7: frozen Steinberger / EBA / M2-A reference viscosity

## Selection and scope

Set `rheol = 7` in `[CitcomS.solver.visc]`. This option computes the reference even with `TDEPV=off`. With `TDEPV=on`,
it multiplies nuref by exp(s*Az*(1/T_K-1/Tref_K)), where s=cold_scale (default0.5) on the cold
side and 1 on the hot side. Actual nodal temperature is interpolated to the
same quadrature point as Tref; conversion to Kelvin uses Ttop and ref_temperature. Existing rheol options retain their behavior.
The paired runs-cmbhf_EBA/cmbhf_EBA.cfg now selects rheol=7 and TDEPV=on.

The base viscosity is evaluated at viscosity quadrature points using interpolated
radius and the existing reference-state temperature:

    Tref_K = Ttop + ref_temperature * Tref_nd
    nuref = (1e21 Pa s / refvisc) * f_M2A(z) * exp(Araw(z)/Tref_K + B(z) - C)

`nuref` denotes dimensionless dynamic viscosity eta/refvisc, not kinematic
viscosity eta/rho. No instantaneous temperature, visc0 multiplier, viscE/viscT,
activation-volume or strain-rate law is introduced in this base calculation.
Az equals Araw from the accepted Figure1 fits, excluding the B/C normalization
offsets. Cold scaling affects only the anomaly exponent and leaves nuref unchanged.
No input temperature clipping or empirical viscT offset is used: positive absolute
temperature is required. Nonfinite or nonrepresentable float viscosities fail
explicitly; the usual physical min/max caps are applied downstream.

## Frozen stages in lib/Steinberger_nuref.h

1. Figure 1 upper mantle: fifth-degree H polynomial in x=z/660, in kJ/mol.
2. Figure 1 lower mantle: quadratic melting-temperature polynomial in
   q=(z-660)/2231, in K. Upper H is truncated at660; lower fit extends to660.
3. Raw exponent: H*1000/(3.5*8.3144*Tref_K) above660, 12*Tm/Tref_K below660.
4. Lower-mantle log multiplier B=-7.466790490109812, upper B=0;
   log normalization C=13.163586986394256. These are the frozen accepted
   Figure4 reconstruction convention, not a new additive correction to H.
5. M2-A factor: quadratic in u=(z-segment_top)/(segment_bottom-segment_top)
   on 0–410 and660–2891 km, constant on410–520 and520–660 km.
   The lithosphere and basal thermal-boundary-layer anchors are excluded.
   Polynomial continuation reaches segment endpoints; deeper side owns a boundary.

All polynomial coefficients are hardcoded in the header. Sources in the parent
workspace: output/steinberger_adiabatic/fit_and_validation.json and
output/steinberger_M2_polynomial_factors/validation.json (results.A).
Original paper: Steinberger & Calderwood (2006), Figures1,4,9.
M2-A is the previously named A vector curve; its individual lmax is unassigned.
Changing Tref changes the computed profile; frozen constants are calibrated to
the current EBA background, not automatically refitted to a different background.

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

Compiles the actual production header with warnings-as-errors, compares all65
nodes with prior Python profile results, checks refvisc scaling, phase jumps and
invalid inputs. Viscosity_structures.c also passed MPI compiler syntax checking
(legacy non-prototype warnings remain). No full MPI convection run performed.

## Initial-temperature validation (70 Ma)

Run `python3 tests/validate_rheol7_initial_profile.py WORKSPACE_ROOT`. It uses the
current cfg and refstate_EBA.txt, the same erfc anomaly construction as
Lith_age.c/init_validate.py, and the compiled production C kernel. Surface
age is70 Ma with full lith_age_depth (127.42 km); the trench-specific half-depth
branch and initial slab anomalies are outside this representative radial column.
Bottom TBL thickness is approximately500 km at the1% amplitude height, with
Tbottom=3700 K. Geological start_age=249.9 Ma is distinct from plate thermal age.

The diagnostic evaluates65 nodes and128 radial Gauss points, checks reference
identity, cold half exponent, hot weakening and the Kelvin formula, and applies
the current depth-dependent viscosity limits. Actual production uses Gauss-point
viscosity; nodal limited values in the plot are illustrative. Outputs are in
output/rheol7_initial_validation. This is not a full3D initialization or MPI solve.

`cold_scale` is declared in the Pyre Visc inventory, transferred to the C state
by `pyCitcom_Visc_set_properties`, and read by the standalone C parser. It must be
finite and nonnegative; 0 disables cold strengthening, 1 restores full Az.
The hot-side exponent is independent of cold_scale. Rebuild the EBA solver including
the Python component before using the new cfg/lsf; jobs are not auto-submitted.

## Runtime binding and failure diagnostics

Both input paths validate cold_scale and print rheol=7, TDEPV, cold_scale and
hot_scale at startup. The property bridge also records cold_scale in its normal
parameter output. This fixes an omitted Pyre-to-C assignment that left the
malloc-allocated cold_scale field uninitialized.

The existing limits remain after temperature/composition/strain-rate/plasticity/
plate-boundary/channel corrections, before the two GP-to-node-to-GP filtering
cycles. The maximum is visc_max where element-center radius > 0.89641 and
5*visc_max below; the minimum is visc_min. These physical limits are not relocated.
For rheol7, finite nodal temperatures used for viscosity are now clipped to [0,1]
before Gauss-point interpolation, matching rheol3. E->T and Tref are unchanged.
With T_K = Ttop + ref_temperature*T_nd, these bounds are 300 and 3700 K for
the current experiments, matching their top and bottom boundary temperatures.
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
boundary types/values, effective clip bounds and individual failed-check flags
(1 means failed). Each node reports its original and clipped temperature plus
the Gauss interpolation weight. A nonfinite interpolated rheology temperature
is identified before calling the viscosity kernel. Finite out-of-range inputs
are clipped; invalid settings/nonfinite inputs return NaN and abort, without a
silent fallback to [0,1]. These diagnostics do not change the clipping policy.

Run `python3 tests/test_rheol7_runtime.py` for the production bridge and diagnostic
helpers with mocked Python property APIs and real two-process MPI abort tests.
This does not replace a build and end-to-end run with the cluster's Python2/Pyre.

Run `python3 tests/test_rheol7_temperature_clip.py` to check the production nodal
interpolation block with finite overshoots, unchanged in-range inputs, and
nonfinite temperatures, including the nodal temperatures from the step3 failure.
