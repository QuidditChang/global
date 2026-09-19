# EBA CBF native output

Runtime now writes the original Q1 boundary data, without GRD interpolation or
NetCDF. The Appendix C Galerkin residual and GLL flux calculation are unchanged.

Set `monitoringFrequency_cmbhf_CBF=-1` to inherit monitoringFrequency, `0` to
disable, or a positive interval. Enable cbf_output_shflux/cbf_output_bhflux.

Files relative to the solver working directory:

```
PostProc/HF_CBF/cmbhf_CBF_<step>.rank<rank:06d>.dat
PostProc/HF_CBF/eshf_CBF_<step>.rank<rank:06d>.dat
```

Only ranks touching that physical boundary write its file. Files are text with
17 significant digits, preserving double values through decimal round-trip.
Header metadata records time, scales, sign, global native power and area.

- `N cap_global node_local x_nd y_nd z_nd theta_rad phi_rad r_nd q_W_m2 rhs_nd mass_nd`
- `F cap_global element_local node1 node2 node3 node4 weight1_nd weight2_nd weight3_nd weight4_nd`

Node IDs are scoped by **rank and cap**. Shared MPI/cap nodes are duplicated;
do not sum assembled nodal masses across all files. Faces are owned uniquely:
sum `weight_nd * q_W_m2 * length_scale_m^2` over each face corner to reproduce
the native power. Face node IDs refer to N records in the same rank/cap file.
Cartesian and spherical coordinates are original solver coordinates.

RHS is the assembled nondimensional **outward** energy residual. With
`scale=k0*deltaT/length`, top q=scale*rhs/mass and bottom q=-scale*rhs/mass.
These are total residuals, not a separate decomposition of every source term.
The q sign is mantle-to-surface at top and core-to-mantle at bottom.

Writes use temporary files followed by rename and collective error checks.
A `CBF_NATIVE_COMPLETE step=... max_wall_seconds=...` solver log entry is
printed only after both enabled boundaries finish; partial files from a failed
job do not establish a complete step. GRD generation is a separate future
postprocessing action, not part of this runtime output.

State caveat: output T with solver Tdot excludes assimilation/filter increments
from Tdot. The initial_or_restart_state header identifies launch-state outputs.

## Benchmark 13600 -> 13601

Use runs/cmbhf_EBA_Q0_30_rheol7_scold1.0.CBF.benchmark.cfg and .lsf.
Rebuild/reinstall the **cmbhf_EBA C library and Python Controller.py** into the
HPC build before submitting. The fixed Controller early exit permits exactly
one remaining step. No NetCDF flags/library are required by CBF anymore.

The LSF script keeps the existing model directory as cwd, checks all 384 binary
Restart/global.chkpt.<rank>.13600 headers for step 13600 and time_nd 0.647491,
and runs 384 ranks with the original 4x4x2 decomposition over 12 caps. It reads
checkpoint time rather than inventing a cfg time override. Config intervals are
13601, so no CBF/regular/profile/checkpoint output is requested at 13600;
13601 is the terminal step. Ordinary output uses CBF_benchmark/DATA/%RANK.
Each boundary produces 192 rank files. Restart inputs are not modified.

If the prepared Restart contains ASCII temperature files rather than binary
checkpoints, this full-state restart is not applicable: it needs the original
binary checkpoints including temperature derivative, velocity and tracers.

## Local verification

```
python3 tests/cbf/test_controller_one_step.py
python3 tests/cbf/test_cbf_kernel.py
python3 tests/cbf/build_validation.py --build-dir /tmp/cbf-native-build
python3 tests/cbf/verify_native_outputs.py RUN_DIR --step 1 --boundary-ranks 12
```

For the HPC benchmark use `--step 13601 --boundary-ranks 192`. The verifier
checks node equations, finite values, positive masses, expected rank-file count
and reconstructed face power. It requires Python 3, independently of Pyre's
Python 2 runtime. Historical GRD validation is in VALIDATION_REPORT.md; those
GRD-only tests were retired with the runtime GRD writer.
