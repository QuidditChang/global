# EBA CBF native output

Runtime now writes the original Q1 boundary data, without GRD interpolation or
NetCDF. The Appendix C Galerkin residual and GLL flux calculation are unchanged.

Set `monitoringFrequency_CBF=-1` to inherit monitoringFrequency, `0` to
disable, or a positive interval. Enable output_q_surf_CBF/output_q_botm_CBF.

Files in each rank's **expanded cfg datadir** (no global prefix):

```
<datadir>/q.botm.<rank>.<step>
<datadir>/q.surf.<rank>.<step>
```

For this benchmark: `CBF_benchmark/DATA/30/q.botm.30.13601` and
`CBF_benchmark/DATA/29/q.surf.29.13601` under the existing model directory.
The step is the actual solver step (13601, not the example typo 13061).

Canonical configuration keys: `output_q_surf_CBF`, `output_q_botm_CBF`,
`CBF_use_advection`, `monitoringFrequency_CBF`. Standalone uses `CBF_frequency`.
Old key names have been removed: update both C library and Python bindings and
use the migrated cfg files. Internal legacy arrays are q_surf/q_botm;
CBF arrays are q_surf_CBF/q_botm_CBF.

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
Restart/global.chkpt.<rank>.13600 headers for step 13600 and consistent finite timing across ranks,
and runs 384 ranks with the original 4x4x2 decomposition over 12 caps. It reads
checkpoint time rather than inventing a cfg time override. The user-reported
0.647491 is geological age in Ma; the actual checkpoint elapsed time is
approximately 0.000193788626348 nondimensional. Keep start_age=249.9 in cfg;
the full checkpoint restores both elapsed_time and its original start_age. Config intervals are
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

## Existing global.surf / global.botm files

These legacy ASCII products remain unchanged in layout. Each cap starts with
`cap_local nsf`; each following node line contains four nondimensional solver
values: (1) surface dynamic topography (or pseudo-free-surface topography when
enabled), or bottom dynamic topography; (2) legacy q_surf / q_botm; (3) V_theta;
(4) V_phi. Coordinates are in the normal coordinate output, with the same
surface-node ordering. Standard topography is currently a normal-stress-derived
diagnostic with topo_scaling=1, not a value already converted to metres. The values are printed with %.4e.

The legacy heat_flux() projects u_r*T - dT/dr to nodes and extrapolates the
boundary value. Its source explicitly notes missing conductivity, heat capacity
and unit conversion. It is not the CBF W/m2 result and not the variable-k dTdr
postprocessor. Renaming does not change that numerical definition.

## Restart clock repair (2026-09-20)

Full binary checkpoint time and original start_age are now preloaded before
lith_age_init and velocity boundary initialization. solution_cycles remains at
its initialization value so boundary arrays still allocate; the later full
checkpoint read restores the actual step. zero_elapsed_time/reset_startage do
not override the clock during full restart. At tracer checkpoint loading the
age-crossing cache is set to the current integer geological age, so restoring
the same checkpoint is not interpreted as a new age crossing.

Keep benchmark start_age=249.9; the authoritative start_age and elapsed time
come from the checkpoint. 0.647491 Ma is the current geological age, not the
initial age or nondimensional elapsed time. Check the first boundary messages:
they should now read the 0/1 Ma files, preceded by RESTART_CLOCK_PRELOAD.

Local binary-header tests and a 12-rank full checkpoint restart passed, including
an intentionally wrong cfg start_age and zero_elapsed_time=on. The synthetic
small-model checkpoint was relabelled to 13600; exactly one step and native CBF
output at 13601 were observed. This is not a rerun of the HPC production model.

The prior production data have complete files but a surface power of -54.12 TW,
strong heating/assimilation and the previously identified stale plate boundary.
The time fix does not by itself prove the CBF physical result correct. Rerun the
HPC case; output T / solver Tdot and velocity/source time-level consistency still
require physical budget analysis. No sign flip or value clipping was applied.

Native postprocessing is now provided in the scripts repository:
merge_CBF_native.py and cmbhf_EBA_Q0_30_rheol7_scold1.0.CBF.lsf.
It writes lossless float64 native-mesh NPZ archives; see CBF_POSTPROCESS.md.
