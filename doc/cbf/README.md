# EBA CBF heat flux

Read [method/audit/plan](APPENDIX_C_AUDIT_PLAN.md) and [actual validation](VALIDATION_REPORT.md).

## Runtime

In the Pyre cfg:

```ini
[CitcomS.controller]
monitoringFrequency = 50
monitoringFrequency_cmbhf_CBF = -1

[CitcomS.solver]
cbf_output_shflux = on
cbf_output_bhflux = on
cbf_use_advection = on
```

`-1` inherits the ordinary monitoring interval; `0` disables CBF; a positive
integer is an independent interval. Disabling advection is rejected when CBF
runs. The main `runs-cmbhf_EBA/cmbhf_EBA.cfg` now inherits its 50-step interval.

Output is relative to the solver's launch working directory:

```
PostProc/HF_CBF/cmbhf_CBF_50.grd
PostProc/HF_CBF/eshf_CBF_50.grd
```

Both are W/m². Bottom positive means core to mantle; top positive means mantle
to surface. Each file is classic NetCDF, `z[y,x]`, 0.5-degree gridline nodes,
longitude 0..360 and latitude -90..90. Step 0 is explicitly marked as an initial
diagnostic. Read `state` and assimilation/filter metadata before interpreting
transient energy balance. Native integrated heat is stored as a global
attribute; point sampling to GRD is not conservative remapping.

## Build

NetCDF C is required for enabled CBF output. Supply the flags reported by
`nc-config --cflags` and `nc-config --libs` when configuring the usual build,
while preserving existing MPI/HDF5/Pyre options. Regenerate `configure` from
`configure.ac` with the project's autotools workflow so the NetCDF detection
runs. The configure check defines `USE_CBF_NETCDF` on success. Only having an
old configure script or only installing Python netCDF packages is insufficient.

The pre-existing startup pointer errors in `Parsing.c` are fixed. Validation
builds now compile repository sources directly, without a temporary parser patch.

## Tests

From the solver worktree:

```sh
python3 tests/cbf/test_parser_control.py
python3 tests/cbf/test_cbf_kernel.py
python3 tests/phase_energy/test_phase_energy_geometry.py
MPIRUN=/path/to/matching/mpirun python3 tests/cbf/test_grid_output.py
python3 tests/cbf/build_validation.py --build-dir /tmp/cbf-validation
```

The grid test requires NumPy/SciPy only for independent test reading, not for
runtime output. Use an MPI launcher matching `mpicc`; local sockets must be
allowed. The standalone driver uses `cmbhf_CBF_freq=-1` to inherit
`storage_spacing`, or a positive explicit interval.
