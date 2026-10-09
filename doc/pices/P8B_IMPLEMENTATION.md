# P8b: flow, heating and automatic steps

Implementation dates: 2026-10-08 to 2026-10-09. Target branch: cmbhf_EBA_PICES.
Fresh synthetic gates only; production approval remains pending P8c and the
P7 temperature-error investigation. No PG-to-PICES conversion is introduced.

## Accepted-state ordering

An accepted Stokes state supplies velocity and viscosity to particle transport
and the following thermal stage. The existing EBA source implementation is
shared by CBF and PICES; qvis_mode=2 applies the existing pressure-dependent
thermal conversion cap. Thermal substeps rebuild temperature-dependent thermal
coefficients, without rebuilding Stokes viscosity. PICES_VISCOUS records raw
and applied integrated heating; the original six-term heat ledger is unchanged.

After thermal advancement, reclassification and TA, elapsed_time is at arrival.
PICES updates the existing plate-file reader before assembling the next Stokes
system. Rheol7 then evaluates the new temperature and composition. The old
late plate refresh remains unchanged for PG. Prescribed PICES plates require
remove_rigid_rotation=off, matching the target cfg; subtracting a rigid rotation
would otherwise change the prescribed boundary velocities after solving.

## Automatic outer step

fixed_timestep=0 retains the shared advective/diffusive candidate and adds the
minimum of particle CFL=0.24, pices_max_timestep_Ma (default 0.1 Ma), the next
integer geological forcing age and the present day. Fixed positive timesteps
retain their existing behavior and CFL rejection. Heat subcycling retains its
independent stability check. New parameters are wired through both parsers.

Legacy elapsed_time and timestep are float. The selected step is rounded down,
with a representability check; age-knot resolution is
4*FLT_EPSILON*max(1,abs(start_age)) Ma. This is approximately 119 years at
249.9 Ma and is a clock precision bound, not a requested physical timestep.
The standalone driver and Pyre endTimestep share present-day termination.
The Python bindings are compile checked; local runtime tests use CitcomSFull,
not a complete Pythia installation.

## Test matrix and limits

Seven sequential 12-rank cases, one cap/rank, 9^3 nodes/cap, 64 particles/element
(393216 total), 25 flavors, primordial24, kC=0.8 and relaxed TA:

- coupled: rheol7, qvis2 (1e7 Pa, 0.085 rad), changing plate files, auto dt;
- uncapped: same initial state, qvis0;
- cap_probe: low synthetic yield threshold to exercise active limiting;
- static_plate: stationary surface control;
- present_day: 0.03 Ma start, automatic stop after one step;
- fixed_step: four equal physical steps, bypassing the new automatic limiter;
- multigrid: initial solve plus one accepted step, three levels on the small mesh;
  50 coarse sweeps instead of production vlowstep=2000 to avoid oversolving
  the 3^3-node coarse grid. The same 1e-6 inner/outer tolerances remain.

Four-step cases start at 2.13 Ma. Auto dt reaches the 2 Ma knot with
0.05, 0.05, 0.03 Ma, then 0.05 Ma. Synthetic constant tangential plate components
probe file interpolation and boundary timing; they are not Earth reconstructions.
The six-column reference state uses a synthetic monotone pressure profile.
Rayleigh/Di/Q0 remain small-test values. The actual 384-rank target mesh and
production Rayleigh/forcing still require the P8c pilot; small three-level MG
is not validation of production resolution or cost.

The audit independently checks output surface velocity against the input
interpolation at the accepted age, finite outputs, particle/primordial counts,
viscosity bounds, heat balance, raw/applied source, one TA per step, projection,
Stokes convergence and native-face CBF integrals. CG and three-level MG must
agree on their initial velocity field within relative L2=1e-4 (solver
consistency gate, not a production temperature-error criterion). Negative cases reject new
checkpoint combinations, rotation with prescribed plates and invalid limits.

P8b checkpoint/restart with plates, qvis, automatic dt or rheol7 remains blocked.
P8c must establish accepted-state continuation equivalence before lifting this
guard. P8a constant-viscosity composition restart stays supported.

## Multigrid correction discovered by the gate

The existing V-cycle retained a residual as the next cycle's original RHS,
while retaining the accumulated solution. That subtracts K*u twice and can
report convergence for an incorrect solution. Restore the original level RHS
at the start of each V-cycle (including intermediate levels). Also normalize
its returned residual with the same global equation count as solve_del2_u.
The correction is in the shared MG implementation; it does not replace MG with CG.
PICES independently evaluates ||K*u - rhs|| after MG and rejects a false pass.
The CG path is untouched by this correction.

## Prescribed-velocity assembly consistency

The coupled CG/MG comparison also exposed duplicate boundary lifting in the
legacy element-assembled CG path. The nodal matrix masks Dirichlet columns,
so assemble_forces must subtract their contribution there. The element matrix
retains those columns; initial_vel_residual already subtracts K*V, including
prescribed values. Select lifting from the existing assembly mode, rather than
applying it twice. Evaluate viscosity before force assembly so nodal lifting
and the subsequently built matrix use the same rheol7 state. These shared
linear-Stokes corrections are covered by a PG CG/MG comparison, independently
of PICES thermal advancement. Pseudo-free-surface and nonlinear rheology loops
are not extended by P8b.

For PICES with VISC_UPDATE=off, plate values still conform to U before solving;
the lack of a viscosity rebuild must not freeze a time-dependent boundary.

## Local validation result (2026-10-09)

PASS: seven P8b cases, rejection guards, frozen-viscosity plate updates and
an additional high-velocity case that activates CFL=0.23999999951257159.
PG and PICES CG/MG initial velocity comparisons both give relative L2
2.472685098304224e-6. The P8a regression preserves 480 decoded restart
outputs, 168 CBF outputs and identical checkpoint live-state payloads.
Standalone O2 compilation and all three changed Python C binding sources
compile successfully. Shell syntax and git diff checks pass.

HPC audit is pending. The runs worktree stores PICES_P8B_LOCAL_AUDIT.json
(with source/input hashes), PICES_P8B_LOCAL_AUDIT.md and PICES_P8B_RUNBOOK.md.
