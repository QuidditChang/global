# PICES route B — P1 implementation

Status (2026-09-29): local standalone MPI validation passed; HPC/Pyre installation and HPC run remain to be validated. Solver branch: `cmbhf_EBA_PICES`; runs branch: `cmbhf_EBA`. P0 contract remains in `P0_CONTRACT.md`.

## Scope and sequence

`energy_solver=pg` remains the default. `energy_solver=pices` selects the new kernel. P1 supports full sphere, one cap per rank, fixed positive timestep, fixed radial temperature boundaries, constant positive rho/Cp/conductivity, and no physical heat sources or phase changes. Restart, strict ALA, composition feedback, temperature filtering, rheo.dat temperature edits, legacy CBF and checkpoint publication are excluded. Unsupported configurations terminate with `PICES_ERROR` rather than silently using PG.

Tp is a persistent double tracer extra quantity, registered after optional flavor data before allocation. Initialization interpolates the initialized nodal temperature to particles. Existing predictor/corrector, relocation and MPI migration carry Tp. The thermal step advances particles once; the subsequent legacy tracer call is suppressed only for PICES.

Both particle-to-grid P and grid-to-particle Q use the existing gnomonic wedge/radial interpolation geometry. P exchanges weighted numerators and denominators before normalization. Free nodes without coverage abort; absolute grid projection imposes fixed boundary values. Increment projection does not overwrite boundary increments. Empty elements are counted in the log.

After transport, project Tp to g once. With particle positions frozen, each explicit lumped-FE diffusion substep updates g, then applies:

```
sub = (Q(g_new) - Tp) * (-expm1(-ds/tau))
Tp += sub + Q(g_new - g_old - P(sub))
tau = Le**2 / kappa
```

The original particle subgrid increment is retained. There is no repeated projection of Tp back to g between thermal substeps. Le is the minimum of the three directional RMS Cartesian edge lengths and is evaluated at initialization. Physical transport coefficients come from the same exported kernel used by PG; PG arithmetic is unchanged.

The heat-step limit uses 0.8 divided by the largest free-node accumulated element absolute row sum / lumped capacity. This is a conservative upper bound on the assembled absolute row sum, because cancellation between elements can reduce the latter. Particle CFL uses global maximum nodal speed and minimum Cartesian edge length, and must be <= 0.25 before movement. The maximum heat-substep count is also checked before movement.

## Output semantics

`PICES_INIT` records the Tp slot and diffusion bound. `PICES_STEP` records accepted time, particle count, substeps, CFL, temperature range, free-node grid/particle mismatch, maximum subgrid increment, and separate remapping/heat energy changes. Energy uses local element capacity integrals with one global reduction; shared nodes are not double-counted. These are grid energy ledgers, not a particle-mass conservation proof or CBF certification.

Legacy Tdot is zero and legacy PG thermal-budget/heat-flux output is suppressed. No P1 restart checkpoint is written even if a checkpoint frequency is configured. The tracer output has theta, phi, radius, Tp for the supplied zero-flavor smoke case; Tp is nondimensional, with Kelvin = 300 + 3400*Tp. Velocity output temperatures are Kelvin.

The early log stream is initialized to stderr so tracer input can safely log before rank-specific files open. Pyre inventory and C property bindings expose the same options, but this phase is exercised through standalone CitcomSFull; Pyre runtime is not locally certified.

## Local evidence

See `P1_LOCAL_RESULT.json`. The actual 12-rank, two-step Stokes/PICES smoke has 98,304 particles at all three outputs, finite 300–3700 K grid values and maximum CFL 0.2024548855. The prescribed-field driver links production objects and tests constant preservation, pure rotation with exact persistent-Tp multiset preservation, six-substep diffusion against an independently assembled NumPy reference, flavor-slot registration, and rejection of empty coverage, sources, excessive CFL, excessive substeps and restart. Diffusion particle error is 4.44e-16; the raw subgrid survival test detects accidental cancellation.

The matrix checks include element symmetry, row-sum zero, positive semidefiniteness, and the assembled free-node spectral stability bound. They do not establish spatial convergence order. The two-step PG regression matches all 36 decoded velocity/temperature outputs from the P0 local baseline exactly. Existing CBF kernel (6), P0 auditor (7) and temperature-audit (3) tests also pass.

Reproduce locally (working MPI and NumPy required; use fresh output directories):

```bash
MPICC=/usr/local/bin/mpicc python3 tests/cbf/build_validation.py --build-dir /tmp/pices-p1-build-new
python3 tests/pices/build_manufactured.py /tmp/pices-p1-build-new
MPIEXEC=/usr/local/bin/mpiexec python3 tests/pices/run_local.py /tmp/pices-p1-build-new ../runs-cmbhf_EBA/cmbhf_EBA_PICES_P1.cfg --output /tmp/pices-p1-tests-new
```

The test build does not replace the HPC Autotools/Python 2 installation. After the HPC run, use `tests/pices/verify_p1.py JOB_DIRECTORY --summary result.json` without `--local` to check provenance as well as numerical output. Systematic cap-seam transport, resolution convergence, restart equivalence, variable coefficients and physical heat-source extensions remain later work. Do not advance to P2 before the downloaded P1 HPC evidence is audited.
