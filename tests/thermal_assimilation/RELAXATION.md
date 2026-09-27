# TA physical-time relaxation

For `lith_age_asml=1`, `lith_age_asml_tau_Ma` is positive and finite, in Myr
(default 1.0). The TA experiment explicitly selects 0.1 at the user's request,
not a calibrated value or a guarantee of removing the observed thermal layer.

After an accepted forward thermal step:

    alpha = -expm1(-dt_nd * scalet_Myr * w(z) / tau_Myr)
    T += alpha * (Ttarget - T)

The existing Tref/HSC target, plate-age cap, exponential spatial rate weight,
trench support, and initial full HSC anomaly are retained. The effective local
relaxation time is tau/w. For tau=0.1 Myr, the full-thickness column's last interior
node at 247.59 km has w=0.0046684, hence tau_eff≈21.42 Myr. This intentionally
changes the previous fixed-fraction, multiple-callback constraint.

`temperatures_conform_bcs` no longer performs TA interior assimilation. It
preserves radial boundary flags and restores prescribed top/bottom values for
Dirichlet boundaries. TA interior nodes have their TB/FB flags cleared: they
must participate in the physical heat equation, rather than have their
`DTdot` suppressed in `pg_solver`. Removing only repeated calls would leave
that stronger constraint in place. The current TA experiment is a full sphere.

The accepted-step operator uses geometric support directly, independently of
TB flags. It skips both radial boundaries, records only the actual stored
interior temperature increment in `assim_delta_T`, and runs before
`measure_temperature_assimilation`. Corrector callbacks cannot add relaxation
heat. Rejected attempts reset increments and skip the final TA operator when
`iredo` remains set. Legacy mode 0 retains its original callbacks/weights.
The unused backward timestep implementation explicitly rejects TA rather than
silently use an inappropriate accumulated time increment.

The global fields in `global_defs.h`, C parser, Pyre inventory and C extension
bridge all include the new parameter. Rebuild/install the executable, extension,
and Python components together. The online launcher checks for the new Pyre
parameter; that check does not establish binary/Python ABI consistency.
Existing running jobs do not acquire these changes. To compare experiments,
use a fresh run, or treat restart into this formulation as a changed experiment.

Local validation:
- Compiled production C target, taper, relaxation, initialization and boundary
  functions in a column harness; one accepted dt vs ten dt/10 operators agree
  within float temperature precision, and accumulated temperature increments
  match the actual change.
- Repeated boundary callbacks do not modify interior temperatures, while top
  and bottom fixed values remain exact; stale interior TB/FB flags are removed.
- Legacy target, zero-age handling, trench geometry and startup reference guards
  remain covered; source schedule check places the sole forward operator after
  retry handling and before heat-source measurement.
- Syntax checks for Lith_age.c, BC_util.c, Advection_diffusion.c. These do not
  replace a full linked HPC build or a production evolution comparison.
