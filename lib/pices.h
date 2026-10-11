#ifndef CITCOMS_PICES_H
#define CITCOMS_PICES_H
struct All_variables;
enum { PD_INTERNAL, PD_ADI, PD_VISC, PD_PHASE, PD_RAW, PD_TA, PD_NQ };
enum { PD_DT, PD_TIME, PD_BEFORE, PD_AFTER, PD_REMAP, PD_REBASE,
 PD_HEAT, PD_BOUNDARY, PD_TASTORE, PD_RECON, PD_RESIDUAL, PD_HEAT_RESIDUAL,
 PD_PARTICLES, PD_TPMIN, PD_TPMAX, PD_TPMEAN, PD_GRIDMIN, PD_GRIDMAX,
 PD_ENERGY_SCALE, PD_POWER_DENSITY_SCALE, PD_QINTERNAL, PD_QADI, PD_QVISC, PD_QPHASE, PD_QRAW, PD_NVALUES };
struct Pices_diagnostics {
 int valid, step;
 double values[PD_NVALUES], *q, *trial, *volume;
};
const struct Pices_diagnostics *pices_diagnostics(struct All_variables *);
extern const char *pices_diagnostic_names[PD_NVALUES];
double pices_profile_source(struct All_variables *,int,int);

void pices_consistent_project(struct All_variables *, const double *, double *, int);
void pices_parameters(struct All_variables *);
void pices_validate(struct All_variables *);
void pices_initialize(struct All_variables *);
void pices_advance(struct All_variables *);
void pices_limit_timestep(struct All_variables *);
int pices_time_finished(struct All_variables *);
void pices_fail(struct All_variables *, const char *);
void tracer_temperature_weights(struct All_variables *, int, int, int *, double *);
void thermal_transport_at_gp(struct All_variables *, int, int, const double *,
                            double, double *, double *, double *, double *);
void pices_restore(struct All_variables *);
int pices_checkpoint_coupled(struct All_variables *);
void pices_checkpoint_preflight(struct All_variables *, const char *);
void pices_checkpoint_publish(struct All_variables *, const char *, const char *);
void pices_checkpoint_check_state(struct All_variables *, const char *);
void pices_checkpoint_restore_state(struct All_variables *, const char *);
#endif
