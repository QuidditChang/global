#ifndef CITCOMS_PICES_H
#define CITCOMS_PICES_H
struct All_variables;
void pices_parameters(struct All_variables *);
void pices_validate(struct All_variables *);
void pices_initialize(struct All_variables *);
void pices_advance(struct All_variables *);
void pices_fail(struct All_variables *, const char *);
void tracer_temperature_weights(struct All_variables *, int, int, int *, double *);
void thermal_transport_at_gp(struct All_variables *, int, int, const double *,
                            double, double *, double *, double *, double *);
void pices_restore(struct All_variables *);
void pices_checkpoint_preflight(struct All_variables *, const char *);
void pices_checkpoint_publish(struct All_variables *, const char *, const char *);
void pices_checkpoint_check_state(struct All_variables *, const char *);
void pices_checkpoint_restore_velocity(struct All_variables *, const char *);
#endif
