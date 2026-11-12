#ifndef CITCOMS_PICES_VOLUME_CONTROL_H
#define CITCOMS_PICES_VOLUME_CONTROL_H
struct All_variables;
int pices_volume_control_enabled(struct All_variables *);
void pices_volume_control_register(struct All_variables *);
void pices_volume_control_run(struct All_variables *);
void pices_volume_control_check_velocity(struct All_variables *,const double *);
void pices_volume_control_check_bounds(struct All_variables *,const double [3],double,double,double);
/* Exact transport uses the same host search and MPI exchange as advection. */
void tracer_relocate_particles(struct All_variables *);
#endif
