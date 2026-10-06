#ifndef PICES_BENCHMARK_H
#define PICES_BENCHMARK_H
struct All_variables;
void p5_parameters(struct All_variables *);
void p5_initial(struct All_variables *);
void p5_particle_initial(struct All_variables *);
int p5_coupled(void);
int p5_velocity(struct All_variables *);
double p5_length_scale(void);
void p5_tick(int);
void p5_tock(int);
void p5_report(struct All_variables *);
#endif
