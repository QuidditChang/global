#ifndef CITCOMS_TRACER_SEED_H
#define CITCOMS_TRACER_SEED_H
struct All_variables;
void tracer_seed_configure(struct All_variables *,const char *,int);
void tracer_seed_check_ic(struct All_variables *);
void tracer_seed_generate(struct All_variables *);
void tracer_seed_report(struct All_variables *);
int pices_seed_host_valid(struct All_variables *,int,int,double,double,double);
#endif
