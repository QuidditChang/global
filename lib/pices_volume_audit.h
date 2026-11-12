#ifndef CITCOMS_PICES_VOLUME_AUDIT_H
#define CITCOMS_PICES_VOLUME_AUDIT_H
struct All_variables;
/* Observational hooks only. No particle, composition or checkpoint writes. */
void pices_volume_audit_initialize(struct All_variables *);
void pices_volume_audit_begin(struct All_variables *);
void pices_volume_audit_moved(struct All_variables *);
void pices_volume_audit_transition(struct All_variables *,int,int,int,int);
void pices_volume_audit_rebuilt(struct All_variables *);
void pices_volume_audit_commit(struct All_variables *);
#endif
