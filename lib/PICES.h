/*
 * PICES.h — Particle-In-Cell Energy Solver declarations
 *
 * Dong et al. (2024) Earth and Planetary Physics, doi:10.26464/epp2024021
 *
 * Functions that require gnomonic shape-function helpers (static in
 * Full_tracer_advection.c) are defined in that file.  Functions that
 * do not need those helpers are defined in PICES.c.  All external
 * declarations appear here so both translation units share a single
 * prototype source.
 */

#ifndef PICES_H
#define PICES_H

struct All_variables;   /* forward declaration */

/* ---- defined in Full_tracer_advection.c -------------------------------- */

/* One-time init: interpolate E->T → extraq[itracer_pices_temp] (Eq. 8 t=0) */
void PICES_init_particle_T(struct All_variables *E);                   //PICES

/* Eq. 13: accumulate particle T → nodal numerator/denominator arrays.
 * Caller must pre-zero, MPI-exchange, then normalise to build T_node. */
void PICES_particles_to_nodes(struct All_variables *E,                 //PICES
                               double **T_num,                         //PICES
                               double **W_den);                        //PICES

/* Eq. 14: interpolate a nodal field → particle extraq slot dest_idx. */
void PICES_nodes_to_particles(struct All_variables *E,                 //PICES
                               double **T_node,                        //PICES
                               int dest_extraq_idx);                   //PICES

/* Eq. 11: single-pass subgrid diffusion.
 * Caller must pre-zero arrays; function accumulates and returns.
 * MPI exchange + normalisation done by caller before PICES_update_particle_T. */
void PICES_subgrid_diffusion(struct All_variables *E,                  //PICES
                              double **dT_sub_node_num,                //PICES
                              double **dT_sub_node_den,                //PICES
                              double dt);                              //PICES

/* ---- defined in PICES.c ----------------------------------------------- */

/* Eq. 8+12: update particle T using diffusion-solve result and subgrid δT.
 * Uses extraq[pices_dtsub_idx] as internal workspace (overwritten). */
void PICES_update_particle_T(struct All_variables *E,                  //PICES
                              double **T_new_node,                     //PICES
                              double **dT_sub_node_num,                //PICES
                              double **dT_sub_node_den);               //PICES

#endif /* PICES_H */
