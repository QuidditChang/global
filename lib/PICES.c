/*
 * PICES.c — Particle-In-Cell Energy Solver (non-shape-function functions)
 *
 * Dong et al. (2024) Earth and Planetary Physics, doi:10.26464/epp2024021
 *
 * Functions that require the gnomonic shape-function helpers
 * (spherical_to_uv, get_2dshape, get_radial_shape — all static in
 * Full_tracer_advection.c) are defined there.  Only functions that do
 * NOT need those helpers live here.
 */

#include <math.h>
#include <stdlib.h>
#include "global_defs.h"
#include "PICES.h"


/****** PICES_update_particle_T ********************************************
 * Update each particle's temperature after the Eulerian diffusion solve.
 * Implements paper Eq. 8 and 12:
 *
 *   δT_sub_norm[node]  = dT_sub_node_num / dT_sub_node_den       (Eq. 11 nodal)
 *   δT_remain[node]    = δT_diff_node[node] - δT_sub_norm[node]  (Eq. 12)
 *   T_p^{n+1} = T_p^n + δT_sub_p + δT_remain_p                  (Eq. 8)
 *
 * Execution order:
 *   1. Build δT_sub_norm nodal field (all caps).
 *   2. PICES_nodes_to_particles → extraq[pDS]: δT_sub per particle.
 *   3. Per-cap: malloc dT_sub_p, save extraq[pDS], overwrite nodal
 *      array with dT_remain = δT_diff - δT_sub_norm.
 *   4. PICES_nodes_to_particles → extraq[pDS]: δT_remain per particle.
 *   5. Per-cap: T_p += dT_sub_p + extraq[pDS];  free dT_sub_p.
 *   6. Free temporary nodal scratch.
 *
 * extraq[pDS] is workspace; its content after this call is undefined.
 *
 * Inputs (all already MPI-exchanged by caller):
 *   dT_diff_node     — nodal temperature change caused by diffusion
 *   dT_sub_node_num  — Σ(δT_sub * w) from PICES_subgrid_diffusion
 *   dT_sub_node_den  — Σ(w)           from PICES_subgrid_diffusion
 **************************************************************************/

void PICES_update_particle_T(struct All_variables *E,    //PICES
                              double **dT_diff_node,     //PICES
                              double **dT_sub_node_num,  //PICES
                              double **dT_sub_node_den)  //PICES
{
    int j, n, kk;                                        //PICES
    int pT  = E->trace.itracer_pices_temp;               //PICES: particle T slot
    int pDS = E->trace.pices_dtsub_idx;                  //PICES: δT workspace slot
    int caps = E->sphere.caps_per_proc;                  //PICES
    int nno  = E->lmesh.nno;                             //PICES

    /* temporary nodal scratch (all caps) for δT_sub_norm, then dT_remain */
    double **dT_scratch;                                 //PICES
    dT_scratch = (double **)malloc((caps + 1) * sizeof(double *));  //PICES
    for (j = 1; j <= caps; j++)                                     //PICES
        dT_scratch[j] = (double *)malloc((nno + 1) * sizeof(double)); //PICES

    /* per-cap pointer array for δT_sub_p buffers (buffers alloc'd below) */
    double **dT_sub_p;                                   //PICES
    dT_sub_p = (double **)malloc((caps + 1) * sizeof(double *));  //PICES

    /* --- Step 1: normalise δT_sub to nodal field (all caps) ------------ */
    for (j = 1; j <= caps; j++) {                                   //PICES
        for (n = 1; n <= nno; n++) {                                //PICES
            if (dT_sub_node_den[j][n] > 0.0)                       //PICES
                dT_scratch[j][n] = dT_sub_node_num[j][n]           //PICES
                                 / dT_sub_node_den[j][n];          //PICES
            else                                                    //PICES
                dT_scratch[j][n] = 0.0;  /* no tracer coverage */  //PICES
        }                                                           //PICES
    }                                                               //PICES

    /* --- Step 2: interpolate δT_sub_norm → particles (all caps, → pDS) */
    PICES_nodes_to_particles(E, dT_scratch, pDS);                  //PICES

    /* --- Step 3 (per-cap): save δT_sub_p; overwrite nodal with dT_remain */
    for (j = 1; j <= caps; j++) {                                   //PICES
        int ntracers = E->trace.ntracers[j];                        //PICES

        dT_sub_p[j] = (double *)malloc((ntracers + 1) * sizeof(double)); //PICES

        for (kk = 1; kk <= ntracers; kk++)                         //PICES
            dT_sub_p[j][kk] = E->trace.extraq[j][pDS][kk];        //PICES: save δT_sub_p

        /* Eq. 12: δT_remain = δT_diff - δT_sub_norm (reuse dT_scratch) */
        for (n = 1; n <= nno; n++)                                  //PICES
            dT_scratch[j][n] = dT_diff_node[j][n] - dT_scratch[j][n]; //PICES
    }                                                               //PICES

    /* --- Step 4: interpolate dT_remain → particles (all caps, → pDS) -- */
    PICES_nodes_to_particles(E, dT_scratch, pDS);                  //PICES

    /* --- Step 5 (per-cap): T_p += δT_sub_p + δT_remain_p; free buffer - */
    for (j = 1; j <= caps; j++) {                                   //PICES
        int ntracers = E->trace.ntracers[j];                        //PICES
        for (kk = 1; kk <= ntracers; kk++)                         //PICES
            E->trace.extraq[j][pT][kk] +=                          //PICES
                dT_sub_p[j][kk]                                    //PICES: δT_sub_p (Eq. 11)
              + E->trace.extraq[j][pDS][kk];                       //PICES: δT_remain_p (Eq. 12)
        free(dT_sub_p[j]);                                          //PICES
    }                                                               //PICES

    /* --- Step 6: free temporary nodal scratch -------------------------- */
    free(dT_sub_p);                                                 //PICES
    for (j = 1; j <= caps; j++) free(dT_scratch[j]);               //PICES
    free(dT_scratch);                                               //PICES

    return;                                                         //PICES
}  /* end PICES_update_particle_T */                                //PICES
