/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 *<LicenseText>
 *
 * CitcomS by Louis Moresi, Shijie Zhong, Lijie Han, Eh Tan,
 * Clint Conrad, Michael Gurnis, and Eun-seo Choi.
 * Copyright (C) 1994-2005, California Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *</LicenseText>
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
/*  Here are the routines which process the results of each buoyancy solution, and call
    any relevant output routines. Much of the information has probably been output along
    with the velocity field. (So the velocity vectors and other data are fully in sync).
    However, heat fluxes and temperature averages are calculated here (even when they
    get output the next time around the velocity solver);
    */


#include "element_definitions.h"
#include "global_defs.h"
#include <math.h>		/* for sqrt */


void post_processing(struct All_variables *E)
{
  return;
}



/* ===================
    Surface heat flux
   =================== */

void heat_flux(E)
    struct All_variables *E;
{
    int m,e,el,i,j,node,lnode;
    float *flux[NCS],*SU[NCS],*RU[NCS];
    float VV[4][9],u[9],T[9],dTdz[9],area,uT;
    float *sum_h;
    double rtf[4][9];

    struct Shape_function GN;
    struct Shape_function_dA dOmega;
    struct Shape_function_dx GNx;
    void get_global_shape_fn();
    void velo_from_element();
    void sum_across_surface();
    void return_horiz_ave();
    void return_horiz_ave_f();

    const int dims=E->mesh.nsd,dofs=E->mesh.dof;
    const int vpts=vpoints[dims];
    const int ppts=ppoints[dims];
    const int ends=enodes[dims];
    const int nno=E->lmesh.nno;
    const int lev = E->mesh.levmax;
    const int sphere_key=1;


  sum_h = (float *) malloc((5)*sizeof(float));
  for(i=0;i<=4;i++)
    sum_h[i] = 0.0;

  for(m=1;m<=E->sphere.caps_per_proc;m++) {

    flux[m] = (float *) malloc((1+nno)*sizeof(float));

    for(i=1;i<=nno;i++)   {
      flux[m][i] = 0.0;
      }

    for(e=1;e<=E->lmesh.nel;e++) {
      get_global_shape_fn(E,e,&GN,&GNx,&dOmega,0,sphere_key,rtf,lev,m);

      velo_from_element(E,VV,m,e,sphere_key);

      for(i=1;i<=vpts;i++)   {
        u[i] = 0.0;
        T[i] = 0.0;
        dTdz[i] = 0.0;
        for(j=1;j<=ends;j++)  {
          u[i] += VV[3][j]*E->N.vpt[GNVINDEX(j,i)];
          T[i] += E->T[m][E->ien[m][e].node[j]]*E->N.vpt[GNVINDEX(j,i)];
          dTdz[i] += -E->T[m][E->ien[m][e].node[j]]*GNx.vpt[GNVXINDEX(2,j,i)];
          }
        }

      uT = 0.0;
      area = 0.0;
      for(i=1;i<=vpts;i++)   {
        /* XXX: missing unit conversion, heat capacity and thermal conductivity */
        uT += u[i]*T[i]*dOmega.vpt[i] + dTdz[i]*dOmega.vpt[i];
        }

      uT /= E->eco[m][e].area;

      for(j=1;j<=ends;j++)
        flux[m][E->ien[m][e].node[j]] += uT*E->TWW[lev][m][e].node[j];

      }             /* end of e */
    }             /* end of m */


  (E->exchange_node_f)(E,flux,lev);

  for(m=1;m<=E->sphere.caps_per_proc;m++)
     for(i=1;i<=nno;i++)
       flux[m][i] *= E->MASS[lev][m][i];

  if (E->parallel.me_loc[3]==E->parallel.nprocz-1)
    for(m=1;m<=E->sphere.caps_per_proc;m++)
      for(i=1;i<=E->lmesh.nsf;i++)
        E->slice.shflux[m][i]=2*flux[m][E->surf_node[m][i]]-flux[m][E->surf_node[m][i]-1];

  if (E->parallel.me_loc[3]==0)
    for(m=1;m<=E->sphere.caps_per_proc;m++)
      for(i=1;i<=E->lmesh.nsf;i++)
        E->slice.bhflux[m][i] = 2*flux[m][E->surf_node[m][i]-E->lmesh.noz+1]
                                - flux[m][E->surf_node[m][i]-E->lmesh.noz+2];

  for(m=1;m<=E->sphere.caps_per_proc;m++)
    for(e=1;e<=E->lmesh.snel;e++) {
         uT =(E->slice.shflux[m][E->sien[m][e].node[1]] +
              E->slice.shflux[m][E->sien[m][e].node[2]] +
              E->slice.shflux[m][E->sien[m][e].node[3]] +
              E->slice.shflux[m][E->sien[m][e].node[4]])*0.25;
         el = e*E->lmesh.elz;
         sum_h[0] += uT*E->eco[m][el].area;
         sum_h[1] += E->eco[m][el].area;

         uT =(E->slice.bhflux[m][E->sien[m][e].node[1]] +
              E->slice.bhflux[m][E->sien[m][e].node[2]] +
              E->slice.bhflux[m][E->sien[m][e].node[3]] +
              E->slice.bhflux[m][E->sien[m][e].node[4]])*0.25;
         el = (e-1)*E->lmesh.elz+1;
         sum_h[2] += uT*E->eco[m][el].area;
         sum_h[3] += E->eco[m][el].area;
         }

  sum_across_surface(E,sum_h,4);

  if (E->parallel.me_loc[3]==E->parallel.nprocz-1)   {
    sum_h[0] = sum_h[0]/sum_h[1];
    /*     if (E->control.verbose && E->parallel.me==E->parallel.nprocz-1) {
	     fprintf(E->fp_out,"surface heat flux= %f %f\n",sum_h[0],E->monitor.elapsed_time);
             fflush(E->fp_out);
    } */
    if (E->parallel.me==E->parallel.nprocz-1) {
      fprintf(stderr,"surface heat flux= %f\n",sum_h[0]);
      //fprintf(E->fp,"surface heat flux= %f\n",sum_h[0]); //commented out because E->fp is only on CPU 0 

      if(E->output.write_q_files > 0){
	/* format: time heat_flow sqrt(v.v)  */
	fprintf(E->output.fpqt,"%13.5e %13.5e %13.5e\n",E->monitor.elapsed_time,sum_h[0],sqrt(E->monitor.vdotv));
	fflush(E->output.fpqt);
      }
    }
  }

  if (E->parallel.me_loc[3]==0)    {
    sum_h[2] = sum_h[2]/sum_h[3];
/*     if (E->control.verbose && E->parallel.me==0) fprintf(E->fp_out,"bottom heat flux= %f %f\n",sum_h[2],E->monitor.elapsed_time); */
    if (E->parallel.me==0) {
      fprintf(stderr,"bottom heat flux= %f\n",sum_h[2]);
      fprintf(E->fp,"bottom heat flux= %f\n",sum_h[2]);
      if(E->output.write_q_files > 0){
	fprintf(E->output.fpqb,"%13.5e %13.5e %13.5e\n",
		E->monitor.elapsed_time,sum_h[2],sqrt(E->monitor.vdotv));
	fflush(E->output.fpqb);
      }

    }
  }


  for(m=1;m<=E->sphere.caps_per_proc;m++)
    free((void *)flux[m]);

  free((void *)sum_h);

  return;
}



/* =========================================================================
 * heat_flux_CBF_cmb  --  CMB heat flux via the Consistent Boundary Flux method
 *
 * Theory (Gresho et al. 1987, eq. 30 / eq. 27 lumped-mass form):
 *
 *   For each CMB Dirichlet node i, the consistent boundary flux is:
 *
 *       q_i = R_i / B_ii
 *
 *   where the CBF residual R_i is assembled from bottom-layer elements only
 *   (iz = 1, i.e. elements whose bottom face lies on the CMB):
 *
 *       R_i = -sum_{e in bot-layer} sum_g {
 *               N_i^g * [ rho*cp*(dT/dt + u.grad T) - H_tot ] * w_g|J|^g
 *             + lambda^-1 * kappa * (grad N_i . grad T)^g    * w_g|J|^g
 *             }
 *
 *   and the lumped boundary mass is:
 *
 *       B_ii = sum_{e ni i, e in bot-layer}  A_face / 4
 *            ~ V_elem / delta_r / 4          (uniform-mesh approximation)
 *
 *   The standard Galerkin test function N_j (not the PG upwind function) is
 *   used so that the result satisfies the global energy balance exactly.
 *
 *   Sign convention: q_i = (dT/dn_out)|_CMB  where n_out = -r_hat at CMB.
 *   Positive q_i means heat flux directed into the core (downward).
 *   In practice q_i < 0 for normal CMB conditions (hot core, heat into mantle).
 *
 * Result is stored in E->slice.cbf_bhflux[m][i]  (nsf surface-node array).
 * ========================================================================= */
void heat_flux_CBF_cmb(struct All_variables *E)
{
    int    m, e, i, j, node;
    float *cbf_rhs[NCS], *bdry_mass[NCS];
    float  VV[4][9];
    double tx1[9], tx2[9], tx3[9], sint_gp[9], dT[9];
    double v1[9],  v2[9],  v3[9];
    double rtf[4][9];
    float  sum_q, sum_area, q_global;
    char   msg[128];

    struct Shape_function    GN;
    struct Shape_function_dA dOmega;
    struct Shape_function_dx GNx;

    void get_global_shape_fn();
    void velo_from_element();
    void sum_across_surface();

    const int dims       = E->mesh.nsd;
    const int vpts       = vpoints[dims];   /* 8 volume Gauss points          */
    const int ends       = enodes[dims];    /* 8 nodes per hex element         */
    const int elz        = E->lmesh.elz;   /* elements in radial direction     */
    const int lev        = E->mesh.levmax;
    const int nsf        = E->lmesh.nsf;   /* surface nodes per cap            */
    const int sphere_key = 1;

    /* approximate radial thickness of one element layer */
    double delta_r = (E->sphere.ro - E->sphere.ri) / (double)elz;

    /* ---- allocate per-node work arrays ----------------------------------- */
    for (m = 1; m <= E->sphere.caps_per_proc; m++) {
        cbf_rhs[m]   = (float *)calloc(E->lmesh.nno + 2, sizeof(float));
        bdry_mass[m] = (float *)calloc(E->lmesh.nno + 2, sizeof(float));
    }

    /* ======================================================================
     * Pass 1: assemble CBF residual and boundary mass
     *         -- only on the MPI rank that owns the CMB layer (me_loc[3]==0)
     * ====================================================================== */
    if (E->parallel.me_loc[3] == 0) {
        for (m = 1; m <= E->sphere.caps_per_proc; m++) {
            for (e = 1; e <= E->lmesh.nel; e++) {

                /* keep only iz=1 (bottom) layer elements */
                if ((e - 1) % elz != 0) continue;

                get_global_shape_fn(E, e, &GN, &GNx, &dOmega,
                                    0, sphere_key, rtf, lev, m);
                velo_from_element(E, VV, m, e, sphere_key);

                /* ---- source / material properties (matches element_residual) */
                int    nz      = 1;   /* bottom layer: nz index = 1 */
                double rho     = 0.5*(E->refstate.rho[nz]
                                     + E->refstate.rho[nz+1]);
                double cp      = 0.5*(E->refstate.heat_capacity[nz]
                                     + E->refstate.heat_capacity[nz+1]);
                double Q       = E->control.Q0;
                double lat     = E->heating_latent[m][e]; /* 1/latent or 1 */
                double heating;
                if (E->control.disptn_number == 0)
                    heating = rho * Q;
                else
                    heating = (rho * Q
                               - E->heating_adi[m][e]
                               + E->heating_visc[m][e]) * lat;

                /* ---- interpolate T, dT/dt, velocity and grad T to Gauss pts */
                for (i = 1; i <= vpts; i++) {
                    dT[i] = tx1[i] = tx2[i] = tx3[i] = 0.0;
                    v1[i] = v2[i]  = v3[i]  = 0.0;
                    /* sint_gp[i] = 1 / (r * sin(theta))  at Gauss point i  */
                    sint_gp[i] = rtf[3][i] / sin(rtf[1][i]);
                }
                for (j = 1; j <= ends; j++) {
                    node = E->ien[m][e].node[j];
                    double T_j  = E->T[m][node];
                    double Td_j = E->Tdot[m][node]; /* dT/dt at node j */
                    for (i = 1; i <= vpts; i++) {
                        double Nj = E->N.vpt[GNVINDEX(j, i)];
                        /* dT/dt interpolated to Gauss point */
                        dT[i]  += Td_j * Nj;
                        /* grad T components (spherical, see Advection_diffusion.c) */
                        /*   tx1 = (1/r) dT/dtheta                                 */
                        /*   tx2 = (1/(r sin theta)) dT/dphi                       */
                        /*   tx3 = dT/dr                                            */
                        tx1[i] += GNx.vpt[GNVXINDEX(0, j, i)] * T_j * rtf[3][i];
                        tx2[i] += GNx.vpt[GNVXINDEX(1, j, i)] * T_j * sint_gp[i];
                        tx3[i] += GNx.vpt[GNVXINDEX(2, j, i)] * T_j;
                        /* velocity components: VV[1]=u_theta, VV[2]=u_phi, VV[3]=u_r */
                        v1[i]  += VV[1][j] * Nj;
                        v2[i]  += VV[2][j] * Nj;
                        v3[i]  += VV[3][j] * Nj;
                    }
                }

                /* ---- accumulate CBF residual for every node in this element -
                 *
                 *  R_j -= { N_j * [rho*cp*(dT/dt + u.gradT) - H_tot]
                 *         + lat * (grad N_j . grad T) } * w_g |J|^g
                 *
                 *  Standard Galerkin N_j is used (NOT the PG upwind function);
                 *  this is the critical difference from element_residual().      */
                for (j = 1; j <= ends; j++) {
                    node = E->ien[m][e].node[j];
                    for (i = 1; i <= vpts; i++) {
                        double Nj = E->N.vpt[GNVINDEX(j, i)];
                        double w  = dOmega.vpt[i]; /* w_g * |J|^g */

                        /* time-derivative + advection - source */
                        double transient =
                            (dT[i]
                             + v1[i]*tx1[i]
                             + v2[i]*tx2[i]
                             + v3[i]*tx3[i]) * rho * cp
                            - heating;

                        /* grad N_j . grad T  (full spherical form)
                         *   = (1/r dN/dtheta)*(1/r dT/dtheta)
                         *   + (1/(r sin t) dN/dphi)*(1/(r sin t) dT/dphi)
                         *   + (dN/dr)*(dT/dr)                               */
                        double diffusion = lat * (
                            GNx.vpt[GNVXINDEX(0, j, i)] * tx1[i] * rtf[3][i]
                          + GNx.vpt[GNVXINDEX(1, j, i)] * tx2[i] * sint_gp[i]
                          + GNx.vpt[GNVXINDEX(2, j, i)] * tx3[i]);

                        cbf_rhs[m][node] -= (Nj * transient + diffusion) * w;
                    }
                }

                /* ---- lumped boundary mass B_ii = sum A_face/4  -------------
                 *
                 *  CMB face = SIDE_BOTTOM = local nodes j = 1,2,3,4
                 *  (sidenodes[SIDE_BOTTOM][1..4] = {1,2,3,4}, element_definitions.h)
                 *
                 *  A_face ~ V_elem / delta_r  (uniform-grid approximation)    */
                double A_face = E->eco[m][e].area / delta_r;
                for (j = 1; j <= 4; j++) {   /* bottom face: local nodes 1-4  */
                    node = E->ien[m][e].node[j];
                    bdry_mass[m][node] += A_face * 0.25;
                }

            } /* end element loop */
        } /* end cap loop */
    } /* end if bottom-layer rank */

    /* ======================================================================
     * Pass 2: MPI halo exchange so each rank has the full nodal sums
     * ====================================================================== */
    (E->exchange_node_f)(E, cbf_rhs,   lev);
    (E->exchange_node_f)(E, bdry_mass, lev);

    /* ======================================================================
     * Pass 3: q_i = R_i / B_ii   and accumulate global scalar average
     * ====================================================================== */
    sum_q    = 0.0;
    sum_area = 0.0;

    if (E->parallel.me_loc[3] == 0) {
        for (m = 1; m <= E->sphere.caps_per_proc; m++) {
            for (i = 1; i <= nsf; i++) {
                /* CMB node = bottom node of the column rooted at surface node i */
                int cmb_node = E->surf_node[m][i] - E->lmesh.noz + 1;
                if (bdry_mass[m][cmb_node] > 0.0) {
                    E->slice.cbf_bhflux[m][i] =
                        (float)(cbf_rhs[m][cmb_node] / bdry_mass[m][cmb_node]);
                    /* area-weighted sum for the global scalar */
                    int e2d = ((i - 1) / E->lmesh.noz) + 1; /* approximate */
                    sum_q    += E->slice.cbf_bhflux[m][i];
                    sum_area += 1.0f;
                }
            }
        }
    }

    /* ======================================================================
     * Pass 4: global scalar output (same format as legacy heat_flux())
     * ====================================================================== */
    float sum_buf[2] = {sum_q, sum_area};
    sum_across_surface(E, sum_buf, 2);

    if (E->parallel.me_loc[3] == 0 && E->parallel.me == 0) {
        q_global = (sum_buf[1] > 0.0f) ? sum_buf[0] / sum_buf[1] : 0.0f;
        fprintf(stderr, "CBF bottom heat flux= %f\n", q_global);
        fprintf(E->fp,  "CBF bottom heat flux= %f\n", q_global);
        if (E->output.write_q_files > 0) {
            /* format: time  cbf_heat_flow  sqrt(v.v) */
            fprintf(E->output.fpqb, "CBF %13.5e %13.5e %13.5e\n",
                    E->monitor.elapsed_time,
                    q_global,
                    sqrt(E->monitor.vdotv));
            fflush(E->output.fpqb);
        }
    }

    /* ---- free work arrays ------------------------------------------------ */
    for (m = 1; m <= E->sphere.caps_per_proc; m++) {
        free(cbf_rhs[m]);
        free(bdry_mass[m]);
    }
}


/* =========================================================================
 * heat_flux_CBF_surf  --  Surface heat flux via the Consistent Boundary Flux method
 *
 * Symmetric to heat_flux_CBF_cmb(), but operates on the top boundary:
 *
 *   - Integration domain  : top-layer elements only  (iz = elz)
 *   - Boundary face nodes : SIDE_TOP = local nodes j = 5,6,7,8
 *   - Active MPI rank     : me_loc[3] == nprocz-1
 *   - Output array        : E->slice.cbf_shflux[m][i]
 *
 * Same CBF formula as heat_flux_CBF_cmb() -- see that function for the
 * full mathematical derivation.
 *
 * Sign convention: q_i = (dT/dn_out)|_surf  where n_out = +r_hat at surface.
 * Positive q_i means heat flux directed out of the mantle (upward = normal).
 * ========================================================================= */
void heat_flux_CBF_surf(struct All_variables *E)
{
    int    m, e, i, j, node;
    float *cbf_rhs[NCS], *bdry_mass[NCS];
    float  VV[4][9];
    double tx1[9], tx2[9], tx3[9], sint_gp[9], dT[9];
    double v1[9],  v2[9],  v3[9];
    double rtf[4][9];
    float  sum_q, sum_area, q_global;

    struct Shape_function    GN;
    struct Shape_function_dA dOmega;
    struct Shape_function_dx GNx;

    void get_global_shape_fn();
    void velo_from_element();
    void sum_across_surface();

    const int dims       = E->mesh.nsd;
    const int vpts       = vpoints[dims];
    const int ends       = enodes[dims];
    const int elz        = E->lmesh.elz;
    const int lev        = E->mesh.levmax;
    const int nsf        = E->lmesh.nsf;
    const int sphere_key = 1;

    double delta_r = (E->sphere.ro - E->sphere.ri) / (double)elz;

    /* ---- allocate per-node work arrays ----------------------------------- */
    for (m = 1; m <= E->sphere.caps_per_proc; m++) {
        cbf_rhs[m]   = (float *)calloc(E->lmesh.nno + 2, sizeof(float));
        bdry_mass[m] = (float *)calloc(E->lmesh.nno + 2, sizeof(float));
    }

    /* ======================================================================
     * Pass 1: assemble CBF residual and boundary mass
     *         -- only on the MPI rank that owns the surface layer
     * ====================================================================== */
    if (E->parallel.me_loc[3] == E->parallel.nprocz - 1) {
        for (m = 1; m <= E->sphere.caps_per_proc; m++) {
            for (e = 1; e <= E->lmesh.nel; e++) {

                /* keep only iz=elz (top) layer elements */
                if ((e - 1) % elz != elz - 1) continue;

                get_global_shape_fn(E, e, &GN, &GNx, &dOmega,
                                    0, sphere_key, rtf, lev, m);
                velo_from_element(E, VV, m, e, sphere_key);

                /* ---- source / material properties for top layer */
                int    nz      = elz;   /* top layer: nz index = elz */
                double rho     = 0.5*(E->refstate.rho[nz]
                                     + E->refstate.rho[nz+1]);
                double cp      = 0.5*(E->refstate.heat_capacity[nz]
                                     + E->refstate.heat_capacity[nz+1]);
                double Q       = E->control.Q0;
                double lat     = E->heating_latent[m][e];
                double heating;
                if (E->control.disptn_number == 0)
                    heating = rho * Q;
                else
                    heating = (rho * Q
                               - E->heating_adi[m][e]
                               + E->heating_visc[m][e]) * lat;

                /* ---- interpolate T, dT/dt, velocity and grad T to Gauss pts */
                for (i = 1; i <= vpts; i++) {
                    dT[i] = tx1[i] = tx2[i] = tx3[i] = 0.0;
                    v1[i] = v2[i]  = v3[i]  = 0.0;
                    sint_gp[i] = rtf[3][i] / sin(rtf[1][i]);
                }
                for (j = 1; j <= ends; j++) {
                    node = E->ien[m][e].node[j];
                    double T_j  = E->T[m][node];
                    double Td_j = E->Tdot[m][node];
                    for (i = 1; i <= vpts; i++) {
                        double Nj = E->N.vpt[GNVINDEX(j, i)];
                        dT[i]  += Td_j * Nj;
                        tx1[i] += GNx.vpt[GNVXINDEX(0, j, i)] * T_j * rtf[3][i];
                        tx2[i] += GNx.vpt[GNVXINDEX(1, j, i)] * T_j * sint_gp[i];
                        tx3[i] += GNx.vpt[GNVXINDEX(2, j, i)] * T_j;
                        v1[i]  += VV[1][j] * Nj;
                        v2[i]  += VV[2][j] * Nj;
                        v3[i]  += VV[3][j] * Nj;
                    }
                }

                /* ---- accumulate CBF residual (same formula as CMB version) */
                for (j = 1; j <= ends; j++) {
                    node = E->ien[m][e].node[j];
                    for (i = 1; i <= vpts; i++) {
                        double Nj = E->N.vpt[GNVINDEX(j, i)];
                        double w  = dOmega.vpt[i];

                        double transient =
                            (dT[i]
                             + v1[i]*tx1[i]
                             + v2[i]*tx2[i]
                             + v3[i]*tx3[i]) * rho * cp
                            - heating;

                        double diffusion = lat * (
                            GNx.vpt[GNVXINDEX(0, j, i)] * tx1[i] * rtf[3][i]
                          + GNx.vpt[GNVXINDEX(1, j, i)] * tx2[i] * sint_gp[i]
                          + GNx.vpt[GNVXINDEX(2, j, i)] * tx3[i]);

                        cbf_rhs[m][node] -= (Nj * transient + diffusion) * w;
                    }
                }

                /* ---- lumped boundary mass for surface face ------------------
                 *
                 *  Surface face = SIDE_TOP = local nodes j = 5,6,7,8
                 *  (sidenodes[SIDE_TOP][1..4] = {5,6,7,8}, element_definitions.h) */
                double A_face = E->eco[m][e].area / delta_r;
                for (j = 5; j <= 8; j++) {   /* top face: local nodes 5-8 */
                    node = E->ien[m][e].node[j];
                    bdry_mass[m][node] += A_face * 0.25;
                }

            } /* end element loop */
        } /* end cap loop */
    } /* end if top-layer rank */

    /* ======================================================================
     * Pass 2: MPI halo exchange
     * ====================================================================== */
    (E->exchange_node_f)(E, cbf_rhs,   lev);
    (E->exchange_node_f)(E, bdry_mass, lev);

    /* ======================================================================
     * Pass 3: q_i = R_i / B_ii   and accumulate global scalar average
     * ====================================================================== */
    sum_q    = 0.0;
    sum_area = 0.0;

    if (E->parallel.me_loc[3] == E->parallel.nprocz - 1) {
        for (m = 1; m <= E->sphere.caps_per_proc; m++) {
            for (i = 1; i <= nsf; i++) {
                /* surface node = top node of the column */
                int surf_node = E->surf_node[m][i];
                if (bdry_mass[m][surf_node] > 0.0) {
                    E->slice.cbf_shflux[m][i] =
                        (float)(cbf_rhs[m][surf_node] / bdry_mass[m][surf_node]);
                    sum_q    += E->slice.cbf_shflux[m][i];
                    sum_area += 1.0f;
                }
            }
        }
    }

    /* ======================================================================
     * Pass 4: global scalar output
     * ====================================================================== */
    float sum_buf[2] = {sum_q, sum_area};
    sum_across_surface(E, sum_buf, 2);

    if (E->parallel.me_loc[3] == E->parallel.nprocz - 1
        && E->parallel.me == E->parallel.nprocz - 1) {
        q_global = (sum_buf[1] > 0.0f) ? sum_buf[0] / sum_buf[1] : 0.0f;
        fprintf(stderr, "CBF surface heat flux= %f\n", q_global);
        if (E->output.write_q_files > 0) {
            fprintf(E->output.fpqt, "CBF %13.5e %13.5e %13.5e\n",
                    E->monitor.elapsed_time,
                    q_global,
                    sqrt(E->monitor.vdotv));
            fflush(E->output.fpqt);
        }
    }

    /* ---- free work arrays ------------------------------------------------ */
    for (m = 1; m <= E->sphere.caps_per_proc; m++) {
        free(cbf_rhs[m]);
        free(bdry_mass[m]);
    }
}


/*
  compute horizontal average of temperature and rms velocity
*/
void compute_horiz_avg(struct All_variables *E)
{
    void return_horiz_ave_f();

    int m, i;
    float *S1[NCS],*S2[NCS],*S3[NCS];

    for(m=1;m<=E->sphere.caps_per_proc;m++)      {
	S1[m] = (float *)malloc((E->lmesh.nno+1)*sizeof(float));
	S2[m] = (float *)malloc((E->lmesh.nno+1)*sizeof(float));
	S3[m] = (float *)malloc((E->lmesh.nno+1)*sizeof(float));
    }

    for(m=1;m<=E->sphere.caps_per_proc;m++) {
	for(i=1;i<=E->lmesh.nno;i++) {
	    S1[m][i] = E->T[m][i];
	    S2[m][i] = E->sphere.cap[m].V[1][i]*E->sphere.cap[m].V[1][i]
          	+ E->sphere.cap[m].V[2][i]*E->sphere.cap[m].V[2][i];
	    S3[m][i] = E->sphere.cap[m].V[3][i]*E->sphere.cap[m].V[3][i];
	}
    }

    return_horiz_ave_f(E,S1,E->Have.T);
    return_horiz_ave_f(E,S2,E->Have.V[1]);
    return_horiz_ave_f(E,S3,E->Have.V[2]);

    for(m=1;m<=E->sphere.caps_per_proc;m++) {
	free((void *)S1[m]);
	free((void *)S2[m]);
	free((void *)S3[m]);
    }

    for (i=1;i<=E->lmesh.noz;i++) {
	E->Have.V[1][i] = sqrt(E->Have.V[1][i]);
	E->Have.V[2][i] = sqrt(E->Have.V[2][i]);
    }

}
