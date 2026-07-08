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


/* ==========================================================================
   Heat flux via row-sum lumped Consistent Boundary Flux (CBF) method.

   Reference: Gresho et al. (1987), Dannberg et al. (2024, GJI).

   For each boundary-layer element, integrates the energy-equation weak form
   over the element volume to obtain the nodal residual b_i, and integrates
   the shape function over the boundary face to obtain the row-sum lumped CBF
   boundary mass A_i.
   The boundary heat flux at node i is:

       q_i = sign * Q_SCALE * b_i / A_i   (W/m^2)

   where b_i mirrors the temperature residual terms in Advection_diffusion.c:
         int_Omega [ kgp * latent * grad(N_i).grad(T)
                     + N_i*((dT/dt + u.grad(T))*rho*cp - heating) ] dOmega
         A_i = int_Gamma N_i dGamma

   This is a lumped CBF-style reconstruction, not a full B_ij boundary-mass
   solve. Stored sign conventions:
     top/surface: positive means heat flows from mantle to surface.
     bottom/CMB:  positive means heat flows from core into mantle.

   Results are stored in E->slice.shflux_CBF or E->slice.bhflux_CBF.
   ========================================================================== */

static void heat_flux_CBF_boundary(struct All_variables *E, int top,
                                   float *slice_flux[NCS],
                                   const char *label,
                                   const char *positive_direction)
{
    const int target_zproc = top ? E->parallel.nprocz - 1 : 0;
    const int face_id = top ? 1 : 0;
    const int side = top ? SIDE_TOP : SIDE_BOTTOM;
    const double sign = top ? -1.0 : 1.0;

    int    m, e, el, i, j, k, node, a, face_node;
    double b_node, A_node;
    double gradT[4];       /* grad(T) at a volume Gauss point, 1-indexed    */
    double Tdot_gp;        /* dT/dt interpolated to Gauss point              */
    double udotgradT;      /* u.grad(T) at Gauss point (advection, optional) */
    double T_gp;           /* T interpolated to Gauss point for k~_T         */
    double DT;
    double gradN_dot_gradT, mass_residual;
    double Q, rho, cp, heating;
    double kd_el, kC, kE, kT, T_dim, kgp;
    double Q_SCALE;        /* dimensional scale: k*DeltaT/R_earth  (W/m^2)  */

    float *b_cbf[NCS];     /* nodal CBF residual accumulator                 */
    float *A_cbf[NCS];     /* nodal face-area accumulator                    */
    float *sum_h;

    struct Shape_function     GN;
    struct Shape_function_dA  dOmega;
    struct Shape_function_dx  GNx;
    struct Shape_function1    GM;
    struct Shape_function1_dA dGamma;
    float  VV[4][9];
    double rtf[4][9];

    void get_global_shape_fn();
    void get_global_1d_shape_fn();
    void velo_from_element();
    void sum_across_surface();
    void parallel_process_termination();

    const int dims       = E->mesh.nsd;
    const int vpts       = vpoints[dims];        /* 8 volume Gauss points    */
    const int ends       = enodes[dims];         /* 8 nodes per element      */
    const int oned       = onedvpoints[dims];    /* 4 face Gauss points      */
    const int nno        = E->lmesh.nno;
    const int nsf        = E->lmesh.nsf;
    const int elz        = E->lmesh.elz;
    const int lev        = E->mesh.levmax;
    const int sphere_key = 1;

    if(E->parallel.me_loc[3] != target_zproc) return;

    /* Q_SCALE = k0 * DeltaT / R_earth  ~2.59e-3 W/m^2
       k0 = kappa * rho0 * Cp0  (baseline thermal conductivity)
       E->data.ref_temperature = DeltaT (K)
       E->data.radius_km * 1e3 = R_earth (m)

       Variable-k multipliers and inputdiffusivity are included once in kgp
       below, matching element_residual(). Q_SCALE only dimensionalizes the
       nondimensional boundary flux.                                      */
    Q_SCALE = (E->data.therm_diff * E->data.density * E->data.Cp)
              * E->data.ref_temperature / (E->data.radius_km * 1.0e3);

    /* allocate per-cap arrays */
    sum_h = (float *)malloc(5 * sizeof(float));
    for(i = 0; i <= 4; i++) sum_h[i] = 0.0;

    for(m = 1; m <= E->sphere.caps_per_proc; m++) {
        b_cbf[m] = (float *)malloc((nno + 2) * sizeof(float));
        A_cbf[m] = (float *)malloc((nno + 2) * sizeof(float));
        for(i = 1; i <= nno; i++) {
            b_cbf[m][i] = 0.0f;
            A_cbf[m][i] = 0.0f;
        }
    }

    /* ------------------------------------------------------------------
       Loop over the selected boundary-layer elements.
       ------------------------------------------------------------------ */
    for(m = 1; m <= E->sphere.caps_per_proc; m++) {
        for(e = 1; e <= E->lmesh.nel; e++) {

            if(top) {
                if(e % elz != 0) continue;          /* skip non-top elements */
            }
            else {
                if((e - 1) % elz != 0) continue;    /* skip non-bottom elements */
            }

            /* volume shape functions and Jacobian (sphere_key=1: Cartesian GNx) */
            get_global_shape_fn(E, e, &GN, &GNx, &dOmega, 0, sphere_key, rtf, lev, m);

            /* face shape functions and Jacobian. top=1 fills both bottom and
               top entries; top=0 fills the bottom entry only. */
            get_global_1d_shape_fn(E, e, &GM, &dGamma, top, m);

            if(E->output.cbf_use_advection)
                velo_from_element(E, VV, m, e, sphere_key);

            /* Keep this in sync with Advection_diffusion.c:element_residual().
               This is a local mirror, not a shared helper, to keep the CBF
               correction tightly scoped. */
            Q = E->control.Q0;
            if(E->control.tracer_enriched) {
                Q *= (1.0 - E->composition.comp_el[m][0][e]);
                Q += E->composition.comp_el[m][0][e] * E->control.Q0ER;
            }

            {
                int nz = ((e - 1) % E->lmesh.elz) + 1;
                rho = 0.5 * (E->refstate.rho[nz] + E->refstate.rho[nz+1]);
                cp = 0.5 * (E->refstate.heat_capacity[nz] + E->refstate.heat_capacity[nz+1]);
                kd_el = 0.5 * (E->refstate.thermal_conductivity[nz]
                              + E->refstate.thermal_conductivity[nz+1]);
            }

            /* k~_C is currently disabled in element_residual(); keep CBF in
               lockstep until composition conductivity is wired there. */
            kC = 1.0;
            kE = kd_el * kC;

            if(E->control.disptn_number == 0)
                heating = rho * Q;
            else
                heating = (rho * Q - E->heating_adi[m][e] + E->heating_visc[m][e])
                    * E->heating_latent[m][e];

            /* ---- volume integral: b_i ---- */
            for(a = 1; a <= ends; a++) {
                node    = E->ien[m][e].node[a];
                b_node  = 0.0;

                for(i = 1; i <= vpts; i++) {

                    /* Cartesian grad(T) at volume Gauss point i               */
                    gradT[1] = 0.0;  gradT[2] = 0.0;  gradT[3] = 0.0;
                    Tdot_gp  = 0.0;
                    T_gp     = 0.0;
                    for(j = 1; j <= ends; j++) {
                        int nj = E->ien[m][e].node[j];
                        double sfn = E->N.vpt[GNVINDEX(j,i)];
                        gradT[1] += GNx.vpt[GNVXINDEX(0,j,i)] * E->T[m][nj];
                        gradT[2] += GNx.vpt[GNVXINDEX(1,j,i)] * E->T[m][nj];
                        gradT[3] += GNx.vpt[GNVXINDEX(2,j,i)] * E->T[m][nj];
                        T_gp     += sfn * E->T[m][nj];

                        if(E->node[m][nj] & (TBX | TBY | TBZ))
                            DT = 0.0;
                        else
                            DT = E->Tdot[m][nj];
                        Tdot_gp += sfn * DT;
                    }

                    if(E->control.kT_exponent == 0.0) {
                        kT = 1.0;
                    }
                    else {
                        T_dim = (T_gp + E->control.surface_temp)
                                * E->data.ref_temperature;
                        if(T_dim <= 0.0) {
                            fprintf(stderr,
                                    "Invalid dimensional temperature in heat_flux_CBF: "
                                    "el=%d cap=%d gauss=%d T_gp=%e surface_temp=%e "
                                    "ref_temperature=%e T_dim=%e\n",
                                    e, m, i, T_gp, E->control.surface_temp,
                                    E->data.ref_temperature, T_dim);
                            parallel_process_termination();
                        }
                        kT = pow(300.0 / T_dim, (double)E->control.kT_exponent);
                    }
                    kgp = E->control.inputdiff * kE * kT;

                    gradN_dot_gradT =
                        GNx.vpt[GNVXINDEX(0,a,i)] * gradT[1]
                      + GNx.vpt[GNVXINDEX(1,a,i)] * gradT[2]
                      + GNx.vpt[GNVXINDEX(2,a,i)] * gradT[3];

                    b_node += kgp * E->heating_latent[m][e]
                        * gradN_dot_gradT * dOmega.vpt[i];

                    mass_residual = Tdot_gp;
                    if(E->output.cbf_use_advection) {
                        udotgradT = 0.0;
                        for(k = 1; k <= dims; k++) {
                            double u_k = 0.0;
                            for(j = 1; j <= ends; j++)
                                u_k += VV[k][j] * E->N.vpt[GNVINDEX(j,i)];
                            udotgradT += u_k * gradT[k];
                        }
                        mass_residual += udotgradT;
                    }

                    b_node += E->N.vpt[GNVINDEX(a,i)]
                        * (mass_residual * rho * cp - heating) * dOmega.vpt[i];

                }  /* end volume Gauss loop */

                b_cbf[m][node] += (float)b_node;

            }  /* end volume nodes loop */

            /* ---- face integral: A_i for boundary-face nodes ----
               Uses E->M.vpt[GMVINDEX(a,k)]  (face shape fn, node a, Gauss pt k)
               and  dGamma.vpt[GMVGAMMA(face_id,k)] (face Jacobian at Gauss pt k)
               Gauss weights are 1.0 for 2-pt quadrature on [-1,1], so omitted. */
            for(a = 1; a <= oned; a++) {
                face_node = sidenodes[side][a];
                node = E->ien[m][e].node[face_node];
                A_node = 0.0;
                for(k = 1; k <= oned; k++)
                    A_node += E->M.vpt[GMVINDEX(a,k)] * dGamma.vpt[GMVGAMMA(face_id,k)];
                A_cbf[m][node] += (float)A_node;
            }

        }  /* end element loop */
    }  /* end cap loop */

    /* ---- MPI exchange: sum shared boundary nodes ---- */
    (E->exchange_node_f)(E, b_cbf, lev);
    (E->exchange_node_f)(E, A_cbf, lev);

    /* ---- compute q = sign * Q_SCALE * b / A and store in the target slice ---- */
    for(m = 1; m <= E->sphere.caps_per_proc; m++) {
        for(i = 1; i <= nsf; i++) {
            node = top ? E->surf_node[m][i] : E->surf_node[m][i] - E->lmesh.noz + 1;
            if(A_cbf[m][node] > 0.0f)
                slice_flux[m][i] = (float)(sign * Q_SCALE * b_cbf[m][node] / A_cbf[m][node]);
            else
                slice_flux[m][i] = 0.0f;
        }
    }

    /* ---- global average for monitoring ---- */
    for(m = 1; m <= E->sphere.caps_per_proc; m++) {
        for(e = 1; e <= E->lmesh.snel; e++) {
            float qavg = ( slice_flux[m][E->sien[m][e].node[1]]
                         + slice_flux[m][E->sien[m][e].node[2]]
                         + slice_flux[m][E->sien[m][e].node[3]]
                         + slice_flux[m][E->sien[m][e].node[4]] ) * 0.25f;
            el = top ? e * elz : (e - 1) * elz + 1;
            sum_h[0] += qavg * E->eco[m][el].area;
            sum_h[1] += E->eco[m][el].area;
        }
    }

    sum_across_surface(E, sum_h, 2);

    if((top && E->parallel.me == E->parallel.nprocz - 1) ||
       (!top && E->parallel.me == 0)) {
        if(sum_h[1] > 0.0f) {
            float q_ave = sum_h[0] / sum_h[1];
            fprintf(stderr,
                    "%s heat flux (lumped CBF, positive %s, advection %s) = %g W/m^2\n",
                    label, positive_direction,
                    E->output.cbf_use_advection ? "on" : "off", q_ave);
            fprintf(E->fp,
                    "%s heat flux (lumped CBF, positive %s, advection %s) = %g W/m^2\n",
                    label, positive_direction,
                    E->output.cbf_use_advection ? "on" : "off", q_ave);
        }
    }

    for(m = 1; m <= E->sphere.caps_per_proc; m++) {
        free((void *)b_cbf[m]);
        free((void *)A_cbf[m]);
    }
    free((void *)sum_h);

    return;
}

void heat_flux_CBF(struct All_variables *E)
{
    if(E->output.cbf_output_shflux)
        heat_flux_CBF_boundary(E, 1, E->slice.shflux_CBF,
                               "surface", "mantle->surface");

    if(E->output.cbf_output_bhflux)
        heat_flux_CBF_boundary(E, 0, E->slice.bhflux_CBF,
                               "bottom/CMB", "core->mantle");

    return;
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
