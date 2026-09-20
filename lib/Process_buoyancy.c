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
#include "material_properties.h"
#include "advection_diffusion.h"
#include "CBF_face_geometry.h"
#include "CBF_native_output.h"
#include <float.h>
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
        E->slice.q_surf[m][i]=2*flux[m][E->surf_node[m][i]]-flux[m][E->surf_node[m][i]-1];

  if (E->parallel.me_loc[3]==0)
    for(m=1;m<=E->sphere.caps_per_proc;m++)
      for(i=1;i<=E->lmesh.nsf;i++)
        E->slice.q_botm[m][i] = 2*flux[m][E->surf_node[m][i]-E->lmesh.noz+1]
                                - flux[m][E->surf_node[m][i]-E->lmesh.noz+2];

  for(m=1;m<=E->sphere.caps_per_proc;m++)
    for(e=1;e<=E->lmesh.snel;e++) {
         uT =(E->slice.q_surf[m][E->sien[m][e].node[1]] +
              E->slice.q_surf[m][E->sien[m][e].node[2]] +
              E->slice.q_surf[m][E->sien[m][e].node[3]] +
              E->slice.q_surf[m][E->sien[m][e].node[4]])*0.25;
         el = e*E->lmesh.elz;
         sum_h[0] += uT*E->eco[m][el].area;
         sum_h[1] += E->eco[m][el].area;

         uT =(E->slice.q_botm[m][E->sien[m][e].node[1]] +
              E->slice.q_botm[m][E->sien[m][e].node[2]] +
              E->slice.q_botm[m][E->sien[m][e].node[3]] +
              E->slice.q_botm[m][E->sien[m][e].node[4]])*0.25;
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


/* Dannberg 2024 Appendix C: Galerkin energy RHS and support-point GLL
 * boundary mass. All ranks participate in the existing nodal exchange, which
 * includes radial messages. Only physical-boundary rows are assembled. */
static void heat_flux_CBF_boundary(struct All_variables *E, int top,
                                   double *slice_flux[NCS], int write_native, double stats[4])
{
    int m,e,a,i,d,node,bad=0,global_bad;
    const int lev=E->mesh.levmax, elz=E->lmesh.elz;
    const int target=top ? E->parallel.nprocz-1 : 0;
    const int active=E->parallel.me_loc[3]==target;
    const int side=top ? SIDE_TOP : SIDE_BOTTOM;
    const double sign=top ? 1.0 : -1.0;
    const double scale=(double)E->data.k0*E->data.ref_temperature
                         /(E->data.radius_km*1000.0);
    double *rhs[NCS], *mass[NCS], *qnodal[NCS];
    double er[9],x[4][3],dm[4],totals[2]={0,0},global_totals[2];
    double extrema[2]={-DBL_MAX,DBL_MAX}, global_extrema[2];
    void parallel_process_termination();

    for(m=1;m<=E->sphere.caps_per_proc;++m) {
        rhs[m]=(double *)calloc(E->lmesh.nno+2,sizeof(double));
        mass[m]=(double *)calloc(E->lmesh.nno+2,sizeof(double));
        qnodal[m]=(double *)calloc(E->lmesh.nno+2,sizeof(double));
        if(!rhs[m] || !mass[m] || !qnodal[m]) bad=1;
    }
    MPI_Allreduce(&bad,&global_bad,1,MPI_INT,MPI_MAX,E->parallel.world);
    if(global_bad) parallel_process_termination();
    if(active) for(m=1;m<=E->sphere.caps_per_proc;++m)
        for(e=top ? elz : 1;e<=E->lmesh.nel;e+=elz) {
            CBF_element_thermal_residual(E,m,e,er);
            for(a=0;a<4;++a) {
                node=E->ien[m][e].node[sidenodes[side][a+1]];
                if(!(E->node[m][node] & TBZ)) bad=1;
                for(d=0;d<3;++d) x[a][d]=E->X[lev][m][d+1][node];
            }
            CBF_face_gll_mass(x,dm);
            for(a=0;a<4;++a) {
                node=E->ien[m][e].node[sidenodes[side][a+1]];
                if(!(dm[a]>0) || !isfinite(dm[a]) ||
                   !isfinite(er[sidenodes[side][a+1]])) bad=1;
                rhs[m][node]+=er[sidenodes[side][a+1]];
                mass[m][node]+=dm[a];
            }
        }
    MPI_Allreduce(&bad,&global_bad,1,MPI_INT,MPI_MAX,E->parallel.world);
    if(global_bad) {
        if(E->parallel.me==0)
            fprintf(stderr,"CBF requires finite RHS, positive faces and Dirichlet radial boundaries\n");
        parallel_process_termination();
    }
    E->exchange_node_d(E,rhs,lev);
    E->exchange_node_d(E,mass,lev);
    if(active) for(m=1;m<=E->sphere.caps_per_proc;++m)
        for(i=1;i<=E->lmesh.nsf;++i) {
            node=top ? E->surf_node[m][i] : E->surf_node[m][i]-E->lmesh.noz+1;
            if(!(mass[m][node]>0)) bad=1;
            else {
                qnodal[m][node]=sign*scale*rhs[m][node]/mass[m][node];
                slice_flux[m][i]=qnodal[m][node];
                extrema[0]=fmax(extrema[0],qnodal[m][node]);
                extrema[1]=fmin(extrema[1],qnodal[m][node]);
                if(!isfinite(qnodal[m][node])) bad=1;
            }
        }
    /* Each physical face occurs once; use its LOCAL GLL weights, not the
     * globally assembled duplicate-node mass, for the integral. */
    if(active) for(m=1;m<=E->sphere.caps_per_proc;++m)
        for(e=top ? elz : 1;e<=E->lmesh.nel;e+=elz) {
            for(a=0;a<4;++a) {
                node=E->ien[m][e].node[sidenodes[side][a+1]];
                for(d=0;d<3;++d) x[a][d]=E->X[lev][m][d+1][node];
            }
            CBF_face_gll_mass(x,dm);
            for(a=0;a<4;++a) {
                node=E->ien[m][e].node[sidenodes[side][a+1]];
                totals[0]+=dm[a]*qnodal[m][node];
                totals[1]+=dm[a];
            }
        }
    MPI_Allreduce(&bad,&global_bad,1,MPI_INT,MPI_MAX,E->parallel.world);
    if(global_bad) parallel_process_termination();
    MPI_Allreduce(totals,global_totals,2,MPI_DOUBLE,MPI_SUM,E->parallel.world);
    MPI_Allreduce(&extrema[0],&global_extrema[0],1,MPI_DOUBLE,MPI_MAX,E->parallel.world);
    MPI_Allreduce(&extrema[1],&global_extrema[1],1,MPI_DOUBLE,MPI_MIN,E->parallel.world);
    if(stats) {
        double length=E->data.radius_km*1000.0;
        stats[0]=global_totals[0]*length*length;
        stats[1]=global_extrema[0]; stats[2]=global_extrema[1];
        stats[3]=global_totals[1]*length*length;
    }
    if(write_native && E->parallel.me==0) {
        fprintf(E->fp,"CBF_GLL_Q1 boundary=%s mean_W_m2=%.16e area_nd=%.16e\n",
                top ? "top" : "bottom",global_totals[0]/global_totals[1],global_totals[1]);
        fflush(E->fp);
    }
    if(write_native) CBF_native_boundary(E,top,rhs,mass,qnodal,global_totals);
    for(m=1;m<=E->sphere.caps_per_proc;++m) { free(rhs[m]); free(mass[m]); free(qnodal[m]); }
}

static void evaluate_heat_flux_CBF(struct All_variables *E, int write_native,
                                   double stats[2][4])
{
    int m,bad=0,global_bad;
    double *saved_adi[NCS],*saved_visc[NCS],*adi[NCS],*visc[NCS];
    struct CC saved_cc=E->element_Cc;
    struct CCX saved_ccx=E->element_Ccx;
    void parallel_process_termination();
    if(write_native && !E->output.CBF_use_advection) {
        if(E->parallel.me==0) fprintf(stderr,"CBF requires CBF_use_advection=on\n");
        parallel_process_termination();
    }
    for(m=1;m<=E->sphere.caps_per_proc;++m) {
        adi[m]=(double *)calloc(E->lmesh.nel+1,sizeof(double));
        visc[m]=(double *)calloc(E->lmesh.nel+1,sizeof(double));
        if(!adi[m] || !visc[m])bad=1;
    }
    MPI_Allreduce(&bad,&global_bad,1,MPI_INT,MPI_MAX,E->parallel.world);
    if(global_bad)parallel_process_termination();
    for(m=1;m<=E->sphere.caps_per_proc;++m) {
        CBF_heat_sources(E,m,adi[m],visc[m]);
        saved_adi[m]=E->heating_adi[m];saved_visc[m]=E->heating_visc[m];
        E->heating_adi[m]=adi[m];E->heating_visc[m]=visc[m];
    }
    if(write_native ? E->output.output_q_surf_CBF : E->mesh.toptbc==1)
        heat_flux_CBF_boundary(E,1,E->slice.q_surf_CBF,write_native,stats ? stats[0]:NULL);
    if(write_native ? E->output.output_q_botm_CBF : E->mesh.bottbc==1)
        heat_flux_CBF_boundary(E,0,E->slice.q_botm_CBF,write_native,stats ? stats[1]:NULL);
    for(m=1;m<=E->sphere.caps_per_proc;++m) {
        E->heating_adi[m]=saved_adi[m];E->heating_visc[m]=saved_visc[m];
        free(adi[m]);free(visc[m]);
    }
    E->element_Cc=saved_cc;E->element_Ccx=saved_ccx;
}


void heat_flux_CBF(struct All_variables *E)
{
    evaluate_heat_flux_CBF(E,1,NULL);
}

/* Every-step diagnostics use the full physical residual, independently of
 * native-file frequency/toggles. Non-Dirichlet boundaries are unavailable. */
void CBF_boundary_stats(struct All_variables *E, double stats[2][4])
{
    int i,j;
    for(i=0;i<2;i++) for(j=0;j<4;j++) stats[i][j]=NAN;
    evaluate_heat_flux_CBF(E,0,stats);
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
