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
#include <math.h>
#include <float.h>
#include <string.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "drive_solvers.h"

double global_vdot();
double vnorm_nonnewt();

static void write_stokes_diagnostics(struct All_variables *E);


/************************************************************/

void general_stokes_solver_setup(struct All_variables *E)
{
  int i, m;
  void construct_node_maps();

  if (E->control.NMULTIGRID || E->control.NASSEMBLE)
    construct_node_maps(E);
  else
    for (i=E->mesh.gridmin;i<=E->mesh.gridmax;i++)
      for (m=1;m<=E->sphere.caps_per_proc;m++)
	E->elt_k[i][m]=(struct EK *)malloc((E->lmesh.NEL[i]+1)*sizeof(struct EK));


  return;
}




void general_stokes_solver(struct All_variables *E)
{
  void solve_constrained_flow_iterative();
  void construct_stiffness_B_matrix();
  void velocities_conform_bcs();
  void assemble_forces();
  void sphere_harmonics_layer();
  void get_system_viscosity();
  void remove_rigid_rot();

  float vmag;

  double Udot_mag, dUdot_mag,omega[3];
  int m,count,i,j,k;

  double *oldU[NCS], *delta_U[NCS];

  const int nno = E->lmesh.nno;
  const int nel = E->lmesh.nel;
  const int nnov = E->lmesh.nnov;
  const int neq = E->lmesh.neq;
  const int vpts = vpoints[E->mesh.nsd];
  const int dims = E->mesh.nsd;
  const int addi_dof = additional_dof[dims];

  //velocities_conform_bcs(E,E->U);

  assemble_forces(E,0);

  if(E->monitor.solution_cycles==0 || E->viscosity.update_allowed) {
    get_system_viscosity(E,1,E->EVI[E->mesh.levmax],E->VI[E->mesh.levmax]);
    velocities_conform_bcs(E,E->U);
    construct_stiffness_B_matrix(E);
  }

  solve_constrained_flow_iterative(E);

  if (E->viscosity.SDEPV || E->viscosity.PDEPV) {

    for (m=1;m<=E->sphere.caps_per_proc;m++)  {
      delta_U[m] = (double *)malloc(neq*sizeof(double));
      oldU[m] = (double *)malloc(neq*sizeof(double));
      for(i=0;i<neq;i++)
	oldU[m][i]=0.0;
    }

    Udot_mag=dUdot_mag=0.0;
    count=1;

    while (1) {    
     

      for (m=1;m<=E->sphere.caps_per_proc;m++)
	for (i=0;i<neq;i++) {
	  delta_U[m][i] = E->U[m][i] - oldU[m][i];
	  oldU[m][i] = E->U[m][i];
	}

      Udot_mag  = sqrt(global_vdot(E,oldU,oldU,E->mesh.levmax));
      dUdot_mag = vnorm_nonnewt(E,delta_U,oldU,E->mesh.levmax);


      if(E->parallel.me==0){
	fprintf(stderr,"Stress dep. visc./plast.: DUdot = %.4e (%.4e) for iteration %d\n",
		dUdot_mag,Udot_mag,count);
	fprintf(E->fp,"Stress dep. visc./plast.: DUdot = %.4e (%.4e) for iteration %d\n",
		dUdot_mag,Udot_mag,count);
        //fprintf(E->fp,"sdev_misfit=%f\n",E->viscosity.sdepv_misfit);
	fflush(E->fp);
      }
      if ((count>1000) || (dUdot_mag < E->viscosity.sdepv_misfit))
	break;
      
      get_system_viscosity(E,0,E->EVI[E->mesh.levmax],E->VI[E->mesh.levmax]);
      velocities_conform_bcs(E,E->U);
      construct_stiffness_B_matrix(E);
      solve_constrained_flow_iterative(E);
      
      count++;

    } /*end while*/

    for (m=1;m<=E->sphere.caps_per_proc;m++)  {
      free((void *) oldU[m]);
      free((void *) delta_U[m]);
    }

  } /*end if SDEPV or PDEPV */

  write_stokes_diagnostics(E);

  /* remove the rigid rotation component from the velocity solution */
  if(E->sphere.caps == 12 && E->control.remove_rigid_rotation) {
      remove_rigid_rot(E);
  }

  return;
}

void general_stokes_solver_pseudo_surf(struct All_variables *E)
{
  void solve_constrained_flow_iterative_pseudo_surf();
  void construct_stiffness_B_matrix();
  void velocities_conform_bcs();
  void assemble_forces_pseudo_surf();
  void get_system_viscosity();
  void std_timestep();
  void remove_rigid_rot();
  void get_STD_freesurf(struct All_variables *, float**);

  float vmag;

  double Udot_mag, dUdot_mag;
  int m,count,i,j,k,topo_loop;

  double *oldU[NCS], *delta_U[NCS];

  const int nno = E->lmesh.nno;
  const int nel = E->lmesh.nel;
  const int nnov = E->lmesh.nnov;
  const int neq = E->lmesh.neq;
  const int vpts = vpoints[E->mesh.nsd];
  const int dims = E->mesh.nsd;
  const int addi_dof = additional_dof[dims];

  velocities_conform_bcs(E,E->U);

  E->monitor.stop_topo_loop = 0;
  E->monitor.topo_loop = 0;
  if(E->monitor.solution_cycles==0) std_timestep(E);
  while(E->monitor.stop_topo_loop == 0) {

	  assemble_forces_pseudo_surf(E,0);
	  if(E->monitor.solution_cycles==0 || E->viscosity.update_allowed) {
		  get_system_viscosity(E,1,E->EVI[E->mesh.levmax],E->VI[E->mesh.levmax]);
		  construct_stiffness_B_matrix(E);
	  }
	  solve_constrained_flow_iterative_pseudo_surf(E);

	  if (E->viscosity.SDEPV || E->viscosity.PDEPV) {

		  for (m=1;m<=E->sphere.caps_per_proc;m++)  {
			  delta_U[m] = (double *)malloc(neq*sizeof(double));
			  oldU[m] = (double *)malloc(neq*sizeof(double));
			  for(i=0;i<neq;i++)
				  oldU[m][i]=0.0;
		  }

		  Udot_mag=dUdot_mag=0.0;
		  count=1;

		  while (1) {

			  for (m=1;m<=E->sphere.caps_per_proc;m++)
				  for (i=0;i<neq;i++) {
					  delta_U[m][i] = E->U[m][i] - oldU[m][i];
					  oldU[m][i] = E->U[m][i];
				  }

			  Udot_mag  = sqrt(global_vdot(E,oldU,oldU,E->mesh.levmax));
			  dUdot_mag = vnorm_nonnewt(E,delta_U,oldU,E->mesh.levmax);

			  if(E->parallel.me==0){
				  fprintf(stderr,"Stress dependent viscosity: DUdot = %.4e (%.4e) for iteration %d\n",dUdot_mag,Udot_mag,count);
				  fprintf(E->fp,"Stress dependent viscosity: DUdot = %.4e (%.4e) for iteration %d\n",dUdot_mag,Udot_mag,count);
				  fflush(E->fp);
			  }

			  if (count>50 || dUdot_mag<E->viscosity.sdepv_misfit)
				  break;

			  get_system_viscosity(E,1,E->EVI[E->mesh.levmax],E->VI[E->mesh.levmax]);
			  construct_stiffness_B_matrix(E);
			  solve_constrained_flow_iterative_pseudo_surf(E);

			  count++;

		  } /*end while */
		  for (m=1;m<=E->sphere.caps_per_proc;m++)  {
			  free((void *) oldU[m]);
			  free((void *) delta_U[m]);
		  }

	  } /*end if SDEPV or PDEPV */
	  E->monitor.topo_loop++;
  }

  write_stokes_diagnostics(E);

  /* remove the rigid rotation component from the velocity solution */
  if(E->sphere.caps == 12 && E->control.remove_rigid_rotation)
      remove_rigid_rot(E);

  get_STD_freesurf(E,E->slice.freesurf);

  return;
}


struct power_stats {
    double total;
    double maximum;
    double minimum;
};


static double **allocate_nodal_field(struct All_variables *E)
{
    int m;
    double **field;

    field = (double **)malloc(NCS * sizeof(double *));
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        field[m] = (double *)calloc(E->lmesh.nno+1, sizeof(double));
    return field;
}


static double **allocate_element_field(struct All_variables *E)
{
    int m, field_size;
    double **field;

    field_size = max(E->lmesh.nel, E->lmesh.npno) + 1;
    field = (double **)malloc(NCS * sizeof(double *));
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        field[m] = (double *)calloc(field_size, sizeof(double));
    return field;
}


static double **allocate_equation_field(struct All_variables *E)
{
    int m;
    double **field;

    field = (double **)malloc(NCS * sizeof(double *));
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        field[m] = (double *)calloc(E->lmesh.neq+1, sizeof(double));
    return field;
}


static void free_nodal_field(struct All_variables *E, double **field)
{
    int m;

    for(m=1; m<=E->sphere.caps_per_proc; m++)
        free((void *)field[m]);
    free((void *)field);
}


static void free_equation_field(struct All_variables *E, double **field)
{
    free_nodal_field(E, field);
}


static void free_element_field(struct All_variables *E, double **field)
{
    free_nodal_field(E, field);
}


static struct power_stats reduce_power_stats(struct All_variables *E,
                                              double local_total,
                                              double local_max,
                                              double local_min)
{
    struct power_stats result;

    MPI_Allreduce(&local_total, &result.total, 1, MPI_DOUBLE, MPI_SUM,
                  E->parallel.world);
    MPI_Allreduce(&local_max, &result.maximum, 1, MPI_DOUBLE, MPI_MAX,
                  E->parallel.world);
    MPI_Allreduce(&local_min, &result.minimum, 1, MPI_DOUBLE, MPI_MIN,
                  E->parallel.world);
    return result;
}


static struct power_stats nodal_radial_power(struct All_variables *E,
                                             double **force,
                                             double scale)
{
    int m, e, i, a, node;
    double force_gp, velocity_gp, density, weight;
    double local_total, local_max, local_min;
    double rtf[4][9];
    float VV[4][9];
    struct Shape_function GN;
    struct Shape_function_dA dOmega;
    struct Shape_function_dx GNx;
    void get_global_shape_fn();
    void velo_from_element();

    const int ends = enodes[E->mesh.nsd];
    const int vpts = vpoints[E->mesh.nsd];
    const int lev = E->mesh.levmax;

    local_total = 0.0;
    local_max = -DBL_MAX;
    local_min = DBL_MAX;
    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        for(e=1; e<=E->lmesh.nel; e++) {
            get_global_shape_fn(E, e, &GN, &GNx, &dOmega, 0, 1,
                                rtf, lev, m);
            velo_from_element(E, VV, m, e, 1);
            for(i=1; i<=vpts; i++) {
                force_gp = 0.0;
                velocity_gp = 0.0;
                for(a=1; a<=ends; a++) {
                    node = E->ien[m][e].node[a];
                    weight = E->N.vpt[GNVINDEX(a,i)];
                    force_gp += force[m][node] * weight;
                    velocity_gp += VV[3][a] * weight;
                }
                density = scale * force_gp * velocity_gp;
                weight = dOmega.vpt[i]
                    * g_point[i].weight[E->mesh.nsd-1];
                local_total += density * weight;
                local_max = max(local_max, density);
                local_min = min(local_min, density);
            }
        }
    }
    return reduce_power_stats(E, local_total, local_max, local_min);
}


static struct power_stats pressure_power(struct All_variables *E,
                                         double scale)
{
    int m, e, i, a;
    double velocity_gp, density, weight, beta, volume;
    double gradient_power_density;
    double local_total, local_max, local_min;
    double **div_u;
    double rtf[4][9];
    float VV[4][9];
    struct Shape_function GN;
    struct Shape_function_dA dOmega;
    struct Shape_function_dx GNx;
    void get_global_shape_fn();
    void velo_from_element();
    void assemble_div_u();

    const int ends = enodes[E->mesh.nsd];
    const int vpts = vpoints[E->mesh.nsd];
    const int lev = E->mesh.levmax;

    local_total = 0.0;
    local_max = -DBL_MAX;
    local_min = DBL_MAX;
    div_u = allocate_nodal_field(E);
    assemble_div_u(E, E->U, div_u, lev);
    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        for(e=1; e<=E->lmesh.nel; e++) {
            beta = (E->control.inv_gruneisen != 0.0)
                ? E->refstate.ala_beta[((e-1) % E->lmesh.elz)+1]
                : 0.0;
            get_global_shape_fn(E, e, &GN, &GNx, &dOmega, 0, 1,
                                rtf, lev, m);
            velo_from_element(E, VV, m, e, 1);
            volume = E->eco[m][e].area;
            gradient_power_density = (volume > 0.0)
                ? scale * E->P[m][e] * div_u[m][e] / volume
                : 0.0;
            for(i=1; i<=vpts; i++) {
                velocity_gp = 0.0;
                for(a=1; a<=ends; a++)
                    velocity_gp += VV[3][a]
                        * E->N.vpt[GNVINDEX(a,i)];
                /* Net pressure contribution: ALA pressure body force minus
                 * the standard discrete u.grad(p) work. */
                density = scale * (-beta * E->P[m][e]) * velocity_gp
                    - gradient_power_density;
                weight = dOmega.vpt[i]
                    * g_point[i].weight[E->mesh.nsd-1];
                local_total += density * weight;
                local_max = max(local_max, density);
                local_min = min(local_min, density);
            }
        }
    }
    free_nodal_field(E, div_u);
    return reduce_power_stats(E, local_total, local_max, local_min);
}


static struct power_stats viscous_power(struct All_variables *E)
{
    int m, e, i;
    double value, viscosity, scale;
    double local_total, local_max, local_min;
    float *strain_sqr;
    void strain_rate_2_inv();

    const int vpts = vpoints[E->mesh.nsd];

    strain_sqr = (float *)malloc((E->lmesh.nel+1)*sizeof(float));
    scale = E->control.disptn_number / E->control.Atemp / vpts;
    local_total = 0.0;
    local_max = -DBL_MAX;
    local_min = DBL_MAX;
    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        strain_rate_2_inv(E, m, strain_sqr, 0);
        for(e=1; e<=E->lmesh.nel; e++) {
            viscosity = 0.0;
            for(i=1; i<=vpts; i++)
                viscosity += E->EVi[m][(e-1)*vpts+i];
            value = scale * viscosity * strain_sqr[e];
            local_total += value * E->eco[m][e].area;
            local_max = max(local_max, value);
            local_min = min(local_min, value);
        }
    }
    free((void *)strain_sqr);
    return reduce_power_stats(E, local_total, local_max, local_min);
}


static void assemble_grad_p_unstripped(struct All_variables *E,
                                       double **grad_p)
{
    int m, e, a, p, node, j1, j2, j3;
    const int lev = E->mesh.levmax;
    const int ends = enodes[E->mesh.nsd];
    const int dims = E->mesh.nsd;

    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        memset(grad_p[m], 0, E->lmesh.neq*sizeof(double));
        for(e=1; e<=E->lmesh.nel; e++) {
            for(a=1; a<=ends; a++) {
                p = (a-1)*dims;
                node = E->IEN[lev][m][e].node[a];
                j1 = E->ID[lev][m][node].doff[1];
                j2 = E->ID[lev][m][node].doff[2];
                j3 = E->ID[lev][m][node].doff[3];
                grad_p[m][j1] += E->elt_del[lev][m][e].g[p][0]
                    * E->P[m][e];
                grad_p[m][j2] += E->elt_del[lev][m][e].g[p+1][0]
                    * E->P[m][e];
                grad_p[m][j3] += E->elt_del[lev][m][e].g[p+2][0]
                    * E->P[m][e];
            }
        }
    }
    (E->solver.exchange_id_d)(E, grad_p, lev);
}


static void add_element_force_to(struct All_variables *E, int element,
                                 const double *element_force, int cap,
                                 double **force)
{
    int a, node, p;
    const int dims = E->mesh.nsd;
    const int ends = enodes[dims];

    for(a=1; a<=ends; a++) {
        node = E->ien[cap][element].node[a];
        p = (a-1)*dims;
        force[cap][E->id[cap][node].doff[1]] += element_force[p];
        force[cap][E->id[cap][node].doff[2]] += element_force[p+1];
        force[cap][E->id[cap][node].doff[3]] += element_force[p+2];
    }
}


static struct power_stats plate_power(struct All_variables *E, double scale)
{
    int m, e, i, a, node, d, eq;
    double elt_f[24], contribution;
    double local_max, local_min;
    double **ku, **grad_p, **external_force, **plate_velocity;
    struct power_stats result;
    void assemble_del2_u();
    void get_elt_f();
    void get_elt_tr();

    const int lev = E->mesh.levmax;
    const int neq = E->lmesh.neq;
    const unsigned int flags[4] = {0, VBX, VBY, VBZ};

    ku = allocate_equation_field(E);
    grad_p = allocate_equation_field(E);
    external_force = allocate_equation_field(E);
    plate_velocity = allocate_equation_field(E);

    assemble_del2_u(E, E->U, ku, lev, 0);
    assemble_grad_p_unstripped(E, grad_p);

    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        memset(external_force[m], 0, neq*sizeof(double));
        for(e=1; e<=E->lmesh.nel; e++) {
            get_elt_f(E, e, elt_f, 0, m);
            add_element_force_to(E, e, elt_f, m, external_force);
        }
        for(i=1; i<=E->boundary.nel; i++) {
            e = E->boundary.element[m][i];
            for(a=0; a<24; a++) elt_f[a] = 0.0;
            for(a=SIDE_BEGIN; a<=SIDE_END; a++)
                get_elt_tr(E, i, a, elt_f, m);
            add_element_force_to(E, e, elt_f, m, external_force);
        }
    }
    (E->solver.exchange_id_d)(E, external_force, lev);

    /* The constrained-DOF reaction is K*u + grad(p) - external force. */
    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        for(i=0; i<neq; i++)
            external_force[m][i] = ku[m][i] + grad_p[m][i]
                - external_force[m][i];
        if(E->parallel.me_loc[3] == E->parallel.nprocz-1) {
            for(node=E->lmesh.noz; node<=E->lmesh.nno;
                node+=E->lmesh.noz) {
                for(d=1; d<=E->mesh.nsd; d++) {
                    if(E->node[m][node] & flags[d]) {
                        eq = E->id[m][node].doff[d];
                        plate_velocity[m][eq] = E->U[m][eq];
                    }
                }
            }
        }
    }

    result.total = scale * global_vdot(E, plate_velocity, external_force, lev);
    local_max = -DBL_MAX;
    local_min = DBL_MAX;
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(i=0; i<neq; i++) {
            if(plate_velocity[m][i] != 0.0) {
                contribution = scale * plate_velocity[m][i]
                    * external_force[m][i];
                local_max = max(local_max, contribution);
                local_min = min(local_min, contribution);
            }
        }
    if(local_max == -DBL_MAX) local_max = 0.0;
    if(local_min == DBL_MAX) local_min = 0.0;
    MPI_Allreduce(&local_max, &result.maximum, 1, MPI_DOUBLE, MPI_MAX,
                  E->parallel.world);
    MPI_Allreduce(&local_min, &result.minimum, 1, MPI_DOUBLE, MPI_MIN,
                  E->parallel.world);

    free_equation_field(E, ku);
    free_equation_field(E, grad_p);
    free_equation_field(E, external_force);
    free_equation_field(E, plate_velocity);
    return result;
}


static void print_power_row(FILE *fp, const char *name,
                            struct power_stats value)
{
    fprintf(fp, "%-20s  %+16.8e  %+16.8e  %+16.8e\n",
            name, value.total, value.maximum, value.minimum);
}


static void write_mechanical_power(struct All_variables *E)
{
    int m, i, j, nz, phase_index;
    double scale;
    double **thermal, **chemical, **phase[PHASE_TRANSITIONS], **phase_total;
    struct power_stats pplate, qvisc, wthermal, wchemical;
    struct power_stats wphase[PHASE_TRANSITIONS], wphase_total, wpressure;
    struct power_stats residual;
    static const char *phase_names[PHASE_TRANSITIONS] = {
        "Wphase_410", "Wphase_520", "Wphase_660"
    };
    static const char separator[] =
        "==============================================================================\n";
    static const char rule[] =
        "------------------------------------------------------------------------------\n";
    void remove_horiz_ave2();

    scale = E->control.disptn_number / E->control.Atemp;
    thermal = allocate_nodal_field(E);
    chemical = allocate_nodal_field(E);
    phase_total = allocate_nodal_field(E);
    for(phase_index=0; phase_index<PHASE_TRANSITIONS; phase_index++)
        phase[phase_index] = allocate_nodal_field(E);

    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        for(i=1; i<=E->lmesh.nno; i++) {
            nz = ((i-1) % E->lmesh.noz) + 1;
            thermal[m][i] = E->control.Atemp * E->refstate.rho[nz]
                * E->refstate.thermal_expansivity[nz] * E->T[m][i];
            if(E->control.tracer && E->composition.ichemical_buoyancy)
                for(j=0; j<E->composition.ncomp; j++)
                    chemical[m][i] -= E->control.Atemp
                        * E->composition.buoyancy_ratio[j]
                        * E->composition.comp_node[m][j][i];
            for(phase_index=0; phase_index<PHASE_TRANSITIONS; phase_index++) {
                if(E->control.phase[phase_index].Ra != 0.0)
                    phase[phase_index][m][i] =
                        -E->control.phase[phase_index].Ra
                        * E->phase_B[phase_index][m][i];
                else
                    phase[phase_index][m][i] = 0.0;
            }
        }
    }

    remove_horiz_ave2(E, thermal);
    remove_horiz_ave2(E, chemical);
    for(phase_index=0; phase_index<PHASE_TRANSITIONS; phase_index++)
        remove_horiz_ave2(E, phase[phase_index]);

    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(i=1; i<=E->lmesh.nno; i++) {
            nz = ((i-1) % E->lmesh.noz) + 1;
            thermal[m][i] *= E->refstate.gravity[nz];
            chemical[m][i] *= E->refstate.gravity[nz];
            for(phase_index=0; phase_index<PHASE_TRANSITIONS; phase_index++) {
                phase[phase_index][m][i] *= E->refstate.gravity[nz];
                phase_total[m][i] += phase[phase_index][m][i];
            }
        }

    pplate = plate_power(E, scale);
    qvisc = viscous_power(E);
    wthermal = nodal_radial_power(E, thermal, scale);
    wchemical = nodal_radial_power(E, chemical, scale);
    for(phase_index=0; phase_index<PHASE_TRANSITIONS; phase_index++)
        wphase[phase_index] = nodal_radial_power(E, phase[phase_index], scale);
    wphase_total = nodal_radial_power(E, phase_total, scale);
    wpressure = pressure_power(E, scale);

    residual.total = pplate.total + wthermal.total + wchemical.total
        + wphase_total.total + wpressure.total - qvisc.total;
    residual.maximum = NAN;
    residual.minimum = NAN;

    if(E->parallel.me == 0) {
        fputs(separator, E->fp);
        fprintf(E->fp, "MECHANICAL_POWER  step=%d  scale=Di/Atemp\n",
                E->monitor.solution_cycles);
        fputs(rule, E->fp);
        fprintf(E->fp, "%-20s  %16s  %16s  %16s\n",
                "TERM", "TOTAL", "MAX", "MIN");
        fputs(rule, E->fp);
        print_power_row(E->fp, "Pplate", pplate);
        print_power_row(E->fp, "Qvisc", qvisc);
        print_power_row(E->fp, "Wthermal", wthermal);
        print_power_row(E->fp, "Wchemical", wchemical);
        for(phase_index=0; phase_index<PHASE_TRANSITIONS; phase_index++)
            print_power_row(E->fp, phase_names[phase_index],
                            wphase[phase_index]);
        print_power_row(E->fp, "Wphase_total", wphase_total);
        print_power_row(E->fp, "Wpressure", wpressure);
        print_power_row(E->fp, "Rmechanical", residual);
        fputs(separator, E->fp);
        fflush(E->fp);
    }

    free_nodal_field(E, thermal);
    free_nodal_field(E, chemical);
    free_nodal_field(E, phase_total);
    for(phase_index=0; phase_index<PHASE_TRANSITIONS; phase_index++)
        free_nodal_field(E, phase[phase_index]);
}


static void write_ala_residual(struct All_variables *E)
{
    int m, e;
    double volume, rdiv, rbeta, rala, relative, denominator;
    double scale_density, characteristic_scale, active_scale_cutoff;
    double global_cancellation_l2, elementwise_relative_l2;
    double active_relative_l2, active_volume_fraction;
    double active_fraction_above_threshold;
    double local_sum[6], global_sum[5], global_linf;
    double local_active[3], global_active[3];
    double **div_u, **beta_u;
    static const double active_scale_fraction = 1.0e-6;
    static const double relative_alert_threshold = 1.0e-1;
    static const char separator[] =
        "==============================================================================\n";
    static const char rule[] =
        "------------------------------------------------------------------------------\n";
    void assemble_div_u();
    void assemble_c_u();

    div_u = allocate_element_field(E);
    beta_u = allocate_element_field(E);
    assemble_div_u(E, E->U, div_u, E->mesh.levmax);
    assemble_c_u(E, E->U, beta_u, E->mesh.levmax);

    for(e=0; e<6; e++) local_sum[e] = 0.0;
    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        for(e=1; e<=E->lmesh.nel; e++) {
            volume = E->eco[m][e].area;
            if(volume <= 0.0) continue;
            rdiv = div_u[m][e];
            rbeta = beta_u[m][e];
            rala = (rdiv + rbeta) / volume;
            denominator = fabs(rdiv) + fabs(rbeta);
            relative = fabs(rdiv + rbeta)
                / max(denominator, 1.0e-30 * volume);
            scale_density = denominator / volume;
            local_sum[0] += volume * rala * rala;
            local_sum[1] += volume * fabs(rala);
            local_sum[2] += volume * relative * relative;
            local_sum[3] += volume;
            local_sum[4] += volume * scale_density * scale_density;
            local_sum[5] = max(local_sum[5], fabs(rala));
        }
    }
    MPI_Allreduce(local_sum, global_sum, 5, MPI_DOUBLE, MPI_SUM,
                  E->parallel.world);
    MPI_Allreduce(&local_sum[5], &global_linf, 1, MPI_DOUBLE, MPI_MAX,
                  E->parallel.world);

    characteristic_scale = sqrt(global_sum[4]
                                  / max(global_sum[3], 1.0e-300));
    active_scale_cutoff = active_scale_fraction * characteristic_scale;
    global_cancellation_l2 = sqrt(global_sum[0]
                                   / max(global_sum[4], 1.0e-300));
    elementwise_relative_l2 = sqrt(global_sum[2]
                                    / max(global_sum[3], 1.0e-300));

    for(e=0; e<3; e++) local_active[e] = 0.0;
    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        for(e=1; e<=E->lmesh.nel; e++) {
            volume = E->eco[m][e].area;
            if(volume <= 0.0) continue;
            rdiv = div_u[m][e];
            rbeta = beta_u[m][e];
            denominator = fabs(rdiv) + fabs(rbeta);
            scale_density = denominator / volume;
            if(scale_density < active_scale_cutoff) continue;
            relative = fabs(rdiv + rbeta)
                / max(denominator, 1.0e-30 * volume);
            local_active[0] += volume;
            local_active[1] += volume * relative * relative;
            if(relative > relative_alert_threshold)
                local_active[2] += volume;
        }
    }
    MPI_Allreduce(local_active, global_active, 3, MPI_DOUBLE, MPI_SUM,
                  E->parallel.world);
    active_volume_fraction = global_active[0]
        / max(global_sum[3], 1.0e-300);
    active_relative_l2 = sqrt(global_active[1]
                               / max(global_active[0], 1.0e-300));
    active_fraction_above_threshold = global_active[2]
        / max(global_active[0], 1.0e-300);

    if(E->parallel.me == 0) {
        fputs(separator, E->fp);
        fprintf(E->fp, "ALA_CONTINUITY_RESIDUAL  step=%d\n",
                E->monitor.solution_cycles);
        fputs(rule, E->fp);
        fprintf(E->fp, "%-30s  %16s\n", "METRIC", "VALUE");
        fputs(rule, E->fp);
        fprintf(E->fp, "%-30s  %+16.8e\n", "L2",
                sqrt(global_sum[0] / global_sum[3]));
        fprintf(E->fp, "%-30s  %+16.8e\n", "Linf", global_linf);
        fprintf(E->fp, "%-30s  %+16.8e\n",
                "volume_weighted_mean_abs",
                global_sum[1] / global_sum[3]);
        fprintf(E->fp, "%-30s  %+16.8e\n", "solver_div_v",
                E->monitor.incompressibility);
        fprintf(E->fp, "%-30s  %+16.8e\n", "global_cancellation_L2",
                global_cancellation_l2);
        fprintf(E->fp, "%-30s  %+16.8e\n",
                "elementwise_relative_L2", elementwise_relative_l2);
        fprintf(E->fp, "%-30s  %+16.8e\n", "active_scale_cutoff",
                active_scale_cutoff);
        fprintf(E->fp, "%-30s  %+16.8e\n", "active_volume_fraction",
                active_volume_fraction);
        fprintf(E->fp, "%-30s  %+16.8e\n", "active_relative_L2",
                active_relative_l2);
        fprintf(E->fp, "%-30s  %+16.8e\n",
                "active_fraction_rel_gt_0.1",
                active_fraction_above_threshold);
        fputs(separator, E->fp);
        fflush(E->fp);
    }

    free_element_field(E, div_u);
    free_element_field(E, beta_u);
}


static void write_stokes_diagnostics(struct All_variables *E)
{
    if(!E->control.ala_pressure_buoyancy)
        return;
    write_mechanical_power(E);
    write_ala_residual(E);
}
