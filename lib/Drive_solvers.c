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
#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "drive_solvers.h"
#include "phase_change.h"

double global_vdot();
double global_pdot();
double vnorm_nonnewt();
void myerror(struct All_variables *,char *);

static void write_stokes_diagnostics(struct All_variables *E);

/* VC1 is deliberately process-level continuation: every viscosity cap starts
 * a new CitcomS process, so no Krylov vector or numerical cache can cross a
 * stage boundary.  Only the accepted equation-point U and P arrays are moved
 * between processes.  Temperature, tracers/composition, phase state, time and
 * viscosity are reconstructed from the unchanged Tnew inputs in each process.
 */
#define STRICT_ALA_VC1_WARM_MAGIC "ALAVC1UP"
#define STRICT_ALA_VC1_WARM_VERSION 1U

static unsigned long long strict_ala_vc1_hash_bytes(
    unsigned long long hash,const void *data,size_t size)
{
  const unsigned char *bytes=(const unsigned char *)data;
  size_t i;
  for(i=0;i<size;i++) {
    hash^=(unsigned long long)bytes[i];
    hash*=1099511628211ULL;
  }
  return hash;
}

static void strict_ala_vc1_warm_path(struct All_variables *E,
    const char *directory,char *path,size_t size,const char *suffix)
{
  int written=snprintf(path,size,"%s/warm_state.rank%06d.bin%s",directory,
                       E->parallel.me,suffix);
  if(written<0 || (size_t)written>=size)
    myerror(E,"VC1 warm-state path is too long");
}

static void strict_ala_vc1_transfer_header(struct All_variables *E,FILE *fp,
    int writing,unsigned long long *hash)
{
  char magic[8];
  uint32_t header[9],expected[9];
  size_t count;
  memcpy(magic,STRICT_ALA_VC1_WARM_MAGIC,sizeof(magic));
  expected[0]=STRICT_ALA_VC1_WARM_VERSION;
  expected[1]=(uint32_t)E->parallel.me;
  expected[2]=(uint32_t)E->parallel.nproc;
  expected[3]=(uint32_t)E->sphere.caps_per_proc;
  expected[4]=(uint32_t)E->lmesh.neq;
  expected[5]=(uint32_t)E->lmesh.npno;
  expected[6]=(uint32_t)E->lmesh.nox;
  expected[7]=(uint32_t)E->lmesh.noz;
  expected[8]=(uint32_t)E->lmesh.noy;
  if(writing) {
    if(fwrite(magic,1,sizeof(magic),fp)!=sizeof(magic) ||
       fwrite(expected,sizeof(uint32_t),9,fp)!=9)
      myerror(E,"Unable to write VC1 warm-state header");
    *hash=strict_ala_vc1_hash_bytes(*hash,magic,sizeof(magic));
    *hash=strict_ala_vc1_hash_bytes(*hash,expected,sizeof(expected));
    return;
  }
  count=fread(magic,1,sizeof(magic),fp);
  count+=fread(header,sizeof(uint32_t),9,fp)*sizeof(uint32_t);
  if(count!=sizeof(magic)+sizeof(header) ||
     memcmp(magic,STRICT_ALA_VC1_WARM_MAGIC,sizeof(magic))!=0 ||
     memcmp(header,expected,sizeof(header))!=0)
    myerror(E,"VC1 warm-state header/mesh/decomposition mismatch");
  *hash=strict_ala_vc1_hash_bytes(*hash,magic,sizeof(magic));
  *hash=strict_ala_vc1_hash_bytes(*hash,header,sizeof(header));
}

static void strict_ala_vc1_transfer_warm_state(struct All_variables *E,
    const char *directory,int writing)
{
  void v_from_vector();
  void p_to_nodes();
  char path[1024],temporary[1056];
  FILE *fp;
  int m;
  unsigned long long hash=1469598103934665603ULL,stored_hash=0,global_hash=0;
  unsigned long long p_hash=1469598103934665603ULL,u_hash=1469598103934665603ULL;
  unsigned long long global_p_hash=0,global_u_hash=0;
  strict_ala_vc1_warm_path(E,directory,path,sizeof(path),"");
  if(writing) {
    strict_ala_vc1_warm_path(E,directory,temporary,sizeof(temporary),".tmp");
    fp=fopen(temporary,"wb");
  }
  else
    fp=fopen(path,"rb");
  if(fp==NULL) {
    fprintf(stderr,"VC1 warm-state open failed path=%s errno=%d\n",
            writing ? temporary : path,errno);
    myerror(E,"Unable to open VC1 warm state");
  }
  strict_ala_vc1_transfer_header(E,fp,writing,&hash);
  for(m=1;m<=E->sphere.caps_per_proc;m++) {
    if(writing) {
      if(fwrite(&E->P[m][1],sizeof(double),E->lmesh.npno,fp)
             !=(size_t)E->lmesh.npno ||
         fwrite(E->U[m],sizeof(double),E->lmesh.neq,fp)
             !=(size_t)E->lmesh.neq)
        myerror(E,"Unable to write VC1 U/P warm state");
    }
    else if(fread(&E->P[m][1],sizeof(double),E->lmesh.npno,fp)
                 !=(size_t)E->lmesh.npno ||
            fread(E->U[m],sizeof(double),E->lmesh.neq,fp)
                 !=(size_t)E->lmesh.neq)
      myerror(E,"Unable to read VC1 U/P warm state");
    E->P[m][0]=0.0;
    hash=strict_ala_vc1_hash_bytes(hash,&E->P[m][1],
                                  E->lmesh.npno*sizeof(double));
    hash=strict_ala_vc1_hash_bytes(hash,E->U[m],
                                  E->lmesh.neq*sizeof(double));
    p_hash=strict_ala_vc1_hash_bytes(p_hash,&E->P[m][1],
                                    E->lmesh.npno*sizeof(double));
    u_hash=strict_ala_vc1_hash_bytes(u_hash,E->U[m],
                                    E->lmesh.neq*sizeof(double));
  }
  if(writing) {
    if(fwrite(&hash,sizeof(hash),1,fp)!=1 || fclose(fp)!=0 ||
       rename(temporary,path)!=0)
      myerror(E,"Unable to finalize VC1 U/P warm state");
  }
  else {
    if(fread(&stored_hash,sizeof(stored_hash),1,fp)!=1 ||
       fgetc(fp)!=EOF || fclose(fp)!=0 || stored_hash!=hash)
      myerror(E,"VC1 warm-state checksum/length mismatch");
    v_from_vector(E);
    p_to_nodes(E,E->P,E->NP,E->mesh.levmax);
  }
  MPI_Allreduce(&hash,&global_hash,1,MPI_UNSIGNED_LONG_LONG,MPI_BXOR,
                E->parallel.world);
  MPI_Allreduce(&p_hash,&global_p_hash,1,MPI_UNSIGNED_LONG_LONG,MPI_BXOR,
                E->parallel.world);
  MPI_Allreduce(&u_hash,&global_u_hash,1,MPI_UNSIGNED_LONG_LONG,MPI_BXOR,
                E->parallel.world);
  if(E->parallel.me==0) {
    fprintf(E->fp,"STRICT_ALA_VC1_WARM_STATE action=%s directory=%s "
            "global_xor_checksum=%016llx u_checksum=%016llx p_checksum=%016llx "
            "pressure_gauge=exact_no_regauge\n",
            writing ? "write" : "read",directory,global_hash,global_u_hash,
            global_p_hash);
    fflush(E->fp);
  }
}

static void strict_ala_vc1_viscosity_audit(struct All_variables *E)
{
  int m,e,j;
  long long local_count=0,local_lower=0,local_upper=0;
  long long global_count,global_lower,global_upper;
  double value,upper_cap,rad,local_min=DBL_MAX,local_max=-DBL_MAX;
  double global_min,global_max;
  unsigned long long local_hash=1469598103934665603ULL,global_hash;
  const int vpts=vpoints[E->mesh.nsd];
  for(m=1;m<=E->sphere.caps_per_proc;m++)
    for(e=1;e<=E->lmesh.nel;e++)
      for(j=1;j<=vpts;j++) {
        value=E->EVI[E->mesh.levmax][m][(e-1)*vpts+j];
        rad=0.5*(E->sx[m][3][E->ien[m][e].node[1]]+
                 E->sx[m][3][E->ien[m][e].node[8]]);
        upper_cap=(rad>0.89641) ? E->viscosity.max_value
                               : 5.0*E->viscosity.max_value;
        local_min=min(local_min,value); local_max=max(local_max,value);
        if(E->viscosity.MIN && value==E->viscosity.min_value) local_lower++;
        if(E->viscosity.MAX && value==upper_cap) local_upper++;
        local_count++;
        local_hash=strict_ala_vc1_hash_bytes(local_hash,&value,sizeof(value));
      }
  MPI_Allreduce(&local_min,&global_min,1,MPI_DOUBLE,MPI_MIN,E->parallel.world);
  MPI_Allreduce(&local_max,&global_max,1,MPI_DOUBLE,MPI_MAX,E->parallel.world);
  MPI_Allreduce(&local_count,&global_count,1,MPI_LONG_LONG,MPI_SUM,
                E->parallel.world);
  MPI_Allreduce(&local_lower,&global_lower,1,MPI_LONG_LONG,MPI_SUM,
                E->parallel.world);
  MPI_Allreduce(&local_upper,&global_upper,1,MPI_LONG_LONG,MPI_SUM,
                E->parallel.world);
  MPI_Allreduce(&local_hash,&global_hash,1,MPI_UNSIGNED_LONG_LONG,MPI_BXOR,
                E->parallel.world);
  if(E->parallel.me==0) {
    fprintf(E->fp,"STRICT_ALA_VC1_VISCOSITY configured_visc_max=%.17e "
            "eta_min=%.17e eta_max=%.17e eta_ratio=%.17e "
            "lower_clamp_fraction=%.17e upper_clamp_fraction=%.17e "
            "sample_count=%lld global_xor_checksum=%016llx\n",
            E->viscosity.max_value,global_min,global_max,
            global_max/max(global_min,DBL_MIN),
            (double)global_lower/max((double)global_count,1.0),
            (double)global_upper/max((double)global_count,1.0),
            global_count,global_hash);
    fflush(E->fp);
  }
}

static unsigned long long strict_ala_vc1_physical_state_hash(
    struct All_variables *E,unsigned long long *local_hash)
{
  int m,i,q;
  unsigned long long local=1469598103934665603ULL,global;
  local=strict_ala_vc1_hash_bytes(local,&E->monitor.solution_cycles,
                                  sizeof(E->monitor.solution_cycles));
  local=strict_ala_vc1_hash_bytes(local,&E->monitor.elapsed_time,
                                  sizeof(E->monitor.elapsed_time));
  local=strict_ala_vc1_hash_bytes(local,&E->advection.timestep,
                                  sizeof(E->advection.timestep));
  for(m=1;m<=E->sphere.caps_per_proc;m++) {
    local=strict_ala_vc1_hash_bytes(local,&E->T[m][1],
                                   E->lmesh.nno*sizeof(double));
    local=strict_ala_vc1_hash_bytes(local,&E->Tdot[m][1],
                                   E->lmesh.nno*sizeof(double));
    if(E->composition.on)
      for(q=0;q<E->composition.ncomp;q++) {
        local=strict_ala_vc1_hash_bytes(local,&E->composition.comp_el[m][q][1],
                                       E->lmesh.nel*sizeof(double));
        local=strict_ala_vc1_hash_bytes(local,&E->composition.comp_node[m][q][1],
                                       E->lmesh.nno*sizeof(double));
      }
    if(E->control.tracer) {
      local=strict_ala_vc1_hash_bytes(local,&E->trace.ntracers[m],
                                     sizeof(E->trace.ntracers[m]));
      for(q=0;q<E->trace.number_of_basic_quantities;q++)
        local=strict_ala_vc1_hash_bytes(local,&E->trace.basicq[m][q][1],
          E->trace.ntracers[m]*sizeof(double));
      for(q=0;q<E->trace.number_of_extra_quantities;q++)
        local=strict_ala_vc1_hash_bytes(local,&E->trace.extraq[m][q][1],
          E->trace.ntracers[m]*sizeof(double));
      local=strict_ala_vc1_hash_bytes(local,&E->trace.ielement[m][1],
          E->trace.ntracers[m]*sizeof(int));
      for(i=0;i<E->trace.nflavors;i++)
        local=strict_ala_vc1_hash_bytes(local,&E->trace.ntracer_flavor[m][i][1],
                                       E->lmesh.nel*sizeof(int));
    }
  }
  MPI_Allreduce(&local,&global,1,MPI_UNSIGNED_LONG_LONG,MPI_BXOR,
                E->parallel.world);
  if(local_hash!=NULL) *local_hash=local;
  return global;
}

static void strict_ala_vc1_compare_states(struct All_variables *E,
    const char *left,const char *right)
{
  double *left_u[NCS],*left_p[NCS];
  double u_difference,u_reference,p_difference,p_reference;
  int m,i;
  strict_ala_vc1_transfer_warm_state(E,left,0);
  for(m=1;m<=E->sphere.caps_per_proc;m++) {
    left_u[m]=(double *)malloc((E->lmesh.neq+1)*sizeof(double));
    left_p[m]=(double *)malloc((E->lmesh.npno+1)*sizeof(double));
    if(left_u[m]==NULL || left_p[m]==NULL)
      myerror(E,"Unable to allocate VC1 comparison state");
    memcpy(left_u[m],E->U[m],E->lmesh.neq*sizeof(double));
    left_u[m][E->lmesh.neq]=0.0;
    memcpy(left_p[m],E->P[m],(E->lmesh.npno+1)*sizeof(double));
  }
  strict_ala_vc1_transfer_warm_state(E,right,0);
  for(m=1;m<=E->sphere.caps_per_proc;m++) {
    for(i=0;i<E->lmesh.neq;i++) E->U[m][i]-=left_u[m][i];
    E->U[m][E->lmesh.neq]=0.0;
    for(i=1;i<=E->lmesh.npno;i++) E->P[m][i]-=left_p[m][i];
  }
  u_difference=sqrt(max(global_vdot(E,E->U,E->U,E->mesh.levmax),0.0));
  u_reference=sqrt(max(global_vdot(E,left_u,left_u,E->mesh.levmax),0.0));
  p_difference=sqrt(max(global_pdot(E,E->P,E->P,E->mesh.levmax),0.0));
  p_reference=sqrt(max(global_pdot(E,left_p,left_p,E->mesh.levmax),0.0));
  if(E->parallel.me==0) {
    fprintf(E->fp,"STRICT_ALA_VC1_STATE_COMPARISON left=%s right=%s "
            "relative_velocity_difference=%.17e pressure_difference=%.17e "
            "relative_pressure_difference=%.17e "
            "pressure_gauge=exact_no_regauge_global_pdot\n",left,right,
            u_difference/max(u_reference,DBL_MIN),p_difference,
            p_difference/max(p_reference,DBL_MIN));
    fflush(E->fp);
  }
  for(m=1;m<=E->sphere.caps_per_proc;m++) {
    free(left_u[m]); free(left_p[m]);
  }
}


/************************************************************/

void general_stokes_solver_setup(struct All_variables *E)
{
  int i, m;
  void construct_node_maps();

  if ((E->control.NMULTIGRID || E->control.EMULTIGRID) &&
      E->parallel.me == 0) {
    fprintf(E->fp,
            "Velocity MG hierarchy levels=%d range=%d:%d cycle=%d "
            "smooth=(down:%d,up:%d,fine:%d,coarse_max:%d) "
            "smoother=%s coarse_solver=%s vanka_velocity_damping=%g "
            "vanka_pressure_damping=%g "
            "vanka_regularization=%g\n",
            E->mesh.levmax - E->mesh.levmin + 1,
            E->mesh.levmin, E->mesh.levmax, E->control.mg_cycle,
            E->control.down_heavy, E->control.up_heavy,
            E->control.v_steps_high, E->control.v_steps_low,
            E->control.ala_element_vanka_smoother
              ? "full_element_vanka" : "point_gauss_seidel",
            E->control.ala_element_vanka_smoother
              ? "point_gauss_seidel" : "same_as_smoother",
            E->control.ala_element_vanka_damping,
            E->control.ala_element_vanka_pressure_damping,
            E->control.ala_element_vanka_regularization);
    for (i=E->mesh.levmin;i<=E->mesh.levmax;i++)
      fprintf(E->fp,
              "Velocity MG level=%d global_nodes=%dx%dx%d "
              "local_nodes=%dx%dx%d local_elements=%dx%dx%d\n",
              i, E->mesh.NOX[i], E->mesh.NOY[i], E->mesh.NOZ[i],
              E->lmesh.NOX[i], E->lmesh.NOY[i], E->lmesh.NOZ[i],
              E->lmesh.ELX[i], E->lmesh.ELY[i], E->lmesh.ELZ[i]);
    fflush(E->fp);
  }

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
  void strict_ala_schur_diagnostic();
  void strict_ala_stage_B_diagnostic();
  void remove_rigid_rot();

  float vmag;

  double Udot_mag, dUdot_mag,omega[3];
  double vc1_operator_seconds=0.0,vc1_solver_seconds=0.0,vc1_clock;
  unsigned long long vc1_physical_hash_before=0,vc1_physical_hash_after=0;
  unsigned long long vc1_local_physical_hash_before=0,vc1_local_physical_hash_after=0;
  int vc1_local_physical_mismatch=0,vc1_global_physical_mismatch=0;
  int m,count,i,j,k;
  const char *vc1_stage=getenv("STRICT_ALA_VC1_STAGE");
  const char *vc1_warm_input=getenv("STRICT_ALA_VC1_WARM_INPUT");
  const char *vc1_warm_output=getenv("STRICT_ALA_VC1_WARM_OUTPUT");
  const char *vc1_compare_left=getenv("STRICT_ALA_VC1_COMPARE_LEFT");
  const char *vc1_compare_right=getenv("STRICT_ALA_VC1_COMPARE_RIGHT");

  double *oldU[NCS], *delta_U[NCS];

  const int nno = E->lmesh.nno;
  const int nel = E->lmesh.nel;
  const int nnov = E->lmesh.nnov;
  const int neq = E->lmesh.neq;
  const int vpts = vpoints[E->mesh.nsd];
  const int dims = E->mesh.nsd;
  const int addi_dof = additional_dof[dims];

  //velocities_conform_bcs(E,E->U);

  if(vc1_stage!=NULL && E->monitor.solution_cycles==0 &&
     vc1_compare_left!=NULL && vc1_compare_right!=NULL) {
    strict_ala_vc1_compare_states(E,vc1_compare_left,vc1_compare_right);
    MPI_Barrier(E->parallel.world);
    MPI_Finalize();
    exit(EXIT_SUCCESS);
  }
  if(vc1_stage!=NULL && E->monitor.solution_cycles==0 &&
     vc1_warm_input!=NULL && vc1_warm_input[0]!='\0')
    strict_ala_vc1_transfer_warm_state(E,vc1_warm_input,0);
  if(vc1_stage!=NULL && E->monitor.solution_cycles==0)
    vc1_physical_hash_before=strict_ala_vc1_physical_state_hash(
        E,&vc1_local_physical_hash_before);

  assemble_forces(E,0);

  if(E->monitor.solution_cycles==0 || E->viscosity.update_allowed) {
    if(E->control.ala_schur_diagnostic &&
       E->monitor.solution_cycles==0)
      strict_ala_schur_diagnostic(E);
    else {
      vc1_clock=MPI_Wtime();
      get_system_viscosity(E,1,E->EVI[E->mesh.levmax],E->VI[E->mesh.levmax]);
      velocities_conform_bcs(E,E->U);
      construct_stiffness_B_matrix(E);
      vc1_operator_seconds+=MPI_Wtime()-vc1_clock;
    }

    /* assemble_forces() includes nonzero velocity-boundary terms formed
       with the element stiffness.  A viscosity update changes that
       stiffness, so the force assembled above is stale for K_gamma.  Refresh
       only the active AL strict-ALA path; gamma=0 and legacy paths retain
       their historical call sequence. */
    if(E->control.ala_pressure_buoyancy &&
       E->control.ala_augmented_lagrangian_gamma > 0.0) {
      assemble_forces(E,0);
      if(E->parallel.me==0) {
        fprintf(E->fp,"ALA strict force reassembled after stiffness update "
                "cycle=%d gamma=%e\n",E->monitor.solution_cycles,
                E->control.ala_augmented_lagrangian_gamma);
        fflush(E->fp);
      }
    }
  }

  if(E->monitor.solution_cycles==0 && E->parallel.me==0) {
    fprintf(E->fp,"STRICT_ALA_STAGE_ABC_CONTROLS stage_B=%s stage_C=%s "
            "operator=B_Kgamma_inverse_BT gamma=%e outer_solver=%s "
            "outer_iterations=%d inner_max_cycles=%d "
            "continuity_tolerance=%e momentum_tolerance=%e\n",
            E->control.ala_stage_abc_adjoint_diagnostic ? "on" : "off",
            E->control.ala_stage_abc_production_logging ? "on" : "off",
            E->control.ala_augmented_lagrangian_gamma,
            E->control.ala_outer_solver,E->control.p_iterations,
            E->control.v_steps_low,E->control.tole_comp,
            E->control.ala_unaugmented_momentum_tolerance);
    fflush(E->fp);
  }

  if(vc1_stage!=NULL && E->monitor.solution_cycles==0) {
    strict_ala_vc1_viscosity_audit(E);
    if(E->parallel.me==0) {
      fprintf(E->fp,"STRICT_ALA_VC1_FROZEN_STATE_GUARD before=%016llx "
              "solver_mutable_scope=U_P_only phase_state=derived_from_T_and_cfg\n",
              vc1_physical_hash_before);
      fprintf(E->fp,"STRICT_ALA_VC1_TIMING operator_rebuild_seconds=%.17e "
              "fgmres_and_preconditioner_seconds=0.00000000000000000e+00\n",
              vc1_operator_seconds);
      fflush(E->fp);
    }
  }

  if(E->control.ala_stage_abc_adjoint_diagnostic &&
     E->monitor.solution_cycles==0)
    strict_ala_stage_B_diagnostic(E);

  vc1_clock=MPI_Wtime();
  solve_constrained_flow_iterative(E);
  vc1_solver_seconds+=MPI_Wtime()-vc1_clock;

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
      
      vc1_clock=MPI_Wtime();
      get_system_viscosity(E,0,E->EVI[E->mesh.levmax],E->VI[E->mesh.levmax]);
      velocities_conform_bcs(E,E->U);
      construct_stiffness_B_matrix(E);
      vc1_operator_seconds+=MPI_Wtime()-vc1_clock;
      vc1_clock=MPI_Wtime();
      solve_constrained_flow_iterative(E);
      vc1_solver_seconds+=MPI_Wtime()-vc1_clock;
      
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

  if(vc1_stage!=NULL && E->monitor.solution_cycles==0) {
    vc1_physical_hash_after=strict_ala_vc1_physical_state_hash(
        E,&vc1_local_physical_hash_after);
    vc1_local_physical_mismatch=
        vc1_local_physical_hash_after!=vc1_local_physical_hash_before;
    MPI_Allreduce(&vc1_local_physical_mismatch,&vc1_global_physical_mismatch,
                  1,MPI_INT,MPI_MAX,E->parallel.world);
    if(vc1_global_physical_mismatch)
      myerror(E,"VC1 frozen T/C/tracer/time state changed during Stokes solve");
    strict_ala_vc1_viscosity_audit(E);
    if(E->parallel.me==0) {
      fprintf(E->fp,"STRICT_ALA_VC1_FROZEN_STATE before=%016llx after=%016llx "
              "bitwise_equal=true phase_state=derived_from_frozen_T_and_cfg\n",
              vc1_physical_hash_before,vc1_physical_hash_after);
      fprintf(E->fp,"STRICT_ALA_VC1_TIMING operator_rebuild_seconds=%.17e "
              "fgmres_and_preconditioner_seconds=%.17e\n",
              vc1_operator_seconds,vc1_solver_seconds);
      fflush(E->fp);
    }
    if(vc1_warm_output!=NULL && vc1_warm_output[0]!='\0')
      strict_ala_vc1_transfer_warm_state(E,vc1_warm_output,1);
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
    void assemble_unaugmented_del2_u();
    void get_elt_f();
    void get_elt_tr();

    const int lev = E->mesh.levmax;
    const int neq = E->lmesh.neq;
    const unsigned int flags[4] = {0, VBX, VBY, VBZ};

    ku = allocate_equation_field(E);
    grad_p = allocate_equation_field(E);
    external_force = allocate_equation_field(E);
    plate_velocity = allocate_equation_field(E);

    assemble_unaugmented_del2_u(E, E->U, ku, lev, 0);
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
    if(isfinite(value.maximum) && isfinite(value.minimum))
        fprintf(fp, "%-20s  %+16.8e  %+16.8e  %+16.8e\n",
                name, value.total, value.maximum, value.minimum);
    else
        fprintf(fp, "%-20s  %+16.8e  %16s  %16s\n",
                name, value.total, "N/A", "N/A");
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
                        * (E->phase_B[phase_index][m][i]
                           - phase_change_reference_fraction(
                               E, phase_index, m, i));
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
    int m, e, local_ez;
    int local_linf_cap, local_linf_element, local_linf_global_ez;
    int global_linf_location[3];
    double volume, rdiv, rbeta, rala, relative, denominator;
    double scale_density, characteristic_scale, active_scale_cutoff;
    double global_cancellation_l2, elementwise_relative_l2;
    double active_relative_l2, active_volume_fraction;
    double active_fraction_above_threshold;
    double local_sum[6], global_sum[5], global_linf;
    double local_linf_depth, global_linf_depth;
    double local_linf_signed, global_linf_signed;
    double local_active[3], global_active[3];
    double **div_u, **beta_u;
    struct { double value; int rank; } local_maxloc, global_maxloc;
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
    local_maxloc.value = -1.0;
    local_maxloc.rank = E->parallel.me;
    local_linf_cap = 0;
    local_linf_element = 0;
    local_linf_global_ez = 0;
    local_linf_depth = 0.0;
    local_linf_signed = 0.0;
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
            if(fabs(rala) > local_maxloc.value) {
                local_ez = (e-1) % E->lmesh.ELZ[E->mesh.levmax] + 1;
                local_maxloc.value = fabs(rala);
                local_linf_signed = rala;
                local_linf_cap = m;
                local_linf_element = e;
                local_linf_global_ez =
                    E->lmesh.EZS[E->mesh.levmax] + local_ez;
                local_linf_depth = (1.0 - 0.5 *
                    (E->sx[m][3][local_ez] +
                     E->sx[m][3][local_ez+1])) * E->data.radius_km;
            }
        }
    }
    MPI_Allreduce(local_sum, global_sum, 5, MPI_DOUBLE, MPI_SUM,
                  E->parallel.world);
    MPI_Allreduce(&local_maxloc, &global_maxloc, 1, MPI_DOUBLE_INT,
                  MPI_MAXLOC, E->parallel.world);
    global_linf = global_maxloc.value;
    if(E->parallel.me == global_maxloc.rank) {
        global_linf_location[0] = local_linf_cap;
        global_linf_location[1] = local_linf_element;
        global_linf_location[2] = local_linf_global_ez;
        global_linf_depth = local_linf_depth;
        global_linf_signed = local_linf_signed;
    }
    MPI_Bcast(global_linf_location, 3, MPI_INT, global_maxloc.rank,
              E->parallel.world);
    MPI_Bcast(&global_linf_depth, 1, MPI_DOUBLE, global_maxloc.rank,
              E->parallel.world);
    MPI_Bcast(&global_linf_signed, 1, MPI_DOUBLE, global_maxloc.rank,
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
        fprintf(E->fp, "%-30s  %+16.8e\n", "Linf_signed",
                global_linf_signed);
        fprintf(E->fp, "%-30s  %16d\n", "Linf_owner_rank",
                global_maxloc.rank);
        fprintf(E->fp, "%-30s  %16d\n", "Linf_owner_cap",
                global_linf_location[0]);
        fprintf(E->fp, "%-30s  %16d\n", "Linf_owner_element",
                global_linf_location[1]);
        fprintf(E->fp, "%-30s  %16d\n", "Linf_global_radial_element",
                global_linf_location[2]);
        fprintf(E->fp, "%-30s  %+16.8e\n", "Linf_depth_km",
                global_linf_depth);
        fprintf(E->fp, "%-30s  %+16.8e\n",
                "volume_weighted_mean_abs",
                global_sum[1] / global_sum[3]);
        if(strcmp(E->control.ala_outer_solver,"coupled_fgmres") == 0) {
            fprintf(E->fp, "%-30s  %16s\n", "solver_div_v", "N/A");
            fprintf(E->fp, "%-30s  %16s\n", "solver_div_v_status",
                    "not_applicable");
        }
        else {
            fprintf(E->fp, "%-30s  %+16.8e\n", "solver_div_v",
                    E->monitor.incompressibility);
            fprintf(E->fp, "%-30s  %16s\n", "solver_div_v_status",
                    "legacy_populated");
        }
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


static void write_ala_unaugmented_momentum_residual(struct All_variables *E)
{
    int m,i;
    double force_norm,residual_norm,rms,relative;
    double **ku,**grad_p,**force,**residual;
    const int lev=E->mesh.levmax;
    const int neq=E->lmesh.neq;
    const int gneq=E->mesh.neq;
    void assemble_forces();
    void assemble_grad_p();
    void assemble_unaugmented_del2_u();
    void strip_bcs_from_residual();

    /* The strict-ALA force contains the C^T p pressure-buoyancy term;
       assemble_grad_p supplies D^T p.  Their combination audits
       f-(D+C)^T p-Ku using the original K. */
    assemble_forces(E,0);
    ku=allocate_equation_field(E);
    grad_p=allocate_equation_field(E);
    force=allocate_equation_field(E);
    residual=allocate_equation_field(E);
    assemble_unaugmented_del2_u(E,E->U,ku,lev,1);
    assemble_grad_p(E,E->P,grad_p,lev);

    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        for(i=0;i<neq;i++) {
            force[m][i]=E->F[m][i];
            residual[m][i]=E->F[m][i]-grad_p[m][i]-ku[m][i];
        }
        strip_bcs_from_residual(E,force,lev);
        strip_bcs_from_residual(E,residual,lev);
    }
    force_norm=sqrt(global_vdot(E,force,force,lev));
    residual_norm=sqrt(global_vdot(E,residual,residual,lev));
    rms=residual_norm/sqrt((double)gneq);
    relative=residual_norm/(force_norm+1.0e-32);

    if(E->parallel.me==0) {
        fprintf(stderr,"ALA unaugmented momentum residual: gamma=%e "
                "rms=%e relative=%e\n",
                E->control.ala_augmented_lagrangian_gamma,rms,relative);
        fprintf(E->fp,"ALA unaugmented momentum residual: gamma=%e "
                "rms=%e relative=%e\n",
                E->control.ala_augmented_lagrangian_gamma,rms,relative);
        fflush(E->fp);
    }

    free_equation_field(E,ku);
    free_equation_field(E,grad_p);
    free_equation_field(E,force);
    free_equation_field(E,residual);
}


static void write_stokes_diagnostics(struct All_variables *E)
{
    if(!E->control.ala_pressure_buoyancy)
        return;
    write_ala_unaugmented_momentum_residual(E);
    write_mechanical_power(E);
    write_ala_residual(E);
}
