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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"

void myerror(struct All_variables *,char *);

#ifdef _UNICOS
#include <fortran.h>
#endif

int epsilon[4][4] = {   /* Levi-Cita epsilon */
  {0, 0, 0, 0},
  {0, 1,-1, 1},
  {0,-1, 1,-1},
  {0, 1,-1, 1} };


/*=====================================================================
  Variable dimension matrix allocation  function from numerical recipes
  Note: ANSII consistency requires some additional features !
  =====================================================================  */

double **dmatrix(nrl,nrh,ncl,nch)
     int nrl,nrh,ncl,nch;
{
  int i,nrow = nrh-nrl+1,ncol=nch-ncl+1;
  double **m;

  /* allocate pointer to rows  */
  m=(double **) malloc((nrow+1)* sizeof(double *));
  m+=1;
  m-= nrl;

  /*  allocate rows and set the pointers accordingly   */
  m[nrl] = (double *) malloc((nrow*ncol+1)* sizeof(double));
  m[nrl] += 1;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++)
     m[i] = m[i-1] + ncol;

  return(m);		}


float **fmatrix(nrl,nrh,ncl,nch)
     int nrl,nrh,ncl,nch;
{
  int i,nrow = nrh-nrl+1,ncol=nch-ncl+1;
  float **m;

  /* allocate pointer to rows  */
  m=(float **) malloc((unsigned)((nrow+1)* sizeof(float *)));
  m+=1;
  m-= nrl;

  /*  allocate rows and set the pointers accordingly   */
  m[nrl] = (float *) malloc((unsigned)((nrow*ncol+1)* sizeof(float)));
  m[nrl] += 1;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++)
     m[i] = m[i-1] + ncol;

  return(m);		}


void dfree_matrix(m,nrl,nrh,ncl,nch)
     double **m;
     int nrl,nrh,ncl,nch;
{
  int i;
  for(i=nrh;i>=nrl;i--)
    free((void *)(m[i] + ncl));
  free((void *) (m+nrl));
  return;
}

void ffree_matrix(m,nrl,nrh,ncl,nch)
     float **m;
     int nrl,nrh,ncl,nch;
{
  int i;
  for(i=nrh;i>=nrl;i--)
    free((void *)(m[i] + ncl));
  free((void *) (m+nrl));
  return;
}

/*=============================================================
  Functions to allocate/remove space for variable sized vector.
  =============================================================  */

double *dvector(nl,nh)
     int nl,nh;
{
  double *v;
  v=(double *) malloc((unsigned) ( nh - nl +1)* sizeof(double));
  return( v-nl );  }

float *fvector(nl,nh)
     int nl,nh;
{
  float *v;
  v=(float *) malloc((unsigned) ( nh - nl +1)* sizeof(float));
  return( v-nl );  }

void dfree_vector(v,nl,nh)
     double *v;
     int nl,nh;
{
  free((char*) (v+nl));	}

void ffree_vector(v,nl,nh)
     float *v;
     int nl,nh;
{
  free((char*) (v+nl));	}

int *sivector(nl,nh)
     int nl,nh;
{
  int *v;
  v=(int*) malloc((unsigned)(nh-nl +1) * sizeof(int));
  return (v-nl);
}

void sifree_vector(v,nl,nh)
     int *v;
     int nl,nh;
{ free((char *) (v+nl));    }



void dvcopy(E,A,B,a,b)
     struct All_variables *E;
     double **A,**B;
     int a,b;

{   int i,m;

  for (m=1;m<=E->sphere.caps_per_proc;m++)
    for(i=a;i<=b;i++)
      A[m][i] = B[m][i];

    return; }

void vcopy(A,B,a,b)
     float *A,*B;
     int a,b;

{   int i;

    for(i=a;i<=b;i++)
      A[i] = B[i];

    return; }



/*  ===========================================================
    Iterative solver also using multigrid  ........
    ===========================================================  */

static int solve_del2_u_internal(struct All_variables *E, double **d0,
                                 double **F, double acc, int high_lev,
                                 int max_cycles, int progress_interval)
{
  void assemble_del2_u();
  void e_assemble_del2_u();
  void n_assemble_del2_u();
  void strip_bcs_from_residual();
  void gauss_seidel();

  double conj_grad();
  double multi_grid();
  double global_vdot();
  void record();
  void report();

  int count,counts,cycles,convergent,valid;
  int i, neq, gneq,m;

  char message[200];

  double CPU_time0(),initial_time,time;
  double residual,prior_residual,r0,elapsed,requested_relative;
  double *D1[NCS], *r[NCS], *Au[NCS];

  gneq  = E->mesh.NEQ[high_lev];
  gneq  = E->mesh.neq;
  neq  = E->lmesh.NEQ[high_lev];

  for (m=1;m<=E->sphere.caps_per_proc;m++)
      for(i=0;i<neq;i++)  {
	d0[m][i] = 0.0;
      }

  r0=residual=sqrt(global_vdot(E,F,F,high_lev)/gneq);

  prior_residual=2*residual;
  count = 0;
  initial_time=CPU_time0();

  if (!(E->control.NMULTIGRID || E->control.EMULTIGRID)) {
    cycles = E->control.v_steps_low;
    if(max_cycles>0)
      cycles=min(cycles,max_cycles);
    residual = conj_grad(E,d0,F,acc,&cycles,high_lev);
    valid = (residual < acc)? 1:0;
  }

  else  {

    counts =0;
    if(E->parallel.me==0){	/* output */
      snprintf(message,200,"resi = %.6e for iter %d acc %.6e",residual,counts,acc);
      record(E,message);
      report(E,message);
    }

    do {
      residual=multi_grid(E,d0,F,acc,high_lev);
      valid = (residual < acc)?1:0;
      counts ++;
      if(E->parallel.me==0){	/* output  */
	snprintf(message,200,"resi = %.6e for iter %d acc %.6e",residual,counts,acc);
	record(E,message);
	report(E,message);
      }
      if(max_cycles>0 && progress_interval>0 &&
         (counts%progress_interval)==0 && E->parallel.me==0) {
        fprintf(E->fp,"ALA COUPLED INNER VELOCITY progress cycles=%d "
                "residual=%e target=%e relative=%e\n",counts,residual,acc,
                residual/max(r0,1.0e-300));
        fprintf(stderr,"ALA COUPLED INNER VELOCITY progress cycles=%d "
                "residual=%e target=%e relative=%e\n",counts,residual,acc,
                residual/max(r0,1.0e-300));
        fflush(E->fp);
        fflush(stderr);
      }
    }  while (!valid && (max_cycles<=0 || counts<max_cycles));

    cycles = counts;
  }


  /* Convergence check .....
     We should give it a chance to recover if it briefly diverges initially, and
     don't worry about slower convergence if it is close to the answer   */

  if((count > 0) &&
     (residual > r0*2.0)  ||
     (fabs(residual-prior_residual) < acc*0.1 && (residual > acc * 10.0))   )
    convergent=0;
  else {
    convergent=1;
    prior_residual=residual;
  }

  if(E->control.print_convergence&&E->parallel.me==0)   {
    fprintf(E->fp,"%s residual (%03d) = %.3e from %.3e to %.3e in %5.2f secs \n",
	    (convergent ? " * ":"!!!"),cycles,residual,r0,acc,CPU_time0()-initial_time);
    fflush(E->fp);
  }

  elapsed=CPU_time0()-initial_time;
  if(max_cycles>0 && E->control.ala_stage_abc_production_logging) {
    /* Stage C treats reaching the limit as exhaustion even when the last
     * cycle happens to cross the target.  This makes the max-cycle contract
     * explicit and leaves one cycle of headroom for a valid solve. */
    valid=(isfinite(residual) && residual<=acc && cycles<max_cycles);
    E->control.ala_stage_abc_inner_call_count++;
    E->control.ala_stage_abc_inner_cycle_count+=cycles;
    E->control.ala_stage_abc_inner_seconds+=elapsed;
    requested_relative=acc/max(r0,1.0e-300);
    if(E->parallel.me==0) {
      FILE *inner_fp;
      char inner_path[512];
      int write_header;
      snprintf(inner_path,sizeof(inner_path),
               "%s.strict_ala_stage_C_inner_solves.csv",
               E->control.data_file);
      inner_fp=fopen(inner_path,"a+");
      if(inner_fp==NULL)
        myerror(E,"Unable to open strict-ALA Stage-C inner-solve log");
      fseek(inner_fp,0,SEEK_END);
      write_header=(ftell(inner_fp)==0);
      if(write_header)
        fprintf(inner_fp,"case,call_id,outer_iteration,role,rhs_rms,requested_relative_tolerance,target_absolute,achieved_absolute,achieved_relative,cycles,max_cycles,seconds,status\n");
      fprintf(inner_fp,
              "%s,%lld,%d,%s,%.17e,%.17e,%.17e,%.17e,%.17e,%d,%d,%.17e,%s\n",
              getenv("STRICT_ALA_CASE") ? getenv("STRICT_ALA_CASE") : "UNKNOWN",
              E->control.ala_stage_abc_inner_call_count,
              E->control.ala_stage_abc_outer_iteration,
              E->control.ala_stage_abc_inner_role,r0,
              requested_relative,acc,residual,
              residual/max(r0,1.0e-300),cycles,max_cycles,elapsed,
              valid ? "CONVERGED" : "INVALID");
      fclose(inner_fp);
    }
  }

  if(max_cycles>0 && E->parallel.me==0) {
    fprintf(E->fp,"ALA COUPLED INNER VELOCITY summary status=%s "
            "cycles=%d max_cycles=%d residual=%e initial=%e target=%e "
            "relative=%e seconds=%e\n",valid ? "converged" : "cycle_limit",
            cycles,max_cycles,residual,r0,acc,
            residual/max(r0,1.0e-300),elapsed);
    fprintf(stderr,"ALA COUPLED INNER VELOCITY summary status=%s "
            "cycles=%d max_cycles=%d residual=%e target=%e\n",
            valid ? "converged" : "cycle_limit",cycles,max_cycles,
            residual,acc);
    fflush(E->fp);
    fflush(stderr);
  }

  count++;


  E->control.total_iteration_cycles += count;
  E->control.total_v_solver_calls += 1;

  if(max_cycles>0 && E->control.ala_stage_abc_production_logging && !valid)
    myerror(E,"Strict-ALA Stage-C K_gamma inverse contract failed");

  return(valid);
}


/* Preserve the historical unbounded velocity solve for every legacy caller.
 * The coupled strict-ALA prototype uses the bounded entry point below so a
 * difficult Arnoldi block cannot remain forever inside one MG application. */
int solve_del2_u(E,d0,F,acc,high_lev)
     struct All_variables *E;
     double **d0;
     double **F;
     double acc;
     int high_lev;
{
  return solve_del2_u_internal(E,d0,F,acc,high_lev,0,0);
}


int solve_del2_u_bounded(struct All_variables *E, double **d0, double **F,
                         double acc, int high_lev, int max_cycles,
                         int progress_interval)
{
  if(max_cycles<1 || progress_interval<1)
    myerror(E,"Invalid coupled strict-ALA velocity MG limits");
  return solve_del2_u_internal(E,d0,F,acc,high_lev,max_cycles,
                               progress_interval);
}

/* =================================
   recursive multigrid function ....
   ================================= */

/* Stage F2b observes the existing velocity multigrid without changing its
 * arithmetic.  The caller enables selected cycles through a process-local
 * environment gate; every norm remains a collective over the normal world
 * communicator. */
static int ala_f2b_observer_enabled(void)
{
  const char *value=getenv("STRICT_ALA_STAGE_F2B_OBSERVE");
  return value!=NULL && strcmp(value,"1")==0;
}

static int ala_f2c_active(void)
{
  const char *value=getenv("STRICT_ALA_STAGE_F2C_ACTIVE");
  return value!=NULL && strcmp(value,"1")==0;
}

static int ala_f2c_finest_sweeps(struct All_variables *E,int configured,
                                int level)
{
  const char *value;
  char *end;
  long result;
  if(!ala_f2c_active() || level!=E->mesh.levmax) return configured;
  value=getenv("STRICT_ALA_STAGE_F2C_FINEST_SWEEPS");
  if(value==NULL) myerror(E,"Incomplete strict-ALA Stage-F2c sweep contract");
  result=strtol(value,&end,10);
  if(*value=='\0' || *end!='\0' || result<1 || result>1000)
    myerror(E,"Invalid strict-ALA Stage-F2c finest sweep count");
  return (int)result;
}

static double ala_f2c_finest_damping(struct All_variables *E,double configured,
                                     int level)
{
  const char *value;
  char *end;
  double result;
  if(!ala_f2c_active() || level!=E->mesh.levmax) return configured;
  value=getenv("STRICT_ALA_STAGE_F2C_FINEST_DAMPING");
  if(value==NULL) myerror(E,"Incomplete strict-ALA Stage-F2c damping contract");
  result=strtod(value,&end);
  if(*value=='\0' || *end!='\0' || !isfinite(result) || result<=0.0 || result>1.0)
    myerror(E,"Invalid strict-ALA Stage-F2c finest damping");
  return result;
}

/* Stage F2d changes the composition of the existing local K_gamma blocks,
 * never their factors or the assembled operator.  Eight global-parity colors
 * make same-color element corrections additive and introduce a true residual
 * update between colors.  Mode 2 reverses the color order on alternate
 * sweeps, testing order sensitivity without a second cache. */
static int ala_f2d_finest_coloring(struct All_variables *E,int level)
{
  const char *active=getenv("STRICT_ALA_STAGE_F2D_ACTIVE");
  const char *mode;
  if(active==NULL || strcmp(active,"1")!=0 || level!=E->mesh.levmax)
    return 0;
  mode=getenv("STRICT_ALA_STAGE_F2D_L4_MODE");
  if(mode==NULL)
    myerror(E,"Incomplete strict-ALA Stage-F2d coloring contract");
  if(strcmp(mode,"configured")==0) return 0;
  if(strcmp(mode,"colored_forward")==0) return 1;
  if(strcmp(mode,"colored_alternating")==0) return 2;
  myerror(E,"Invalid strict-ALA Stage-F2d L4 mode");
  return 0;
}

/* Stage F2e replaces only the finest-level diagnostic smoother block.  A
 * candidate joins two face-neighbouring Q1 elements, giving twelve unique
 * nodes and a 36x36 K_gamma patch.  Exactly one orientation is cached at a
 * time, so the experiment never retains the three candidate caches
 * simultaneously. */
#define ALA_F2E_FACE_NODES 12
#define ALA_F2E_FACE_DOF 36
#define ALA_F2E_FACE_CHOL_SIZE 666

struct ala_f2e_face_cache {
  int orientation;
  int level;
  int patch_count[NCS];
  int *equations[NCS];
  higher_precision *chol[NCS];
  double *overlap_inverse[NCS];
};

static struct ala_f2e_face_cache ala_f2e_cache;

static int ala_f2e_finest_orientation(struct All_variables *E,int level)
{
  const char *active=getenv("STRICT_ALA_STAGE_F2E_ACTIVE");
  const char *mode;
  if(active==NULL || strcmp(active,"1")!=0 || level!=E->mesh.levmax)
    return 0;
  mode=getenv("STRICT_ALA_STAGE_F2E_L4_MODE");
  if(mode==NULL)
    myerror(E,"Incomplete strict-ALA Stage-F2e face-patch contract");
  if(strcmp(mode,"configured")==0) return 0;
  if(strcmp(mode,"face_z")==0) return 1;
  if(strcmp(mode,"face_x")==0) return 2;
  if(strcmp(mode,"face_y")==0) return 3;
  myerror(E,"Invalid strict-ALA Stage-F2e face-patch mode");
  return 0;
}

static const char *ala_f2e_orientation_name(int orientation)
{
  if(orientation==1) return "face_z";
  if(orientation==2) return "face_x";
  if(orientation==3) return "face_y";
  return "invalid";
}

static void ala_f2e_free_face_cache(void)
{
  int m;
  for(m=0;m<NCS;m++) {
    if(ala_f2e_cache.equations[m]) free(ala_f2e_cache.equations[m]);
    if(ala_f2e_cache.chol[m]) free(ala_f2e_cache.chol[m]);
    if(ala_f2e_cache.overlap_inverse[m])
      free(ala_f2e_cache.overlap_inverse[m]);
  }
  memset(&ala_f2e_cache,0,sizeof(ala_f2e_cache));
}

static int ala_f2e_element_number(struct All_variables *E,int level,
                                  int x,int y,int z)
{
  return z+(x-1)*E->lmesh.ELZ[level]
      +(y-1)*E->lmesh.ELZ[level]*E->lmesh.ELX[level];
}

static void ala_f2e_patch_elements(struct All_variables *E,int level,
                                   int orientation,int patch,
                                   int *first,int *second)
{
  int nx=E->lmesh.ELX[level],ny=E->lmesh.ELY[level];
  int nz=E->lmesh.ELZ[level],x,y,z,q;
  if(patch<0 || patch>=E->lmesh.NEL[level]/2)
    myerror(E,"Stage-F2e face-patch index is invalid");
  if(orientation==1) {
    z=2*(patch%(nz/2))+1;
    q=patch/(nz/2);
    x=q%nx+1;
    y=q/nx+1;
  }
  else if(orientation==2) {
    z=patch%nz+1;
    q=patch/nz;
    x=2*(q%(nx/2))+1;
    y=q/(nx/2)+1;
  }
  else {
    z=patch%nz+1;
    q=patch/nz;
    x=q%nx+1;
    y=2*(q/nx)+1;
  }
  if(x<1 || x>nx || y<1 || y>ny || z<1 || z>nz)
    myerror(E,"Stage-F2e face-patch coordinate is invalid");
  *first=ala_f2e_element_number(E,level,x,y,z);
  if(orientation==1) z++;
  else if(orientation==2) x++;
  else y++;
  *second=ala_f2e_element_number(E,level,x,y,z);
}

static void ala_f2e_build_face_cache(struct All_variables *E,int level,
                                     int orientation)
{
  int m,p,e,which,a,da,i,j,k,node,node_count,match,eq;
  int elements[2],nodes[ALA_F2E_FACE_NODES],fixed[ALA_F2E_FACE_DOF];
  int map[ALA_VANKA_DOF],patches,paired_dimension,neq;
  double element_matrix[ALA_VANKA_DOF*ALA_VANKA_DOF];
  double matrix[ALA_F2E_FACE_DOF*ALA_F2E_FACE_DOF];
  double L[ALA_F2E_FACE_DOF*ALA_F2E_FACE_DOF];
  double sum,pivot,maxdiag,global_diag,external,shift;
  double local_mb,global_mb,start,seconds,seconds_max;
  double local_min_ratio=1.0e300,global_min_ratio;
  double *overlap;
  higher_precision *factor;
  FILE *fp;
  char path[1024];
  const char *output_dir;
  void get_elt_k(),get_ala_aug_k();
  double CPU_time0();
  const int dims=E->mesh.nsd;
  const int ends=enodes[dims];

  if(dims!=3 || loc_mat_size[dims]!=ALA_VANKA_DOF)
    myerror(E,"Stage-F2e requires 3-D Q1 velocity elements");
  if(orientation<1 || orientation>3 || level!=E->mesh.levmax)
    myerror(E,"Stage-F2e face-cache request is invalid");
  paired_dimension=orientation==1 ? E->lmesh.ELZ[level] :
                   (orientation==2 ? E->lmesh.ELX[level] :
                                      E->lmesh.ELY[level]);
  if(paired_dimension<2 || (paired_dimension&1))
    myerror(E,"Stage-F2e local paired element dimension must be positive and even");
  patches=E->lmesh.NEL[level]/2;
  neq=E->lmesh.NEQ[level];
  if(patches<=0 || 2*patches!=E->lmesh.NEL[level])
    myerror(E,"Stage-F2e face-patch coverage is incomplete");

  ala_f2e_free_face_cache();
  ala_f2e_cache.orientation=orientation;
  ala_f2e_cache.level=level;
  start=CPU_time0();
  local_mb=0.0;
  for(m=1;m<=E->sphere.caps_per_proc;m++) {
    ala_f2e_cache.patch_count[m]=patches;
    ala_f2e_cache.equations[m]=(int *)malloc(
        patches*ALA_F2E_FACE_DOF*sizeof(int));
    ala_f2e_cache.chol[m]=(higher_precision *)malloc(
        patches*ALA_F2E_FACE_CHOL_SIZE*sizeof(higher_precision));
    ala_f2e_cache.overlap_inverse[m]=(double *)malloc(neq*sizeof(double));
    if(!ala_f2e_cache.equations[m] || !ala_f2e_cache.chol[m] ||
       !ala_f2e_cache.overlap_inverse[m])
      myerror(E,"Unable to allocate Stage-F2e face-patch cache");
    local_mb += (patches*(ALA_F2E_FACE_DOF*sizeof(int)
        +ALA_F2E_FACE_CHOL_SIZE*sizeof(higher_precision))
        +neq*sizeof(double))/(1024.0*1024.0);
    overlap=ala_f2e_cache.overlap_inverse[m];
    for(eq=0;eq<neq;eq++) overlap[eq]=0.0;
    for(p=0;p<patches;p++) {
      ala_f2e_patch_elements(E,level,orientation,p,&elements[0],&elements[1]);
      node_count=0;
      for(which=0;which<2;which++) {
        e=elements[which];
        for(a=1;a<=ends;a++) {
          node=E->IEN[level][m][e].node[a];
          match=-1;
          for(i=0;i<node_count;i++) if(nodes[i]==node) { match=i; break; }
          if(match<0) {
            if(node_count>=ALA_F2E_FACE_NODES)
              myerror(E,"Stage-F2e face patch has too many nodes");
            nodes[node_count]=node;
            match=node_count++;
          }
          for(da=0;da<dims;da++) map[(a-1)*dims+da]=match*dims+da;
        }
        get_elt_k(E,e,element_matrix,level,m,0);
        get_ala_aug_k(E,e,element_matrix,level,m);
        if(which==0)
          for(i=0;i<ALA_F2E_FACE_DOF*ALA_F2E_FACE_DOF;i++) matrix[i]=0.0;
        for(i=0;i<ALA_VANKA_DOF;i++)
          for(j=0;j<ALA_VANKA_DOF;j++)
            matrix[map[i]*ALA_F2E_FACE_DOF+map[j]] +=
                element_matrix[i*ALA_VANKA_DOF+j];
      }
      if(node_count!=ALA_F2E_FACE_NODES)
        myerror(E,"Stage-F2e face patch does not have twelve unique nodes");
      maxdiag=0.0;
      for(a=0;a<ALA_F2E_FACE_NODES;a++)
        for(da=0;da<dims;da++) {
          i=a*dims+da;
          node=nodes[a];
          eq=E->ID[level][m][node].doff[da+1];
          fixed[i]=(da==0 && (E->NODE[level][m][node]&VBX)) ||
                   (da==1 && (E->NODE[level][m][node]&VBY)) ||
                   (da==2 && (E->NODE[level][m][node]&VBZ));
          ala_f2e_cache.equations[m][p*ALA_F2E_FACE_DOF+i]=fixed[i] ? -1 : eq;
          if(!fixed[i]) {
            if(eq<0 || eq>=neq || !isfinite(E->BI[level][m][eq]) ||
               E->BI[level][m][eq]<=0.0)
              myerror(E,"Stage-F2e assembled diagonal is invalid");
            global_diag=1.0/E->BI[level][m][eq];
            external=global_diag-matrix[i*ALA_F2E_FACE_DOF+i];
            if(external < -1.0e-8*max(global_diag,1.0e-300))
              myerror(E,"Stage-F2e face-patch external diagonal is negative");
            matrix[i*ALA_F2E_FACE_DOF+i] +=
              E->control.ala_element_vanka_external_diagonal_weight*
              max(external,0.0);
            maxdiag=max(maxdiag,matrix[i*ALA_F2E_FACE_DOF+i]);
            overlap[eq] += 1.0;
          }
        }
      if(!isfinite(maxdiag) || maxdiag<=0.0)
        myerror(E,"Stage-F2e face-patch diagonal is invalid");
      shift=E->control.ala_element_vanka_regularization*maxdiag;
      for(i=0;i<ALA_F2E_FACE_DOF;i++) {
        if(fixed[i]) {
          for(j=0;j<ALA_F2E_FACE_DOF;j++)
            matrix[i*ALA_F2E_FACE_DOF+j]=
              matrix[j*ALA_F2E_FACE_DOF+i]=0.0;
          matrix[i*ALA_F2E_FACE_DOF+i]=maxdiag;
        }
        else {
          matrix[i*ALA_F2E_FACE_DOF+i] += shift;
          for(j=0;j<i;j++) {
            sum=0.5*(matrix[i*ALA_F2E_FACE_DOF+j]
                    +matrix[j*ALA_F2E_FACE_DOF+i]);
            matrix[i*ALA_F2E_FACE_DOF+j]=
              matrix[j*ALA_F2E_FACE_DOF+i]=sum;
          }
        }
      }
      for(i=0;i<ALA_F2E_FACE_DOF*ALA_F2E_FACE_DOF;i++) L[i]=0.0;
      for(i=0;i<ALA_F2E_FACE_DOF;i++)
        for(j=0;j<=i;j++) {
          sum=matrix[i*ALA_F2E_FACE_DOF+j];
          for(k=0;k<j;k++)
            sum-=L[i*ALA_F2E_FACE_DOF+k]*L[j*ALA_F2E_FACE_DOF+k];
          if(i==j) {
            pivot=sum;
            if(!isfinite(pivot) || pivot<=1.0e-14*maxdiag)
              myerror(E,"Stage-F2e face-patch Cholesky failed");
            L[i*ALA_F2E_FACE_DOF+j]=sqrt(pivot);
            local_min_ratio=min(local_min_ratio,pivot/maxdiag);
          }
          else L[i*ALA_F2E_FACE_DOF+j]=sum/L[j*ALA_F2E_FACE_DOF+j];
        }
      factor=ala_f2e_cache.chol[m]+p*ALA_F2E_FACE_CHOL_SIZE;
      k=0;
      for(i=0;i<ALA_F2E_FACE_DOF;i++)
        for(j=0;j<=i;j++) factor[k++]=(higher_precision)L[i*ALA_F2E_FACE_DOF+j];
    }
  }
  (E->solver.exchange_id_d)(E,ala_f2e_cache.overlap_inverse,level);
  for(m=1;m<=E->sphere.caps_per_proc;m++)
    for(eq=0;eq<neq;eq++) {
      sum=ala_f2e_cache.overlap_inverse[m][eq];
      if(!isfinite(sum) || sum<=0.0)
        myerror(E,"Stage-F2e face-patch overlap is incomplete");
      ala_f2e_cache.overlap_inverse[m][eq]=1.0/sum;
    }
  seconds=CPU_time0()-start;
  MPI_Allreduce(&seconds,&seconds_max,1,MPI_DOUBLE,MPI_MAX,E->parallel.world);
  MPI_Allreduce(&local_mb,&global_mb,1,MPI_DOUBLE,MPI_MAX,E->parallel.world);
  MPI_Allreduce(&local_min_ratio,&global_min_ratio,1,MPI_DOUBLE,MPI_MIN,
                E->parallel.world);
  if(E->parallel.me==0) {
    fprintf(E->fp,"STRICT_ALA_STAGE_F2E_FACE_CACHE orientation=%s "
        "patches_per_cap=%d dof=36 cache_mb_max=%.17e "
        "build_seconds_max=%.17e min_pivot_ratio=%.17e\n",
        ala_f2e_orientation_name(orientation),patches,global_mb,seconds_max,
        global_min_ratio);
    fflush(E->fp);
    output_dir=getenv("STRICT_ALA_STAGE_F2E_OUTPUT_DIR");
    if(output_dir==NULL)
      myerror(E,"Missing Stage-F2e output directory");
    snprintf(path,sizeof(path),"%s/strict_ala_stage_F2e_cache.csv",output_dir);
    fp=fopen(path,"a");
    if(fp==NULL) myerror(E,"Unable to append Stage-F2e cache output");
    fprintf(fp,"%s,%d,36,%.17e,%.17e,%.17e,1\n",
        ala_f2e_orientation_name(orientation),patches,global_mb,seconds_max,
        global_min_ratio);
    fclose(fp);
  }
}

static void ala_f2b_write_level_row(struct All_variables *E,
                                    const char *phase,int top_level,int level,
                                    double input2,double output2,double alpha)
{
  const char *directory=getenv("STRICT_ALA_STAGE_F2B_OUTPUT_DIR");
  const char *mode=getenv("STRICT_ALA_STAGE_F2B_MODE");
  const char *weight=getenv("STRICT_ALA_STAGE_F2B_WEIGHT");
  const char *cycle=getenv("STRICT_ALA_STAGE_F2B_CYCLE");
  const char *candidate=getenv("STRICT_ALA_STAGE_F2E_CANDIDATE");
  int f2e=candidate!=NULL;
  int f2d;
  double input_rms,output_rms,ratio;
  char path[1024];
  FILE *fp;

  if(directory==NULL || mode==NULL || weight==NULL || cycle==NULL)
    myerror(E,"Incomplete strict-ALA Stage-F2b observer environment");
  if(!isfinite(input2) || !isfinite(output2) || input2<0.0 || output2<0.0)
    myerror(E,"Invalid strict-ALA Stage-F2b level norm");
  input_rms=sqrt(input2/max((double)E->mesh.NEQ[level],1.0));
  output_rms=sqrt(output2/max((double)E->mesh.NEQ[level],1.0));
  ratio=output_rms/max(input_rms,1.0e-300);
  if(E->parallel.me==0) {
    if(candidate==NULL) candidate=getenv("STRICT_ALA_STAGE_F2D_CANDIDATE");
    f2d=candidate!=NULL && !f2e;
    if(candidate==NULL) candidate=getenv("STRICT_ALA_STAGE_F2C_CANDIDATE");
    snprintf(path,sizeof(path),"%s/%s",directory,f2e ?
             "strict_ala_stage_F2e_level_transfer.csv" : (f2d ?
             "strict_ala_stage_F2d_level_transfer.csv" :
             (candidate!=NULL ? "strict_ala_stage_F2c_level_transfer.csv" :
              "strict_ala_stage_F2b_level_transfer.csv")));
    fp=fopen(path,"a");
    if(fp==NULL)
      myerror(E,"Unable to append strict-ALA Stage-F2b level output");
    if(candidate!=NULL) fprintf(fp,"%s,",candidate);
    fprintf(fp,"%s,%s,%s,%d,%d,%s,%.17e,%.17e,%.17e,%.17e,1\n",
            mode,weight,cycle,top_level,level,phase,input_rms,output_rms,
            ratio,alpha);
    fclose(fp);
  }
}

static void ala_f2b_record_same_level(struct All_variables *E,
    const char *phase,int top_level,int level,double **input,double **action,
    double **scratch,double scale)
{
  int m,i;
  double input2,output2;
  double global_vdot();
  if(!ala_f2b_observer_enabled()) return;
  input2=global_vdot(E,input,input,level);
  for(m=1;m<=E->sphere.caps_per_proc;m++)
    for(i=0;i<E->lmesh.NEQ[level];i++)
      scratch[m][i]=input[m][i]-scale*action[m][i];
  output2=global_vdot(E,scratch,scratch,level);
  ala_f2b_write_level_row(E,phase,top_level,level,input2,output2,scale);
}

static void ala_f2b_record_restriction(struct All_variables *E,int top_level,
    int fine_level,double **fine,int coarse_level,double **coarse)
{
  double fine2,coarse2,fine_rms,coarse_rms;
  double global_vdot();
  const char *directory,*mode,*weight,*cycle;
  const char *candidate;
  int f2d,f2e;
  char path[1024];
  FILE *fp;
  if(!ala_f2b_observer_enabled()) return;
  fine2=global_vdot(E,fine,fine,fine_level);
  coarse2=global_vdot(E,coarse,coarse,coarse_level);
  if(!isfinite(fine2) || !isfinite(coarse2) || fine2<0.0 || coarse2<0.0)
    myerror(E,"Invalid strict-ALA Stage-F2b restriction norm");
  fine_rms=sqrt(fine2/max((double)E->mesh.NEQ[fine_level],1.0));
  coarse_rms=sqrt(coarse2/max((double)E->mesh.NEQ[coarse_level],1.0));
  directory=getenv("STRICT_ALA_STAGE_F2B_OUTPUT_DIR");
  mode=getenv("STRICT_ALA_STAGE_F2B_MODE");
  weight=getenv("STRICT_ALA_STAGE_F2B_WEIGHT");
  cycle=getenv("STRICT_ALA_STAGE_F2B_CYCLE");
  candidate=getenv("STRICT_ALA_STAGE_F2E_CANDIDATE");
  f2e=candidate!=NULL;
  if(candidate==NULL) candidate=getenv("STRICT_ALA_STAGE_F2D_CANDIDATE");
  f2d=candidate!=NULL && !f2e;
  if(candidate==NULL) candidate=getenv("STRICT_ALA_STAGE_F2C_CANDIDATE");
  if(directory==NULL || mode==NULL || weight==NULL || cycle==NULL)
    myerror(E,"Incomplete strict-ALA Stage-F2b restriction environment");
  if(E->parallel.me==0) {
    snprintf(path,sizeof(path),"%s/%s",directory,f2e ?
             "strict_ala_stage_F2e_level_transfer.csv" : (f2d ?
             "strict_ala_stage_F2d_level_transfer.csv" :
             (candidate!=NULL ? "strict_ala_stage_F2c_level_transfer.csv" :
              "strict_ala_stage_F2b_level_transfer.csv")));
    fp=fopen(path,"a");
    if(fp==NULL)
      myerror(E,"Unable to append strict-ALA Stage-F2b restriction output");
    if(candidate!=NULL) fprintf(fp,"%s,",candidate);
    fprintf(fp,"%s,%s,%s,%d,%d,RESTRICTION,%.17e,%.17e,%.17e,0,1\n",
            mode,weight,cycle,top_level,fine_level,fine_rms,coarse_rms,
            coarse_rms/max(fine_rms,1.0e-300));
    fclose(fp);
  }
}

double multi_grid(E,d1,F,acc,hl)
     struct All_variables *E;
     double **d1;
     double **F;
     double acc;
     int hl;  /* higher level of two */
{
    double residual,AudotAu,Audotres;
    void interp_vector();
    void project_vector();
    int m,i,j,Vn,Vnmax,cycles,f2b_observe;
    double alpha,beta;
    void gauss_seidel();
    void element_gauss_seidel();
    void e_assemble_del2_u();
    void strip_bcs_from_residual();
    void n_assemble_del2_u();

    double conj_grad(),global_vdot();

    FILE *fp;
    char filename[1000];
    int lev,ic,ulev,dlev;

    const int levmin = E->mesh.levmin;
    const int levmax = E->mesh.levmax;

    double time1,time,CPU_time0();
    double *res[MAX_LEVELS][NCS],*AU[MAX_LEVELS][NCS];
    double *vel[MAX_LEVELS][NCS],*del_vel[MAX_LEVELS][NCS];
    double *rhs[MAX_LEVELS][NCS],*fl[MAX_LEVELS][NCS];
				/* because it's recursive, need a copy at
				    each level */

    f2b_observe=ala_f2b_observer_enabled();

    for(i=E->mesh.levmin;i<=E->mesh.levmax;i++)
      for(m=1;m<=E->sphere.caps_per_proc;m++)    {
	del_vel[i][m]=(double *)malloc((E->lmesh.NEQ[i]+1)*sizeof(double));
	AU[i][m] = (double *)malloc((E->lmesh.NEQ[i]+1)*sizeof(double));
	vel[i][m]=(double *)malloc((E->lmesh.NEQ[i]+1)*sizeof(double));
	res[i][m]=(double *)malloc((E->lmesh.NEQ[i])*sizeof(double));
	if(f2b_observe)
	  rhs[i][m]=(double *)malloc((E->lmesh.NEQ[i]+1)*sizeof(double));
	if (i<E->mesh.levmax)
	  fl[i][m]=(double *)malloc((E->lmesh.NEQ[i])*sizeof(double));
      }

    Vnmax = E->control.mg_cycle;

        /* Project residual onto all the lower levels */

    project_vector(E,levmax,F,fl[levmax-1],1);
    strip_bcs_from_residual(E,fl[levmax-1],levmax-1);
    for(lev=levmax-1;lev>levmin;lev--) {
        project_vector(E,lev,fl[lev],fl[lev-1],1);
        strip_bcs_from_residual(E,fl[lev-1],lev-1);
       }

        /* Solve for the lowest level */

/*    time=CPU_time0(); */
    cycles = E->control.v_steps_low;

    gauss_seidel(E,vel[levmin],fl[levmin],AU[levmin],acc*0.01,&cycles,levmin,0);
    if(f2b_observe)
      ala_f2b_record_same_level(E,"INITIAL_BOTTOM_SOLVE",levmin,levmin,
          fl[levmin],AU[levmin],rhs[levmin],1.0);

    for(lev=levmin+1;lev<=levmax;lev++) {
      time=CPU_time0();

                         /* Utilize coarse solution and smooth at this level */
      interp_vector(E,lev-1,vel[lev-1],vel[lev]);
      strip_bcs_from_residual(E,vel[lev],lev);

      if (lev==levmax)
        for(m=1;m<=E->sphere.caps_per_proc;m++)
          for(j=0;j<E->lmesh.NEQ[lev];j++)
             res[lev][m][j]=F[m][j];
      else
        for(m=1;m<=E->sphere.caps_per_proc;m++)
          for(j=0;j<E->lmesh.NEQ[lev];j++)
             res[lev][m][j]=fl[lev][m][j];

      for(Vn=1;Vn<=Vnmax;Vn++)   {
                                        /*    Downward stoke of the V    */
        for (dlev=lev;dlev>=levmin+1;dlev--)   {

                                      /* Pre-smoothing  */
          cycles=((dlev==levmax)?E->control.v_steps_high:E->control.down_heavy);
          cycles=ala_f2c_finest_sweeps(E,cycles,dlev);
          ic = ((dlev==lev)?1:0);
          gauss_seidel(E,vel[dlev],res[dlev],AU[dlev],0.01,&cycles,dlev,ic);

          if(f2b_observe)
            ala_f2b_record_same_level(E,"DOWN_SMOOTH",lev,dlev,
                res[dlev],AU[dlev],rhs[dlev],1.0);

          for(m=1;m<=E->sphere.caps_per_proc;m++)
             for(i=0;i<E->lmesh.NEQ[dlev];i++)  {
               res[dlev][m][i]  = res[dlev][m][i] - AU[dlev][m][i];
	       }

          project_vector(E,dlev,res[dlev],res[dlev-1],1);
          strip_bcs_from_residual(E,res[dlev-1],dlev-1);
          if(f2b_observe)
            ala_f2b_record_restriction(E,lev,dlev,res[dlev],dlev-1,
                                       res[dlev-1]);
          }

                                        /*    Bottom of the V    */
       cycles = E->control.v_steps_low;
       gauss_seidel(E,vel[levmin],res[levmin],AU[levmin],acc*0.01,&cycles,levmin,0);
       if(f2b_observe)
         ala_f2b_record_same_level(E,"BOTTOM_SOLVE",lev,levmin,
             res[levmin],AU[levmin],rhs[levmin],1.0);
                                        /*    Upward stoke of the V    */
        for (ulev=levmin+1;ulev<=lev;ulev++)   {
            cycles=((ulev==levmax)?E->control.v_steps_high:E->control.up_heavy);
            cycles=ala_f2c_finest_sweeps(E,cycles,ulev);

            interp_vector(E,ulev-1,vel[ulev-1],del_vel[ulev]);
            strip_bcs_from_residual(E,del_vel[ulev],ulev);
            gauss_seidel(E,del_vel[ulev],res[ulev],AU[ulev],0.01,&cycles,ulev,1);

            AudotAu = global_vdot(E,AU[ulev],AU[ulev],ulev);
            if(!isfinite(AudotAu) || AudotAu<0.0)
              myerror(E,"Invalid multigrid upward action norm");
            if(AudotAu<=1.0e-300)
              alpha=0.0;
            else {
              Audotres=global_vdot(E,AU[ulev],res[ulev],ulev);
              if(!isfinite(Audotres))
                myerror(E,"Invalid multigrid upward action product");
              alpha=Audotres/AudotAu;
              if(!isfinite(alpha))
                myerror(E,"Invalid multigrid upward correction scale");
            }

            if(f2b_observe)
              ala_f2b_record_same_level(E,"UP_CORRECTION",lev,ulev,
                  res[ulev],AU[ulev],rhs[ulev],alpha);

            for(m=1;m<=E->sphere.caps_per_proc;m++)
              for(i=0;i<E->lmesh.NEQ[ulev];i++)   {
               vel[ulev][m][i] += alpha*del_vel[ulev][m][i];
               }

            if (ulev ==levmax)
              for(m=1;m<=E->sphere.caps_per_proc;m++)
                for(i=0;i<E->lmesh.NEQ[ulev];i++)   {
                  res[ulev][m][i] -= alpha*AU[ulev][m][i];
                  }

            }
        }
      }

   for(m=1;m<=E->sphere.caps_per_proc;m++)
     for(j=0;j<E->lmesh.NEQ[levmax];j++)   {
        F[m][j]=res[levmax][m][j];
        d1[m][j]+=vel[levmax][m][j];
        }

     residual = sqrt(global_vdot(E,F,F,hl)/E->mesh.NEQ[hl]);

      for(i=E->mesh.levmin;i<=E->mesh.levmax;i++)
        for(m=1;m<=E->sphere.caps_per_proc;m++) {
	  free((double*) del_vel[i][m]);
	  free((double*) AU[i][m]);
	  free((double*) vel[i][m]);
	  free((double*) res[i][m]);
	  if(f2b_observe) free((double*) rhs[i][m]);
	  if (i<E->mesh.levmax)
	    free((double*) fl[i][m]);
	  }


    return(residual);
}


/*  ===========================================================
    Conjugate gradient relaxation for the matrix equation Kd = f
    Returns the residual reduction after itn iterations ...
    ===========================================================  */


double conj_grad(E,d0,F,acc,cycles,level)
     struct All_variables *E;
     double **d0;
     double **F;
     double acc;
     int *cycles;
     int level;
{
    double *r0[NCS],*r1[NCS],*r2[NCS];
    double *z0[NCS],*z1[NCS],*z2[NCS];
    double *p1[NCS],*p2[NCS];
    double *Ap[NCS];
    double *shuffle[NCS];

    int m,count,i,steps;
    double residual;
    double alpha,beta,dotprod,dotr1z1,dotr0z0;

    double CPU_time0(),time;

    void parallel_process_termination();
    void assemble_del2_u();
    void strip_bcs_from_residual();
    double global_vdot();

    const int mem_lev=E->mesh.levmax;
    const int high_neq = E->lmesh.NEQ[level];

    steps = *cycles;

    for(m=1;m<=E->sphere.caps_per_proc;m++)    {
      r0[m] = (double *)malloc(E->lmesh.NEQ[mem_lev]*sizeof(double));
      r1[m] = (double *)malloc(E->lmesh.NEQ[mem_lev]*sizeof(double));
      r2[m] = (double *)malloc(E->lmesh.NEQ[mem_lev]*sizeof(double));
      z0[m] = (double *)malloc(E->lmesh.NEQ[mem_lev]*sizeof(double));
      z1[m] = (double *)malloc(E->lmesh.NEQ[mem_lev]*sizeof(double));
      p1[m] = (double *)malloc((1+E->lmesh.NEQ[mem_lev])*sizeof(double));
      p2[m] = (double *)malloc((1+E->lmesh.NEQ[mem_lev])*sizeof(double));
      Ap[m] = (double *)malloc((1+E->lmesh.NEQ[mem_lev])*sizeof(double));
    }

    for(m=1;m<=E->sphere.caps_per_proc;m++)
      for(i=0;i<high_neq;i++) {
	r1[m][i] = F[m][i];
        d0[m][i] = 0.0;
	}

    residual = sqrt(global_vdot(E,r1,r1,level)/E->mesh.NEQ[level]);

    assert(residual != 0.0  /* initial residual for CG = 0.0 */);
    count = 0;

    while (((residual > acc) && (count < steps)) || count == 0)  {

    for(m=1;m<=E->sphere.caps_per_proc;m++)
      for(i=0;i<high_neq;i++)
         z1[m][i] = E->BI[level][m][i] * r1[m][i];

	dotr1z1 = global_vdot(E,r1,z1,level);

	if (0==count)
          for(m=1;m<=E->sphere.caps_per_proc;m++)
	    for(i=0;i<high_neq;i++)
	      p2[m][i] = z1[m][i];
	else {
	    assert(dotr0z0 != 0.0 /* in head of conj_grad */);
	    beta = dotr1z1/dotr0z0;
            for(m=1;m<=E->sphere.caps_per_proc;m++)
	      for(i=0;i<high_neq;i++)
		p2[m][i] = z1[m][i] + beta * p1[m][i];
	}

    dotr0z0 = dotr1z1;

	assemble_del2_u(E,p2,Ap,level,1);

	dotprod=global_vdot(E,p2,Ap,level);

	if(0.0==dotprod)
	    alpha=1.0e-3;
	else
	    alpha = dotr1z1/dotprod;

        for(m=1;m<=E->sphere.caps_per_proc;m++)
          for(i=0;i<high_neq;i++) {
	    d0[m][i] += alpha * p2[m][i];
	    r2[m][i] = r1[m][i] - alpha * Ap[m][i];
	    }

	residual = sqrt(global_vdot(E,r2,r2,level)/E->mesh.NEQ[level]);

        for(m=1;m<=E->sphere.caps_per_proc;m++)    {
	  shuffle[m] = r0[m]; r0[m] = r1[m]; r1[m] = r2[m]; r2[m] = shuffle[m];
	  shuffle[m] = z0[m]; z0[m] = z1[m]; z1[m] = shuffle[m];
	  shuffle[m] = p1[m]; p1[m] = p2[m]; p2[m] = shuffle[m];
          }

	count++;
	/* end of while-loop */

    }

    *cycles=count;

    strip_bcs_from_residual(E,d0,level);

    for(m=1;m<=E->sphere.caps_per_proc;m++)    {
      free((double*) r0[m]);
      free((double*) r1[m]);
      free((double*) r2[m]);
      free((double*) z0[m]);
      free((double*) z1[m]);
      free((double*) p1[m]);
      free((double*) p2[m]);
      free((double*) Ap[m]);
    }

    return(residual);   }


/* ========================================================================================
   An element by element version of the gauss-seidel routine. Initially this is a test
   platform, we want to know if it handles discontinuities any better than the node/equation
   versions
   =========================================================================================*/

void element_gauss_seidel(E,d0,F,Ad,acc,cycles,level,guess)
    struct All_variables *E;
    double **d0;
    double **F,**Ad;
    double acc;
    int *cycles;
    int level;
    int guess;
{
    int count,i,j,k,l,m,ns,nc,d,steps,loc;
    int p1,p2,p3,q1,q2,q3;
    int e,eq,node,node1;
    int element,eqn1,eqn2,eqn3,eqn11,eqn12,eqn13;

    void e_assemble_del2_u();
    void n_assemble_del2_u();
    void strip_bcs_from_residual();

    double U1[24],AD1[24],F1[24];
    double w1,w2,w3;
    double w11,w12,w13;
    double w[24];

    double *dd[NCS],*elt_k;
    int *vis[NCS];

    const int dims=E->mesh.nsd;
    const int ends=enodes[dims];
    const int n=loc_mat_size[E->mesh.nsd];
    const int neq=E->lmesh.NEQ[level];
    const int nel=E->lmesh.NEL[level];
    const int nno=E->lmesh.NNO[level];


    steps=*cycles;

    for (m=1;m<=E->sphere.caps_per_proc;m++) {
      dd[m] = (double *)malloc(neq*sizeof(double));
      vis[m] = (int *)malloc((nno+1)*sizeof(int));
    }
    elt_k=(double *)malloc((24*24)*sizeof(double));

    if(guess){
	e_assemble_del2_u(E,d0,Ad,level,1);
    }
    else {
      for (m=1;m<=E->sphere.caps_per_proc;m++)
	for(i=0;i<neq;i++)
	     Ad[m][i]=d0[m][i]=0.0;
    }

   count=0;
   while (count <= steps)      {
     for (m=1;m<=E->sphere.caps_per_proc;m++) {
	for(i=1;i<=nno;i++)
	    vis[m][i]=0;

	for(e=1;e<=nel;e++) {

	    elt_k = E->elt_k[level][m][e].k;

	    for(i=1;i<=ends;i++)  {
		node=E->IEN[level][m][e].node[i];
		p1=(i-1)*dims;
		w[p1] = w[p1+1] = w[p1+2] = 1.0;
		if(E->NODE[level][m][node] & VBX)
		    w[p1] = 0.0;
		if(E->NODE[level][m][node] & VBY)
		    w[p1+1] = 0.0;
		if(E->NODE[level][m][node] & VBZ)
		    w[p1+2] = 0.0;

	        }


	    for(i=1;i<=ends;i++)   {
		node=E->IEN[level][m][e].node[i];
		if(!vis[m][node])
		    continue;

		eqn1=E->ID[level][m][node].doff[1];
		eqn2=E->ID[level][m][node].doff[2];
                eqn3=E->ID[level][m][node].doff[3];
		p1=(i-1)*dims*n;
		p2=p1+n;
		p3=p2+n;


		/* update Au */
		for(j=1;j<=ends;j++) {
		    node1=E->IEN[level][m][e].node[j];

	            eqn11=E->ID[level][m][node1].doff[1];
		    eqn12=E->ID[level][m][node1].doff[2];
		    eqn13=E->ID[level][m][node1].doff[3];
		    q1=(j-1)*3;

                    Ad[m][eqn11] += w[q1]*(elt_k[p1+q1] * dd[m][eqn1] + elt_k[p2+q1] * dd[m][eqn2] + elt_k[p3+q1] * dd[m][eqn3]);
                    Ad[m][eqn12] += w[q1+1]*(elt_k[p1+q1+1] * dd[m][eqn1] + elt_k[p2+q1+1] * dd[m][eqn2] + elt_k[p3+q1+1] * dd[m][eqn3]);
	            Ad[m][eqn13] += w[q1+2]*(elt_k[p1+q1+2] * dd[m][eqn1] + elt_k[p2+q1+2] * dd[m][eqn2] + elt_k[p3+q1+2] * dd[m][eqn3]);

		   }
	        }


	    for(i=1;i<=ends;i++)   {
		node=E->IEN[level][m][e].node[i];
		if(vis[m][node])
		    continue;

		eqn1=E->ID[level][m][node].doff[1];
		eqn2=E->ID[level][m][node].doff[2];
		eqn3=E->ID[level][m][node].doff[3];
		p1=(i-1)*dims*n;
		p2=p1+n;
		p3=p2+n;

		/* update dd, d0 */
		d0[m][eqn1] += (dd[m][eqn1] = w[(i-1)*dims]*(F[m][eqn1]-Ad[m][eqn1])*E->BI[level][m][eqn1]);
		d0[m][eqn2] += (dd[m][eqn2] = w[(i-1)*dims+1]*(F[m][eqn2]-Ad[m][eqn2])*E->BI[level][m][eqn2]);
		d0[m][eqn3] += (dd[m][eqn3] = w[(i-1)*dims+2]*(F[m][eqn3]-Ad[m][eqn3])*E->BI[level][m][eqn3]);

		vis[m][node]=1;

		/* update Au */
		for(j=1;j<=ends;j++) {
		   node1=E->IEN[level][m][e].node[j];

		   eqn11=E->ID[level][m][node1].doff[1];
		   eqn12=E->ID[level][m][node1].doff[2];
		   eqn13=E->ID[level][m][node1].doff[3];
		   q1=(j-1)*3;
		   q2=q1+1;
		   q3=q1+2;

		   Ad[m][eqn11] += w[q1]*(elt_k[p1+q1] * dd[m][eqn1] + elt_k[p2+q1] * dd[m][eqn2] + elt_k[p3+q1] * dd[m][eqn3]);
		   Ad[m][eqn12] += w[q2]*(elt_k[p1+q2] * dd[m][eqn1] + elt_k[p2+q2] * dd[m][eqn2] + elt_k[p3+q2] * dd[m][eqn3]);
		   Ad[m][eqn13] += w[q3]*(elt_k[p1+q3] * dd[m][eqn1] + elt_k[p2+q3] * dd[m][eqn2] + elt_k[p3+q3] * dd[m][eqn3]);
		   }
	       }

	    }          /* end for el  */
	 }               /* end for m */

        (E->solver.exchange_id_d)(E, Ad, level);
        (E->solver.exchange_id_d)(E, d0, level);

	/* completed cycle */

	    count++;

    }

    for (m=1;m<=E->sphere.caps_per_proc;m++) {
      free((double*) dd[m]);
      free((int*) vis[m]);
    }
    free((double*) elt_k);

    return;
}


static void ala_f2e_face_patch_smooth(struct All_variables *E,double **d0,
                                      double **F,double **Ad,int *cycles,
                                      int level,int guess,int orientation)
{
  int m,p,i,j,eq,count,steps,fixed;
  double rhs[ALA_F2E_FACE_DOF],forward[ALA_F2E_FACE_DOF];
  double solution[ALA_F2E_FACE_DOF],sum,diagonal,correction,damping;
  double *delta[NCS];
  higher_precision *chol;
  void n_assemble_del2_u();
  const int neq=E->lmesh.NEQ[level];

  if(ala_f2e_cache.orientation!=orientation ||
     ala_f2e_cache.level!=level)
    ala_f2e_build_face_cache(E,level,orientation);
  steps=*cycles;
  damping=ala_f2c_finest_damping(E,
      E->control.ala_element_vanka_damping,level);
  for(m=1;m<=E->sphere.caps_per_proc;m++) {
    delta[m]=(double *)malloc(neq*sizeof(double));
    if(delta[m]==NULL)
      myerror(E,"Unable to allocate Stage-F2e face-patch correction");
  }
  if(guess) n_assemble_del2_u(E,d0,Ad,level,1);
  else
    for(m=1;m<=E->sphere.caps_per_proc;m++)
      for(eq=0;eq<neq;eq++) d0[m][eq]=Ad[m][eq]=0.0;

  for(count=0;count<steps;count++) {
    for(m=1;m<=E->sphere.caps_per_proc;m++)
      for(eq=0;eq<neq;eq++) delta[m][eq]=0.0;
    for(m=1;m<=E->sphere.caps_per_proc;m++)
      for(p=0;p<ala_f2e_cache.patch_count[m];p++) {
        chol=ala_f2e_cache.chol[m]+p*ALA_F2E_FACE_CHOL_SIZE;
        for(i=0;i<ALA_F2E_FACE_DOF;i++) {
          eq=ala_f2e_cache.equations[m][p*ALA_F2E_FACE_DOF+i];
          fixed=eq<0 || eq>=neq;
          rhs[i]=fixed ? 0.0 : F[m][eq]-Ad[m][eq];
        }
        for(i=0;i<ALA_F2E_FACE_DOF;i++) {
          sum=rhs[i];
          for(j=0;j<i;j++) sum-=chol[i*(i+1)/2+j]*forward[j];
          diagonal=chol[i*(i+1)/2+i];
          if(!isfinite(diagonal) || diagonal<=0.0)
            myerror(E,"Stage-F2e face-patch factor diagonal is invalid");
          forward[i]=sum/diagonal;
        }
        for(i=ALA_F2E_FACE_DOF-1;i>=0;i--) {
          sum=forward[i];
          for(j=i+1;j<ALA_F2E_FACE_DOF;j++)
            sum-=chol[j*(j+1)/2+i]*solution[j];
          solution[i]=sum/chol[i*(i+1)/2+i];
          if(!isfinite(solution[i]))
            myerror(E,"Non-finite Stage-F2e face-patch solve");
        }
        for(i=0;i<ALA_F2E_FACE_DOF;i++) {
          eq=ala_f2e_cache.equations[m][p*ALA_F2E_FACE_DOF+i];
          if(eq>=0 && eq<neq) delta[m][eq]+=solution[i];
        }
      }
    (E->solver.exchange_id_d)(E,delta,level);
    for(m=1;m<=E->sphere.caps_per_proc;m++)
      for(eq=0;eq<neq;eq++) {
        correction=damping*ala_f2e_cache.overlap_inverse[m][eq]*delta[m][eq];
        if(!isfinite(correction))
          myerror(E,"Non-finite Stage-F2e face-patch correction");
        d0[m][eq]+=correction;
      }
    n_assemble_del2_u(E,d0,Ad,level,1);
  }
  for(m=1;m<=E->sphere.caps_per_proc;m++) free(delta[m]);
  *cycles=steps;
}

/* Apply an overlapping full-element correction for the strict-ALA augmented
   velocity operator.  Each cached factor contains the complete 24x24 local
   K_gamma block plus the assembled diagonal contribution from neighboring
   elements, so local rigid modes are controlled without discarding component
   coupling. */
void ala_element_vanka_smooth(struct All_variables *E, double **d0,
                              double **F, double **Ad, double acc,
                              int *cycles, int level, int guess)
{
    int m,e,a,da,ia,j,eq,count,steps,node,fixed;
    int coloring,colors,color,active_color,ex,ey,ez,element_color;
    int eqs[24];
    double rhs[24],forward[24],solution[24],sum,diagonal;
    double correction,damping,residual;
    double *delta[NCS];
    higher_precision *chol;
    static int announced[MAX_LEVELS];
    const int dims=E->mesh.nsd;
    const int ends=enodes[dims];
    const int n=loc_mat_size[dims];
    const int neq=E->lmesh.NEQ[level];
    const int nel=E->lmesh.NEL[level];
    void n_assemble_del2_u();
    void myerror();

    (void)acc;
    steps=*cycles;
    coloring=ala_f2d_finest_coloring(E,level);
    colors=coloring ? 8 : 1;
    damping=ala_f2c_finest_damping(E,
        E->control.ala_element_vanka_damping,level);
    if(E->parallel.me==0 && !announced[level]) {
        fprintf(E->fp,
                "ALA full element-Vanka active level=%d local_elements=%d "
                "sweeps=%d damping=%g\n",
                level,nel,steps,damping);
        fprintf(stderr,
                "ALA full element-Vanka active level=%d local_elements=%d "
                "sweeps=%d damping=%g\n",
                level,nel,steps,damping);
        fflush(E->fp);
        fflush(stderr);
        announced[level]=1;
    }
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        delta[m]=(double *)malloc(neq*sizeof(double));
        if(delta[m]==NULL)
            myerror(E,"Unable to allocate ALA element-Vanka correction");
    }

    if(guess)
        n_assemble_del2_u(E,d0,Ad,level,1);
    else
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(eq=0;eq<neq;eq++)
                d0[m][eq]=Ad[m][eq]=0.0;

    for(count=0;count<steps;count++)
      for(color=0;color<colors;color++) {
        active_color=(coloring==2 && (count&1)) ? 7-color : color;
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(eq=0;eq<neq;eq++)
                delta[m][eq]=0.0;

        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=1;e<=nel;e++) {
                if(coloring) {
                    ez=(e-1)%E->lmesh.ELZ[level];
                    ex=((e-1)/E->lmesh.ELZ[level])%
                        E->lmesh.ELX[level];
                    ey=(e-1)/(E->lmesh.ELZ[level]*E->lmesh.ELX[level]);
                    element_color=((E->lmesh.EXS[level]+ex)&1)
                        +2*((E->lmesh.EYS[level]+ey)&1)
                        +4*((E->lmesh.EZS[level]+ez)&1);
                    if(element_color!=active_color) continue;
                }
                if(!E->ALA_vanka_valid[level][m][e])
                    myerror(E,"ALA element-Vanka factor cache is invalid");
                chol=E->ALA_vanka_chol[level][m]
                    +e*ALA_VANKA_CHOL_SIZE;
                for(a=1;a<=ends;a++) {
                    node=E->IEN[level][m][e].node[a];
                    for(da=0;da<dims;da++) {
                        ia=(a-1)*dims+da;
                        eq=E->ID[level][m][node].doff[da+1];
                        eqs[ia]=eq;
                        fixed=(da==0 && (E->NODE[level][m][node]&VBX)) ||
                              (da==1 && (E->NODE[level][m][node]&VBY)) ||
                              (da==2 && (E->NODE[level][m][node]&VBZ));
                        residual=F[m][eq]-Ad[m][eq];
                        rhs[ia]=fixed ? 0.0 : residual;
                    }
                }
                for(ia=0;ia<n;ia++) {
                    sum=rhs[ia];
                    for(j=0;j<ia;j++)
                        sum -= chol[ia*(ia+1)/2+j]*forward[j];
                    diagonal=chol[ia*(ia+1)/2+ia];
                    if(!isfinite(diagonal) || diagonal<=0.0)
                        myerror(E,"ALA element-Vanka factor diagonal is invalid");
                    forward[ia]=sum/diagonal;
                }
                for(ia=n-1;ia>=0;ia--) {
                    sum=forward[ia];
                    for(j=ia+1;j<n;j++)
                        sum -= chol[j*(j+1)/2+ia]*solution[j];
                    diagonal=chol[ia*(ia+1)/2+ia];
                    solution[ia]=sum/diagonal;
                    if(!isfinite(solution[ia])) {
                        fprintf(stderr,"rank=%d level=%d sweep=%d cap=%d "
                                "element=%d local_dof=%d rhs=%e forward=%e "
                                "diagonal=%e\n",E->parallel.me,level,count,m,e,
                                ia,rhs[ia],forward[ia],diagonal);
                        myerror(E,"Non-finite full element-Vanka block solve");
                    }
                }
                for(ia=0;ia<n;ia++)
                    delta[m][eqs[ia]] += solution[ia];
            }

        (E->solver.exchange_id_d)(E,delta,level);
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(eq=0;eq<neq;eq++) {
                correction=damping
                    *E->ALA_vanka_overlap_BI[level][m][eq]*delta[m][eq];
                if(!isfinite(correction)) {
                    fprintf(stderr,"rank=%d level=%d sweep=%d cap=%d eq=%d "
                            "overlap_inverse=%e delta=%e damping=%e\n",
                            E->parallel.me,level,count,m,eq,
                            E->ALA_vanka_overlap_BI[level][m][eq],
                            delta[m][eq],damping);
                    myerror(E,"Non-finite full element-Vanka correction");
                }
                d0[m][eq] += correction;
            }
        n_assemble_del2_u(E,d0,Ad,level,1);
      }

    for(m=1;m<=E->sphere.caps_per_proc;m++)
        free(delta[m]);
    *cycles=steps;
}


/* ============================================================================
   Multigrid Gauss-Seidel relaxation scheme which requires the storage of local
   information, otherwise some other method is required. NOTE this is a bit worse
   than real gauss-seidel because it relaxes all the equations for a node at one
   time (Jacobi at a node). It does the job though.
   ============================================================================ */

void gauss_seidel(E,d0,F,Ad,acc,cycles,level,guess)
     struct All_variables *E;
     double **d0;
     double **F,**Ad;
     double acc;
     int *cycles;
     int level;
     int guess;
{
    int count,i,j,k,l,m,ns,steps,f2e_orientation;
    int *C;
    int eqn1,eqn2,eqn3;

    void parallel_process_termination();
    void n_assemble_del2_u();

    double U1,U2,U3,UU;
    double sor,residual,global_vdot();

    higher_precision *B1,*B2,*B3;


    const int dims=E->mesh.nsd;
    const int ends=enodes[dims];
    const int n=loc_mat_size[E->mesh.nsd];
    const int neq=E->lmesh.NEQ[level];
    const int num_nodes=E->lmesh.NNO[level];
    const int nox=E->lmesh.NOX[level];
    const int noz=E->lmesh.NOY[level];
    const int noy=E->lmesh.NOZ[level];
    const int gneq  = E->mesh.neq;
    const int max_eqn=14*dims;

    const double zeroo = 0.0;

    /* A logical sweep count is identical on every rank.  It is intentionally
       not process-summed in Stage-E output. */
    if(E->control.ala_stage_e_diagnostic)
      E->control.ala_stage_e_velocity_smoother_sweeps += *cycles;

    /* Vanka is a multigrid smoother, not the coarse-grid solver.  In
       particular, v_steps_low can be O(1000), which would turn the additive
       element update and its halo exchange into thousands of global sweeps. */
    if(E->control.ala_element_vanka_smoother &&
       E->control.ala_augmented_lagrangian_gamma>0.0 &&
       level>E->mesh.levmin) {
      f2e_orientation=ala_f2e_finest_orientation(E,level);
      if(f2e_orientation) {
        ala_f2e_face_patch_smooth(E,d0,F,Ad,cycles,level,guess,
                                  f2e_orientation);
        return;
      }
      ala_element_vanka_smooth(E,d0,F,Ad,acc,cycles,level,guess);
      return;
    }

    steps=*cycles;
    sor = 1.3;

    if(guess) {
      n_assemble_del2_u(E,d0,Ad,level,1);
    }
    else
      for (m=1;m<=E->sphere.caps_per_proc;m++)
	for(i=0;i<neq;i++) {
	    d0[m][i]=Ad[m][i]=zeroo;
	}

    count = 0;


    while (count < steps) {
      for (m=1;m<=E->sphere.caps_per_proc;m++)
 	for(j=0;j<=E->lmesh.NEQ[level];j++)
          E->temp[m][j] = zeroo;

      for (m=1;m<=E->sphere.caps_per_proc;m++)
          Ad[m][neq] = zeroo;

      for (m=1;m<=E->sphere.caps_per_proc;m++)
 	for(i=1;i<=E->lmesh.NNO[level];i++)
          if(E->NODE[level][m][i] & OFFSIDE)   {

	    eqn1=E->ID[level][m][i].doff[1];
	    eqn2=E->ID[level][m][i].doff[2];
	    eqn3=E->ID[level][m][i].doff[3];
	    E->temp[m][eqn1] = (F[m][eqn1] - Ad[m][eqn1])*E->BI[level][m][eqn1];
	    E->temp[m][eqn2] = (F[m][eqn2] - Ad[m][eqn2])*E->BI[level][m][eqn2];
	    E->temp[m][eqn3] = (F[m][eqn3] - Ad[m][eqn3])*E->BI[level][m][eqn3];
	    E->temp1[m][eqn1] = Ad[m][eqn1];
	    E->temp1[m][eqn2] = Ad[m][eqn2];
	    E->temp1[m][eqn3] = Ad[m][eqn3];
            }

      for (m=1;m<=E->sphere.caps_per_proc;m++)
 	for(i=1;i<=E->lmesh.NNO[level];i++)     {

	    eqn1=E->ID[level][m][i].doff[1];
	    eqn2=E->ID[level][m][i].doff[2];
	    eqn3=E->ID[level][m][i].doff[3];
            C=E->Node_map[level][m]+(i-1)*max_eqn;
	    B1=E->Eqn_k1[level][m]+(i-1)*max_eqn;
	    B2=E->Eqn_k2[level][m]+(i-1)*max_eqn;
 	    B3=E->Eqn_k3[level][m]+(i-1)*max_eqn;

                 /* Ad on boundaries differs after the following operation, but
                  no communications are needed yet, because boundary Ad will
                  not be used for the G-S iterations for interior nodes */

            for(j=3;j<max_eqn;j++)  {
                 UU = E->temp[m][C[j]];
                 Ad[m][eqn1] += B1[j]*UU;
                 Ad[m][eqn2] += B2[j]*UU;
                 Ad[m][eqn3] += B3[j]*UU;
                 }

            if (!(E->NODE[level][m][i]&OFFSIDE))   {
               E->temp[m][eqn1] = (F[m][eqn1] - Ad[m][eqn1])*E->BI[level][m][eqn1];
               E->temp[m][eqn2] = (F[m][eqn2] - Ad[m][eqn2])*E->BI[level][m][eqn2];
               E->temp[m][eqn3] = (F[m][eqn3] - Ad[m][eqn3])*E->BI[level][m][eqn3];
	       }

                 /* Ad on boundaries differs after the following operation */
	    for(j=0;j<max_eqn;j++)
		    Ad[m][C[j]]  += B1[j]*E->temp[m][eqn1]
                                 +  B2[j]*E->temp[m][eqn2]
                                 +  B3[j]*E->temp[m][eqn3];

	    d0[m][eqn1] += E->temp[m][eqn1];
	    d0[m][eqn2] += E->temp[m][eqn2];
	    d0[m][eqn3] += E->temp[m][eqn3];
  	    }

      for (m=1;m<=E->sphere.caps_per_proc;m++)
 	for(i=1;i<=E->lmesh.NNO[level];i++)
          if(E->NODE[level][m][i] & OFFSIDE)   {
	    eqn1=E->ID[level][m][i].doff[1];
	    eqn2=E->ID[level][m][i].doff[2];
	    eqn3=E->ID[level][m][i].doff[3];
	    Ad[m][eqn1] -= E->temp1[m][eqn1];
	    Ad[m][eqn2] -= E->temp1[m][eqn2];
	    Ad[m][eqn3] -= E->temp1[m][eqn3];
	    }

      (E->solver.exchange_id_d)(E, Ad, level);

      for (m=1;m<=E->sphere.caps_per_proc;m++)
 	for(i=1;i<=E->lmesh.NNO[level];i++)
          if(E->NODE[level][m][i] & OFFSIDE)   {
	    eqn1=E->ID[level][m][i].doff[1];
	    eqn2=E->ID[level][m][i].doff[2];
	    eqn3=E->ID[level][m][i].doff[3];
	    Ad[m][eqn1] += E->temp1[m][eqn1];
	    Ad[m][eqn2] += E->temp1[m][eqn2];
	    Ad[m][eqn3] += E->temp1[m][eqn3];
	    }


	count++;

/*     for (m=1;m<=E->sphere.caps_per_proc;m++)
	for(i=0;i<neq;i++)          {
	   F[m][i] -= Ad[m][i];
	   Ad[m][i] = 0.0;
	  }
*/
      }

    *cycles=count;
    return;

}

void print_elt_k(E,a)
     struct All_variables *E;
     double a[24*24];

{ int l,ll,n;

  printf("elt k is ...\n");


  n = loc_mat_size[E->mesh.nsd];

  for(l=0;l<n;l++)
    { fprintf(stderr,"\n");fflush(stderr);
      for(ll=0;ll<n;ll++)
	{ fprintf(stderr,"%s%.3e ",a[ll*n+l] >= 0.0 ? "+" : "",a[ll*n+l]);
	  fflush(stderr);
	}
    }
  fprintf(stderr,"\n"); fflush(stderr);

  return; }


double cofactor(A,i,j,n)
     double A[4][4];
     int i,j,n;

{ int k,l,p,q;
  double determinant();
  double **dmatrix();

  double B[4][4]; /* because of recursive behaviour of det/cofac, need to use
			       new copy of B at each 'n' level of this routine */

  if (n>3) printf("Error, no cofactors for matrix more than 3x3\n");

  p=q=1;

  for(k=1;k<=n;k++)    {
     if(k==i) continue;
     for(l=1;l<=n;l++)      {
	   if (l==j) continue;
           B[p][q]=A[k][l];
	   q++ ;
	   }
     q=1;p++;
     }


  return(epsilon[i][j]*determinant(B,n-1));


}


/* Fast (conditional) determinant for 3x3 or 2x2 ... otherwise calls general routine */

double determinant(A,n)
     double A[4][4];
     int n;

{ double gen_determinant();

  switch (n)
    { case 1:
	return(A[1][1]);
      case 2:
	return(A[1][1]*A[2][2]-A[1][2]*A[2][1]);
      case 3:
	return(A[1][1]*(A[2][2]*A[3][3]-A[2][3]*A[3][2])-
	       A[1][2]*(A[2][1]*A[3][3]-A[2][3]*A[3][1])+
	       A[1][3]*(A[2][1]*A[3][2]-A[2][2]*A[3][1]));
      default:
	return(1);
/*	return(gen_determinant(A,n)); */
      }
}


/* recursive function to determine matrix determinant */

double gen_determinant(A,n)
     double **A;
     int n;

{ double det;
  double cofactor();

  int i;


  if(n==1) return(A[1][1]); /* need a way to break the recursion */

  det=0.0;
  for(i=1;i<=n;i++)
    det += A[1][i]*cofactor(A,1,i,n);

  return(det);
}


 float area_of_4node(x1,y1,x2,y2,x3,y3,x4,y4)
 float x1,y1,x2,y2,x3,y3,x4,y4;

 {
 float area;

 area = fabs(0.5*(x1*(y2-y4)+x2*(y4-y1)+x4*(y1-y2)))
      + fabs(0.5*(x2*(y3-y4)+x3*(y4-y2)+x4*(y2-y3)));

 return area;
 }

/* =====================================*/
 double sphere_h(l,m,t,f,ic)
 int l,m,ic;
 double t,f;
 {

 double plgndr_a(),sphere_hamonics;

 sphere_hamonics = 0.0;
 if (ic==0)
    sphere_hamonics = cos(m*f)*plgndr_a(l,m,t);
 else if (m)
    sphere_hamonics = sin(m*f)*plgndr_a(l,m,t);

 return sphere_hamonics;
 }

/* =====================================*/
 double plgndr_a(l,m,t)
 int l,m;
 double t;
 {

  int i,ll;
  double x,fact,pll,pmm,pmmp1,somx2,plgndr;
  const double two=2.0;
  const double one=1.0;

  x = cos(t);
  pmm=one;
  if(m>0) {
    somx2=sqrt((one-x)*(one+x));
    fact = one;
    for (i=1;i<=m;i++)   {
      pmm = -pmm*fact*somx2;
      fact = fact + two;
      }
    }

  if (l==m)
     plgndr = pmm;
  else  {
     pmmp1 = x*(2*m+1)*pmm;
     if(l==m+1)
       plgndr = pmmp1;
     else   {
       for (ll=m+2;ll<=l;ll++)  {
         pll = (x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
         pmm = pmmp1;
         pmmp1 = pll;
         }
       plgndr = pll;
       }
     }

 return plgndr;
 }

/* =====================================
 =====================================*/
 double modified_plgndr_a(l,m,t)
 int l,m;
 double t;
 {

  int i,ll;
  double x,fact1,fact2,fact,pll,pmm,pmmp1,somx2,plgndr;
  const double three=3.0;
  const double two=2.0;
  const double one=1.0;

  x = cos(t);
  pmm=one;
  if(m>0) {
    somx2=sqrt((one-x)*(one+x));
    fact1= three;
    fact2= two;
    for (i=1;i<=m;i++)   {
      fact=sqrt(fact1/fact2);
      pmm = -pmm*fact*somx2;
      fact1+=  two;
      fact2+=  two;
      }
    }

  if (l==m)
     plgndr = pmm;
  else  {
     pmmp1 = x*sqrt(two*m+three)*pmm;
     if(l==m+1)
       plgndr = pmmp1;
     else   {
       for (ll=m+2;ll<=l;ll++)  {
	 fact1= sqrt((4.0*ll*ll-one)*(double)(ll-m)/(double)(ll+m));
	 fact2= sqrt((2.0*ll+one)*(ll-m)*(ll+m-one)*(ll-m-one)
		     /(double)((two*ll-three)*(ll+m)));
         pll = ( x*fact1*pmmp1-fact2*pmm)/(ll-m);
         pmm = pmmp1;
         pmmp1 = pll;
         }
       plgndr = pll;
       }
     }

 plgndr /= sqrt(4.0*M_PI);

 if (m!=0) plgndr *= sqrt(two);

 return plgndr;
 }

 /* ===================================  */
  double sqrt_multis(jj,ii)
  int ii,jj;
 {
  int i;
  double sqrt_multisa;

  sqrt_multisa = 1.0;
  if(jj>ii)
    for (i=jj;i>ii;i--)
      sqrt_multisa *= 1.0/sqrt((double)i);

  return sqrt_multisa;
  }

 /* ===================================  */
  double multis(ii)
  int ii;
 {
  int i;
  double multisa;

  multisa = 1.0;
  if (ii)
    for (i=2;i<=ii;i++)
      multisa *= (double)i;

  return multisa;
  }


 /* ===================================  */
 int int_multis(ii)
 int ii;
 {
 int i,multisa;

 multisa = 1;
 if (ii)
   for (i=2;i<=ii;i++)
     multisa *= i;

 return multisa;
 }
