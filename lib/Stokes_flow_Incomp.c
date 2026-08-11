/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


/*   Functions which solve for the velocity and pressure fields using Uzawa-type iteration loop.  */

#include <math.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "ala_block_vector.h"
#include "ala_coupled_operator.h"
#include <stdlib.h>
#include <string.h>

void assemble_grad_rho_p(struct All_variables *,double **,double **,int);
int solve_del2_u_bounded(struct All_variables *,double **,double **,double,
                         int,int,int);

void myerror(struct All_variables *,char *);
double assemble_Ahatp_jacobi_entry(struct All_variables *,int,int,int,int);

static float solve_Ahat_p_fhat(struct All_variables *E,
                               double **V, double **P, double **F,
                               double imp, int *steps_max);
static float solve_Ahat_p_fhat_CG(struct All_variables *E,
                                  double **V, double **P, double **F,
                                  double imp, int *steps_max);
static float solve_Ahat_p_fhat_BiCG(struct All_variables *E,
                                    double **V, double **P, double **F,
                                    double imp, int *steps_max);
static float solve_Ahat_p_fhat_ALA_PCG(struct All_variables *E,
                                       double **V, double **P, double **F,
                                       double imp, int *steps_max);
struct ala_pressure_preconditioner_cache;
static void strict_ala_continuity_metrics(struct All_variables *E,
                                          double **V, double **r,
                                          double **div_u, int lev,
                                          double *mass_norm,
                                          double *cancellation_l2);
static double strict_ala_continuity_term_strength(struct All_variables *E,
                                                  double **r,
                                                  double **div_u, int lev);
static void strict_ala_momentum_residual_audit(struct All_variables *E,
                                               double **V, double **P,
                                               double **grad_work,
                                               double **residual_work,
                                               int lev, double *rms,
                                               double *relative);
static void strict_ala_momentum_decomposition_audit(
    struct All_variables *E, double **V, double **P,
    double **work_a, double **work_b, int lev, const char *stage,
    int iteration, double *raw_rms, double *raw_relative);
static double strict_ala_inner_accuracy(struct All_variables *E,
                                        double **F, int lev,
                                        double relative_accuracy);
static void apply_ala_pressure_preconditioner(struct All_variables *E,
                                              double **r, double **z,
                                              double **work, int lev,
                                              int iteration,
                                              struct ala_pressure_preconditioner_cache
                                              *cache);
static void strict_ala_depth_diagnostics(struct All_variables *E,
                                         double **r, double **div_u,
                                         int lev, int iteration);
static void strict_ala_beta_causal_diagnostics(
    struct All_variables *E, double **V, double **active_r,
    double **alternate_r, double **div_u, int lev, int iteration);
static void strict_ala_coarse_residual_diagnostics(
    struct All_variables *E, double **r, int lev, int iteration);


static void assemble_ala_coupled_residual(
    struct All_variables *E, double **V, double **P,
    const struct ala_block_vector *rhs, struct ala_block_vector *residual,
    struct ala_block_vector *action, struct ala_block_vector *operator_work,
    int lev)
{
    int m,i,e,neq,npno;

    neq=E->lmesh.NEQ[lev];
    npno=E->lmesh.NPNO[lev];
    apply_ala_coupled_operator(E,V,P,action->velocity,action->pressure,
                               operator_work->velocity,lev);
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        for(i=0;i<neq;i++)
            residual->velocity[m][i]=rhs->velocity[m][i]
                                      -action->velocity[m][i];
        residual->velocity[m][neq]=0.0;
        residual->pressure[m][0]=0.0;
        for(e=1;e<=npno;e++)
            residual->pressure[m][e]=-action->pressure[m][e];
    }
}


static void apply_ala_coupled_triangular_once(
    struct All_variables *E, const struct ala_block_vector *residual,
    struct ala_block_vector *correction, double **pressure_work,
    double **velocity_work, int lev, int iteration,
    struct ala_pressure_preconditioner_cache *cache)
{
    int m,e,i,valid,npno,neq;
    double inner_accuracy;
    void parallel_process_termination();

    npno=E->lmesh.NPNO[lev];
    neq=E->lmesh.NEQ[lev];
    ala_block_vector_zero(E,correction);

    /* Upper block-triangular inverse approximation:
     *
     *   dp = -Sapprox^-1 rp,
     *   du = K_gamma^-1 (ru-G^T dp).
     *
     * This is the measured current-rheology baseline: one bounded velocity
     * solve per Arnoldi application.  The complete LDU experiment exposed an
     * unusable diagonal Schur approximation and increased inner-solve time by
     * more than four times, so it is not retained in the active path. */
    apply_ala_pressure_preconditioner(
        E,residual->pressure,correction->pressure,pressure_work,lev,
        iteration,cache);
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(e=1;e<=npno;e++)
            correction->pressure[m][e]=-correction->pressure[m][e];

    assemble_grad_rho_p(E,correction->pressure,velocity_work,lev);
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(i=0;i<neq;i++)
            velocity_work[m][i]=residual->velocity[m][i]
                                  -velocity_work[m][i];
    inner_accuracy=strict_ala_inner_accuracy(
        E,velocity_work,lev,
        E->control.ala_coupled_inner_relative_tolerance);
    valid=solve_del2_u_bounded(
        E,correction->velocity,velocity_work,inner_accuracy,lev,
        E->control.ala_coupled_inner_max_cycles,
        E->control.ala_coupled_inner_progress_interval);
    if(!valid)
        parallel_process_termination();
}


static int ala_solve_cached_element_k(const higher_precision *chol,
                                      const double rhs[ALA_VANKA_DOF],
                                      double solution[ALA_VANKA_DOF])
{
    int i,j;
    double forward[ALA_VANKA_DOF],sum,diagonal;

    for(i=0;i<ALA_VANKA_DOF;i++) {
        sum=rhs[i];
        for(j=0;j<i;j++)
            sum -= chol[i*(i+1)/2+j]*forward[j];
        diagonal=chol[i*(i+1)/2+i];
        if(!isfinite(diagonal) || diagonal<=0.0)
            return(0);
        forward[i]=sum/diagonal;
    }
    for(i=ALA_VANKA_DOF-1;i>=0;i--) {
        sum=forward[i];
        for(j=i+1;j<ALA_VANKA_DOF;j++)
            sum -= chol[j*(j+1)/2+i]*solution[j];
        solution[i]=sum/chol[i*(i+1)/2+i];
        if(!isfinite(solution[i]))
            return(0);
    }
    return(1);
}


/* One genuine mixed element-Vanka application for Stage 9f.1.  The cached
 * 24x24 K_gamma Cholesky factors define a stable local velocity inverse.  A
 * scalar Schur complement then solves the element's P0 pressure together
 * with its 24 velocity dofs:
 *
 *   du = K_P^-1 (ru-G_P^T dp),
 *   dp = (G_P K_P^-1 ru-rp)/(G_P K_P^-1 G_P^T+epsilon M_P).
 *
 * Velocity overlap is assembled with the same partition of unity as the
 * existing full-element smoother; P0 pressure is element-local. */
static void apply_ala_coupled_element_vanka_region(
    struct All_variables *E, const struct ala_block_vector *residual,
    struct ala_block_vector *correction, struct ala_block_vector *delta,
    int lev, int shallow_layers, int core_layers, int band_only)
{
    int m,e,a,d,i,node,eq,fixed,local_patches,global_patches;
    int radial_element,global_radial_element,selected_layers,selected_core;
    int depth_layer,shallow;
    int report_diagnostics;
    static int reported_cycle[MAX_LEVELS];
    static unsigned char reported_valid[MAX_LEVELS];
    static int reported_shallow_cycle[MAX_LEVELS];
    static unsigned char reported_shallow_valid[MAX_LEVELS];
    static int reported_band_cycle[MAX_LEVELS];
    static unsigned char reported_band_valid[MAX_LEVELS];
    double velocity_rhs[ALA_VANKA_DOF],gradient[ALA_VANKA_DOF];
    double local_weight[ALA_VANKA_DOF];
    double velocity_base[ALA_VANKA_DOF],velocity_pressure[ALA_VANKA_DOF];
    double pressure_rhs,pressure_solution,schur;
    double damping,weight,region_weight,taper_position;
    double local_min,local_max,global_min,global_max;
    higher_precision *chol;
    const int dims=E->mesh.nsd;
    const int ends=enodes[dims];
    const int neq=E->lmesh.NEQ[lev];
    const int npno=E->lmesh.NPNO[lev];
    const int elz=E->lmesh.ELZ[lev];

    if(loc_mat_size[dims]!=ALA_VANKA_DOF)
        myerror(E,"Coupled element-Vanka requires 3-D 24-dof elements");
    damping=E->control.ala_element_vanka_damping;
    shallow=(shallow_layers>0);
    selected_layers=min(shallow_layers,E->mesh.ELZ[lev]);
    selected_core=core_layers<0 ? -1 : min(core_layers,selected_layers);
    ala_block_vector_zero(E,correction);
    ala_block_vector_zero(E,delta);
    report_diagnostics=band_only
        ? (!reported_band_valid[lev] ||
           reported_band_cycle[lev]!=E->monitor.solution_cycles)
        : shallow
        ? (!reported_shallow_valid[lev] ||
           reported_shallow_cycle[lev]!=E->monitor.solution_cycles)
        : (!reported_valid[lev] ||
           reported_cycle[lev]!=E->monitor.solution_cycles);
    local_patches=0;
    local_min=1.0e300;
    local_max=0.0;

    /* Build the regional velocity partition before any local solve.  The
     * same inverse-square-root multiplicity is used to restrict the velocity
     * residual and G^T column and to prolong the velocity correction.  This
     * makes the local mixed Schur denominator belong to the actual weighted
     * patch map instead of combining an unweighted Schur with a one-sided
     * overlap average. */
    if(shallow) {
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=1;e<=E->lmesh.NEL[lev];e++) {
                radial_element=(e-1)%elz;
                global_radial_element=E->lmesh.EZS[lev]+radial_element;
                if(global_radial_element<
                   E->mesh.ELZ[lev]-selected_layers)
                    continue;
                depth_layer=E->mesh.ELZ[lev]-1-global_radial_element;
                if(band_only && depth_layer<selected_core)
                    continue;
                for(a=1;a<=ends;a++) {
                    node=E->IEN[lev][m][e].node[a];
                    for(d=0;d<dims;d++) {
                        eq=E->ID[lev][m][node].doff[d+1];
                        fixed=(d==0 && (E->NODE[lev][m][node]&VBX)) ||
                              (d==1 && (E->NODE[lev][m][node]&VBY)) ||
                              (d==2 && (E->NODE[lev][m][node]&VBZ));
                        if(!fixed)
                            correction->velocity[m][eq] += 1.0;
                    }
                }
            }
        (E->solver.exchange_id_d)(E,correction->velocity,lev);
        for(m=1;m<=E->sphere.caps_per_proc;m++) {
            for(eq=0;eq<neq;eq++)
                correction->velocity[m][eq]=
                    correction->velocity[m][eq]>0.0
                    ? 1.0/sqrt(correction->velocity[m][eq]) : 0.0;
            correction->velocity[m][neq]=0.0;
        }
    }

    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(e=1;e<=E->lmesh.NEL[lev];e++) {
            radial_element=(e-1)%elz;
            global_radial_element=E->lmesh.EZS[lev]+radial_element;
            if(shallow && global_radial_element<
                          E->mesh.ELZ[lev]-selected_layers)
                continue;
            depth_layer=E->mesh.ELZ[lev]-1-global_radial_element;
            if(band_only && depth_layer<selected_core)
                continue;
            region_weight=1.0;
            if(band_only) {
                taper_position=(depth_layer-selected_core+0.5)
                    /(double)(selected_layers-selected_core);
                region_weight=sin(3.14159265358979323846*taper_position);
            }
            if(!E->ALA_vanka_valid[lev][m][e])
                myerror(E,"Coupled element-Vanka velocity factor is invalid");
            chol=E->ALA_vanka_chol[lev][m]+e*ALA_VANKA_CHOL_SIZE;
            for(a=1;a<=ends;a++) {
                node=E->IEN[lev][m][e].node[a];
                for(d=0;d<dims;d++) {
                    i=(a-1)*dims+d;
                    eq=E->ID[lev][m][node].doff[d+1];
                    fixed=(d==0 && (E->NODE[lev][m][node]&VBX)) ||
                          (d==1 && (E->NODE[lev][m][node]&VBY)) ||
                          (d==2 && (E->NODE[lev][m][node]&VBZ));
                    weight=shallow
                        ? correction->velocity[m][eq]
                        : sqrt(E->ALA_vanka_overlap_BI[lev][m][eq]);
                    local_weight[i]=fixed ? 0.0 : weight;
                    velocity_rhs[i]=local_weight[i]
                        *residual->velocity[m][eq];
                    gradient[i]=local_weight[i]
                        *(E->elt_del[lev][m][e].g[i][0]
                          +E->elt_c[lev][m][e].c[i][0]);
                }
            }
            for(i=0;i<ALA_VANKA_DOF;i++)
                velocity_base[i]=velocity_pressure[i]=0.0;
            if(!ala_solve_cached_element_k(
                   chol,velocity_rhs,velocity_base) ||
               !ala_solve_cached_element_k(
                   chol,gradient,velocity_pressure))
                myerror(E,"Coupled element-Vanka velocity solve failed");
            pressure_rhs=residual->pressure[m][e];
            schur=E->control.ala_element_vanka_regularization
                *E->ECO[lev][m][e].area;
            pressure_solution=-pressure_rhs;
            for(i=0;i<ALA_VANKA_DOF;i++) {
                pressure_solution += gradient[i]*velocity_base[i];
                schur += gradient[i]*velocity_pressure[i];
            }
            if(!isfinite(schur) || schur<=1.0e-300)
                myerror(E,"Coupled element-Vanka pressure Schur is invalid");
            pressure_solution /= schur;
            /* Pressure and velocity corrections have different global
             * scalings in the mixed block. Keep pressure damping independent
             * so a small local Schur value cannot dominate the global step. */
            delta->pressure[m][e]=region_weight
                *E->control.ala_element_vanka_pressure_damping
                *pressure_solution;
            for(a=1;a<=ends;a++) {
                node=E->IEN[lev][m][e].node[a];
                for(d=0;d<dims;d++) {
                    i=(a-1)*dims+d;
                    eq=E->ID[lev][m][node].doff[d+1];
                    delta->velocity[m][eq] += region_weight*local_weight[i]
                        *(velocity_base[i]
                          -velocity_pressure[i]*pressure_solution);
                }
            }
            if(report_diagnostics) {
                local_min=min(local_min,schur);
                local_max=max(local_max,schur);
                local_patches++;
            }
        }

    (E->solver.exchange_id_d)(E,delta->velocity,lev);
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        for(eq=0;eq<neq;eq++) {
            correction->velocity[m][eq]=damping*delta->velocity[m][eq];
        }
        correction->velocity[m][neq]=0.0;
        correction->pressure[m][0]=0.0;
        for(e=1;e<=npno;e++)
            correction->pressure[m][e]=damping*delta->pressure[m][e];
    }
    if(report_diagnostics) {
        MPI_Allreduce(&local_patches,&global_patches,1,MPI_INT,MPI_SUM,
                      E->parallel.world);
        MPI_Allreduce(&local_min,&global_min,1,MPI_DOUBLE,MPI_MIN,
                      E->parallel.world);
        MPI_Allreduce(&local_max,&global_max,1,MPI_DOUBLE,MPI_MAX,
                      E->parallel.world);
        if(band_only) {
            reported_band_cycle[lev]=E->monitor.solution_cycles;
            reported_band_valid[lev]=1;
        }
        else if(shallow) {
            reported_shallow_cycle[lev]=E->monitor.solution_cycles;
            reported_shallow_valid[lev]=1;
        }
        else {
            reported_cycle[lev]=E->monitor.solution_cycles;
            reported_valid[lev]=1;
        }
        if(E->parallel.me==0) {
            fprintf(E->fp,"ALA COUPLED ELEMENT VANKA APPLICATION level=%d "
                    "solution_cycle=%d global_patches=%d velocity_damping=%e "
                    "pressure_damping=%e "
                    "pressure_regularization=%e "
                    "local_schur_range=[%e,%e] fallback_count=0 "
                    "region=%s radial_layers=%d core_layers=%d window=%s "
                    "overlap_contract=sqrt_partition "
                    "diagnostic_scope=first_application_per_cycle\n",
                    lev,E->monitor.solution_cycles,global_patches,damping,
                    E->control.ala_element_vanka_pressure_damping,
                    E->control.ala_element_vanka_regularization,
                    global_min,global_max,
                    band_only ? "shallow_band" :
                        (shallow ? "shallow_core" : "full"),selected_layers,
                    selected_core,band_only ? "sine" : "off");
            fprintf(stderr,"ALA COUPLED ELEMENT VANKA APPLICATION level=%d "
                    "solution_cycle=%d global_patches=%d velocity_damping=%e "
                    "pressure_damping=%e "
                    "pressure_regularization=%e "
                    "local_schur_range=[%e,%e] fallback_count=0 "
                    "region=%s radial_layers=%d core_layers=%d window=%s "
                    "overlap_contract=sqrt_partition "
                    "diagnostic_scope=first_application_per_cycle\n",
                    lev,E->monitor.solution_cycles,global_patches,damping,
                    E->control.ala_element_vanka_pressure_damping,
                    E->control.ala_element_vanka_regularization,
                    global_min,global_max,
                    band_only ? "shallow_band" :
                        (shallow ? "shallow_core" : "full"),selected_layers,
                    selected_core,band_only ? "sine" : "off");
            fflush(E->fp);
            fflush(stderr);
        }
    }
}


static void apply_ala_coupled_element_vanka_once(
    struct All_variables *E, const struct ala_block_vector *residual,
    struct ala_block_vector *correction, struct ala_block_vector *delta,
    int lev)
{
    apply_ala_coupled_element_vanka_region(
        E,residual,correction,delta,lev,0,-1,0);
}


static void assemble_ala_coupled_block_defect(
    struct All_variables *E, const struct ala_block_vector *residual,
    const struct ala_block_vector *correction,
    struct ala_block_vector *defect, struct ala_block_vector *action,
    struct ala_block_vector *velocity_work, int lev)
{
    int m,i,e;
    int neq=E->lmesh.NEQ[lev];
    int npno=E->lmesh.NPNO[lev];

    apply_ala_coupled_operator(
        E,correction->velocity,correction->pressure,
        action->velocity,action->pressure,velocity_work->velocity,lev);
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        for(i=0;i<neq;i++)
            defect->velocity[m][i]=residual->velocity[m][i]
                                      -action->velocity[m][i];
        defect->velocity[m][neq]=0.0;
        defect->pressure[m][0]=0.0;
        for(e=1;e<=npno;e++)
            defect->pressure[m][e]=residual->pressure[m][e]
                                      -action->pressure[m][e];
    }
}


/* Observe, but never modify, the coefficient consistency of the current
 * rediscretized hierarchy.  A Galerkin-consistent constant-P0 aggregate has
 * S_c=sum(children S_f).  The ratio and relative mismatch below decide
 * whether a new hierarchy is required; they are not used as a correction. */
static void ala_audit_local_schur_coarsening(
    struct All_variables *E, int fine_level)
{
    int m,cx,cy,cz,ce,dx,dy,dz,ex,ey,ez,e;
    int local_count,global_count;
    int report_diagnostics;
    static int reported_cycle[MAX_LEVELS];
    static unsigned char reported_valid[MAX_LEVELS];
    double fine_sum,coarse_schur,ratio,mismatch;
    double local_min,global_min,local_max,global_max,local_sum,global_sum;
    double local_mismatch_sum,global_mismatch_sum;
    double local_mismatch_max,global_mismatch_max;
    const int coarse_level=fine_level-1;
    const int felz=E->lmesh.ELZ[fine_level];
    const int felx=E->lmesh.ELX[fine_level];
    const int celz=E->lmesh.ELZ[coarse_level];
    const int celx=E->lmesh.ELX[coarse_level];
    const int cely=E->lmesh.ELY[coarse_level];

    if(!E->control.ala_viscosity_spectrum_diagnostics)
        return;
    report_diagnostics=!reported_valid[fine_level] ||
        reported_cycle[fine_level]!=E->monitor.solution_cycles;
    local_count=0;
    local_min=1.0e300;
    local_max=local_sum=0.0;
    local_mismatch_sum=local_mismatch_max=0.0;
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(cy=1;cy<=cely;cy++)
            for(cx=1;cx<=celx;cx++)
                for(cz=1;cz<=celz;cz++) {
                    ce=cz+(cx-1)*celz+(cy-1)*celz*celx;
                    coarse_schur=E->ALA_vanka_schur[coarse_level][m][ce];
                    fine_sum=0.0;
                    for(dy=0;dy<2;dy++)
                        for(dx=0;dx<2;dx++)
                            for(dz=0;dz<2;dz++) {
                                ex=2*cx-1+dx;
                                ey=2*cy-1+dy;
                                ez=2*cz-1+dz;
                                e=ez+(ex-1)*felz+(ey-1)*felz*felx;
                                fine_sum +=
                                    E->ALA_vanka_schur[fine_level][m][e];
                            }
                    if(!isfinite(coarse_schur) || coarse_schur<=0.0 ||
                       !isfinite(fine_sum) || fine_sum<=0.0)
                        myerror(E,"ALA Schur coarsening audit "
                                "has an invalid local Schur scale");
                    ratio=coarse_schur/fine_sum;
                    mismatch=fabs(coarse_schur-fine_sum)
                        /(fabs(coarse_schur)+fabs(fine_sum));
                    if(report_diagnostics) {
                        local_min=min(local_min,ratio);
                        local_max=max(local_max,ratio);
                        local_sum += ratio;
                        local_mismatch_sum += mismatch;
                        local_mismatch_max=max(local_mismatch_max,mismatch);
                        local_count++;
                    }
                }
    if(report_diagnostics) {
        MPI_Allreduce(&local_count,&global_count,1,MPI_INT,MPI_SUM,
                      E->parallel.world);
        MPI_Allreduce(&local_min,&global_min,1,MPI_DOUBLE,MPI_MIN,
                      E->parallel.world);
        MPI_Allreduce(&local_max,&global_max,1,MPI_DOUBLE,MPI_MAX,
                      E->parallel.world);
        MPI_Allreduce(&local_sum,&global_sum,1,MPI_DOUBLE,MPI_SUM,
                      E->parallel.world);
        MPI_Allreduce(&local_mismatch_sum,&global_mismatch_sum,1,
                      MPI_DOUBLE,MPI_SUM,E->parallel.world);
        MPI_Allreduce(&local_mismatch_max,&global_mismatch_max,1,
                      MPI_DOUBLE,MPI_MAX,E->parallel.world);
        if(global_count<=0)
            myerror(E,"ALA Schur coarsening audit is empty");
        if(E->parallel.me==0) {
            fprintf(E->fp,"ALA SCHUR COARSENING AUDIT fine_level=%d "
                    "coarse_level=%d ratio_range=[%e,%e] ratio_mean=%e "
                    "relative_mismatch_mean=%e relative_mismatch_max=%e "
                    "elements=%d cycle=%d action=observe_only\n",
                    fine_level,coarse_level,global_min,global_max,
                    global_sum/global_count,
                    global_mismatch_sum/global_count,global_mismatch_max,
                    global_count,
                    E->monitor.solution_cycles);
            fprintf(stderr,"ALA SCHUR COARSENING AUDIT fine_level=%d "
                    "coarse_level=%d ratio_range=[%e,%e] ratio_mean=%e "
                    "relative_mismatch_mean=%e relative_mismatch_max=%e "
                    "elements=%d cycle=%d action=observe_only\n",
                    fine_level,coarse_level,global_min,global_max,
                    global_sum/global_count,
                    global_mismatch_sum/global_count,global_mismatch_max,
                    global_count,
                    E->monitor.solution_cycles);
            fflush(E->fp);
            fflush(stderr);
        }
        reported_cycle[fine_level]=E->monitor.solution_cycles;
        reported_valid[fine_level]=1;
    }
}


static void apply_ala_coupled_multilevel_vcycle(
    struct All_variables *E, const struct ala_block_vector *residual,
    struct ala_block_vector *correction, int lev)
{
    struct ala_block_vector *smooth,*delta,*defect,*action,*velocity_work;
    struct ala_block_vector *coarse_residual=NULL,*coarse_correction=NULL;
    int band_sweep,band_sweeps,coarse_sweep,shallow_sweep,shallow_sweeps;
    double coarse_weight;

    smooth=ala_block_vector_create(E,lev);
    delta=ala_block_vector_create(E,lev);
    defect=ala_block_vector_create(E,lev);
    action=ala_block_vector_create(E,lev);
    velocity_work=ala_block_vector_create(E,lev);
    ala_block_vector_zero(E,correction);
    shallow_sweeps=E->control.ala_coupled_shallow_vanka_sweeps;
    band_sweeps=E->control.ala_coupled_shallow_vanka_band_sweeps;
    if(E->monitor.solution_cycles>0 &&
       E->control.ala_coupled_shallow_vanka_warm_sweeps>=0)
        shallow_sweeps=E->control.ala_coupled_shallow_vanka_warm_sweeps;

    /* One mixed pre-sweep. */
    apply_ala_coupled_element_vanka_once(
        E,residual,smooth,delta,lev);
    ala_block_vector_axpy(E,1.0,smooth,correction);

    if(lev>E->mesh.levmin) {
        assemble_ala_coupled_block_defect(
            E,residual,correction,defect,action,velocity_work,lev);
        coarse_residual=ala_block_vector_create(E,lev-1);
        coarse_correction=ala_block_vector_create(E,lev-1);
        ala_coupled_restrict_velocity(
            E,lev,defect->velocity,coarse_residual->velocity);
        ala_coupled_restrict_pressure_p0_transpose(
            E,lev,defect->pressure,coarse_residual->pressure);
        ala_audit_local_schur_coarsening(E,lev);
        apply_ala_coupled_multilevel_vcycle(
            E,coarse_residual,coarse_correction,lev-1);
        ala_block_vector_zero(E,smooth);
        ala_coupled_prolong_velocity(
            E,lev-1,coarse_correction->velocity,smooth->velocity);
        ala_coupled_prolong_pressure_p0(
            E,lev-1,coarse_correction->pressure,smooth->pressure);
        /* The configured damping stabilizes the one coarse correction seen
         * by the finest operator.  Applying it at every recursive interface
         * attenuates the globally coarsest mode by weight^nlevels (1/16 for
         * the production 0.5/five-level hierarchy), precisely where the
         * residual diagnostics show the remaining energy. */
        coarse_weight=(lev==E->mesh.levmax)
            ? E->control.ala_coupled_multilevel_coarse_weight : 1.0;
        ala_block_vector_axpy(
            E,coarse_weight,smooth,correction);
        ala_block_vector_destroy(E,coarse_residual);
        ala_block_vector_destroy(E,coarse_correction);
    }

    /* One mixed post-sweep; on the coarsest level this is the second local
     * solve and serves as the deliberately bounded coarse solver. */
    assemble_ala_coupled_block_defect(
        E,residual,correction,defect,action,velocity_work,lev);
    apply_ala_coupled_element_vanka_once(E,defect,smooth,delta,lev);
    ala_block_vector_axpy(E,1.0,smooth,correction);
    if(lev==E->mesh.levmin)
        for(coarse_sweep=2;
            coarse_sweep<E->control.ala_coupled_multilevel_coarse_sweeps;
            coarse_sweep++) {
            assemble_ala_coupled_block_defect(
                E,residual,correction,defect,action,velocity_work,lev);
            apply_ala_coupled_element_vanka_once(
                E,defect,smooth,delta,lev);
            ala_block_vector_axpy(E,1.0,smooth,correction);
        }

    /* Prescribed plate velocities excite a shallow high-frequency
     * continuity mode that barely restricts to the next pressure level.
     * Apply defect-based mixed smoothing only in the configured top radial
     * layers of the finest grid.  Dynamic selected-patch multiplicities keep
     * the velocity partition of unity exact at the lower region interface. */
    if(lev==E->mesh.levmax)
        for(shallow_sweep=0;
            shallow_sweep<shallow_sweeps;
            shallow_sweep++) {
            assemble_ala_coupled_block_defect(
                E,residual,correction,defect,action,velocity_work,lev);
            apply_ala_coupled_element_vanka_region(
                E,defect,smooth,delta,lev,
                E->control.ala_coupled_shallow_vanka_core_layers>=0
                    ? E->control.ala_coupled_shallow_vanka_core_layers
                    : E->control.ala_coupled_shallow_vanka_layers,
                -1,0);
            ala_block_vector_axpy(E,1.0,smooth,correction);
        }

    /* Apply the upper-mantle extension independently from the proven core.
     * Its sine window vanishes at both radial boundaries, while the bounded
     * sweep count prevents the extension from being repeated with every
     * strong surface correction. */
    if(lev==E->mesh.levmax)
        for(band_sweep=0;band_sweep<band_sweeps;band_sweep++) {
            assemble_ala_coupled_block_defect(
                E,residual,correction,defect,action,velocity_work,lev);
            apply_ala_coupled_element_vanka_region(
                E,defect,smooth,delta,lev,
                E->control.ala_coupled_shallow_vanka_layers,
                E->control.ala_coupled_shallow_vanka_core_layers,1);
            ala_block_vector_axpy(E,1.0,smooth,correction);
        }

    ala_block_vector_destroy(E,smooth);
    ala_block_vector_destroy(E,delta);
    ala_block_vector_destroy(E,defect);
    ala_block_vector_destroy(E,action);
    ala_block_vector_destroy(E,velocity_work);
}


static void apply_ala_coupled_block_preconditioner(
    struct All_variables *E, const struct ala_block_vector *residual,
    struct ala_block_vector *correction,
    struct ala_block_vector *correction_work,
    struct ala_block_vector *action_work,
    struct ala_block_vector *vanka_delta, double **pressure_work,
    struct ala_block_vector *operator_work, int lev, int iteration,
    struct ala_pressure_preconditioner_cache *cache)
{
    int m,i,e,pass,neq,npno,ex,ey,ez,fx,fy,fz,elx,ely,elz;
    int first_ez,projected_layers;
    double mean;

    neq=E->lmesh.NEQ[lev];
    npno=E->lmesh.NPNO[lev];
    if(E->control.ala_coupled_multilevel_vcycle) {
        static int reported_cycle;
        static int report_valid;

        apply_ala_coupled_multilevel_vcycle(E,residual,correction,lev);
        if(E->parallel.me==0 &&
           (!report_valid || reported_cycle!=E->monitor.solution_cycles)) {
            fprintf(E->fp,"ALA COUPLED MULTILEVEL VCYCLE APPLICATION "
                    "levels=%d range=%d:%d pre_sweeps=1 post_sweeps=1 "
                    "coarse_sweeps=%d coarse_weight=%e "
                    "nested_coarse_weight=1.000000e+00 "
                    "shallow_layers=%d shallow_core_layers=%d "
                    "shallow_band_window=sine shallow_band_sweeps=%d "
                    "shallow_sweeps=%d "
                    "shallow_cold_sweeps=%d shallow_warm_sweeps=%d "
                    "solution_start=%s\n",
                    E->mesh.levmax-E->mesh.levmin+1,
                    E->mesh.levmin,E->mesh.levmax,
                    E->control.ala_coupled_multilevel_coarse_sweeps,
                    E->control.ala_coupled_multilevel_coarse_weight,
                    E->control.ala_coupled_shallow_vanka_layers,
                    E->control.ala_coupled_shallow_vanka_core_layers,
                    E->control.ala_coupled_shallow_vanka_band_sweeps,
                    (E->monitor.solution_cycles>0 &&
                     E->control.ala_coupled_shallow_vanka_warm_sweeps>=0)
                        ? E->control.ala_coupled_shallow_vanka_warm_sweeps
                        : E->control.ala_coupled_shallow_vanka_sweeps,
                    E->control.ala_coupled_shallow_vanka_sweeps,
                    E->control.ala_coupled_shallow_vanka_warm_sweeps,
                    E->monitor.solution_cycles>0 ? "warm" : "cold");
            reported_cycle=E->monitor.solution_cycles;
            report_valid=1;
            fprintf(stderr,"ALA COUPLED MULTILEVEL VCYCLE APPLICATION "
                    "levels=%d range=%d:%d pre_sweeps=1 post_sweeps=1 "
                    "coarse_sweeps=%d coarse_weight=%e "
                    "nested_coarse_weight=1.000000e+00 "
                    "shallow_layers=%d shallow_core_layers=%d "
                    "shallow_band_window=sine shallow_band_sweeps=%d "
                    "shallow_sweeps=%d "
                    "shallow_cold_sweeps=%d shallow_warm_sweeps=%d "
                    "solution_start=%s\n",
                    E->mesh.levmax-E->mesh.levmin+1,
                    E->mesh.levmin,E->mesh.levmax,
                    E->control.ala_coupled_multilevel_coarse_sweeps,
                    E->control.ala_coupled_multilevel_coarse_weight,
                    E->control.ala_coupled_shallow_vanka_layers,
                    E->control.ala_coupled_shallow_vanka_core_layers,
                    E->control.ala_coupled_shallow_vanka_band_sweeps,
                    (E->monitor.solution_cycles>0 &&
                     E->control.ala_coupled_shallow_vanka_warm_sweeps>=0)
                        ? E->control.ala_coupled_shallow_vanka_warm_sweeps
                        : E->control.ala_coupled_shallow_vanka_sweeps,
                    E->control.ala_coupled_shallow_vanka_sweeps,
                    E->control.ala_coupled_shallow_vanka_warm_sweeps,
                    E->monitor.solution_cycles>0 ? "warm" : "cold");
            fflush(E->fp);
            fflush(stderr);
        }
        return;
    }
    if(E->control.ala_coupled_element_vanka) {
        apply_ala_coupled_element_vanka_once(
            E,residual,correction,correction_work,lev);
        return;
    }
    apply_ala_coupled_triangular_once(
        E,residual,correction,pressure_work,operator_work->velocity,lev,
        iteration,cache);

    /* Current-rheology shallow defect correction:
     *
     *   z_0 = P_tri^-1 r,
     *   d   = r-A z_0,
     *   z_1 = z_0+P_Vanka,shallow^-1 d.
     *
     * The triangular map retains the globally useful pressure direction;
     * the true element K_gamma inverse is allowed to act only on the shallow
     * fine-scale defect that dominated the completed 220-iteration audit.
     * This avoids letting element-local coarse modes replace the global BPI
     * direction, which the standalone gamma=1 and gamma=10 audits rejected. */
    if(E->control.ala_coupled_shallow_vanka_sweeps>0) {
        for(pass=0;pass<E->control.ala_coupled_shallow_vanka_sweeps;pass++) {
            assemble_ala_coupled_block_defect(
                E,residual,correction,action_work,vanka_delta,operator_work,
                lev);
            apply_ala_coupled_element_vanka_region(
                E,action_work,correction_work,vanka_delta,lev,
                E->control.ala_coupled_shallow_vanka_layers,-1,0);
            ala_block_vector_axpy(E,1.0,correction_work,correction);
        }
    }

    /* Multiplicative shallow factor-two coarse correction.
     *
     * Form the true coupled defect after the proven triangular plus shallow
     * Vanka map, project only its pressure component onto piecewise constants
     * over aligned 2x2x2 pressure cells, and apply one complete triangular
     * block correction.  Because this acts on r-Az rather than r, the coarse
     * component is not counted once by BPI and then added a second time.  The
     * arithmetic block mean is the Euclidean orthogonal projector Q=Q^T=Q^2;
     * rounding the configured shallow extent to a full radial pair preserves
     * both the factor-two hierarchy and the constant-pressure mode. */
    if(E->control.ala_coupled_factor2_coarse_correction) {
        static int reported_cycle;
        static int report_valid;

        elx=E->lmesh.ELX[lev];
        ely=E->lmesh.ELY[lev];
        elz=E->lmesh.ELZ[lev];
        if(elx%2 || ely%2 || elz%2)
            myerror(E,"ALA factor-two coarse defect requires even mesh");
        projected_layers=2*((min(E->control.ala_coupled_shallow_vanka_layers,
                                 elz)+1)/2);
        projected_layers=min(projected_layers,elz);
        first_ez=elz-projected_layers+1;
        assemble_ala_coupled_block_defect(
            E,residual,correction,action_work,vanka_delta,operator_work,lev);
        for(m=1;m<=E->sphere.caps_per_proc;m++) {
            for(i=0;i<neq;i++)
                action_work->velocity[m][i]=0.0;
            action_work->velocity[m][neq]=0.0;
            for(e=1;e<=npno;e++) {
                vanka_delta->pressure[m][e]=action_work->pressure[m][e];
                action_work->pressure[m][e]=0.0;
            }
            action_work->pressure[m][0]=0.0;
            for(ey=1;ey<=ely;ey+=2)
                for(ex=1;ex<=elx;ex+=2)
                    for(ez=first_ez;ez<=elz;ez+=2) {
                        mean=0.0;
                        for(fy=ey;fy<=ey+1;fy++)
                            for(fx=ex;fx<=ex+1;fx++)
                                for(fz=ez;fz<=ez+1;fz++) {
                                    e=fz+(fx-1)*elz+(fy-1)*elz*elx;
                                    mean += vanka_delta->pressure[m][e];
                                }
                        mean /= 8.0;
                        for(fy=ey;fy<=ey+1;fy++)
                            for(fx=ex;fx<=ex+1;fx++)
                                for(fz=ez;fz<=ez+1;fz++) {
                                    e=fz+(fx-1)*elz+(fy-1)*elz*elx;
                                    action_work->pressure[m][e]=mean;
                                }
                    }
        }
        apply_ala_coupled_triangular_once(
            E,action_work,correction_work,pressure_work,
            operator_work->velocity,lev,iteration,cache);
        ala_block_vector_axpy(E,1.0,correction_work,correction);
        if(E->parallel.me==0 &&
           (!report_valid || reported_cycle!=E->monitor.solution_cycles)) {
            fprintf(E->fp,"ALA COUPLED FACTOR2 DEFECT CORRECTION "
                    "level=%d projected_layers=%d block=2x2x2 "
                    "projection=orthogonal_mean application=multiplicative "
                    "coarse_solver=triangular\n",lev,projected_layers);
            fprintf(stderr,"ALA COUPLED FACTOR2 DEFECT CORRECTION "
                    "level=%d projected_layers=%d block=2x2x2 "
                    "projection=orthogonal_mean application=multiplicative "
                    "coarse_solver=triangular\n",lev,projected_layers);
            fflush(E->fp);
            fflush(stderr);
            reported_cycle=E->monitor.solution_cycles;
            report_valid=1;
        }
    }

    /* Multiplicative block defect correction:
     *
     *   z_0 = P_tri^-1 r,
     *   d_k = r - A z_k,
     *   z_{k+1} = z_k + P_tri^-1 d_k.
     *
     * This is a systematic refinement of the complete block inverse, not a
     * pressure-only blend.  The explicit A action preserves every D+C cross
     * term, and FGMRES permits the bounded velocity solves to vary between
     * passes.  Zero corrections reproduces the Stage 9d.4 map exactly. */
    for(pass=0;pass<E->control.ala_coupled_defect_corrections;pass++) {
        apply_ala_coupled_operator(
            E,correction->velocity,correction->pressure,
            action_work->velocity,action_work->pressure,
            operator_work->velocity,lev);
        for(m=1;m<=E->sphere.caps_per_proc;m++) {
            for(i=0;i<neq;i++)
                action_work->velocity[m][i]=residual->velocity[m][i]
                                                -action_work->velocity[m][i];
            action_work->velocity[m][neq]=0.0;
            action_work->pressure[m][0]=0.0;
            for(e=1;e<=npno;e++)
                action_work->pressure[m][e]=residual->pressure[m][e]
                                                -action_work->pressure[m][e];
        }
        apply_ala_coupled_triangular_once(
            E,action_work,correction_work,pressure_work,
            operator_work->velocity,lev,iteration,cache);
        ala_block_vector_axpy(E,1.0,correction_work,correction);
    }
}


static void strict_ala_pressure_mode_audit(
    struct All_variables *E, const struct ala_block_vector *residual,
    const struct ala_block_vector *action, int lev, int iteration)
{
    int offset,clev,factor,m,cx,cy,cz,fx,fy,fz,e,base,n_offsets;
    int celx,cely,celz,elx,ely,elz,nchildren,i;
    double mean_r,mean_a,value_r,value_a;
    double local[24],global[24],coarse_cosine,fine_cosine;
    double coarse_scale,fine_scale,total_r2;

    elx=E->lmesh.ELX[lev];
    ely=E->lmesh.ELY[lev];
    elz=E->lmesh.ELZ[lev];
    n_offsets=min(4,lev-E->mesh.levmin);
    for(i=0;i<24;i++)
        local[i]=0.0;
    for(offset=1;offset<=n_offsets;offset++) {
        base=6*(offset-1);
        clev=lev-offset;
        factor=1<<offset;
        celx=E->lmesh.ELX[clev];
        cely=E->lmesh.ELY[clev];
        celz=E->lmesh.ELZ[clev];
        if(elx!=factor*celx || ely!=factor*cely || elz!=factor*celz)
            myerror(E,"ALA pressure mode audit mesh mismatch");
        nchildren=factor*factor*factor;
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(cy=1;cy<=cely;cy++)
                for(cx=1;cx<=celx;cx++)
                    for(cz=1;cz<=celz;cz++) {
                        mean_r=0.0;
                        mean_a=0.0;
                        for(fy=(cy-1)*factor+1;fy<=cy*factor;fy++)
                            for(fx=(cx-1)*factor+1;fx<=cx*factor;fx++)
                                for(fz=(cz-1)*factor+1;fz<=cz*factor;fz++) {
                                    e=fz+(fx-1)*elz+(fy-1)*elz*elx;
                                    mean_r += residual->pressure[m][e];
                                    mean_a += action->pressure[m][e];
                                }
                        mean_r /= (double)nchildren;
                        mean_a /= (double)nchildren;
                        local[base+0] += nchildren*mean_r*mean_r;
                        local[base+1] += nchildren*mean_a*mean_a;
                        local[base+2] += nchildren*mean_r*mean_a;
                        for(fy=(cy-1)*factor+1;fy<=cy*factor;fy++)
                            for(fx=(cx-1)*factor+1;fx<=cx*factor;fx++)
                                for(fz=(cz-1)*factor+1;fz<=cz*factor;fz++) {
                                    e=fz+(fx-1)*elz+(fy-1)*elz*elx;
                                    value_r=residual->pressure[m][e]-mean_r;
                                    value_a=action->pressure[m][e]-mean_a;
                                    local[base+3] += value_r*value_r;
                                    local[base+4] += value_a*value_a;
                                    local[base+5] += value_r*value_a;
                                }
                    }
    }
    MPI_Allreduce(local,global,6*n_offsets,MPI_DOUBLE,MPI_SUM,
                  E->parallel.world);
    for(offset=1;offset<=n_offsets;offset++) {
        base=6*(offset-1);
        factor=1<<offset;
        coarse_cosine=global[base+2]
            /sqrt(max(global[base+0]*global[base+1],1.0e-300));
        fine_cosine=global[base+5]
            /sqrt(max(global[base+3]*global[base+4],1.0e-300));
        coarse_scale=global[base+2]/max(global[base+1],1.0e-300);
        fine_scale=global[base+5]/max(global[base+4],1.0e-300);
        total_r2=global[base+0]+global[base+3];
        if(E->parallel.me==0) {
            fprintf(E->fp,"ALA PRESSURE MODE AUDIT iteration=%d offset=%d "
                    "factor=%d residual_fraction=(coarse:%e,fine:%e) "
                    "cosine=(coarse:%e,fine:%e) "
                    "optimal_scale=(coarse:%e,fine:%e)\n",
                    iteration,offset,factor,
                    global[base+0]/max(total_r2,1.0e-300),
                    global[base+3]/max(total_r2,1.0e-300),
                    coarse_cosine,fine_cosine,coarse_scale,fine_scale);
            fflush(E->fp);
        }
    }
}


/* Resolve the pressure part of A P^-1 r by physical depth before changing
 * the shallow Schwarz weights.  The same residual/action pair used by the
 * block audit is split into fixed mantle bands, so each reported least-squares
 * scale is a direct property of the current Tref rheology and preconditioner
 * action rather than a residual-amplitude guess. */
static void strict_ala_pressure_depth_action_audit(
    struct All_variables *E, const struct ala_block_vector *residual,
    const struct ala_block_vector *action, int lev, int iteration)
{
    int m,e,ez,bin,elz,i;
    double depth_km,total_r2,cosine,optimal_scale,action_to_residual;
    double boundaries[7],local[18],global[18];

    boundaries[0]=0.0;
    boundaries[1]=200.0;
    boundaries[2]=410.0;
    boundaries[3]=660.0;
    boundaries[4]=1000.0;
    boundaries[5]=2000.0;
    boundaries[6]=1.0e30;
    for(i=0;i<18;i++)
        local[i]=0.0;
    elz=E->lmesh.ELZ[lev];
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(e=1;e<=E->lmesh.NPNO[lev];e++) {
            ez=(e-1)%elz+1;
            depth_km=(1.0-0.5*(E->sx[m][3][ez]
                              +E->sx[m][3][ez+1]))*E->data.radius_km;
            bin=5;
            for(i=0;i<5;i++)
                if(depth_km<boundaries[i+1]) {
                    bin=i;
                    break;
                }
            local[3*bin+0] += residual->pressure[m][e]
                                  *residual->pressure[m][e];
            local[3*bin+1] += action->pressure[m][e]
                                  *action->pressure[m][e];
            local[3*bin+2] += residual->pressure[m][e]
                                  *action->pressure[m][e];
        }
    MPI_Allreduce(local,global,18,MPI_DOUBLE,MPI_SUM,E->parallel.world);
    total_r2=0.0;
    for(bin=0;bin<6;bin++)
        total_r2 += global[3*bin+0];
    if(E->parallel.me==0)
        for(bin=0;bin<6;bin++) {
            cosine=global[3*bin+2]
                /sqrt(max(global[3*bin+0]*global[3*bin+1],1.0e-300));
            optimal_scale=global[3*bin+2]
                /max(global[3*bin+1],1.0e-300);
            action_to_residual=sqrt(global[3*bin+1]
                /max(global[3*bin+0],1.0e-300));
            fprintf(E->fp,"ALA PRESSURE DEPTH ACTION AUDIT iteration=%d "
                    "bin=%d depth_km=[%e,%e) residual_fraction=%e "
                    "action_to_residual=%e cosine=%e optimal_scale=%e\n",
                    iteration,bin,boundaries[bin],boundaries[bin+1],
                    global[3*bin+0]/max(total_r2,1.0e-300),
                    action_to_residual,cosine,optimal_scale);
        }
    if(E->parallel.me==0)
        fflush(E->fp);
}


static void strict_ala_coupled_preconditioner_audit(
    struct All_variables *E, const struct ala_block_vector *residual,
    const struct ala_block_vector *correction,
    const struct ala_block_vector *action, double velocity_weight,
    double pressure_weight, int iteration)
{
    double rv,rp,zv,zp,av,ap,rav,rap;
    double velocity_cosine,pressure_cosine,block_cosine;
    double velocity_defect,pressure_defect,block_rnorm,block_anorm;
    double velocity_optimal_scale,pressure_optimal_scale,block_optimal_scale;
    double velocity_projected_defect,pressure_projected_defect;
    double block_projected_defect,block_product;

    ala_block_vector_component_norms(E,residual,&rv,&rp);
    ala_block_vector_component_norms(E,correction,&zv,&zp);
    ala_block_vector_component_norms(E,action,&av,&ap);
    ala_block_vector_component_products(
        E,residual,action,&rav,&rap);
    velocity_cosine=rav/max(rv*av,1.0e-300);
    pressure_cosine=rap/max(rp*ap,1.0e-300);
    block_rnorm=sqrt(max(velocity_weight*rv*rv+
                         pressure_weight*rp*rp,0.0));
    block_anorm=sqrt(max(velocity_weight*av*av+
                         pressure_weight*ap*ap,0.0));
    block_cosine=(velocity_weight*rav+pressure_weight*rap)
        /max(block_rnorm*block_anorm,1.0e-300);
    velocity_defect=sqrt(max(rv*rv+av*av-2.0*rav,0.0))
        /max(rv,1.0e-300);
    pressure_defect=sqrt(max(rp*rp+ap*ap-2.0*rap,0.0))
        /max(rp,1.0e-300);
    velocity_optimal_scale=rav/max(av*av,1.0e-300);
    pressure_optimal_scale=rap/max(ap*ap,1.0e-300);
    block_product=velocity_weight*rav+pressure_weight*rap;
    block_optimal_scale=block_product/max(block_anorm*block_anorm,1.0e-300);
    velocity_projected_defect=sqrt(max(1.0-velocity_cosine*velocity_cosine,0.0));
    pressure_projected_defect=sqrt(max(1.0-pressure_cosine*pressure_cosine,0.0));
    block_projected_defect=sqrt(max(1.0-block_cosine*block_cosine,0.0));
    if(E->parallel.me==0) {
        fprintf(E->fp,"ALA COUPLED PRECONDITIONER AUDIT iteration=%d "
                "r_components=(velocity:%e,pressure:%e) "
                "z_components=(velocity:%e,pressure:%e) "
                "Az_components=(velocity:%e,pressure:%e) "
                "Az_to_r=(velocity:%e,pressure:%e) "
                "defect_to_r=(velocity:%e,pressure:%e) "
                "cosine=(velocity:%e,pressure:%e,block:%e) "
                "optimal_scale=(velocity:%e,pressure:%e,block:%e) "
                "projected_defect=(velocity:%e,pressure:%e,block:%e)\n",
                iteration,rv,rp,zv,zp,av,ap,
                av/max(rv,1.0e-300),ap/max(rp,1.0e-300),
                velocity_defect,pressure_defect,
                velocity_cosine,pressure_cosine,block_cosine,
                velocity_optimal_scale,pressure_optimal_scale,
                block_optimal_scale,velocity_projected_defect,
                pressure_projected_defect,block_projected_defect);
        fflush(E->fp);
    }
    strict_ala_pressure_depth_action_audit(
        E,residual,action,E->mesh.levmax,iteration);
    strict_ala_pressure_mode_audit(
        E,residual,action,E->mesh.levmax,iteration);
}


/* First monolithic strict-ALA prototype.  It solves the complete augmented
 * saddle-point block with right-preconditioned restarted FGMRES.  The block
 * triangular preconditioner reuses one velocity solve and the proven pressure
 * BPI/Schwarz map; coupled Vanka and coupled multigrid are later stages.
 * Explicit block residual replacement is performed at every candidate so the
 * initial A/B prioritizes algebraic correctness over final performance. */
static float solve_ala_coupled_fgmres_core(
    struct All_variables *E, double **V, double **P, int *steps_max, int lev,
    struct ala_pressure_preconditioner_cache *cache,
    double **pressure_work)
{
    int m,i,j,e,count,restart,max_basis,used,breakdown,converged;
    int momentum_gate,continuity_converged,momentum_converged;
    int neq,npno;
    double h[65][64],cs[64],sn[64],g[65],y[64],y_old[64];
    double beta,norm,sum,delta,residual_est,explicit_norm,block_relative;
    double velocity_component,pressure_component;
    double velocity_weight,pressure_weight,initial_block_norm;
    double mass_norm,initial_mass_norm,cancellation_l2,best_cancellation;
    double momentum_rms,momentum_relative,augmented_momentum_relative;
    double force_norm,velocity_norm,initial_velocity_norm;
    double global_vdot();
    void assemble_div_u();
    void parallel_process_termination();
    struct ala_block_vector *rhs,*r,*explicit_r,*w,*operator_work;
    struct ala_block_vector *preconditioner_delta,*preconditioner_vanka_delta;
    struct ala_block_vector **vb,**zb;
    const char *acceptance_status;
    const char *preconditioner_name;

    neq=E->lmesh.NEQ[lev];
    npno=E->lmesh.NPNO[lev];
    audit_ala_coupled_multilevel_contracts(E);
    if(E->control.ala_coupled_multilevel_audit_only) {
        if(E->parallel.me==0 && E->fp!=NULL) {
            fprintf(E->fp,
                    "ALA COUPLED MULTILEVEL AUDIT COMPLETE: "
                    "terminating before FGMRES iteration 1\n");
            fflush(E->fp);
        }
        if(E->parallel.me==0) {
            fprintf(stderr,
                    "ALA COUPLED MULTILEVEL AUDIT COMPLETE: "
                    "terminating before FGMRES iteration 1\n");
            fflush(stderr);
        }
        MPI_Barrier(E->parallel.world);
        MPI_Finalize();
        exit(EXIT_SUCCESS);
    }
    restart=E->control.ala_pcg_restart_interval;
    if(restart<1 || restart>64)
        myerror(E,"ALA coupled FGMRES restart must be between 1 and 64");
    max_basis=restart+1;
    rhs=ala_block_vector_create(E,lev);
    r=ala_block_vector_create(E,lev);
    explicit_r=ala_block_vector_create(E,lev);
    w=ala_block_vector_create(E,lev);
    operator_work=ala_block_vector_create(E,lev);
    preconditioner_delta=ala_block_vector_create(E,lev);
    preconditioner_vanka_delta=ala_block_vector_create(E,lev);
    vb=(struct ala_block_vector **)calloc(
        max_basis,sizeof(struct ala_block_vector *));
    zb=(struct ala_block_vector **)calloc(
        restart,sizeof(struct ala_block_vector *));
    if(vb==NULL || zb==NULL)
        myerror(E,"Unable to allocate ALA coupled FGMRES basis tables");
    for(j=0;j<max_basis;j++)
        vb[j]=ala_block_vector_create(E,lev);
    for(j=0;j<restart;j++)
        zb[j]=ala_block_vector_create(E,lev);

    ala_block_vector_zero(E,rhs);
    assemble_ala_pressure_independent_force(
        E,rhs->velocity,operator_work->velocity,lev);
    force_norm=sqrt(max(global_vdot(
        E,rhs->velocity,rhs->velocity,lev),1.0e-300));
    assemble_ala_coupled_residual(
        E,V,P,rhs,r,w,operator_work,lev);
    ala_block_vector_component_norms(
        E,r,&velocity_component,&pressure_component);
    /* The augmented momentum residual already satisfies the configured
     * velocity solve to about 1e-6 of the external force.  Scaling it by its
     * own tiny initial value would make that hard tail half of the Arnoldi
     * objective forever.  Use the physical force scale instead; normalize
     * the active algebraic continuity residual independently. */
    velocity_weight=1.0/max(force_norm*force_norm,1.0e-300);
    pressure_weight=1.0/max(pressure_component*pressure_component,1.0e-300);
    initial_block_norm=ala_block_vector_norm(
        E,r,velocity_weight,pressure_weight);
    if(initial_block_norm<=1.0e-300)
        initial_block_norm=1.0;

    assemble_div_u(E,V,operator_work->pressure,lev);
    strict_ala_continuity_metrics(
        E,V,w->pressure,operator_work->pressure,lev,
        &initial_mass_norm,&cancellation_l2);
    mass_norm=initial_mass_norm;
    best_cancellation=cancellation_l2;
    initial_velocity_norm=sqrt(max(global_vdot(E,V,V,lev),1.0e-300));
    strict_ala_momentum_residual_audit(
        E,V,P,w->velocity,operator_work->velocity,lev,
        &momentum_rms,&momentum_relative);
    momentum_gate=(E->control.ala_unaugmented_momentum_tolerance>0.0);
    continuity_converged=(cancellation_l2<E->control.tole_comp);
    momentum_converged=(!momentum_gate || momentum_relative<=
        E->control.ala_unaugmented_momentum_tolerance);
    converged=(continuity_converged && momentum_converged);
    count=0;
    breakdown=0;
    residual_est=1.0;
    explicit_norm=initial_block_norm;
    block_relative=1.0;
    augmented_momentum_relative=velocity_component/force_norm;
    preconditioner_name=E->control.ala_coupled_multilevel_vcycle
        ? "coupled_multilevel_vcycle"
        : (E->control.ala_coupled_element_vanka
           ? "coupled_element_vanka"
           : (E->control.ala_coupled_shallow_vanka_sweeps>0
              ? "upper_block_triangular_plus_shallow_vanka"
              : "upper_block_triangular"));
    if(E->parallel.me==0) {
        fprintf(E->fp,"ALA COUPLED FGMRES startup restart=%d "
                "preconditioner=%s "
                "defect_corrections=%d "
                "inner_relative_tolerance=%e inner_max_cycles=%d "
                "inner_progress_interval=%d "
                "metric_weights=(velocity_force:%e,pressure_algebraic:%e) "
                "block_norm=%e momentum_component=%e "
                "continuity_component=%e cancellation=%e "
                "raw_momentum_relative=%e\n",restart,preconditioner_name,
                E->control.ala_coupled_defect_corrections,
                E->control.ala_coupled_inner_relative_tolerance,
                E->control.ala_coupled_inner_max_cycles,
                E->control.ala_coupled_inner_progress_interval,
                velocity_weight,pressure_weight,initial_block_norm,velocity_component,
                pressure_component,cancellation_l2,momentum_relative);
        fprintf(stderr,"ALA COUPLED FGMRES startup restart=%d "
                "preconditioner=%s "
                "defect_corrections=%d "
                "inner_relative_tolerance=%e inner_max_cycles=%d "
                "block_norm=%e cancellation=%e "
                "raw_momentum_relative=%e\n",restart,preconditioner_name,
                E->control.ala_coupled_defect_corrections,
                E->control.ala_coupled_inner_relative_tolerance,
                E->control.ala_coupled_inner_max_cycles,initial_block_norm,
                cancellation_l2,momentum_relative);
        fflush(E->fp);
        fflush(stderr);
    }
    strict_ala_depth_diagnostics(
        E,w->pressure,operator_work->pressure,lev,0);
    strict_ala_coarse_residual_diagnostics(
        E,r->pressure,lev,0);

    while(count<*steps_max && !converged) {
        for(j=0;j<65;j++) {
            g[j]=0.0;
            for(i=0;i<64;i++)
                h[j][i]=0.0;
        }
        for(j=0;j<64;j++) {
            cs[j]=0.0;
            sn[j]=0.0;
            y[j]=0.0;
            y_old[j]=0.0;
        }
        beta=ala_block_vector_norm(E,r,velocity_weight,pressure_weight);
        if(!isfinite(beta) || beta<=1.0e-300)
            break;
        g[0]=beta;
        ala_block_vector_copy(E,r,vb[0]);
        ala_block_vector_scale(E,1.0/beta,vb[0]);
        used=0;
        for(j=0;j<restart && count<*steps_max;j++) {
            apply_ala_coupled_block_preconditioner(
                E,vb[j],zb[j],preconditioner_delta,explicit_r,
                preconditioner_vanka_delta,pressure_work,
                operator_work,lev,count,cache);
            apply_ala_coupled_operator(
                E,zb[j]->velocity,zb[j]->pressure,w->velocity,w->pressure,
                operator_work->velocity,lev);
            if(count==0 || j==restart-1 ||
               (E->control.ala_coupled_debug_stop_iteration>0 &&
                count+1==E->control.ala_coupled_debug_stop_iteration))
                strict_ala_coupled_preconditioner_audit(
                    E,vb[j],zb[j],w,velocity_weight,pressure_weight,count+1);
            if(count==0 &&
               E->control.ala_coupled_first_preconditioner_audit_only) {
                if(E->parallel.me==0 && E->fp!=NULL) {
                    fprintf(E->fp,
                            "ALA COUPLED FIRST PRECONDITIONER AUDIT COMPLETE: "
                            "terminating before Arnoldi orthogonalization\n");
                    fflush(E->fp);
                }
                if(E->parallel.me==0) {
                    fprintf(stderr,
                            "ALA COUPLED FIRST PRECONDITIONER AUDIT COMPLETE: "
                            "terminating before Arnoldi orthogonalization\n");
                    fflush(stderr);
                }
                MPI_Barrier(E->parallel.world);
                MPI_Finalize();
                exit(EXIT_SUCCESS);
            }
            for(i=0;i<=j;i++) {
                h[i][j]=ala_block_vector_dot(
                    E,w,vb[i],velocity_weight,pressure_weight);
                ala_block_vector_axpy(E,-h[i][j],vb[i],w);
            }
            h[j+1][j]=ala_block_vector_norm(
                E,w,velocity_weight,pressure_weight);
            if(h[j+1][j]>1.0e-300) {
                ala_block_vector_copy(E,w,vb[j+1]);
                ala_block_vector_scale(E,1.0/h[j+1][j],vb[j+1]);
            }
            else {
                breakdown=1;
                h[j+1][j]=0.0;
            }
            for(i=0;i<j;i++) {
                sum=cs[i]*h[i][j]+sn[i]*h[i+1][j];
                h[i+1][j]=-sn[i]*h[i][j]+cs[i]*h[i+1][j];
                h[i][j]=sum;
            }
            norm=sqrt(h[j][j]*h[j][j]+h[j+1][j]*h[j+1][j]);
            if(norm<=1.0e-300) {
                cs[j]=1.0;
                sn[j]=0.0;
            }
            else {
                cs[j]=h[j][j]/norm;
                sn[j]=h[j+1][j]/norm;
            }
            h[j][j]=cs[j]*h[j][j]+sn[j]*h[j+1][j];
            h[j+1][j]=0.0;
            sum=cs[j]*g[j]+sn[j]*g[j+1];
            g[j+1]=-sn[j]*g[j]+cs[j]*g[j+1];
            g[j]=sum;
            used=j+1;
            for(i=used-1;i>=0;i--) {
                sum=g[i];
                for(e=i+1;e<used;e++)
                    sum-=h[i][e]*y[e];
                if(fabs(h[i][i])<=1.0e-300) {
                    breakdown=1;
                    y[i]=0.0;
                }
                else
                    y[i]=sum/h[i][i];
            }
            for(i=0;i<used;i++) {
                delta=y[i]-y_old[i];
                if(delta==0.0)
                    continue;
                for(m=1;m<=E->sphere.caps_per_proc;m++) {
                    for(e=0;e<neq;e++)
                        V[m][e] += delta*zb[i]->velocity[m][e];
                    for(e=1;e<=npno;e++)
                        P[m][e] += delta*zb[i]->pressure[m][e];
                }
                y_old[i]=y[i];
            }

            assemble_ala_coupled_residual(
                E,V,P,rhs,explicit_r,w,operator_work,lev);
            explicit_norm=ala_block_vector_norm(
                E,explicit_r,velocity_weight,pressure_weight);
            block_relative=explicit_norm/initial_block_norm;
            residual_est=fabs(g[j+1])/initial_block_norm;
            ala_block_vector_component_norms(
                E,explicit_r,&velocity_component,&pressure_component);
            augmented_momentum_relative=velocity_component/force_norm;
            assemble_div_u(E,V,operator_work->pressure,lev);
            strict_ala_continuity_metrics(
                E,V,w->pressure,operator_work->pressure,lev,
                &mass_norm,&cancellation_l2);
            velocity_norm=sqrt(max(global_vdot(E,V,V,lev),0.0));
            count++;
            if(cancellation_l2<best_cancellation)
                best_cancellation=cancellation_l2;
            continuity_converged=(cancellation_l2<E->control.tole_comp);
            if(continuity_converged && momentum_gate) {
                strict_ala_momentum_residual_audit(
                    E,V,P,w->velocity,operator_work->velocity,lev,
                    &momentum_rms,&momentum_relative);
                momentum_converged=(momentum_relative<=
                    E->control.ala_unaugmented_momentum_tolerance);
            }
            else if(!momentum_gate)
                momentum_converged=1;
            converged=(continuity_converged && momentum_converged);
            if(E->parallel.me==0) {
                fprintf(E->fp,"ALA COUPLED FGMRES iteration=%d "
                        "block_relative=%e arnoldi_relative=%e drift=%e "
                        "augmented_momentum_relative=%e "
                        "continuity_mass=%e continuity_component=%e "
                        "cancellation=%e velocity_relative=%e "
                        "raw_momentum_relative_last_audit=%e status=%s\n",count,
                        block_relative,residual_est,
                        residual_est/max(block_relative,1.0e-300),
                        augmented_momentum_relative,mass_norm,
                        pressure_component,cancellation_l2,
                        velocity_norm/max(initial_velocity_norm,1.0e-300),
                        momentum_relative,
                        converged ? "accepted" : "continue");
                fflush(E->fp);
            }
            strict_ala_depth_diagnostics(
                E,w->pressure,operator_work->pressure,lev,count);
            strict_ala_coarse_residual_diagnostics(
                E,explicit_r->pressure,lev,count);
            if(E->control.ala_coupled_debug_stop_iteration>0 &&
               count>=E->control.ala_coupled_debug_stop_iteration) {
                strict_ala_momentum_decomposition_audit(
                    E,V,P,w->velocity,operator_work->velocity,lev,
                    "coupled_debug_stop",count,
                    &momentum_rms,&momentum_relative);
                if(E->parallel.me==0 && E->fp!=NULL) {
                    fprintf(E->fp,
                            "ALA COUPLED DEBUG ITERATION STOP COMPLETE: "
                            "iteration=%d terminating before iteration=%d\n",
                            count,count+1);
                    fflush(E->fp);
                }
                if(E->parallel.me==0) {
                    fprintf(stderr,
                            "ALA COUPLED DEBUG ITERATION STOP COMPLETE: "
                            "iteration=%d terminating before iteration=%d\n",
                            count,count+1);
                    fflush(stderr);
                }
                MPI_Barrier(E->parallel.world);
                MPI_Finalize();
                exit(EXIT_SUCCESS);
            }
            if(converged || breakdown)
                break;
            for(i=0;i<used;i++)
                y_old[i]=y[i];
        }
        ala_block_vector_copy(E,explicit_r,r);
        strict_ala_momentum_decomposition_audit(
            E,V,P,w->velocity,operator_work->velocity,lev,"coupled_restart",
            count,&momentum_rms,&momentum_relative);
        momentum_converged=(!momentum_gate || momentum_relative<=
            E->control.ala_unaugmented_momentum_tolerance);
        converged=(continuity_converged && momentum_converged);
        if(E->parallel.me==0) {
            fprintf(E->fp,"ALA COUPLED FGMRES restart count=%d "
                    "block_relative=%e cancellation=%e "
                    "raw_momentum_relative=%e breakdown=%d\n",
                    count,explicit_norm/initial_block_norm,
                    cancellation_l2,momentum_relative,breakdown);
            fflush(E->fp);
        }
        if(breakdown)
            break;
    }

    strict_ala_momentum_decomposition_audit(
        E,V,P,w->velocity,operator_work->velocity,lev,"coupled_final",
        count,&momentum_rms,&momentum_relative);
    continuity_converged=(cancellation_l2<E->control.tole_comp);
    momentum_converged=(!momentum_gate || momentum_relative<=
        E->control.ala_unaugmented_momentum_tolerance);
    converged=(continuity_converged && momentum_converged);

    if(converged)
        acceptance_status="joint_target_reached";
    else if(continuity_converged && !momentum_converged)
        acceptance_status="momentum_target_not_reached";
    else
        acceptance_status="iteration_budget_exhausted";
    if(E->parallel.me==0) {
        fprintf(E->fp,"ALA_COUPLED_FEASIBILITY_SUMMARY status=%s "
                "iterations=%d cancellation=%e best=%e target=%e "
                "augmented_momentum_relative=%e raw_momentum_relative=%e "
                "momentum_target=%e arnoldi_relative=%e\n",
                acceptance_status,count,cancellation_l2,best_cancellation,
                E->control.tole_comp,augmented_momentum_relative,
                momentum_relative,E->control.ala_unaugmented_momentum_tolerance,
                residual_est);
        fflush(E->fp);
    }
    if(!converged) {
        if(E->parallel.me==0) {
            fprintf(stderr,"Strict ALA coupled FGMRES failed acceptance: "
                    "cancellation=%e momentum_relative=%e iterations=%d "
                    "breakdown=%d\n",cancellation_l2,momentum_relative,
                    count,breakdown);
            fflush(stderr);
        }
        parallel_process_termination();
    }

    *steps_max=count;
    for(j=0;j<max_basis;j++)
        ala_block_vector_destroy(E,vb[j]);
    for(j=0;j<restart;j++)
        ala_block_vector_destroy(E,zb[j]);
    free(vb);
    free(zb);
    ala_block_vector_destroy(E,rhs);
    ala_block_vector_destroy(E,r);
    ala_block_vector_destroy(E,explicit_r);
    ala_block_vector_destroy(E,w);
    ala_block_vector_destroy(E,operator_work);
    ala_block_vector_destroy(E,preconditioner_delta);
    ala_block_vector_destroy(E,preconditioner_vanka_delta);
    return((float)cancellation_l2);
}


/* Restarted flexible GMRES for the strict-ALA Schur equation.  The basis
 * stores both the preconditioned pressure vectors and their velocity
 * corrections, so the original coupled (P,V) iterate is updated with the
 * same transpose and K^{-1} actions used by the PCG path. */
static float solve_ala_fgmres_core(struct All_variables *E, double **V,
                                   double **P, int *steps_max, int lev,
                                   struct ala_pressure_preconditioner_cache
                                   *cache, double **r, double **explicit_r,
                                   double **div_u, double **preconditioner_work)
{
    void assemble_div_rho_u();
    void assemble_grad_rho_p();
    void strip_bcs_from_residual();
    int solve_del2_u();
    void parallel_process_termination();
    double global_pdot();
    double global_vdot();
    int m,j,i,e,count,valid,levnpno,neq,restart,used;
    int arnoldi_breakdown,converged;
    int momentum_gate,continuity_converged,momentum_converged;
    double beta,norm,inner_accuracy,relative_residual;
    double algebraic_relative,residual_drift;
    double cancellation_l2,mass_norm,initial_mass_norm,mass_relative;
    double term_strength,initial_term_strength,velocity_norm;
    double initial_velocity_norm;
    double momentum_rms,momentum_relative;
    double initial_rnorm,residual_est,residual,delta;
    double audit_best_cancellation;
    double sum,explicit_norm;
    double h[65][64],cs[64],sn[64],g[65],y[64],y_old[64];
    double **w,**tmpF,**tmpU;
    double ***vb,***zb,***ub;
    int max_basis;
    const char *acceptance_status;

    levnpno=E->lmesh.NPNO[lev];
    neq=E->lmesh.NEQ[lev];
    restart=E->control.ala_pcg_restart_interval;
    if(restart<1 || restart>64)
        myerror(E,"ALA FGMRES restart interval must be between 1 and 64");
    max_basis=restart+1;
    vb=(double ***)calloc(max_basis,sizeof(double **));
    zb=(double ***)calloc(restart,sizeof(double **));
    ub=(double ***)calloc(restart,sizeof(double **));
    w=(double **)calloc(NCS,sizeof(double *));
    tmpF=(double **)calloc(NCS,sizeof(double *));
    tmpU=(double **)calloc(NCS,sizeof(double *));
    if(vb==NULL || zb==NULL || ub==NULL || w==NULL || tmpF==NULL ||
       tmpU==NULL)
        myerror(E,"Unable to allocate ALA FGMRES basis tables");
    for(j=0;j<max_basis;j++) {
        vb[j]=(double **)calloc(NCS,sizeof(double *));
        if(vb[j]==NULL)
            myerror(E,"Unable to allocate ALA FGMRES pressure basis");
    }
    for(j=0;j<restart;j++) {
        zb[j]=(double **)calloc(NCS,sizeof(double *));
        ub[j]=(double **)calloc(NCS,sizeof(double *));
        if(zb[j]==NULL || ub[j]==NULL)
            myerror(E,"Unable to allocate ALA FGMRES correction basis");
    }
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        w[m]=(double *)calloc(levnpno+1,sizeof(double));
        /*
         * Nodal velocity assembly clears entries 0..neq and writes the
         * constrained zero sentinel at neq.  FGMRES velocity work vectors
         * therefore follow the legacy neq+1 velocity-vector convention.
         */
        tmpF[m]=(double *)calloc(neq+1,sizeof(double));
        tmpU[m]=(double *)calloc(neq+1,sizeof(double));
        if(w[m]==NULL || tmpF[m]==NULL || tmpU[m]==NULL)
            myerror(E,"Unable to allocate ALA FGMRES operator workspace");
        for(j=0;j<max_basis;j++) {
            vb[j][m]=(double *)calloc(levnpno+1,sizeof(double));
            if(vb[j][m]==NULL)
                myerror(E,"Unable to allocate ALA FGMRES Krylov basis");
        }
        for(j=0;j<restart;j++) {
            zb[j][m]=(double *)calloc(levnpno+1,sizeof(double));
            ub[j][m]=(double *)calloc(neq+1,sizeof(double));
            if(zb[j][m]==NULL || ub[j][m]==NULL)
                myerror(E,"Unable to allocate ALA FGMRES flexible basis");
        }
    }

    count=0;
    initial_mass_norm=0.0;
    assemble_div_rho_u(E,V,r,lev);
    initial_rnorm=sqrt(global_pdot(E,r,r,lev));
    assemble_div_rho_u(E,V,explicit_r,lev);
    strict_ala_continuity_metrics(E,V,r,div_u,lev,&initial_mass_norm,
                                  &cancellation_l2);
    initial_term_strength = strict_ala_continuity_term_strength(
        E,r,div_u,lev);
    term_strength = initial_term_strength;
    initial_velocity_norm = sqrt(max(global_vdot(E,V,V,lev),0.0));
    if(initial_rnorm<=1.0e-300)
        initial_rnorm=1.0;
    mass_norm=initial_mass_norm;
    mass_relative=(initial_mass_norm>0.0) ? 1.0 : 0.0;
    residual=cancellation_l2;
    residual_est=1.0;
    algebraic_relative=1.0;
    residual_drift=1.0;
    audit_best_cancellation=cancellation_l2;
    momentum_gate=(E->control.ala_unaugmented_momentum_tolerance>0.0);
    continuity_converged=(cancellation_l2<E->control.tole_comp);
    momentum_converged=!momentum_gate;
    converged=0;
    arnoldi_breakdown=0;
    if(E->parallel.me==0) {
        fprintf(E->fp,"ALA FGMRES startup restart=%d\n",restart);
        fprintf(stderr,"ALA FGMRES startup restart=%d\n",restart);
        fprintf(E->fp,"ALA FGMRES continuity iteration=0 "
                "cancellation=%e mass_norm=%e mass_relative=1.000000e+00 "
                "algebraic_relative=1.000000e+00 "
                "arnoldi_relative=1.000000e+00 drift=1.000000e+00 "
                "Q=%e Q_relative=1.000000e+00 "
                "term_strength=%e term_strength_relative=1.000000e+00 "
                "velocity_relative=1.000000e+00\n",
                cancellation_l2,mass_norm,term_strength,term_strength);
        fflush(E->fp);
        fflush(stderr);
    }
    strict_ala_momentum_decomposition_audit(
        E,V,P,tmpF,tmpU,lev,"startup",0,
        &momentum_rms,&momentum_relative);
    momentum_converged=(!momentum_gate ||
        momentum_relative<=E->control.ala_unaugmented_momentum_tolerance);
    converged=(continuity_converged && momentum_converged);
    if(E->parallel.me==0)
        fprintf(E->fp,"ALA FGMRES momentum audit restart=0 "
                "rms=%e relative=%e\n",momentum_rms,momentum_relative);
    if(E->parallel.me==0) {
        fprintf(E->fp,"ALA FGMRES joint acceptance "
                "cancellation_target=%e momentum_gate=%s "
                "momentum_target=%e\n",E->control.tole_comp,
                momentum_gate ? "on" : "off",
                E->control.ala_unaugmented_momentum_tolerance);
        fflush(E->fp);
    }

    while(count<*steps_max && !converged) {
        for(j=0;j<65;j++) {
            g[j]=0.0;
            for(i=0;i<64;i++)
                h[j][i]=0.0;
        }
        for(j=0;j<64;j++) {
            cs[j]=0.0;
            sn[j]=0.0;
            y[j]=0.0;
            y_old[j]=0.0;
        }
        beta=sqrt(global_pdot(E,r,r,lev));
        if(!isfinite(beta) || beta<=1.0e-300)
            break;
        g[0]=beta;
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=1;e<=levnpno;e++)
                vb[0][m][e]=r[m][e]/beta;
        used=0;
        for(j=0;j<restart && count<*steps_max;j++) {
            apply_ala_pressure_preconditioner(E,vb[j],zb[j],
                                              preconditioner_work,lev,count,
                                              cache);
            assemble_grad_rho_p(E,zb[j],tmpF,lev);
            inner_accuracy=strict_ala_inner_accuracy(
                E,tmpF,lev,E->control.ala_inner_accuracy_max);
            valid=solve_del2_u(E,tmpU,tmpF,inner_accuracy,lev);
            if(!valid)
                parallel_process_termination();
            strip_bcs_from_residual(E,tmpU,lev);
            assemble_div_rho_u(E,tmpU,w,lev);
            for(m=1;m<=E->sphere.caps_per_proc;m++)
                for(e=0;e<neq;e++)
                    ub[j][m][e]=tmpU[m][e];
            for(i=0;i<=j;i++) {
                h[i][j]=global_pdot(E,w,vb[i],lev);
                for(m=1;m<=E->sphere.caps_per_proc;m++)
                    for(e=1;e<=levnpno;e++)
                        w[m][e]-=h[i][j]*vb[i][m][e];
            }
            h[j+1][j]=sqrt(global_pdot(E,w,w,lev));
            if(h[j+1][j]>1.0e-300) {
                for(m=1;m<=E->sphere.caps_per_proc;m++)
                    for(e=1;e<=levnpno;e++)
                        vb[j+1][m][e]=w[m][e]/h[j+1][j];
            }
            else {
                arnoldi_breakdown=1;
                h[j+1][j]=0.0;
            }
            for(i=0;i<j;i++) {
                sum=cs[i]*h[i][j]+sn[i]*h[i+1][j];
                h[i+1][j]=-sn[i]*h[i][j]+cs[i]*h[i+1][j];
                h[i][j]=sum;
            }
            norm=sqrt(h[j][j]*h[j][j]+h[j+1][j]*h[j+1][j]);
            if(norm<=1.0e-300) {
                cs[j]=1.0;
                sn[j]=0.0;
            }
            else {
                cs[j]=h[j][j]/norm;
                sn[j]=h[j+1][j]/norm;
            }
            h[j][j]=cs[j]*h[j][j]+sn[j]*h[j+1][j];
            h[j+1][j]=0.0;
            sum=cs[j]*g[j]+sn[j]*g[j+1];
            g[j+1]=-sn[j]*g[j]+cs[j]*g[j+1];
            g[j]=sum;
            used=j+1;
            for(i=used-1;i>=0;i--) {
                sum=g[i];
                for(e=i+1;e<used;e++)
                    sum-=h[i][e]*y[e];
                y[i]=sum/h[i][i];
            }
            for(i=0;i<used;i++) {
                delta=y[i]-y_old[i];
                if(delta==0.0)
                    continue;
                for(m=1;m<=E->sphere.caps_per_proc;m++) {
                    for(e=1;e<=levnpno;e++)
                        P[m][e]+=delta*zb[i][m][e];
                    for(e=0;e<neq;e++)
                        V[m][e]-=delta*ub[i][m][e];
                }
                y_old[i]=y[i];
            }
            assemble_div_rho_u(E,V,explicit_r,lev);
            strict_ala_continuity_metrics(E,V,explicit_r,div_u,lev,
                                          &mass_norm,&cancellation_l2);
            term_strength = strict_ala_continuity_term_strength(
                E,explicit_r,div_u,lev);
            velocity_norm = sqrt(max(global_vdot(E,V,V,lev),0.0));
            mass_relative=(initial_mass_norm>0.0)
                ? mass_norm/initial_mass_norm : 0.0;
            residual=cancellation_l2;
            explicit_norm=sqrt(global_pdot(
                E,explicit_r,explicit_r,lev));
            algebraic_relative=explicit_norm/initial_rnorm;
            residual_est=fabs(g[j+1])/initial_rnorm;
            residual_drift=residual_est
                /max(algebraic_relative,1.0e-300);
            count++;
            if(cancellation_l2<audit_best_cancellation)
                audit_best_cancellation=cancellation_l2;
            if(E->parallel.me==0) {
                fprintf(E->fp,"ALA FGMRES continuity iteration=%d "
                        "cancellation=%e mass_norm=%e mass_relative=%e "
                        "algebraic_relative=%e arnoldi_relative=%e "
                        "drift=%e Q=%e Q_relative=%e "
                        "term_strength=%e term_strength_relative=%e "
                        "velocity_relative=%e\n",count,
                        cancellation_l2,mass_norm,mass_relative,
                        algebraic_relative,residual_est,residual_drift,
                        term_strength,
                        term_strength/max(initial_term_strength,1.0e-300),
                        term_strength,
                        term_strength/max(initial_term_strength,1.0e-300),
                        velocity_norm/max(initial_velocity_norm,1.0e-300));
                fflush(E->fp);
            }
            strict_ala_depth_diagnostics(E,explicit_r,div_u,lev,count);
            strict_ala_beta_causal_diagnostics(
                E,V,explicit_r,preconditioner_work,div_u,lev,count);
            strict_ala_coarse_residual_diagnostics(E,explicit_r,lev,count);
            continuity_converged=(cancellation_l2<E->control.tole_comp);
            if(continuity_converged && momentum_gate) {
                strict_ala_momentum_residual_audit(
                    E,V,P,tmpF,tmpU,lev,&momentum_rms,&momentum_relative);
                momentum_converged=(momentum_relative<=
                    E->control.ala_unaugmented_momentum_tolerance);
                if(E->parallel.me==0) {
                    fprintf(E->fp,"ALA FGMRES acceptance audit iteration=%d "
                            "cancellation=%e target=%e momentum_rms=%e "
                            "momentum_relative=%e momentum_target=%e "
                            "status=%s\n",count,cancellation_l2,
                            E->control.tole_comp,momentum_rms,
                            momentum_relative,
                            E->control.ala_unaugmented_momentum_tolerance,
                            momentum_converged ? "accepted"
                                               : "continue");
                    fflush(E->fp);
                }
            }
            else if(!momentum_gate)
                momentum_converged=1;
            converged=(continuity_converged && momentum_converged);
            if(converged || arnoldi_breakdown)
                break;
            for(i=0;i<used;i++)
                y_old[i]=y[i];
        }
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=1;e<=levnpno;e++)
                r[m][e]=explicit_r[m][e];
        strict_ala_momentum_decomposition_audit(
            E,V,P,tmpF,tmpU,lev,"restart",count,
            &momentum_rms,&momentum_relative);
        if(E->parallel.me==0)
            fprintf(E->fp,"ALA FGMRES restart cycle completed count=%d "
                    "residual=%e arnoldi_breakdown=%d "
                    "momentum_rms=%e momentum_relative=%e\n",count,
                    residual,arnoldi_breakdown,momentum_rms,
                    momentum_relative);
        if(arnoldi_breakdown)
            break;
    }
    explicit_norm=sqrt(global_pdot(E,explicit_r,explicit_r,lev));
    relative_residual=explicit_norm/initial_rnorm;
    continuity_converged=(residual<E->control.tole_comp);
    momentum_converged=(!momentum_gate ||
        momentum_relative<=E->control.ala_unaugmented_momentum_tolerance);
    if(converged)
        acceptance_status=momentum_gate ? "joint_target_reached"
                                        : "discrete_target_reached";
    else if(continuity_converged && !momentum_converged)
        acceptance_status="momentum_target_not_reached";
    else
        acceptance_status="iteration_budget_exhausted";
    if(E->parallel.me==0) {
        fprintf(E->fp,"ALA FGMRES operator audit restarts=%d iterations=%d "
                "basis=%d final=%e best=%e algebraic_relative=%e "
                "arnoldi_relative=%e drift=%e\n",
                (count+restart-1)/restart,count,restart,residual,
                audit_best_cancellation,relative_residual,residual_est,
                residual_est/max(relative_residual,1.0e-300));
        fprintf(E->fp,"ALA_FEASIBILITY_SUMMARY status=%s final=%e "
                "best=%e target=%e iterations=%d mass_norm=%e "
                "mass_relative=%e Q=%e Q_relative=%e "
                "momentum_relative=%e momentum_target=%e\n",
                acceptance_status,
                residual,audit_best_cancellation,E->control.tole_comp,count,
                mass_norm,mass_relative,term_strength,
                term_strength/max(initial_term_strength,1.0e-300),
                momentum_relative,
                E->control.ala_unaugmented_momentum_tolerance);
        fflush(E->fp);
    }
    if(!converged) {
        if(E->parallel.me==0) {
            fprintf(E->fp,"Strict ALA FGMRES failed joint acceptance: "
                    "cancellation=%e tolerance=%e momentum_relative=%e "
                    "momentum_tolerance=%e iterations=%d "
                    "arnoldi_breakdown=%d\n",residual,
                    E->control.tole_comp,momentum_relative,
                    E->control.ala_unaugmented_momentum_tolerance,count,
                    arnoldi_breakdown);
            fprintf(stderr,"Strict ALA FGMRES failed joint acceptance: "
                    "cancellation=%e tolerance=%e momentum_relative=%e "
                    "momentum_tolerance=%e iterations=%d "
                    "arnoldi_breakdown=%d\n",residual,
                    E->control.tole_comp,momentum_relative,
                    E->control.ala_unaugmented_momentum_tolerance,count,
                    arnoldi_breakdown);
            fflush(E->fp);
        }
        parallel_process_termination();
    }
    *steps_max=count;
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        free((void *)w[m]);
        free((void *)tmpF[m]);
        free((void *)tmpU[m]);
        for(j=0;j<max_basis;j++)
            free((void *)vb[j][m]);
        for(j=0;j<restart;j++) {
            free((void *)zb[j][m]);
            free((void *)ub[j][m]);
        }
    }
    for(j=0;j<max_basis;j++)
        free((void *)vb[j]);
    for(j=0;j<restart;j++) {
        free((void *)zb[j]);
        free((void *)ub[j]);
    }
    free((void *)vb);
    free((void *)zb);
    free((void *)ub);
    free((void *)w);
    free((void *)tmpF);
    free((void *)tmpU);
    return((float)residual);
}


static float solve_Ahat_p_fhat_iterCG(struct All_variables *E,
                                      double **V, double **P, double **F,
                                      double imp, int *steps_max);
static double initial_vel_residual(struct All_variables *E,
                                   double **V, double **P, double **F,
                                   double imp);
static double incompressibility_residual(struct All_variables *E,
                                         double **V, double **r);
static void strict_ala_continuity_metrics(struct All_variables *E,
                                          double **V, double **r,
                                          double **div_u, int lev,
                                          double *mass_norm,
                                          double *cancellation_l2);
static void strict_ala_depth_diagnostics(struct All_variables *E,
                                         double **r, double **div_u,
                                         int lev, int iteration);
static void strict_ala_beta_causal_diagnostics(
    struct All_variables *E, double **V, double **active_r,
    double **alternate_r, double **div_u, int lev, int iteration);
static void strict_ala_coarse_residual_diagnostics(
    struct All_variables *E, double **r, int lev, int iteration);
static double strict_ala_inner_relative_accuracy(struct All_variables *E,
                                                 double outer_relative_residual);
static double strict_ala_inner_accuracy(struct All_variables *E,
                                        double **F, int lev,
                                        double relative_accuracy);
#define ALA_PATCH_RADIAL_ELEMENTS 2
#define ALA_PATCH_RADIAL_STRIDE 1
#define ALA_PATCH_MPI_FACES 4
#define ALA_PATCH_MAX_HORIZONTAL_ELEMENTS 8
#define ALA_PATCH_MAX_MPI_OVERLAP 4
#define ALA_HALO_ELEMENT_RECORD_BASE_DOUBLES 98
#define ALA_HALO_ELEMENT_RECORD_DOUBLES \
    (ALA_HALO_ELEMENT_RECORD_BASE_DOUBLES+ALA_VANKA_CHOL_SIZE)
#define ALA_COARSE_POWER_ITERATIONS 6
#define ALA_VELOCITY_POWER_ITERATIONS 6
#define ALA_GLOBAL_ANGULAR_BASIS 16
#define ALA_GLOBAL_RADIAL_BASIS 6
#define ALA_GLOBAL_BASIS_COUNT \
    (ALA_GLOBAL_ANGULAR_BASIS*ALA_GLOBAL_RADIAL_BASIS)
#define ALA_GENEO_MAX_BINS 256
#define ALA_PATCH_MAX_ELEMENTS \
    (ALA_PATCH_MAX_HORIZONTAL_ELEMENTS \
     *ALA_PATCH_MAX_HORIZONTAL_ELEMENTS \
     *ALA_PATCH_RADIAL_ELEMENTS)
struct ala_pressure_preconditioner_cache {
    int patch_capacity;
    int interface_capacity;
    int blocks[NCS];
    unsigned char *size[NCS];
    unsigned short *multiplicity[NCS];
    int *elements[NCS];
    double *chol[NCS];
    int interface_blocks[NCS];
    unsigned char *interface_size[NCS];
    unsigned char *interface_face[NCS];
    int *interface_elements[NCS];
    double *interface_chol[NCS];
    int halo_send_count[NCS][ALA_PATCH_MPI_FACES];
    int halo_recv_count[NCS][ALA_PATCH_MPI_FACES];
    int *halo_send_elements[NCS][ALA_PATCH_MPI_FACES];
    double *halo_send_records[NCS][ALA_PATCH_MPI_FACES];
    double *halo_recv_records[NCS][ALA_PATCH_MPI_FACES];
    unsigned short *halo_multiplicity[NCS][ALA_PATCH_MPI_FACES];
    double *velocity_block_chol[NCS];
    int velocity_block_fallbacks;
    double *coarse_bpi[NCS];
    double *global_basis[NCS];
    double *global_matrix;
    double *global_chol;
    int global_basis_count;
    double *geneo_basis[NCS];
    double *geneo_matrix;
    double *geneo_chol;
    double *geneo_eigenvalues;
    int geneo_basis_count;
    int geneo_local_modes;
    int geneo_local_offset;
    int geneo_schur_applications;
    double coarse_eigenvalue_min;
    double coarse_eigenvalue_max;
    double velocity_eigenvalue_min;
    double velocity_eigenvalue_max;
    double *pressure_mg_rhs[MAX_LEVELS][NCS];
    double *pressure_mg_x[MAX_LEVELS][NCS];
    double *pressure_mg_residual[MAX_LEVELS][NCS];
    double *pressure_mg_Ax[MAX_LEVELS][NCS];
    double *pressure_mg_velocity[NCS];
    double *pressure_mg_velocity_rhs[NCS];
    double *pressure_mg_galerkin_p[NCS];
    double *pressure_mg_galerkin_Ap[NCS];
    double *pressure_depth_weighted_r[NCS];
    int pressure_mg_min_level;
    int pressure_mg_max_level;
    int pressure_mg_operator_applications;
};
static void build_ala_shallow_patch_cache(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache, int lev);
static void exchange_ala_shallow_halo_values(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache, int lev, int m,
    double **local, double *ghost[ALA_PATCH_MPI_FACES], int return_to_owner);
static void build_ala_two_level_cache(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache, int lev);
static void calibrate_ala_two_level_spectrum(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache, int lev);
static void calibrate_ala_velocity_spectrum(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache, int lev);
static void build_ala_pressure_multigrid_cache(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache, int lev);
static void build_ala_global_coarse_cache(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache, int lev);
static void build_ala_geneo_coarse_cache(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache, int lev);
static void build_ala_cross_rank_geneo_coarse_cache(
    struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache, int lev);
static void apply_ala_geneo_correction(struct All_variables *E,
    double **r, double **z, int lev, int iteration,
    struct ala_pressure_preconditioner_cache *cache);
static void free_ala_pressure_preconditioner_cache(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache);
static void apply_ala_geneo_correction(struct All_variables *E,
    double **r, double **z, int lev, int iteration,
    struct ala_pressure_preconditioner_cache *cache)
{
    int m,e,ex,ey,ez,cx,cy,cz,ce,i,j,n,clev,factor;
    int celx,celz,cnpno,stride;
    double sum,correction,global_energy;
    double residual2,rhs2;
    double *coarse_rhs[NCS];
    double *local_rhs,*global_rhs,*solution;

    n=cache->geneo_basis_count;
    if(n<=0)
        myerror(E,"ALA GenEO correction has no basis");
    clev=lev-E->control.ala_two_level_offset;
    factor=1 << E->control.ala_two_level_offset;
    celx=E->lmesh.ELX[clev];
    celz=E->lmesh.ELZ[clev];
    cnpno=E->lmesh.NPNO[clev];
    stride=cnpno+1;
    local_rhs=(double *)calloc(n,sizeof(double));
    global_rhs=(double *)calloc(n,sizeof(double));
    solution=(double *)calloc(n,sizeof(double));
    if(local_rhs==NULL || global_rhs==NULL || solution==NULL)
        myerror(E,"Unable to allocate ALA GenEO solve vectors");
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        coarse_rhs[m]=(double *)calloc(cnpno+1,sizeof(double));
        if(coarse_rhs[m]==NULL)
            myerror(E,"Unable to allocate ALA GenEO restricted residual");
    }
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(ey=1;ey<=E->lmesh.ELY[lev];ey++)
            for(ex=1;ex<=E->lmesh.ELX[lev];ex++)
                for(ez=1;ez<=E->lmesh.ELZ[lev];ez++) {
                    e=ez+(ex-1)*E->lmesh.ELZ[lev]
                        +(ey-1)*E->lmesh.ELZ[lev]*E->lmesh.ELX[lev];
                    cx=(ex-1)/factor+1;
                    cy=(ey-1)/factor+1;
                    cz=(ez-1)/factor+1;
                    ce=cz+(cx-1)*celz+(cy-1)*celz*celx;
                    coarse_rhs[m][ce] += r[m][e];
                }
    /* A cross-rank aggregate has the same coarse column on every rank in
       its processor group.  Sum every locally supported column; columns
       outside the group are stored as exact zeros. */
    for(i=0;i<n;i++)
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(ce=1;ce<=cnpno;ce++)
                local_rhs[i] += cache->geneo_basis[m][i*stride+ce]
                    *coarse_rhs[m][ce];
    MPI_Allreduce(local_rhs,global_rhs,n,MPI_DOUBLE,MPI_SUM,
                  E->parallel.world);
    for(i=0;i<n;i++) {
        sum=global_rhs[i];
        for(j=0;j<i;j++)
            sum -= cache->geneo_chol[i*n+j]*solution[j];
        solution[i]=sum/cache->geneo_chol[i*n+i];
    }
    for(i=n-1;i>=0;i--) {
        sum=solution[i];
        for(j=i+1;j<n;j++)
            sum -= cache->geneo_chol[j*n+i]*solution[j];
        solution[i]=sum/cache->geneo_chol[i*n+i];
    }
    global_energy=0.0;
    for(i=0;i<n;i++)
        global_energy += global_rhs[i]*solution[i];
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(ey=1;ey<=E->lmesh.ELY[lev];ey++)
            for(ex=1;ex<=E->lmesh.ELX[lev];ex++)
                for(ez=1;ez<=E->lmesh.ELZ[lev];ez++) {
                    e=ez+(ex-1)*E->lmesh.ELZ[lev]
                        +(ey-1)*E->lmesh.ELZ[lev]*E->lmesh.ELX[lev];
                    cx=(ex-1)/factor+1;
                    cy=(ey-1)/factor+1;
                    cz=(ez-1)/factor+1;
                    ce=cz+(cx-1)*celz+(cy-1)*celz*celx;
                    correction=0.0;
                    for(i=0;i<n;i++)
                        correction += cache->geneo_basis[m][i*stride+ce]
                            *solution[i];
                    z[m][e] += E->control.ala_geneo_weight*correction;
                }
    if(iteration==0 ||
       iteration%E->control.ala_coarse_residual_interval==0) {
        residual2=0.0;
        rhs2=0.0;
        for(i=0;i<n;i++) {
            sum=0.0;
            for(j=0;j<n;j++)
                sum += 0.5*(cache->geneo_matrix[i*n+j]
                            +cache->geneo_matrix[j*n+i])*solution[j];
            residual2 += (global_rhs[i]-sum)*(global_rhs[i]-sum);
            rhs2 += global_rhs[i]*global_rhs[i];
        }
        if(E->parallel.me==0) {
            fprintf(E->fp,"ALA_GENEO_COARSE_SOLVE iteration=%d modes=%d "
                    "basis_type=%s "
                    "projected_rhs_norm=%e residual_reduction=%e "
                    "coarse_energy=%e weight=%e\n",iteration,n,
                    E->control.ala_geneo_basis_type,sqrt(rhs2),
                    sqrt(residual2/max(rhs2,1.0e-300)),global_energy,
                    E->control.ala_geneo_weight);
            fflush(E->fp);
        }
    }
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        free((void *)coarse_rhs[m]);
    free((void *)local_rhs);
    free((void *)global_rhs);
    free((void *)solution);
}


static void apply_ala_pressure_preconditioner(struct All_variables *E,
                                              double **r, double **z,
                                              double **work, int lev,
                                              int iteration,
    struct ala_pressure_preconditioner_cache *cache);
static void audit_ala_shallow_patch_preconditioner(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache, int lev);
static void apply_ala_coarse_fixed_schur(struct All_variables *E,
                                         double **p, double **Ap,
                                         double **velocity,
                                         double **velocity_rhs,
                                         double **velocity_Ax,
                                         double **velocity_direction,
    struct ala_pressure_preconditioner_cache *cache,
                                         int lev,
                                         double *residual_reduction);
static void apply_ala_galerkin_fixed_schur(struct All_variables *E,
                                           double **coarse_p,
                                           double **coarse_Ap,
                                           double **fine_p,
                                           double **fine_Ap,
                                           double **fine_velocity,
                                           double **fine_velocity_rhs,
                                           double **fine_velocity_Ax,
                                           double **fine_velocity_direction,
    struct ala_pressure_preconditioner_cache *cache,
                                           int clev, int lev, int factor,
                                           double *residual_reduction);




/* Build a fixed SPD velocity metric at every node.  The legacy path stores
   sqrt(diag(K)); Stage 6c stores the Cholesky factor of the assembled 3x3
   same-node velocity block.  Applying L^-1 to a pressure-gradient vector
   gives a Gram representation of G Kblock^-1 G^T, so every local and
   cross-rank pressure patch remains symmetric positive semidefinite before
   its configured pressure shift. */
static int ala_build_velocity_node_factor(struct All_variables *E,
    int node, int lev, int m, double factor[9])
{
    int row,col,k,index,found,eq,loc0;
    double matrix[9],sum,pivot,maxdiag,shift;
    const int dims=E->mesh.nsd;
    const int max_eqn=14*dims;

    for(row=0;row<9;row++)
        matrix[row]=factor[row]=0.0;
    for(col=0;col<dims;col++) {
        eq=E->ID[lev][m][node].doff[col+1];
        loc0=(node-1)*max_eqn;
        found=0;
        for(k=0;k<max_eqn;k++)
            if(E->Node_map[lev][m][loc0+k]==eq) {
                index=loc0+k;
                matrix[0*dims+col]=E->Eqn_k1[lev][m][index];
                matrix[1*dims+col]=E->Eqn_k2[lev][m][index];
                matrix[2*dims+col]=E->Eqn_k3[lev][m][index];
                found=1;
                break;
            }
        if(!found)
            return(0);
    }
    for(row=0;row<dims;row++) {
        eq=E->ID[lev][m][node].doff[row+1];
        if(!isfinite(E->ALA_velocity_BI[lev][m][eq]) ||
           E->ALA_velocity_BI[lev][m][eq]<=0.0)
            return(0);
        matrix[row*dims+row]=1.0/E->ALA_velocity_BI[lev][m][eq];
        for(col=0;col<row;col++) {
            sum=0.5*(matrix[row*dims+col]+matrix[col*dims+row]);
            matrix[row*dims+col]=matrix[col*dims+row]=sum;
        }
    }
    maxdiag=max(matrix[0],max(matrix[4],matrix[8]));
    if(!isfinite(maxdiag) || maxdiag<=0.0)
        return(0);
    shift=E->control.ala_element_vanka_regularization*maxdiag;
    for(row=0;row<dims;row++)
        matrix[row*dims+row] += shift;
    for(row=0;row<dims;row++)
        for(col=0;col<=row;col++) {
            sum=matrix[row*dims+col];
            for(k=0;k<col;k++)
                sum -= factor[row*dims+k]*factor[col*dims+k];
            if(row==col) {
                pivot=sum;
                if(!isfinite(pivot) || pivot<=1.0e-14*maxdiag)
                    return(0);
                factor[row*dims+col]=sqrt(pivot);
            }
            else
                factor[row*dims+col]=sum/factor[col*dims+col];
        }
    return(1);
}


static void ala_build_velocity_block_cache(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache, int lev,
    int *local_fallbacks)
{
    int m,node,d,eq,nno;
    double *factor;
    const int dims=E->mesh.nsd;

    *local_fallbacks=0;
    nno=E->lmesh.NNO[lev];
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        cache->velocity_block_chol[m]=(double *)calloc(
            (size_t)(nno+1)*dims*dims,sizeof(double));
        if(cache->velocity_block_chol[m]==NULL)
            myerror(E,"Unable to allocate ALA shallow velocity-block cache");
        for(node=1;node<=nno;node++) {
            factor=cache->velocity_block_chol[m]+node*dims*dims;
            if(strcmp(E->control.ala_shallow_patch_velocity_solver,
                      "node_block")==0 &&
               ala_build_velocity_node_factor(E,node,lev,m,factor))
                continue;
            if(strcmp(E->control.ala_shallow_patch_velocity_solver,
                      "node_block")==0)
                (*local_fallbacks)++;
            for(d=0;d<dims;d++) {
                eq=E->ID[lev][m][node].doff[d+1];
                if(!isfinite(E->ALA_velocity_BI[lev][m][eq]) ||
                   E->ALA_velocity_BI[lev][m][eq]<=0.0)
                    myerror(E,"Invalid ALA velocity diagonal in patch cache");
                factor[d*dims+d]
                    =sqrt(1.0/E->ALA_velocity_BI[lev][m][eq]);
            }
        }
    }
}


static void ala_velocity_node_feature(
    struct ala_pressure_preconditioner_cache *cache,
    int m, int node, int dims, const double gradient[3], double feature[3])
{
    int row,col;
    double sum,*factor;

    factor=cache->velocity_block_chol[m]+node*dims*dims;
    for(row=0;row<dims;row++) {
        sum=gradient[row];
        for(col=0;col<row;col++)
            sum -= factor[row*dims+col]*feature[col];
        feature[row]=sum/factor[row*dims+row];
    }
}


static void ala_fill_halo_element_record(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache,
    int e, int lev, int m, int normal_layer, double *record)
{
    int a,d,node,p,eq,k;
    double gradient[3],feature[3];
    higher_precision *chol;
    const int ends=enodes[E->mesh.nsd];
    const int dims=E->mesh.nsd;

    for(k=0;k<ALA_HALO_ELEMENT_RECORD_DOUBLES;k++)
        record[k]=0.0;
    record[0]=(double)normal_layer;
    record[1]=1.0/E->BPI[lev][m][e];
    for(a=1;a<=ends;a++) {
        node=E->IEN[lev][m][e].node[a];
        for(d=1;d<=dims;d++) {
            p=(a-1)*dims+d-1;
            record[2+p]=E->X[lev][m][d][node];
            gradient[d-1]=E->elt_del[lev][m][e].g[p][0];
            if(E->control.ala_pressure_buoyancy)
                gradient[d-1] += E->elt_c[lev][m][e].c[p][0];
            eq=E->ID[lev][m][node].doff[d];
            record[26+p]=gradient[d-1];
            record[50+p]=
                strcmp(E->control.ala_shallow_patch_velocity_solver,
                       "element_vanka")==0
                ? E->ALA_vanka_overlap_BI[lev][m][eq]
                : E->ALA_velocity_BI[lev][m][eq];
        }
        ala_velocity_node_feature(cache,m,node,dims,gradient,feature);
        for(d=0;d<dims;d++)
            record[74+(a-1)*dims+d]=feature[d];
    }
    if(strcmp(E->control.ala_shallow_patch_velocity_solver,
              "element_vanka")==0) {
        if(!E->ALA_vanka_valid[lev][m][e])
            myerror(E,"ALA halo record found invalid Vanka factor");
        chol=E->ALA_vanka_chol[lev][m]+e*ALA_VANKA_CHOL_SIZE;
        for(k=0;k<ALA_VANKA_CHOL_SIZE;k++)
            record[ALA_HALO_ELEMENT_RECORD_BASE_DOUBLES+k]=(double)chol[k];
    }
}


static double ala_halo_record_coupling(struct All_variables *E,
                                        const double *left,
                                        const double *right)
{
    int a,b,d,p1,p2;
    double bi,dx,dy,dz,distance2,value;
    const double coordinate_tolerance2=1.0e-18;

    value=0.0;
    for(a=0;a<8;a++)
        for(b=0;b<8;b++) {
            p1=3*a;
            p2=3*b;
            dx=left[2+p1]-right[2+p2];
            dy=left[3+p1]-right[3+p2];
            dz=left[4+p1]-right[4+p2];
            distance2=dx*dx+dy*dy+dz*dz;
            if(distance2>coordinate_tolerance2)
                continue;
            for(d=0;d<3;d++) {
                if(strcmp(E->control.ala_shallow_patch_velocity_solver,
                          "node_block")==0)
                    value += left[74+p1+d]*right[74+p2+d];
                else {
                    bi=0.5*(left[50+p1+d]+right[50+p2+d]);
                    value += left[26+p1+d]*bi*right[26+p2+d];
                }
            }
        }
    return(value);
}


/* Build an interface pressure Gram entry from the Vanka factors carried by
   every local-plus-ghost pressure element in this interface block.  Each
   pressure element also acts as one additive velocity patch.  Restrict both
   pressure gradients to that patch with the same inverse-square-root
   multiplicity used by the rank-local aggregate, solve its cached K_gamma
   factor, and accumulate the cross product.  Using the identical patch list
   for every matrix entry makes the MPI interface block symmetric positive
   semidefinite before regularization. */
static double ala_halo_element_vanka_coupling(
    const double records[ALA_PATCH_MAX_ELEMENTS]
                        [ALA_HALO_ELEMENT_RECORD_DOUBLES],
    int n, int left_index, int right_index)
{
    int patch,a,b,d,i,j,k;
    double rhs_left[ALA_VANKA_DOF],rhs_right[ALA_VANKA_DOF];
    double forward[ALA_VANKA_DOF],solution[ALA_VANKA_DOF];
    double dx,dy,dz,distance2,weight,sum,value,diagonal;
    const double *left,*right,*support,*chol;
    const double coordinate_tolerance2=1.0e-18;

    left=records[left_index];
    right=records[right_index];
    value=0.0;
    for(patch=0;patch<n;patch++) {
        support=records[patch];
        chol=support+ALA_HALO_ELEMENT_RECORD_BASE_DOUBLES;
        for(i=0;i<ALA_VANKA_DOF;i++)
            rhs_left[i]=rhs_right[i]=forward[i]=solution[i]=0.0;
        for(a=0;a<8;a++)
            for(b=0;b<8;b++) {
                dx=support[2+3*a]-left[2+3*b];
                dy=support[3+3*a]-left[3+3*b];
                dz=support[4+3*a]-left[4+3*b];
                distance2=dx*dx+dy*dy+dz*dz;
                if(distance2<=coordinate_tolerance2)
                    for(d=0;d<3;d++) {
                        i=3*a+d;
                        j=3*b+d;
                        weight=sqrt(max(support[50+i],0.0));
                        rhs_left[i]=weight*left[26+j];
                    }
                dx=support[2+3*a]-right[2+3*b];
                dy=support[3+3*a]-right[3+3*b];
                dz=support[4+3*a]-right[4+3*b];
                distance2=dx*dx+dy*dy+dz*dz;
                if(distance2<=coordinate_tolerance2)
                    for(d=0;d<3;d++) {
                        i=3*a+d;
                        j=3*b+d;
                        weight=sqrt(max(support[50+i],0.0));
                        rhs_right[i]=weight*right[26+j];
                    }
            }
        for(i=0;i<ALA_VANKA_DOF;i++) {
            sum=rhs_right[i];
            for(k=0;k<i;k++)
                sum -= chol[i*(i+1)/2+k]*forward[k];
            diagonal=chol[i*(i+1)/2+i];
            if(!isfinite(diagonal) || diagonal<=0.0)
                return(NAN);
            forward[i]=sum/diagonal;
        }
        for(i=ALA_VANKA_DOF-1;i>=0;i--) {
            sum=forward[i];
            for(k=i+1;k<ALA_VANKA_DOF;k++)
                sum -= chol[k*(k+1)/2+i]*solution[k];
            solution[i]=sum/chol[i*(i+1)/2+i];
        }
        for(i=0;i<ALA_VANKA_DOF;i++)
            value += rhs_left[i]*solution[i];
    }
    return(value);
}


static void ala_halo_record_center(const double *record, double center[3])
{
    int a,d;
    for(d=0;d<3;d++) {
        center[d]=0.0;
        for(a=0;a<8;a++)
            center[d] += record[2+3*a+d];
        center[d] *= 0.125;
    }
}


static int ala_factor_pressure_patch(struct All_variables *E,
    double matrix[ALA_PATCH_MAX_ELEMENTS][ALA_PATCH_MAX_ELEMENTS],
    int n, double *L, int ld, double *minimum_pivot_ratio)
{
    int i,j,k;
    double maxdiag,sum,pivot,ratio;
    const double pivot_tolerance=1.0e-12;

    maxdiag=0.0;
    for(i=0;i<n;i++) {
        matrix[i][i] *=
            1.0+E->control.ala_shallow_patch_regularization;
        maxdiag=max(maxdiag,matrix[i][i]);
    }
    if(!isfinite(maxdiag) || maxdiag<=0.0)
        return(0);
    for(i=0;i<n;i++) {
        for(j=0;j<=i;j++) {
            sum=matrix[i][j];
            for(k=0;k<j;k++)
                sum -= L[i*ld+k]*L[j*ld+k];
            if(i==j) {
                pivot=sum;
                if(!isfinite(pivot) || pivot<=pivot_tolerance*maxdiag)
                    return(0);
                ratio=pivot/maxdiag;
                *minimum_pivot_ratio=min(*minimum_pivot_ratio,ratio);
                L[i*ld+j]=sqrt(pivot);
            }
            else
                L[i*ld+j]=sum/L[j*ld+j];
        }
    }
    return(1);
}


static void ala_exchange_halo_records(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache, int lev, int m)
{
    int face,target,nrequest;
    MPI_Request requests[2*ALA_PATCH_MPI_FACES];
    MPI_Status statuses[2*ALA_PATCH_MPI_FACES];

    nrequest=0;
    for(face=0;face<ALA_PATCH_MPI_FACES;face++) {
        target=E->parallel.PROCESSOR[lev][m].pass[face+1];
        MPI_Irecv(cache->halo_recv_records[m][face],
            cache->halo_recv_count[m][face]
                *ALA_HALO_ELEMENT_RECORD_DOUBLES,
            MPI_DOUBLE,target,31901,E->parallel.world,
            &requests[nrequest++]);
        MPI_Isend(cache->halo_send_records[m][face],
            cache->halo_send_count[m][face]
                *ALA_HALO_ELEMENT_RECORD_DOUBLES,
            MPI_DOUBLE,target,31901,E->parallel.world,
            &requests[nrequest++]);
    }
    MPI_Waitall(nrequest,requests,statuses);
}


static int ala_find_halo_record(const double *records, int count,
    int normal_layer, const double target_center[3])
{
    int i,best;
    double center[3],dx,dy,dz,distance2,best_distance2;

    best=-1;
    best_distance2=1.0e300;
    for(i=0;i<count;i++) {
        if((int)(records[i*ALA_HALO_ELEMENT_RECORD_DOUBLES]+0.5)
           !=normal_layer)
            continue;
        ala_halo_record_center(
            records+i*ALA_HALO_ELEMENT_RECORD_DOUBLES,center);
        dx=center[0]-target_center[0];
        dy=center[1]-target_center[1];
        dz=center[2]-target_center[2];
        distance2=dx*dx+dy*dy+dz*dz;
        if(distance2<best_distance2) {
            best_distance2=distance2;
            best=i;
        }
    }
    if(best_distance2>4.0e-2)
        return(-1);
    return(best);
}


static double ala_velocity_block_schur_entry(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache,
    int e1, int e2, int lev, int m)
{
    int a,b,d,node1,node2,p1,p2;
    double gradient1[3],gradient2[3],feature1[3],feature2[3],value;
    const int ends=enodes[E->mesh.nsd];
    const int dims=E->mesh.nsd;

    value=0.0;
    for(a=1;a<=ends;a++) {
        node1=E->IEN[lev][m][e1].node[a];
        for(b=1;b<=ends;b++) {
            node2=E->IEN[lev][m][e2].node[b];
            if(node1!=node2)
                continue;
            for(d=0;d<dims;d++) {
                p1=(a-1)*dims+d;
                p2=(b-1)*dims+d;
                gradient1[d]=E->elt_del[lev][m][e1].g[p1][0];
                gradient2[d]=E->elt_del[lev][m][e2].g[p2][0];
                if(E->control.ala_pressure_buoyancy) {
                    gradient1[d] += E->elt_c[lev][m][e1].c[p1][0];
                    gradient2[d] += E->elt_c[lev][m][e2].c[p2][0];
                }
            }
            ala_velocity_node_feature(cache,m,node1,dims,gradient1,feature1);
            ala_velocity_node_feature(cache,m,node1,dims,gradient2,feature2);
            for(d=0;d<dims;d++)
                value += feature1[d]*feature2[d];
        }
    }
    return(value);
}


/* Couple two P0 pressure rows through the existing current-rheology
   element-Vanka inverse.  For every velocity element touching e1, restrict
   both pressure gradients to that element's Q1 velocity support and form

       (W G_e1^T)^T Kpatch^-1 (W G_e2^T).

   Summing these Gram products produces a symmetric positive-semidefinite
   pressure metric with genuine neighbour off-diagonals.  Unlike the former
   nodal diagonal and 3x3 approximations, this metric can represent the
   factor-two alternating P0 mode while reusing the already cached 24x24
   K_gamma factors. */
static double ala_element_vanka_schur_entry(struct All_variables *E,
    int e1, int e2, int lev, int m)
{
    int ex1,ey1,ez1,px,py,pz,pe,a,b,d,node,node1,node2,i,j;
    int elx,ely,elz,found1,found2;
    double rhs1[ALA_VANKA_DOF],rhs2[ALA_VANKA_DOF];
    double sol2[ALA_VANKA_DOF],weight,sum,value;
    higher_precision *chol;
    const int dims=E->mesh.nsd;
    const int ends=enodes[dims];

    elx=E->lmesh.ELX[lev];
    ely=E->lmesh.ELY[lev];
    elz=E->lmesh.ELZ[lev];
    ez1=(e1-1)%elz+1;
    ex1=((e1-1)/elz)%elx+1;
    ey1=(e1-1)/(elz*elx)+1;
    value=0.0;
    for(py=max(1,ey1-1);py<=min(ely,ey1+1);py++)
        for(px=max(1,ex1-1);px<=min(elx,ex1+1);px++)
            for(pz=max(1,ez1-1);pz<=min(elz,ez1+1);pz++) {
                pe=pz+(px-1)*elz+(py-1)*elz*elx;
                if(!E->ALA_vanka_valid[lev][m][pe])
                    myerror(E,"ALA pressure aggregate found invalid Vanka factor");
                for(i=0;i<ALA_VANKA_DOF;i++)
                    rhs1[i]=rhs2[i]=sol2[i]=0.0;
                found1=found2=0;
                for(a=1;a<=ends;a++) {
                    node=E->IEN[lev][m][pe].node[a];
                    for(b=1;b<=ends;b++) {
                        node1=E->IEN[lev][m][e1].node[b];
                        node2=E->IEN[lev][m][e2].node[b];
                        for(d=0;d<dims;d++) {
                            i=(a-1)*dims+d;
                            j=(b-1)*dims+d;
                            weight=sqrt(E->ALA_vanka_overlap_BI[lev][m]
                                [E->ID[lev][m][node].doff[d+1]]);
                            if(node==node1) {
                                rhs1[i]=weight*(E->elt_del[lev][m][e1].g[j][0]
                                    +E->elt_c[lev][m][e1].c[j][0]);
                                found1=1;
                            }
                            if(node==node2) {
                                rhs2[i]=weight*(E->elt_del[lev][m][e2].g[j][0]
                                    +E->elt_c[lev][m][e2].c[j][0]);
                                found2=1;
                            }
                        }
                    }
                }
                if(!found1 || !found2)
                    continue;
                chol=E->ALA_vanka_chol[lev][m]+pe*ALA_VANKA_CHOL_SIZE;
                if(!ala_solve_cached_element_k(chol,rhs2,sol2))
                    myerror(E,"ALA pressure aggregate Vanka solve failed");
                sum=0.0;
                for(i=0;i<ALA_VANKA_DOF;i++)
                    sum += rhs1[i]*sol2[i];
                value += sum;
            }
    return(value);
}


/* Build a pressure-space principal block of the complete strict-ALA nodal
   velocity-metric Schur approximation.  Stage 6c uses Gram operators
   G Mv^-1 G^T, so its principal submatrices are symmetric positive
   semidefinite before the configured pressure shift.  The diagonal choice
   deliberately calls the Stage-6b.3 assembly unchanged. */
static void ala_build_pressure_patch_schur(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache,
    const int *elements, int n, int lev, int m, double *schur)
{
    int i,j,e1,e2,ex1,ey1,ez1,ex2,ey2,ez2;
    int elx,elz;
    elx=E->lmesh.ELX[lev];
    elz=E->lmesh.ELZ[lev];
    for(i=0;i<n;i++) {
        e1=elements[i];
        ez1=(e1-1)%elz+1;
        ex1=((e1-1)/elz)%elx+1;
        ey1=(e1-1)/(elz*elx)+1;
        for(j=0;j<=i;j++) {
            e2=elements[j];
            ez2=(e2-1)%elz+1;
            ex2=((e2-1)/elz)%elx+1;
            ey2=(e2-1)/(elz*elx)+1;
            if(abs(ex1-ex2)>1 || abs(ey1-ey2)>1 || abs(ez1-ez2)>1)
                schur[i*ALA_PATCH_MAX_ELEMENTS+j]=0.0;
            else if(strcmp(E->control.ala_shallow_patch_velocity_solver,
                           "element_vanka")==0)
                schur[i*ALA_PATCH_MAX_ELEMENTS+j]
                    =ala_element_vanka_schur_entry(E,e1,e2,lev,m);
            else if(strcmp(E->control.ala_shallow_patch_velocity_solver,
                           "node_block")==0)
                schur[i*ALA_PATCH_MAX_ELEMENTS+j]
                    =ala_velocity_block_schur_entry(E,cache,e1,e2,lev,m);
            else
                schur[i*ALA_PATCH_MAX_ELEMENTS+j]
                    =assemble_Ahatp_jacobi_entry(E,e1,e2,lev,m);
            schur[j*ALA_PATCH_MAX_ELEMENTS+i]
                =schur[i*ALA_PATCH_MAX_ELEMENTS+j];
        }
    }
}


static void build_ala_shallow_patch_cache(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache, int lev)
{
    int m,ex,ey,ez,dx,dy,dz,b,i,j,n,block_count;
    int elx,ely,elz,shallow_min_ez,shallow_layers,e,npno,nno;
    int horizontal_elements,horizontal_stride;
    int patch_capacity,interface_capacity;
    int local_blocks,global_blocks,local_fallback,global_fallback;
    int local_elements,global_elements,local_unique,global_unique;
    int local_overlap_min,global_overlap_min;
    int local_overlap_max,global_overlap_max;
    int face,q,t,tangent,tangent_count,overlap,interface_count;
    int local_interface_blocks,global_interface_blocks;
    int local_interface_entries,global_interface_entries;
    int local_ghost_elements,global_ghost_elements,ghost_index;
    int local_velocity_fallbacks,global_velocity_fallbacks;
    int ref,boundary_coordinate;
    int nrequest,target;
    int received_counts[ALA_PATCH_MPI_FACES];
    unsigned short *ghost_count[ALA_PATCH_MPI_FACES];
    unsigned short *owner_count[ALA_PATCH_MPI_FACES];
    unsigned short *final_count[ALA_PATCH_MPI_FACES];
    double depth_km;
    double build_start,build_seconds,global_build_seconds;
    double local_cache_bytes,global_cache_bytes;
    double matrix[ALA_PATCH_MAX_ELEMENTS][ALA_PATCH_MAX_ELEMENTS];
    double patch_records[ALA_PATCH_MAX_ELEMENTS]
                        [ALA_HALO_ELEMENT_RECORD_DOUBLES];
    double target_record[ALA_HALO_ELEMENT_RECORD_DOUBLES];
    double target_center[3];
    double local_min_pivot_ratio,global_min_pivot_ratio;
    double *L;
    MPI_Request requests[2*ALA_PATCH_MPI_FACES];
    MPI_Status statuses[2*ALA_PATCH_MPI_FACES];
    double CPU_time0();

    memset(cache,0,sizeof(*cache));
    if(E->sphere.caps_per_proc!=1)
        myerror(E,"ALA MPI-overlap Schwarz requires one cap per rank");
    build_start=CPU_time0();
    elx=E->lmesh.ELX[lev];
    ely=E->lmesh.ELY[lev];
    elz=E->lmesh.ELZ[lev];
    npno=E->lmesh.NPNO[lev];
    nno=E->lmesh.NNO[lev];
    horizontal_elements=E->control.ala_shallow_patch_horizontal_elements;
    horizontal_stride=E->control.ala_shallow_patch_horizontal_stride;
    overlap=E->control.ala_shallow_patch_mpi_overlap;
    if(horizontal_elements>elx || horizontal_elements>ely)
        myerror(E,"ALA shallow-patch width exceeds the local horizontal mesh");
    if(overlap>ALA_PATCH_MAX_MPI_OVERLAP)
        myerror(E,"ALA shallow-patch MPI overlap exceeds the compiled maximum");
    if(overlap>elx || overlap>ely)
        myerror(E,"ALA shallow-patch MPI overlap exceeds the local mesh");
    patch_capacity=horizontal_elements*horizontal_elements
        *ALA_PATCH_RADIAL_ELEMENTS;
    interface_capacity=2*overlap*horizontal_elements
        *ALA_PATCH_RADIAL_ELEMENTS;
    if(patch_capacity>ALA_PATCH_MAX_ELEMENTS ||
       interface_capacity>ALA_PATCH_MAX_ELEMENTS)
        myerror(E,"ALA shallow-patch capacity exceeds the compiled maximum");
    cache->patch_capacity=patch_capacity;
    cache->interface_capacity=interface_capacity;
    shallow_min_ez=elz+1;
    for(ez=elz;ez>=1;ez--) {
        depth_km=(1.0-0.5*(E->sx[1][3][ez]+E->sx[1][3][ez+1]))
            *E->data.radius_km;
        if(depth_km>E->control.ala_shallow_patch_depth_km)
            break;
        shallow_min_ez=ez;
    }
    shallow_layers=(shallow_min_ez<=elz)
        ? elz-shallow_min_ez+1 : 0;
    block_count=((elx+horizontal_stride-1)/horizontal_stride)
        *((ely+horizontal_stride-1)/horizontal_stride)
        *((shallow_layers+ALA_PATCH_RADIAL_STRIDE-1)
          /ALA_PATCH_RADIAL_STRIDE);
    local_blocks=0;
    local_fallback=0;
    local_elements=0;
    local_unique=0;
    local_overlap_min=255;
    local_overlap_max=0;
    local_interface_blocks=0;
    local_interface_entries=0;
    local_ghost_elements=0;
    local_min_pivot_ratio=1.0e300;
    local_cache_bytes=0.0;
    ala_build_velocity_block_cache(E,cache,lev,&local_velocity_fallbacks);
    local_cache_bytes += E->sphere.caps_per_proc
        *(double)(nno+1)*E->mesh.nsd*E->mesh.nsd*sizeof(double);

    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        cache->blocks[m]=block_count;
        /* Deep radial ranks legitimately own no cells above the configured
           shallow depth.  The apply path still scans their pressure cells, so
           multiplicity must exist even when this rank has zero patch blocks. */
        cache->multiplicity[m]=(unsigned short *)calloc(npno+1,
                                                       sizeof(unsigned short));
        local_cache_bytes += (npno+1)*sizeof(unsigned short);
        if(cache->multiplicity[m]==NULL)
            myerror(E,"Unable to allocate ALA shallow-patch multiplicity");
        if(block_count==0)
            continue;
        cache->size[m]=(unsigned char *)calloc(block_count,
                                                sizeof(unsigned char));
        cache->elements[m]=(int *)calloc(
            (size_t)block_count*patch_capacity,sizeof(int));
        cache->chol[m]=(double *)calloc(
            (size_t)block_count*patch_capacity*patch_capacity,sizeof(double));
        local_cache_bytes += block_count*sizeof(unsigned char)
            +(double)block_count*patch_capacity*sizeof(int)
            +(double)block_count*patch_capacity*patch_capacity*sizeof(double);
        if(cache->size[m]==NULL || cache->elements[m]==NULL ||
           cache->chol[m]==NULL)
            myerror(E,"Unable to allocate ALA shallow-patch cache");
        b=0;
        for(ey=1;ey<=ely;ey+=horizontal_stride)
            for(ex=1;ex<=elx;ex+=horizontal_stride)
                for(ez=elz;ez>=shallow_min_ez;
                    ez-=ALA_PATCH_RADIAL_STRIDE) {
                    n=0;
                    for(dy=0;dy<horizontal_elements &&
                        ey+dy<=ely;dy++)
                        for(dx=0;dx<horizontal_elements &&
                            ex+dx<=elx;dx++)
                            for(dz=0;dz<ALA_PATCH_RADIAL_ELEMENTS &&
                                ez-dz>=shallow_min_ez;dz++) {
                                e=(ez-dz)+(ex+dx-1)*elz
                                    +(ey+dy-1)*elz*elx;
                                cache->elements[m][b*patch_capacity+n]
                                    =e;
                                n++;
                            }
                    cache->size[m][b]=(unsigned char)n;
                    for(i=0;i<n;i++)
                        for(j=0;j<n;j++)
                            matrix[i][j]=0.0;
                    ala_build_pressure_patch_schur(
                        E,cache,cache->elements[m]+b*patch_capacity,
                        n,lev,m,&matrix[0][0]);
                    L=cache->chol[m]+(size_t)b*patch_capacity*patch_capacity;
                    if(!ala_factor_pressure_patch(E,matrix,n,L,patch_capacity,
                                                  &local_min_pivot_ratio)) {
                        cache->size[m][b]=0;
                        local_fallback++;
                    }
                    else
                        local_blocks++;
                    local_elements += n;
                    b++;
                }
        for(b=0;b<block_count;b++) {
            n=cache->size[m][b];
            for(i=0;i<n;i++) {
                e=cache->elements[m][b*patch_capacity+i];
                if(cache->multiplicity[m][e]==65535)
                    myerror(E,"ALA shallow-patch overlap count overflow");
                cache->multiplicity[m][e]++;
            }
        }

        for(face=0;face<ALA_PATCH_MPI_FACES;face++) {
            tangent_count=(face<2) ? ely : elx;
            cache->halo_send_count[m][face]=
                overlap*tangent_count*shallow_layers;
            cache->halo_recv_count[m][face]=
                cache->halo_send_count[m][face];
            cache->halo_send_elements[m][face]=(int *)calloc(
                max(cache->halo_send_count[m][face],1),sizeof(int));
            cache->halo_send_records[m][face]=(double *)calloc(
                max(cache->halo_send_count[m][face],1)
                    *ALA_HALO_ELEMENT_RECORD_DOUBLES,sizeof(double));
            cache->halo_recv_records[m][face]=(double *)calloc(
                max(cache->halo_recv_count[m][face],1)
                    *ALA_HALO_ELEMENT_RECORD_DOUBLES,sizeof(double));
            cache->halo_multiplicity[m][face]=(unsigned short *)calloc(
                max(cache->halo_recv_count[m][face],1),
                sizeof(unsigned short));
            local_cache_bytes +=
                (double)max(cache->halo_send_count[m][face],1)*sizeof(int)
                +(double)max(cache->halo_send_count[m][face],1)
                    *ALA_HALO_ELEMENT_RECORD_DOUBLES*sizeof(double)
                +(double)max(cache->halo_recv_count[m][face],1)
                    *ALA_HALO_ELEMENT_RECORD_DOUBLES*sizeof(double)
                +(double)max(cache->halo_recv_count[m][face],1)
                    *sizeof(unsigned short);
            if(cache->halo_send_elements[m][face]==NULL ||
               cache->halo_send_records[m][face]==NULL ||
               cache->halo_recv_records[m][face]==NULL ||
               cache->halo_multiplicity[m][face]==NULL)
                myerror(E,"Unable to allocate ALA MPI-overlap halo");
            i=0;
            for(q=0;q<overlap;q++)
                for(t=1;t<=tangent_count;t++)
                    for(ez=shallow_min_ez;ez<=elz;ez++) {
                        if(face==0) {
                            ex=1+q;
                            ey=t;
                        }
                        else if(face==1) {
                            ex=elx-q;
                            ey=t;
                        }
                        else if(face==2) {
                            ex=t;
                            ey=1+q;
                        }
                        else {
                            ex=t;
                            ey=ely-q;
                        }
                        e=ez+(ex-1)*elz+(ey-1)*elz*elx;
                        cache->halo_send_elements[m][face][i]=e;
                        ala_fill_halo_element_record(E,cache,e,lev,m,q,
                            cache->halo_send_records[m][face]
                                +i*ALA_HALO_ELEMENT_RECORD_DOUBLES);
                        i++;
                    }
            local_ghost_elements += cache->halo_recv_count[m][face];
        }
        nrequest=0;
        for(face=0;face<ALA_PATCH_MPI_FACES;face++) {
            target=E->parallel.PROCESSOR[lev][m].pass[face+1];
            MPI_Irecv(&received_counts[face],1,MPI_INT,target,31900,
                E->parallel.world,&requests[nrequest++]);
            MPI_Isend(&cache->halo_send_count[m][face],1,MPI_INT,target,31900,
                E->parallel.world,&requests[nrequest++]);
        }
        MPI_Waitall(nrequest,requests,statuses);
        for(face=0;face<ALA_PATCH_MPI_FACES;face++)
            if(received_counts[face]!=cache->halo_recv_count[m][face]) {
                fprintf(stderr,"rank=%d ALA halo count mismatch face=%d "
                        "target=%d send=%d receive=%d expected=%d\n",
                        E->parallel.me,face+1,
                        E->parallel.PROCESSOR[lev][m].pass[face+1],
                        cache->halo_send_count[m][face],received_counts[face],
                        cache->halo_recv_count[m][face]);
                fflush(stderr);
                myerror(E,"ALA MPI halo face sizes do not match");
            }
        ala_exchange_halo_records(E,cache,lev,m);

        interface_count=0;
        for(face=0;face<ALA_PATCH_MPI_FACES;face++) {
            tangent_count=(face<2) ? ely : elx;
            interface_count +=
                ((tangent_count+horizontal_stride-1)/horizontal_stride)
                *((shallow_layers+ALA_PATCH_RADIAL_STRIDE-1)
                  /ALA_PATCH_RADIAL_STRIDE);
        }
        cache->interface_blocks[m]=interface_count;
        cache->interface_size[m]=(unsigned char *)calloc(
            max(interface_count,1),sizeof(unsigned char));
        cache->interface_face[m]=(unsigned char *)calloc(
            max(interface_count,1),sizeof(unsigned char));
        cache->interface_elements[m]=(int *)calloc(
            (size_t)max(interface_count,1)*interface_capacity,sizeof(int));
        cache->interface_chol[m]=(double *)calloc(
            (size_t)max(interface_count,1)*interface_capacity
                *interface_capacity,sizeof(double));
        local_cache_bytes +=
            2.0*max(interface_count,1)*sizeof(unsigned char)
            +(double)max(interface_count,1)*interface_capacity*sizeof(int)
            +(double)max(interface_count,1)*interface_capacity
                *interface_capacity*sizeof(double);
        if(cache->interface_size[m]==NULL ||
           cache->interface_face[m]==NULL ||
           cache->interface_elements[m]==NULL ||
           cache->interface_chol[m]==NULL)
            myerror(E,"Unable to allocate ALA MPI-overlap patch cache");
        for(face=0;face<ALA_PATCH_MPI_FACES;face++) {
            ghost_count[face]=(unsigned short *)calloc(
                max(cache->halo_recv_count[m][face],1),
                sizeof(unsigned short));
            owner_count[face]=(unsigned short *)calloc(
                max(cache->halo_send_count[m][face],1),
                sizeof(unsigned short));
            final_count[face]=(unsigned short *)calloc(
                max(cache->halo_send_count[m][face],1),
                sizeof(unsigned short));
            if(ghost_count[face]==NULL || owner_count[face]==NULL ||
               final_count[face]==NULL)
                myerror(E,"Unable to allocate ALA halo overlap counts");
        }

        b=0;
        for(face=0;face<ALA_PATCH_MPI_FACES;face++) {
            tangent_count=(face<2) ? ely : elx;
            for(tangent=1;tangent<=tangent_count;
                tangent+=horizontal_stride)
                for(ez=elz;ez>=shallow_min_ez;
                    ez-=ALA_PATCH_RADIAL_STRIDE) {
                    n=0;
                    for(q=0;q<overlap;q++)
                        for(t=0;t<horizontal_elements &&
                            tangent+t<=tangent_count;t++)
                            for(dz=0;dz<ALA_PATCH_RADIAL_ELEMENTS &&
                                ez-dz>=shallow_min_ez;dz++) {
                                if(face==0) {
                                    ex=1+q;
                                    ey=tangent+t;
                                    boundary_coordinate=1;
                                }
                                else if(face==1) {
                                    ex=elx-q;
                                    ey=tangent+t;
                                    boundary_coordinate=elx;
                                }
                                else if(face==2) {
                                    ex=tangent+t;
                                    ey=1+q;
                                    boundary_coordinate=1;
                                }
                                else {
                                    ex=tangent+t;
                                    ey=ely-q;
                                    boundary_coordinate=ely;
                                }
                                e=(ez-dz)+(ex-1)*elz
                                    +(ey-1)*elz*elx;
                                cache->interface_elements[m]
                                    [b*interface_capacity+n]=e;
                                ala_fill_halo_element_record(E,cache,e,lev,m,q,
                                    patch_records[n]);
                                n++;
                            }
                    for(q=0;q<overlap;q++)
                        for(t=0;t<horizontal_elements &&
                            tangent+t<=tangent_count;t++)
                            for(dz=0;dz<ALA_PATCH_RADIAL_ELEMENTS &&
                                ez-dz>=shallow_min_ez;dz++) {
                                if(face<2) {
                                    ex=boundary_coordinate;
                                    ey=tangent+t;
                                }
                                else {
                                    ex=tangent+t;
                                    ey=boundary_coordinate;
                                }
                                e=(ez-dz)+(ex-1)*elz
                                    +(ey-1)*elz*elx;
                                ala_fill_halo_element_record(E,cache,e,lev,m,0,
                                                             target_record);
                                ala_halo_record_center(target_record,
                                                       target_center);
                                ghost_index=ala_find_halo_record(
                                    cache->halo_recv_records[m][face],
                                    cache->halo_recv_count[m][face],q,
                                    target_center);
                                if(ghost_index<0) {
                                    fprintf(stderr,"rank=%d ALA halo match "
                                            "failed face=%d q=%d element=%d "
                                            "target=%d\n",E->parallel.me,
                                            face+1,q,e,
                                            E->parallel.PROCESSOR[lev][m]
                                              .pass[face+1]);
                                    fflush(stderr);
                                    myerror(E,"Unable to match ALA halo element");
                                }
                                cache->interface_elements[m]
                                    [b*interface_capacity+n]
                                    =-(ghost_index+1);
                                memcpy(patch_records[n],
                                    cache->halo_recv_records[m][face]
                                      +ghost_index
                                       *ALA_HALO_ELEMENT_RECORD_DOUBLES,
                                    ALA_HALO_ELEMENT_RECORD_DOUBLES
                                      *sizeof(double));
                                n++;
                            }
                    cache->interface_size[m][b]=(unsigned char)n;
                    cache->interface_face[m][b]=(unsigned char)face;
                    for(i=0;i<n;i++)
                        for(j=0;j<=i;j++) {
                            if(strcmp(E->control
                                      .ala_shallow_patch_velocity_solver,
                                      "element_vanka")==0)
                                matrix[i][j]=ala_halo_element_vanka_coupling(
                                    patch_records,n,i,j);
                            else if(strcmp(E->control
                                      .ala_shallow_patch_velocity_solver,
                                      "node_block")==0)
                                matrix[i][j]=ala_halo_record_coupling(
                                    E,patch_records[i],patch_records[j]);
                            else if(i==j)
                                matrix[i][j]=patch_records[i][1];
                            else
                                matrix[i][j]=ala_halo_record_coupling(
                                    E,patch_records[i],patch_records[j]);
                            matrix[j][i]=matrix[i][j];
                        }
                    L=cache->interface_chol[m]
                        +(size_t)b*interface_capacity*interface_capacity;
                    if(!ala_factor_pressure_patch(E,matrix,n,L,
                                                  interface_capacity,
                                                  &local_min_pivot_ratio)) {
                        cache->interface_size[m][b]=0;
                        local_fallback++;
                    }
                    else {
                        local_interface_blocks++;
                        local_interface_entries += n;
                        for(i=0;i<n;i++) {
                            ref=cache->interface_elements[m]
                                [b*interface_capacity+i];
                            if(ref>0) {
                                if(cache->multiplicity[m][ref]==65535)
                                    myerror(E,"ALA MPI overlap count overflow");
                                cache->multiplicity[m][ref]++;
                            }
                            else {
                                ghost_index=-ref-1;
                                if(ghost_count[face][ghost_index]==65535)
                                    myerror(E,"ALA ghost overlap count overflow");
                                ghost_count[face][ghost_index]++;
                            }
                        }
                    }
                    b++;
                }
        }

        nrequest=0;
        for(face=0;face<ALA_PATCH_MPI_FACES;face++) {
            target=E->parallel.PROCESSOR[lev][m].pass[face+1];
            MPI_Irecv(owner_count[face],cache->halo_send_count[m][face],
                MPI_UNSIGNED_SHORT,target,31902,E->parallel.world,
                &requests[nrequest++]);
            MPI_Isend(ghost_count[face],cache->halo_recv_count[m][face],
                MPI_UNSIGNED_SHORT,target,31902,E->parallel.world,
                &requests[nrequest++]);
        }
        MPI_Waitall(nrequest,requests,statuses);
        for(face=0;face<ALA_PATCH_MPI_FACES;face++)
            for(i=0;i<cache->halo_send_count[m][face];i++) {
                e=cache->halo_send_elements[m][face][i];
                if(65535-cache->multiplicity[m][e]<owner_count[face][i])
                    myerror(E,"ALA global overlap count overflow");
                cache->multiplicity[m][e] += owner_count[face][i];
            }
        for(face=0;face<ALA_PATCH_MPI_FACES;face++)
            for(i=0;i<cache->halo_send_count[m][face];i++) {
                e=cache->halo_send_elements[m][face][i];
                final_count[face][i]=cache->multiplicity[m][e];
            }
        nrequest=0;
        for(face=0;face<ALA_PATCH_MPI_FACES;face++) {
            target=E->parallel.PROCESSOR[lev][m].pass[face+1];
            MPI_Irecv(cache->halo_multiplicity[m][face],
                cache->halo_recv_count[m][face],MPI_UNSIGNED_SHORT,
                target,31903,E->parallel.world,&requests[nrequest++]);
            MPI_Isend(final_count[face],cache->halo_send_count[m][face],
                MPI_UNSIGNED_SHORT,target,31903,E->parallel.world,
                &requests[nrequest++]);
        }
        MPI_Waitall(nrequest,requests,statuses);
        for(face=0;face<ALA_PATCH_MPI_FACES;face++) {
            free((void *)ghost_count[face]);
            free((void *)owner_count[face]);
            free((void *)final_count[face]);
        }

        for(e=1;e<=npno;e++)
            if(cache->multiplicity[m][e]>0) {
                local_unique++;
                local_overlap_min=min(local_overlap_min,
                                      (int)cache->multiplicity[m][e]);
                local_overlap_max=max(local_overlap_max,
                                      (int)cache->multiplicity[m][e]);
            }
    }
    MPI_Allreduce(&local_blocks,&global_blocks,1,MPI_INT,MPI_SUM,
                  E->parallel.world);
    MPI_Allreduce(&local_fallback,&global_fallback,1,MPI_INT,MPI_SUM,
                  E->parallel.world);
    MPI_Allreduce(&local_elements,&global_elements,1,MPI_INT,MPI_SUM,
                  E->parallel.world);
    MPI_Allreduce(&local_unique,&global_unique,1,MPI_INT,MPI_SUM,
                  E->parallel.world);
    MPI_Allreduce(&local_overlap_min,&global_overlap_min,1,MPI_INT,MPI_MIN,
                  E->parallel.world);
    MPI_Allreduce(&local_overlap_max,&global_overlap_max,1,MPI_INT,MPI_MAX,
                  E->parallel.world);
    MPI_Allreduce(&local_interface_blocks,&global_interface_blocks,1,MPI_INT,
                  MPI_SUM,E->parallel.world);
    MPI_Allreduce(&local_interface_entries,&global_interface_entries,1,MPI_INT,
                  MPI_SUM,E->parallel.world);
    MPI_Allreduce(&local_ghost_elements,&global_ghost_elements,1,MPI_INT,
                  MPI_SUM,E->parallel.world);
    MPI_Allreduce(&local_velocity_fallbacks,&global_velocity_fallbacks,1,
                  MPI_INT,MPI_SUM,E->parallel.world);
    cache->velocity_block_fallbacks=global_velocity_fallbacks;
    MPI_Allreduce(&local_min_pivot_ratio,&global_min_pivot_ratio,1,MPI_DOUBLE,
                  MPI_MIN,E->parallel.world);
    MPI_Allreduce(&local_cache_bytes,&global_cache_bytes,1,MPI_DOUBLE,MPI_MAX,
                  E->parallel.world);
    build_seconds=CPU_time0()-build_start;
    MPI_Allreduce(&build_seconds,&global_build_seconds,1,MPI_DOUBLE,MPI_MAX,
                  E->parallel.world);
    if(global_blocks+global_fallback==0)
        myerror(E,"ALA shallow-patch depth selects no global elements");
    if(global_interface_blocks==0)
        myerror(E,"ALA shallow-patch built no MPI interface blocks");
    if(global_unique==0) {
        global_overlap_min=0;
        global_overlap_max=0;
    }
    if(E->parallel.me==0) {
        fprintf(E->fp,"ALA shallow-patch preconditioner depth_km=%e "
                "block=%dx%dx%d stride=%dx%dx1 mpi_halo_overlap=%d "
                "interface_block=%dx%dx%d partition_of_unity=global "
                "weight=%e regularization=%e "
                "velocity_solver=%s "
                "interface_velocity_metric=%s "
                "operator=principal((D+C)Mv^-1(D+C)^T) "
                "velocity_block_fallbacks=%d "
                "patch_entries=%d "
                "unique_elements=%d overlap_range=(%d,%d) "
                "valid_blocks=%d interface_blocks=%d "
                "interface_entries=%d ghost_elements=%d "
                "fallback_blocks=%d min_pivot_ratio=%e "
                "cache_megabytes_max=%e "
                "build_seconds_max=%e\n",
                E->control.ala_shallow_patch_depth_km,
                horizontal_elements,horizontal_elements,
                ALA_PATCH_RADIAL_ELEMENTS,
                horizontal_stride,horizontal_stride,overlap,
                2*overlap,horizontal_elements,ALA_PATCH_RADIAL_ELEMENTS,
                E->control.ala_shallow_patch_weight,
                E->control.ala_shallow_patch_regularization,
                E->control.ala_shallow_patch_velocity_solver,
                strcmp(E->control.ala_shallow_patch_velocity_solver,
                       "element_vanka")==0
                    ? "halo_element_vanka" : "halo_nodal",
                global_velocity_fallbacks,
                global_elements,global_unique,global_overlap_min,
                global_overlap_max,global_blocks,global_interface_blocks,
                global_interface_entries,global_ghost_elements,
                global_fallback,global_min_pivot_ratio,
                global_cache_bytes/(1024.0*1024.0),global_build_seconds);
        fprintf(stderr,"ALA shallow-patch preconditioner depth_km=%e "
                "block=%dx%dx%d stride=%dx%dx1 mpi_halo_overlap=%d "
                "interface_block=%dx%dx%d partition_of_unity=global "
                "weight=%e regularization=%e "
                "velocity_solver=%s "
                "interface_velocity_metric=%s "
                "operator=principal((D+C)Mv^-1(D+C)^T) "
                "velocity_block_fallbacks=%d "
                "patch_entries=%d "
                "unique_elements=%d overlap_range=(%d,%d) "
                "valid_blocks=%d interface_blocks=%d "
                "interface_entries=%d ghost_elements=%d "
                "fallback_blocks=%d min_pivot_ratio=%e "
                "cache_megabytes_max=%e "
                "build_seconds_max=%e\n",
                E->control.ala_shallow_patch_depth_km,
                horizontal_elements,horizontal_elements,
                ALA_PATCH_RADIAL_ELEMENTS,
                horizontal_stride,horizontal_stride,overlap,
                2*overlap,horizontal_elements,ALA_PATCH_RADIAL_ELEMENTS,
                E->control.ala_shallow_patch_weight,
                E->control.ala_shallow_patch_regularization,
                E->control.ala_shallow_patch_velocity_solver,
                strcmp(E->control.ala_shallow_patch_velocity_solver,
                       "element_vanka")==0
                    ? "halo_element_vanka" : "halo_nodal",
                global_velocity_fallbacks,
                global_elements,global_unique,global_overlap_min,
                global_overlap_max,global_blocks,global_interface_blocks,
                global_interface_entries,global_ghost_elements,
                global_fallback,global_min_pivot_ratio,
                global_cache_bytes/(1024.0*1024.0),global_build_seconds);
    }
}


static void build_ala_two_level_cache(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache, int lev)
{
    int m,cx,cy,cz,ce,x1,y1,z1,x2,y2,z2,e1,e2;
    int clev,factor,celx,cely,celz,elx,elz,cnpno;
    int local_invalid,global_invalid;
    double diagonal,local_min,local_max,global_min,global_max;
    clev=lev-E->control.ala_two_level_offset;
    factor=1 << E->control.ala_two_level_offset;
    if(factor>4)
        myerror(E,"ALA Galerkin aggregate factor above four is unsupported");
    celx=E->lmesh.ELX[clev];
    cely=E->lmesh.ELY[clev];
    celz=E->lmesh.ELZ[clev];
    elx=E->lmesh.ELX[lev];
    elz=E->lmesh.ELZ[lev];
    cnpno=E->lmesh.NPNO[clev];
    cache->coarse_eigenvalue_min=
        E->control.ala_two_level_coarse_eigenvalue_min;
    cache->coarse_eigenvalue_max=
        E->control.ala_two_level_coarse_eigenvalue_max;
    cache->velocity_eigenvalue_min=
        E->control.ala_two_level_velocity_eigenvalue_min;
    cache->velocity_eigenvalue_max=
        E->control.ala_two_level_velocity_eigenvalue_max;
    local_invalid=0;
    local_min=1.0e300;
    local_max=0.0;

    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        cache->coarse_bpi[m]=(double *)calloc(cnpno+1,sizeof(double));
        if(cache->coarse_bpi[m]==NULL)
            myerror(E,"Unable to allocate ALA Galerkin coarse diagonal");
        for(cy=1;cy<=cely;cy++)
            for(cx=1;cx<=celx;cx++)
                for(cz=1;cz<=celz;cz++) {
                    ce=cz+(cx-1)*celz+(cy-1)*celz*celx;
                    diagonal=0.0;
                    for(y1=(cy-1)*factor+1;y1<=cy*factor;y1++)
                        for(x1=(cx-1)*factor+1;x1<=cx*factor;x1++)
                            for(z1=(cz-1)*factor+1;z1<=cz*factor;z1++) {
                                e1=z1+(x1-1)*elz+(y1-1)*elz*elx;
                                for(y2=(cy-1)*factor+1;y2<=cy*factor;y2++)
                                    for(x2=(cx-1)*factor+1;x2<=cx*factor;x2++)
                                        for(z2=(cz-1)*factor+1;
                                            z2<=cz*factor;z2++) {
                                            if(abs(x1-x2)>1 ||
                                               abs(y1-y2)>1 ||
                                               abs(z1-z2)>1)
                                                continue;
                                            e2=z2+(x2-1)*elz
                                                +(y2-1)*elz*elx;
                                            diagonal +=
                                                assemble_Ahatp_jacobi_entry(
                                                    E,e1,e2,lev,m);
                                        }
                            }
                    if(!isfinite(diagonal) || diagonal<=0.0) {
                        cache->coarse_bpi[m][ce]=1.0;
                        local_invalid++;
                    }
                    else
                        cache->coarse_bpi[m][ce]=1.0/diagonal;
                    local_min=min(local_min,cache->coarse_bpi[m][ce]);
                    local_max=max(local_max,cache->coarse_bpi[m][ce]);
                }
    }
    MPI_Allreduce(&local_invalid,&global_invalid,1,MPI_INT,MPI_SUM,
                  E->parallel.world);
    MPI_Allreduce(&local_min,&global_min,1,MPI_DOUBLE,MPI_MIN,
                  E->parallel.world);
    MPI_Allreduce(&local_max,&global_max,1,MPI_DOUBLE,MPI_MAX,
                  E->parallel.world);
    if(global_invalid>0 && E->parallel.me==0)
        fprintf(stderr,"ALA Galerkin aggregate diagonal invalid=%d\n",
                global_invalid);
    if(global_invalid>0)
        myerror(E,"ALA Galerkin aggregate diagonal is not positive");
    if(E->parallel.me==0) {
        fprintf(E->fp,"ALA Galerkin aggregate diagonal level=%d factor=%d "
                "BPI_range=(%e,%e) invalid=%d "
                "operator_diag=diag(Pt*((D+C)diag(K)^-1(D+C)^T)*P)\n",
                clev,factor,global_min,global_max,global_invalid);
        fprintf(stderr,"ALA Galerkin aggregate diagonal level=%d factor=%d "
                "BPI_range=(%e,%e) invalid=%d "
                "operator_diag=diag(Pt*((D+C)diag(K)^-1(D+C)^T)*P)\n",
                clev,factor,global_min,global_max,global_invalid);
    }
}


/* Stage 8 pressure multigrid uses the existing distributed pressure grids.
   Its level operator is the complete strict-ALA Gram approximation

       A_l = (D_l+C_l) diag(K_gamma,l)^-1 (D_l+C_l)^T.

   This is deliberately cheaper than the legacy two-level Galerkin
   experiment: no pressure smoothing step invokes a velocity MG solve or
   prolongs back to the finest grid.  FGMRES contains the resulting
   rediscretized, mildly nonsymmetric V-cycle safely. */
static void build_ala_pressure_multigrid_cache(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache, int lev)
{
    int level,m,npno,neq;

    cache->pressure_mg_min_level=
        E->control.ala_pressure_multigrid_min_level;
    cache->pressure_mg_max_level=lev;
    cache->pressure_mg_operator_applications=0;
    for(level=cache->pressure_mg_min_level;level<=lev;level++) {
        npno=E->lmesh.NPNO[level];
        for(m=1;m<=E->sphere.caps_per_proc;m++) {
            cache->pressure_mg_rhs[level][m]=
                (double *)calloc(npno+1,sizeof(double));
            cache->pressure_mg_x[level][m]=
                (double *)calloc(npno+1,sizeof(double));
            cache->pressure_mg_residual[level][m]=
                (double *)calloc(npno+1,sizeof(double));
            cache->pressure_mg_Ax[level][m]=
                (double *)calloc(npno+1,sizeof(double));
            if(cache->pressure_mg_rhs[level][m]==NULL ||
               cache->pressure_mg_x[level][m]==NULL ||
               cache->pressure_mg_residual[level][m]==NULL ||
               cache->pressure_mg_Ax[level][m]==NULL)
                myerror(E,"Unable to allocate ALA pressure multigrid level");
        }
    }
    neq=E->lmesh.NEQ[lev];
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        cache->pressure_mg_velocity[m]=
            (double *)calloc(neq+1,sizeof(double));
        cache->pressure_mg_velocity_rhs[m]=
            (double *)calloc(neq+1,sizeof(double));
        cache->pressure_mg_galerkin_p[m]=
            (double *)calloc(E->lmesh.NPNO[lev]+1,sizeof(double));
        cache->pressure_mg_galerkin_Ap[m]=
            (double *)calloc(E->lmesh.NPNO[lev]+1,sizeof(double));
        if(cache->pressure_mg_velocity[m]==NULL ||
           cache->pressure_mg_velocity_rhs[m]==NULL ||
           cache->pressure_mg_galerkin_p[m]==NULL ||
           cache->pressure_mg_galerkin_Ap[m]==NULL)
            myerror(E,"Unable to allocate ALA pressure multigrid velocity work");
    }
    if(E->parallel.me==0) {
        fprintf(E->fp,"ALA pressure multigrid hierarchy levels=%d range=%d:%d "
                "smooth=(pre:%d,post:%d,coarse:%d) damping=%e weight=%e "
                "operator=%s "
                "transfer=normalized_constant_P/Pt scale=1/sqrt(children)\n",
                lev-cache->pressure_mg_min_level+1,
                cache->pressure_mg_min_level,lev,
                E->control.ala_pressure_multigrid_pre_smooth,
                E->control.ala_pressure_multigrid_post_smooth,
                E->control.ala_pressure_multigrid_coarse_iterations,
                E->control.ala_pressure_multigrid_damping,
                E->control.ala_pressure_multigrid_weight,
                E->control.ala_pressure_multigrid_galerkin
                    ? "Galerkin_Pt*((D+C)diag(Kgamma)^-1(D+C)^T)*P"
                    : "rediscretized_(D+C)diag(Kgamma)^-1(D+C)^T");
        fflush(E->fp);
    }
}


static void apply_ala_pressure_multigrid_operator(
    struct All_variables *E, double **p, double **Ap, int level,
    struct ala_pressure_preconditioner_cache *cache)
{
    int m,j,neq,e,ex,ey,ez,cx,cy,cz,ce,factor;
    int elx,ely,elz,celx,celz,cnpno,finest;
    void assemble_grad_rho_p();
    void assemble_div_rho_u();

    finest=cache->pressure_mg_max_level;
    if(E->control.ala_pressure_multigrid_galerkin && level<finest) {
        factor=1<<(finest-level);
        elx=E->lmesh.ELX[finest];
        ely=E->lmesh.ELY[finest];
        elz=E->lmesh.ELZ[finest];
        celx=E->lmesh.ELX[level];
        celz=E->lmesh.ELZ[level];
        cnpno=E->lmesh.NPNO[level];
        if(elx!=factor*celx || ely!=factor*E->lmesh.ELY[level] ||
           elz!=factor*celz)
            myerror(E,"ALA Galerkin pressure hierarchy mesh mismatch");
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(ey=1;ey<=ely;ey++)
                for(ex=1;ex<=elx;ex++)
                    for(ez=1;ez<=elz;ez++) {
                        e=ez+(ex-1)*elz+(ey-1)*elz*elx;
                        cx=(ex-1)/factor+1;
                        cy=(ey-1)/factor+1;
                        cz=(ez-1)/factor+1;
                        ce=cz+(cx-1)*celz+(cy-1)*celz*celx;
                        cache->pressure_mg_galerkin_p[m][e]=
                            p[m][ce]/sqrt((double)(factor*factor*factor));
                    }
        apply_ala_pressure_multigrid_operator(
            E,cache->pressure_mg_galerkin_p,
            cache->pressure_mg_galerkin_Ap,finest,cache);
        for(m=1;m<=E->sphere.caps_per_proc;m++) {
            for(ce=1;ce<=cnpno;ce++)
                Ap[m][ce]=0.0;
            for(ey=1;ey<=ely;ey++)
                for(ex=1;ex<=elx;ex++)
                    for(ez=1;ez<=elz;ez++) {
                        e=ez+(ex-1)*elz+(ey-1)*elz*elx;
                        cx=(ex-1)/factor+1;
                        cy=(ey-1)/factor+1;
                        cz=(ez-1)/factor+1;
                        ce=cz+(cx-1)*celz+(cy-1)*celz*celx;
                        Ap[m][ce] += cache->pressure_mg_galerkin_Ap[m][e]
                            /sqrt((double)(factor*factor*factor));
                    }
        }
        return;
    }

    neq=E->lmesh.NEQ[level];
    assemble_grad_rho_p(E,p,cache->pressure_mg_velocity_rhs,level);
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        for(j=0;j<neq;j++)
            cache->pressure_mg_velocity[m][j]=
                E->ALA_velocity_BI[level][m][j]
                *cache->pressure_mg_velocity_rhs[m][j];
        cache->pressure_mg_velocity[m][neq]=0.0;
    }
    assemble_div_rho_u(E,cache->pressure_mg_velocity,Ap,level);
    cache->pressure_mg_operator_applications++;
}


static void smooth_ala_pressure_multigrid_level(
    struct All_variables *E, int level, int sweeps,
    struct ala_pressure_preconditioner_cache *cache)
{
    int sweep,m,e,npno;
    double damping;
    double **rhs,**x,**Ax;

    npno=E->lmesh.NPNO[level];
    damping=E->control.ala_pressure_multigrid_damping;
    rhs=cache->pressure_mg_rhs[level];
    x=cache->pressure_mg_x[level];
    Ax=cache->pressure_mg_Ax[level];
    for(sweep=0;sweep<sweeps;sweep++) {
        apply_ala_pressure_multigrid_operator(E,x,Ax,level,cache);
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=1;e<=npno;e++)
                x[m][e] += damping*E->BPI[level][m][e]
                    *(rhs[m][e]-Ax[m][e]);
    }
}


static void ala_pressure_multigrid_vcycle(struct All_variables *E, int level,
    struct ala_pressure_preconditioner_cache *cache)
{
    int m,e,ex,ey,ez,cx,cy,cz,ce;
    int npno,cnpno,elx,ely,elz,celx,celz;
    double **rhs,**x,**residual,**Ax,**coarse_rhs,**coarse_x;

    rhs=cache->pressure_mg_rhs[level];
    x=cache->pressure_mg_x[level];
    residual=cache->pressure_mg_residual[level];
    Ax=cache->pressure_mg_Ax[level];
    npno=E->lmesh.NPNO[level];
    if(level==cache->pressure_mg_min_level) {
        smooth_ala_pressure_multigrid_level(
            E,level,E->control.ala_pressure_multigrid_coarse_iterations,
            cache);
        return;
    }

    smooth_ala_pressure_multigrid_level(
        E,level,E->control.ala_pressure_multigrid_pre_smooth,cache);
    apply_ala_pressure_multigrid_operator(E,x,Ax,level,cache);
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(e=1;e<=npno;e++)
            residual[m][e]=rhs[m][e]-Ax[m][e];

    coarse_rhs=cache->pressure_mg_rhs[level-1];
    coarse_x=cache->pressure_mg_x[level-1];
    cnpno=E->lmesh.NPNO[level-1];
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(ce=1;ce<=cnpno;ce++) {
            coarse_rhs[m][ce]=0.0;
            coarse_x[m][ce]=0.0;
        }
    elx=E->lmesh.ELX[level];
    ely=E->lmesh.ELY[level];
    elz=E->lmesh.ELZ[level];
    celx=E->lmesh.ELX[level-1];
    celz=E->lmesh.ELZ[level-1];
    if(elx!=2*celx || ely!=2*E->lmesh.ELY[level-1] ||
       elz!=2*celz)
        myerror(E,"ALA pressure multigrid requires factor-two mesh levels");
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(ey=1;ey<=ely;ey++)
            for(ex=1;ex<=elx;ex++)
                for(ez=1;ez<=elz;ez++) {
                    e=ez+(ex-1)*elz+(ey-1)*elz*elx;
                    cx=(ex-1)/2+1;
                    cy=(ey-1)/2+1;
                    cz=(ez-1)/2+1;
                    ce=cz+(cx-1)*celz+(cy-1)*celz*celx;
                    coarse_rhs[m][ce] += residual[m][e]/sqrt(8.0);
                }

    ala_pressure_multigrid_vcycle(E,level-1,cache);

    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(ey=1;ey<=ely;ey++)
            for(ex=1;ex<=elx;ex++)
                for(ez=1;ez<=elz;ez++) {
                    e=ez+(ex-1)*elz+(ey-1)*elz*elx;
                    cx=(ex-1)/2+1;
                    cy=(ey-1)/2+1;
                    cz=(ez-1)/2+1;
                    ce=cz+(cx-1)*celz+(cy-1)*celz*celx;
                    x[m][e] += coarse_x[m][ce]/sqrt(8.0);
                }
    smooth_ala_pressure_multigrid_level(
        E,level,E->control.ala_pressure_multigrid_post_smooth,cache);
}


static void apply_ala_pressure_multigrid_correction(
    struct All_variables *E, double **r, double **z, int lev, int iteration,
    struct ala_pressure_preconditioner_cache *cache)
{
    int level,m,e,npno;
    double weight;
    double local[4],global[4];
    double **fine_rhs,**fine_x,**fine_Ax;

    weight=E->control.ala_pressure_multigrid_weight;
    for(level=cache->pressure_mg_min_level;level<=lev;level++) {
        npno=E->lmesh.NPNO[level];
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=1;e<=npno;e++) {
                cache->pressure_mg_rhs[level][m][e]=0.0;
                cache->pressure_mg_x[level][m][e]=0.0;
                cache->pressure_mg_residual[level][m][e]=0.0;
                cache->pressure_mg_Ax[level][m][e]=0.0;
            }
    }
    fine_rhs=cache->pressure_mg_rhs[lev];
    fine_x=cache->pressure_mg_x[lev];
    fine_Ax=cache->pressure_mg_Ax[lev];
    npno=E->lmesh.NPNO[lev];
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(e=1;e<=npno;e++)
            fine_rhs[m][e]=r[m][e];
    ala_pressure_multigrid_vcycle(E,lev,cache);

    local[0]=local[1]=local[2]=local[3]=0.0;
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(e=1;e<=npno;e++) {
            local[0] += r[m][e]*z[m][e];
            local[1] += r[m][e]*fine_x[m][e];
        }
    MPI_Allreduce(local,global,2,MPI_DOUBLE,MPI_SUM,E->parallel.world);
    if(!isfinite(global[1]) || global[1]<=0.0) {
        if(E->parallel.me==0) {
            fprintf(E->fp,"ALA_PRESSURE_MG_ROLLBACK iteration=%d "
                    "base_energy=%e vcycle_energy=%e reason=nonpositive\n",
                    iteration,global[0],global[1]);
            fflush(E->fp);
        }
        return;
    }
    if(iteration==0 ||
       iteration%E->control.ala_coarse_residual_interval==0) {
        apply_ala_pressure_multigrid_operator(E,fine_x,fine_Ax,lev,cache);
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=1;e<=npno;e++) {
                local[2] += r[m][e]*r[m][e];
                local[3] += (r[m][e]-fine_Ax[m][e])
                    *(r[m][e]-fine_Ax[m][e]);
            }
        MPI_Allreduce(local+2,global+2,2,MPI_DOUBLE,MPI_SUM,
                      E->parallel.world);
        if(E->parallel.me==0) {
            fprintf(E->fp,"ALA_PRESSURE_MG iteration=%d levels=%d "
                    "base_energy=%e vcycle_energy=%e "
                    "vcycle_to_base=%e residual_reduction=%e weight=%e "
                    "operator_applications_total=%d\n",iteration,
                    lev-cache->pressure_mg_min_level+1,global[0],global[1],
                    global[1]/max(global[0],1.0e-300),
                    sqrt(global[3]/max(global[2],1.0e-300)),weight,
                    cache->pressure_mg_operator_applications);
            fflush(E->fp);
        }
    }
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(e=1;e<=npno;e++)
            z[m][e]=(1.0-weight)*z[m][e]+weight*fine_x[m][e];
}


static void free_ala_pressure_preconditioner_cache(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache)
{
    int m,face,level;
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        free((void *)cache->size[m]);
        free((void *)cache->multiplicity[m]);
        free((void *)cache->elements[m]);
        free((void *)cache->chol[m]);
        free((void *)cache->interface_size[m]);
        free((void *)cache->interface_face[m]);
        free((void *)cache->interface_elements[m]);
        free((void *)cache->interface_chol[m]);
        free((void *)cache->velocity_block_chol[m]);
        for(face=0;face<ALA_PATCH_MPI_FACES;face++) {
            free((void *)cache->halo_send_elements[m][face]);
            free((void *)cache->halo_send_records[m][face]);
            free((void *)cache->halo_recv_records[m][face]);
            free((void *)cache->halo_multiplicity[m][face]);
        }
        free((void *)cache->coarse_bpi[m]);
        free((void *)cache->global_basis[m]);
        free((void *)cache->geneo_basis[m]);
        free((void *)cache->pressure_mg_velocity[m]);
        free((void *)cache->pressure_mg_velocity_rhs[m]);
        free((void *)cache->pressure_mg_galerkin_p[m]);
        free((void *)cache->pressure_mg_galerkin_Ap[m]);
        free((void *)cache->pressure_depth_weighted_r[m]);
        for(level=0;level<MAX_LEVELS;level++) {
            free((void *)cache->pressure_mg_rhs[level][m]);
            free((void *)cache->pressure_mg_x[level][m]);
            free((void *)cache->pressure_mg_residual[level][m]);
            free((void *)cache->pressure_mg_Ax[level][m]);
        }
    }
    free((void *)cache->global_chol);
    free((void *)cache->global_matrix);
    free((void *)cache->geneo_chol);
    free((void *)cache->geneo_matrix);
    free((void *)cache->geneo_eigenvalues);
    memset(cache,0,sizeof(*cache));
}


static void exchange_ala_shallow_halo_values(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache, int lev, int m,
    double **local, double *ghost[ALA_PATCH_MPI_FACES], int return_to_owner)
{
    int face,i,e,target,nrequest;
    double *send[ALA_PATCH_MPI_FACES];
    double *received[ALA_PATCH_MPI_FACES];
    MPI_Request requests[2*ALA_PATCH_MPI_FACES];
    MPI_Status statuses[2*ALA_PATCH_MPI_FACES];

    for(face=0;face<ALA_PATCH_MPI_FACES;face++) {
        if(return_to_owner) {
            send[face]=ghost[face];
            received[face]=(double *)calloc(
                max(cache->halo_send_count[m][face],1),sizeof(double));
        }
        else {
            send[face]=(double *)calloc(
                max(cache->halo_send_count[m][face],1),sizeof(double));
            received[face]=ghost[face];
            for(i=0;i<cache->halo_send_count[m][face];i++) {
                e=cache->halo_send_elements[m][face][i];
                send[face][i]=local[m][e];
            }
        }
        if(send[face]==NULL || received[face]==NULL)
            myerror(E,"Unable to allocate ALA halo exchange buffer");
    }
    nrequest=0;
    for(face=0;face<ALA_PATCH_MPI_FACES;face++) {
        target=E->parallel.PROCESSOR[lev][m].pass[face+1];
        MPI_Irecv(received[face],
            return_to_owner ? cache->halo_send_count[m][face]
                            : cache->halo_recv_count[m][face],
            MPI_DOUBLE,target,return_to_owner ? 31905 : 31904,
            E->parallel.world,&requests[nrequest++]);
        MPI_Isend(send[face],
            return_to_owner ? cache->halo_recv_count[m][face]
                            : cache->halo_send_count[m][face],
            MPI_DOUBLE,target,return_to_owner ? 31905 : 31904,
            E->parallel.world,&requests[nrequest++]);
    }
    MPI_Waitall(nrequest,requests,statuses);
    if(return_to_owner)
        for(face=0;face<ALA_PATCH_MPI_FACES;face++)
            for(i=0;i<cache->halo_send_count[m][face];i++) {
                e=cache->halo_send_elements[m][face][i];
                local[m][e] += received[face][i];
            }
    for(face=0;face<ALA_PATCH_MPI_FACES;face++) {
        if(return_to_owner)
            free((void *)received[face]);
        else
            free((void *)send[face]);
    }
}


static void apply_ala_coarse_fixed_schur(struct All_variables *E,
                                         double **p, double **Ap,
                                         double **velocity,
                                         double **velocity_rhs,
                                         double **velocity_Ax,
                                         double **velocity_direction,
                                         struct ala_pressure_preconditioner_cache *cache,
                                         int lev,
                                         double *residual_reduction)
{
    int m,j,k,neq;
    double theta,delta,sigma,rho_cheb,rho_new,residual;
    double local_norm[2],global_norm[2];
    void assemble_grad_rho_p();
    void assemble_div_rho_u();
    void assemble_del2_u();

    neq=E->lmesh.NEQ[lev];
    assemble_grad_rho_p(E,p,velocity_rhs,lev);
    if(strcmp(E->control.ala_two_level_velocity_solver,"diagonal")==0) {
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(j=0;j<neq;j++)
                velocity[m][j]=E->ALA_velocity_BI[lev][m][j]
                    *velocity_rhs[m][j];
    }
    else {
        theta=0.5*(cache->velocity_eigenvalue_max+
                   cache->velocity_eigenvalue_min);
        delta=0.5*(cache->velocity_eigenvalue_max-
                   cache->velocity_eigenvalue_min);
        sigma=theta/delta;
        rho_cheb=1.0/sigma;
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(j=0;j<neq;j++) {
                velocity_direction[m][j]=E->ALA_velocity_BI[lev][m][j]
                    *velocity_rhs[m][j]/theta;
                velocity[m][j]=velocity_direction[m][j];
            }
        for(k=1;k<E->control.ala_two_level_velocity_iterations;k++) {
            assemble_del2_u(E,velocity,velocity_Ax,lev,1);
            rho_new=1.0/(2.0*sigma-rho_cheb);
            for(m=1;m<=E->sphere.caps_per_proc;m++)
                for(j=0;j<neq;j++) {
                    residual=velocity_rhs[m][j]-velocity_Ax[m][j];
                    velocity_direction[m][j]=rho_new*rho_cheb
                        *velocity_direction[m][j]
                        +(2.0*rho_new/delta)*E->ALA_velocity_BI[lev][m][j]
                        *residual;
                    velocity[m][j] += velocity_direction[m][j];
                }
            rho_cheb=rho_new;
        }
    }
    if(residual_reduction!=NULL) {
        assemble_del2_u(E,velocity,velocity_Ax,lev,1);
        local_norm[0]=0.0;
        local_norm[1]=0.0;
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(j=0;j<neq;j++) {
                local_norm[0] += velocity_rhs[m][j]*velocity_rhs[m][j];
                residual=velocity_rhs[m][j]-velocity_Ax[m][j];
                local_norm[1] += residual*residual;
            }
        MPI_Allreduce(local_norm,global_norm,2,MPI_DOUBLE,MPI_SUM,
                      E->parallel.world);
        *residual_reduction=sqrt(global_norm[1]
                                 /max(global_norm[0],1.0e-300));
    }
    assemble_div_rho_u(E,velocity,Ap,lev);
}


static void apply_ala_galerkin_fixed_schur(struct All_variables *E,
                                           double **coarse_p,
                                           double **coarse_Ap,
                                           double **fine_p,
                                           double **fine_Ap,
                                           double **fine_velocity,
                                           double **fine_velocity_rhs,
                                           double **fine_velocity_Ax,
                                           double **fine_velocity_direction,
                                           struct ala_pressure_preconditioner_cache *cache,
                                           int clev, int lev, int factor,
                                           double *residual_reduction)
{
    int m,ex,ey,ez,e,cx,cy,cz,ce;
    int elx,ely,elz,celx,celz,cnpno;

    elx=E->lmesh.ELX[lev];
    ely=E->lmesh.ELY[lev];
    elz=E->lmesh.ELZ[lev];
    celx=E->lmesh.ELX[clev];
    celz=E->lmesh.ELZ[clev];
    cnpno=E->lmesh.NPNO[clev];

    /* Constant pressure prolongation P. Every fine pressure cell has one
       geometrically coincident coarse parent. */
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(ey=1;ey<=ely;ey++)
            for(ex=1;ex<=elx;ex++)
                for(ez=1;ez<=elz;ez++) {
                    e=ez+(ex-1)*elz+(ey-1)*elz*elx;
                    cx=(ex-1)/factor+1;
                    cy=(ey-1)/factor+1;
                    cz=(ez-1)/factor+1;
                    ce=cz+(cx-1)*celz+(cy-1)*celz*celx;
                    fine_p[m][e]=coarse_p[m][ce];
                }

    /* Apply the complete fine-grid strict-ALA operator. Its velocity assembly
       carries the pressure correction across MPI and cap boundaries. */
    apply_ala_coarse_fixed_schur(E,fine_p,fine_Ap,fine_velocity,
                                 fine_velocity_rhs,fine_velocity_Ax,
                                 fine_velocity_direction,cache,lev,
                                 residual_reduction);

    /* Exact transpose restriction P^T. Together these operations define the
       matrix-free Galerkin coarse operator S_c=P^T S_f P. */
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        for(ce=1;ce<=cnpno;ce++)
            coarse_Ap[m][ce]=0.0;
        for(ey=1;ey<=ely;ey++)
            for(ex=1;ex<=elx;ex++)
                for(ez=1;ez<=elz;ez++) {
                    e=ez+(ex-1)*elz+(ey-1)*elz*elx;
                    cx=(ex-1)/factor+1;
                    cy=(ey-1)/factor+1;
                    cz=(ez-1)/factor+1;
                    ce=cz+(cx-1)*celz+(cy-1)*celz*celx;
                    coarse_Ap[m][ce] += fine_Ap[m][e];
                }
    }
}


/* Build a replicated lowest-level Galerkin solve.  The basis is the tensor
   product of real Cartesian spherical harmonics through degree three and six
   radial linear hats.  The matrix is assembled with the same frozen
   matrix-free Sc used by the level-one polynomial, then factored once. */
static void build_ala_global_coarse_cache(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache, int lev)
{
    int m,ce,i,j,k,col,row,clev,factor,cnpno,stride,n;
    int npno,neq,radial,angular_index,interval;
    double theta,phi,radius,depth,x,y,z,r2,radial_value;
    double angular[ALA_GLOBAL_ANGULAR_BASIS];
    double radial_hat[ALA_GLOBAL_RADIAL_BASIS];
    double local_norm[ALA_GLOBAL_BASIS_COUNT];
    double global_norm[ALA_GLOBAL_BASIS_COUNT];
    double *local_column,*global_column,*matrix;
    double *coarse_p[NCS],*coarse_Ap[NCS];
    double *fine_p[NCS],*fine_Ap[NCS],*fine_velocity[NCS];
    double *fine_velocity_rhs[NCS],*fine_velocity_Ax[NCS];
    double *fine_velocity_direction[NCS];
    double anti2,total2,symmetric,diagonal_min,diagonal_max,shift;
    double sum,pivot,min_pivot,factor_error2,matrix_norm2;
    static const double radial_knots[ALA_GLOBAL_RADIAL_BASIS] =
        {0.0,200.0,410.0,660.0,1200.0,2891.0};

    clev=lev-E->control.ala_two_level_offset;
    factor=1 << E->control.ala_two_level_offset;
    cnpno=E->lmesh.NPNO[clev];
    stride=cnpno+1;
    npno=E->lmesh.NPNO[lev];
    neq=E->lmesh.NEQ[lev];
    n=ALA_GLOBAL_BASIS_COUNT;
    cache->global_basis_count=n;
    cache->global_chol=(double *)calloc(n*n,sizeof(double));
    cache->global_matrix=(double *)calloc(n*n,sizeof(double));
    matrix=cache->global_matrix;
    local_column=(double *)calloc(n,sizeof(double));
    global_column=(double *)calloc(n,sizeof(double));
    if(cache->global_chol==NULL || matrix==NULL || local_column==NULL ||
       global_column==NULL)
        myerror(E,"Unable to allocate ALA global coarse matrix");

    for(i=0;i<n;i++)
        local_norm[i]=0.0;
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        cache->global_basis[m]=(double *)calloc(n*stride,sizeof(double));
        if(cache->global_basis[m]==NULL)
            myerror(E,"Unable to allocate ALA global pressure basis");
        for(ce=1;ce<=cnpno;ce++) {
            theta=E->ECO[clev][m][ce].centre[1];
            phi=E->ECO[clev][m][ce].centre[2];
            radius=E->ECO[clev][m][ce].centre[3];
            x=sin(theta)*cos(phi);
            y=sin(theta)*sin(phi);
            z=cos(theta);
            r2=x*x+y*y+z*z;
            angular[0]=1.0;
            angular[1]=x;
            angular[2]=y;
            angular[3]=z;
            angular[4]=x*y;
            angular[5]=y*z;
            angular[6]=z*x;
            angular[7]=x*x-y*y;
            angular[8]=3.0*z*z-r2;
            angular[9]=x*(5.0*z*z-r2);
            angular[10]=y*(5.0*z*z-r2);
            angular[11]=z*(5.0*z*z-3.0*r2);
            angular[12]=z*(x*x-y*y);
            angular[13]=x*y*z;
            angular[14]=x*(x*x-3.0*y*y);
            angular[15]=y*(3.0*x*x-y*y);

            for(radial=0;radial<ALA_GLOBAL_RADIAL_BASIS;radial++)
                radial_hat[radial]=0.0;
            depth=(1.0-radius)*E->data.radius_km;
            if(depth<=radial_knots[0])
                radial_hat[0]=1.0;
            else if(depth>=radial_knots[ALA_GLOBAL_RADIAL_BASIS-1])
                radial_hat[ALA_GLOBAL_RADIAL_BASIS-1]=1.0;
            else {
                interval=0;
                while(interval<ALA_GLOBAL_RADIAL_BASIS-1 &&
                      depth>radial_knots[interval+1])
                    interval++;
                radial_value=(depth-radial_knots[interval])
                    /(radial_knots[interval+1]-radial_knots[interval]);
                radial_hat[interval]=1.0-radial_value;
                radial_hat[interval+1]=radial_value;
            }
            for(radial=0;radial<ALA_GLOBAL_RADIAL_BASIS;radial++)
                for(angular_index=0;
                    angular_index<ALA_GLOBAL_ANGULAR_BASIS;
                    angular_index++) {
                    i=radial*ALA_GLOBAL_ANGULAR_BASIS+angular_index;
                    cache->global_basis[m][i*stride+ce]
                        =radial_hat[radial]*angular[angular_index];
                    local_norm[i] += cache->global_basis[m][i*stride+ce]
                        *cache->global_basis[m][i*stride+ce];
                }
        }
    }
    MPI_Allreduce(local_norm,global_norm,n,MPI_DOUBLE,MPI_SUM,
                  E->parallel.world);
    for(i=0;i<n;i++) {
        if(!isfinite(global_norm[i]) || global_norm[i]<=1.0e-20)
            myerror(E,"ALA global pressure basis is rank deficient");
        global_norm[i]=sqrt(global_norm[i]);
    }
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(i=0;i<n;i++)
            for(ce=1;ce<=cnpno;ce++)
                cache->global_basis[m][i*stride+ce] /= global_norm[i];

    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        coarse_p[m]=(double *)calloc(cnpno+1,sizeof(double));
        coarse_Ap[m]=(double *)calloc(cnpno+1,sizeof(double));
        fine_p[m]=(double *)calloc(npno+1,sizeof(double));
        fine_Ap[m]=(double *)calloc(npno+1,sizeof(double));
        fine_velocity[m]=(double *)calloc(neq+1,sizeof(double));
        fine_velocity_rhs[m]=(double *)calloc(neq+1,sizeof(double));
        fine_velocity_Ax[m]=(double *)calloc(neq+1,sizeof(double));
        fine_velocity_direction[m]=(double *)calloc(neq+1,sizeof(double));
        if(coarse_p[m]==NULL || coarse_Ap[m]==NULL || fine_p[m]==NULL ||
           fine_Ap[m]==NULL || fine_velocity[m]==NULL ||
           fine_velocity_rhs[m]==NULL || fine_velocity_Ax[m]==NULL ||
           fine_velocity_direction[m]==NULL)
            myerror(E,"Unable to allocate ALA global coarse workspace");
    }

    for(col=0;col<n;col++) {
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(ce=1;ce<=cnpno;ce++)
                coarse_p[m][ce]=cache->global_basis[m][col*stride+ce];
        apply_ala_galerkin_fixed_schur(
            E,coarse_p,coarse_Ap,fine_p,fine_Ap,fine_velocity,
            fine_velocity_rhs,fine_velocity_Ax,fine_velocity_direction,
            cache,clev,lev,factor,NULL);
        for(row=0;row<n;row++) {
            local_column[row]=0.0;
            for(m=1;m<=E->sphere.caps_per_proc;m++)
                for(ce=1;ce<=cnpno;ce++)
                    local_column[row] +=
                        cache->global_basis[m][row*stride+ce]
                        *coarse_Ap[m][ce];
        }
        MPI_Allreduce(local_column,global_column,n,MPI_DOUBLE,MPI_SUM,
                      E->parallel.world);
        for(row=0;row<n;row++)
            matrix[row*n+col]=global_column[row];
    }

    anti2=0.0;
    total2=0.0;
    diagonal_min=1.0e300;
    diagonal_max=0.0;
    for(i=0;i<n;i++) {
        diagonal_min=min(diagonal_min,matrix[i*n+i]);
        diagonal_max=max(diagonal_max,matrix[i*n+i]);
        for(j=0;j<n;j++) {
            anti2 += (matrix[i*n+j]-matrix[j*n+i])
                *(matrix[i*n+j]-matrix[j*n+i]);
            total2 += matrix[i*n+j]*matrix[i*n+j];
        }
    }
    symmetric=sqrt(anti2/max(total2,1.0e-300));
    if(!isfinite(symmetric) || symmetric>1.0e-8)
        myerror(E,"ALA global Galerkin matrix is not symmetric");
    if(!isfinite(diagonal_min) || diagonal_min<=0.0 ||
       !isfinite(diagonal_max) || diagonal_max<=0.0)
        myerror(E,"ALA global Galerkin diagonal is not positive");
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            cache->global_chol[i*n+j]
                =0.5*(matrix[i*n+j]+matrix[j*n+i]);
    shift=E->control.ala_global_coarse_regularization*diagonal_max;
    for(i=0;i<n;i++)
        cache->global_chol[i*n+i] += shift;

    min_pivot=1.0e300;
    for(i=0;i<n;i++) {
        for(j=0;j<=i;j++) {
            sum=cache->global_chol[i*n+j];
            for(k=0;k<j;k++)
                sum -= cache->global_chol[i*n+k]
                    *cache->global_chol[j*n+k];
            if(i==j) {
                pivot=sum;
                if(!isfinite(pivot) ||
                   pivot<=1.0e-14*max(diagonal_max,1.0e-300))
                    myerror(E,"ALA global Galerkin Cholesky failed");
                cache->global_chol[i*n+j]=sqrt(pivot);
                min_pivot=min(min_pivot,pivot);
            }
            else
                cache->global_chol[i*n+j]
                    =sum/cache->global_chol[j*n+j];
        }
        for(j=i+1;j<n;j++)
            cache->global_chol[i*n+j]=0.0;
    }
    factor_error2=0.0;
    matrix_norm2=0.0;
    for(i=0;i<n;i++)
        for(j=0;j<n;j++) {
            sum=0.0;
            for(k=0;k<=min(i,j);k++)
                sum += cache->global_chol[i*n+k]
                    *cache->global_chol[j*n+k];
            pivot=0.5*(matrix[i*n+j]+matrix[j*n+i]);
            if(i==j)
                pivot += shift;
            factor_error2 += (sum-pivot)*(sum-pivot);
            matrix_norm2 += pivot*pivot;
        }
    if(E->parallel.me==0) {
        fprintf(E->fp,"ALA global coarse basis count=%d angular_degree=3 "
                "radial_knots_km=(0,200,410,660,1200,2891) "
                "symmetry_relative=%e diagonal_range=(%e,%e) "
                "regularization=%e shift=%e min_pivot=%e "
                "factorization_relative=%e factorization=complete\n",
                n,symmetric,diagonal_min,diagonal_max,
                E->control.ala_global_coarse_regularization,shift,min_pivot,
                sqrt(factor_error2/max(matrix_norm2,1.0e-300)));
        fprintf(stderr,"ALA global coarse basis count=%d angular_degree=3 "
                "radial_knots_km=(0,200,410,660,1200,2891) "
                "symmetry_relative=%e diagonal_range=(%e,%e) "
                "regularization=%e shift=%e min_pivot=%e "
                "factorization_relative=%e factorization=complete\n",
                n,symmetric,diagonal_min,diagonal_max,
                E->control.ala_global_coarse_regularization,shift,min_pivot,
                sqrt(factor_error2/max(matrix_norm2,1.0e-300)));
        fflush(E->fp);
        fflush(stderr);
    }

    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        free((void *)coarse_p[m]);
        free((void *)coarse_Ap[m]);
        free((void *)fine_p[m]);
        free((void *)fine_Ap[m]);
        free((void *)fine_velocity[m]);
        free((void *)fine_velocity_rhs[m]);
        free((void *)fine_velocity_Ax[m]);
        free((void *)fine_velocity_direction[m]);
    }
    free((void *)local_column);
    free((void *)global_column);
}


static int ala_geneo_bin(int ex, int ey, int ez, int elx, int ely,
                         int shallow_min_ez, int shallow_layers,
                         int horizontal_bins, int radial_bins)
{
    int bx,by,bz;
    if(ez<shallow_min_ez || shallow_layers<=0)
        return(-1);
    bx=min((ex-1)*horizontal_bins/elx,horizontal_bins-1);
    by=min((ey-1)*horizontal_bins/ely,horizontal_bins-1);
    bz=min((ez-shallow_min_ez)*radial_bins/shallow_layers,radial_bins-1);
    return(bz+radial_bins*(bx+horizontal_bins*by));
}


/* Stage 6d deterministic coarse space.  The radial-bin indicators are
   disjoint on each cross-rank processor aggregate and together span its
   complete shallow constant mode.  Final normalization is performed after
   transfer to the Galerkin pressure level. */
static int ala_build_radial_partition_shapes(const double *aggregate,
    int nbins, int rbins, double *selected_shapes,
    double *selected_values)
{
    int mode,i,j;
    double rayleigh_numerator,rayleigh_denominator;

    for(mode=0;mode<rbins;mode++) {
        for(j=0;j<nbins;j++) {
            selected_shapes[mode*nbins+j]
                =(j%rbins==mode) ? 1.0 : 0.0;
        }
        /* These are geometric indicators, not selected GenEO modes.  Their
           existence must not depend on the selection operator's active map
           or dimensional scale.  The transferred pressure-level support and
           the final Galerkin Cholesky are the authoritative validity tests. */
        rayleigh_numerator=0.0;
        rayleigh_denominator=0.0;
        for(i=0;i<nbins;i++) {
            rayleigh_denominator += aggregate[i*nbins+i]
                *selected_shapes[mode*nbins+i]
                *selected_shapes[mode*nbins+i];
            for(j=0;j<nbins;j++)
                rayleigh_numerator += selected_shapes[mode*nbins+i]
                    *aggregate[i*nbins+j]
                    *selected_shapes[mode*nbins+j];
        }
        selected_values[mode]
            =(isfinite(rayleigh_numerator)
              && isfinite(rayleigh_denominator)
              && fabs(rayleigh_denominator)>1.0e-300)
            ? rayleigh_numerator/rayleigh_denominator : 0.0;
    }
    return(rbins);
}


/* Dense symmetric Jacobi eigensolver for the small local GenEO aggregate
   problem.  Eigenvectors are returned by columns. */
static int ala_geneo_jacobi_eigensolve(double *a, double *vectors,
                                       double *eigenvalues, int n)
{
    int i,j,k,p,q,iteration,max_iterations;
    double offdiag,maximum,app,aqq,apq,tau,t,c,s,akp,akq,vkp,vkq;
    const double tolerance=1.0e-12;

    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            vectors[i*n+j]=(i==j) ? 1.0 : 0.0;
    max_iterations=max(50*n*n,1);
    for(iteration=0;iteration<max_iterations;iteration++) {
        maximum=0.0;
        p=0;
        q=0;
        for(i=0;i<n;i++)
            for(j=i+1;j<n;j++) {
                offdiag=fabs(a[i*n+j]);
                if(offdiag>maximum) {
                    maximum=offdiag;
                    p=i;
                    q=j;
                }
            }
        if(maximum<=tolerance)
            break;
        app=a[p*n+p];
        aqq=a[q*n+q];
        apq=a[p*n+q];
        tau=(aqq-app)/(2.0*apq);
        t=((tau>=0.0) ? 1.0 : -1.0)
            /(fabs(tau)+sqrt(1.0+tau*tau));
        c=1.0/sqrt(1.0+t*t);
        s=t*c;
        for(k=0;k<n;k++)
            if(k!=p && k!=q) {
                akp=a[k*n+p];
                akq=a[k*n+q];
                a[k*n+p]=c*akp-s*akq;
                a[p*n+k]=a[k*n+p];
                a[k*n+q]=s*akp+c*akq;
                a[q*n+k]=a[k*n+q];
            }
        a[p*n+p]=c*c*app-2.0*s*c*apq+s*s*aqq;
        a[q*n+q]=s*s*app+2.0*s*c*apq+c*c*aqq;
        a[p*n+q]=0.0;
        a[q*n+p]=0.0;
        for(k=0;k<n;k++) {
            vkp=vectors[k*n+p];
            vkq=vectors[k*n+q];
            vectors[k*n+p]=c*vkp-s*vkq;
            vectors[k*n+q]=s*vkp+c*vkq;
        }
    }
    for(i=0;i<n;i++)
        eigenvalues[i]=a[i*n+i];
    return(iteration);
}


/* Build a subdomain-adaptive shallow pressure space.  Each rank compresses
   its shallow cells into a small 4x4x2-style aggregate problem and solves
   A_i phi=lambda diag(A_i) phi.  The selected low-energy modes have disjoint
   rank ownership; their complete global coupling is then assembled with the
   same fixed strict-ALA Galerkin Schur map used by the legacy coarse space. */
static void build_ala_geneo_coarse_cache(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache, int lev)
{
    int m,i,j,k,col,row,e1,e2,ex1,ey1,ez1,ex2,ey2,ez2;
    int dx,dy,dz,bin1,bin2,active1,active2,nbins,nactive;
    int hbins,rbins,elx,ely,elz,shallow_min_ez,shallow_layers;
    int clev,factor,celx,cely,celz,cnpno,stride,npno,neq;
    int cx,cy,cz,ce,local_modes,n,iterations,global_iterations;
    int local_active_rank,global_active_ranks;
    int local_mode_min,local_mode_max,global_mode_min,global_mode_max;
    int *active_map,*mode_counts,*mode_offsets;
    double depth,value,norm,local_eigen_min,local_eigen_max;
    double global_eigen_min,global_eigen_max,threshold;
    double local_second_min,local_second_max;
    double global_second_min,global_second_max;
    double *aggregate,*normalized,*vectors,*eigenvalues,*local_selected;
    double *local_column,*global_column,*matrix;
    double *coarse_p[NCS],*coarse_Ap[NCS];
    double *fine_p[NCS],*fine_Ap[NCS],*fine_velocity[NCS];
    double *fine_velocity_rhs[NCS],*fine_velocity_Ax[NCS];
    double *fine_velocity_direction[NCS];
    double anti2,total2,symmetry,diagonal_min,diagonal_max;
    double shift,sum,pivot,min_pivot,setup_start,setup_seconds;
    double global_setup_seconds;
    double CPU_time0();

    setup_start=CPU_time0();
    hbins=E->control.ala_geneo_horizontal_bins;
    rbins=E->control.ala_geneo_radial_bins;
    nbins=hbins*hbins*rbins;
    if(nbins>ALA_GENEO_MAX_BINS)
        myerror(E,"ALA GenEO aggregate dimension exceeds compiled maximum");
    elx=E->lmesh.ELX[lev];
    ely=E->lmesh.ELY[lev];
    elz=E->lmesh.ELZ[lev];
    shallow_min_ez=elz+1;
    for(ez1=elz;ez1>=1;ez1--) {
        depth=(1.0-0.5*(E->sx[1][3][ez1]+E->sx[1][3][ez1+1]))
            *E->data.radius_km;
        if(depth>E->control.ala_shallow_patch_depth_km)
            break;
        shallow_min_ez=ez1;
    }
    shallow_layers=(shallow_min_ez<=elz) ? elz-shallow_min_ez+1 : 0;

    aggregate=(double *)calloc(nbins*nbins,sizeof(double));
    active_map=(int *)calloc(nbins,sizeof(int));
    normalized=(double *)calloc(nbins*nbins,sizeof(double));
    vectors=(double *)calloc(nbins*nbins,sizeof(double));
    eigenvalues=(double *)calloc(nbins,sizeof(double));
    local_selected=(double *)calloc(
        E->control.ala_geneo_max_modes_per_rank,sizeof(double));
    mode_counts=(int *)calloc(E->parallel.nproc,sizeof(int));
    mode_offsets=(int *)calloc(E->parallel.nproc,sizeof(int));
    if(aggregate==NULL || active_map==NULL || normalized==NULL ||
       vectors==NULL || eigenvalues==NULL || local_selected==NULL ||
       mode_counts==NULL || mode_offsets==NULL)
        myerror(E,"Unable to allocate ALA GenEO local eigenproblem");

    if(shallow_layers>0) {
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(ey1=1;ey1<=ely;ey1++)
                for(ex1=1;ex1<=elx;ex1++)
                    for(ez1=shallow_min_ez;ez1<=elz;ez1++) {
                        e1=ez1+(ex1-1)*elz+(ey1-1)*elz*elx;
                        bin1=ala_geneo_bin(ex1,ey1,ez1,elx,ely,
                            shallow_min_ez,shallow_layers,hbins,rbins);
                        for(dy=-1;dy<=1;dy++) {
                            ey2=ey1+dy;
                            if(ey2<1 || ey2>ely)
                                continue;
                            for(dx=-1;dx<=1;dx++) {
                                ex2=ex1+dx;
                                if(ex2<1 || ex2>elx)
                                    continue;
                                for(dz=-1;dz<=1;dz++) {
                                    ez2=ez1+dz;
                                    if(ez2<shallow_min_ez || ez2>elz)
                                        continue;
                                    e2=ez2+(ex2-1)*elz
                                        +(ey2-1)*elz*elx;
                                    bin2=ala_geneo_bin(ex2,ey2,ez2,elx,ely,
                                        shallow_min_ez,shallow_layers,
                                        hbins,rbins);
                                    aggregate[bin1*nbins+bin2] +=
                                        assemble_Ahatp_jacobi_entry(
                                            E,e1,e2,lev,m);
                                }
                            }
                        }
                    }
    }
    nactive=0;
    for(i=0;i<nbins;i++) {
        active_map[i]=-1;
        if(isfinite(aggregate[i*nbins+i]) &&
           aggregate[i*nbins+i]>0.0)
            active_map[i]=nactive++;
    }
    for(i=0;i<nbins;i++)
        if(active_map[i]>=0)
            for(j=0;j<nbins;j++)
                if(active_map[j]>=0) {
                    active1=active_map[i];
                    active2=active_map[j];
                    value=0.5*(aggregate[i*nbins+j]
                               +aggregate[j*nbins+i]);
                    normalized[active1*nactive+active2]=value
                        /sqrt(aggregate[i*nbins+i]
                              *aggregate[j*nbins+j]);
                    if(!isfinite(normalized[active1*nactive+active2]))
                        myerror(E,"ALA GenEO local operator is not finite");
                }
    iterations=0;
    if(nactive>0)
        iterations=ala_geneo_jacobi_eigensolve(
            normalized,vectors,eigenvalues,nactive);
    if(nactive>0 && iterations>=max(50*nactive*nactive,1))
        myerror(E,"ALA GenEO local eigensolve did not converge");
    for(i=0;i<nactive;i++)
        if(!isfinite(eigenvalues[i]))
            myerror(E,"ALA GenEO local eigenvalue is not finite");
    for(i=0;i<nactive;i++)
        for(j=i+1;j<nactive;j++)
            if(eigenvalues[j]<eigenvalues[i]) {
                value=eigenvalues[i];
                eigenvalues[i]=eigenvalues[j];
                eigenvalues[j]=value;
                for(k=0;k<nactive;k++) {
                    value=vectors[k*nactive+i];
                    vectors[k*nactive+i]=vectors[k*nactive+j];
                    vectors[k*nactive+j]=value;
                }
            }
    local_modes=0;
    threshold=E->control.ala_geneo_eigenvalue_threshold;
    while(local_modes<nactive &&
          local_modes<E->control.ala_geneo_max_modes_per_rank &&
          eigenvalues[local_modes]<=threshold)
        local_modes++;
    if(nactive>0)
        local_modes=max(local_modes,
                        E->control.ala_geneo_min_modes_per_rank);
    local_modes=min(local_modes,nactive);
    local_modes=min(local_modes,E->control.ala_geneo_max_modes_per_rank);
    for(i=0;i<local_modes;i++)
        local_selected[i]=eigenvalues[i];

    MPI_Allgather(&local_modes,1,MPI_INT,mode_counts,1,MPI_INT,
                  E->parallel.world);
    n=0;
    for(i=0;i<E->parallel.nproc;i++) {
        mode_offsets[i]=n;
        n += mode_counts[i];
    }
    if(n<1)
        myerror(E,"ALA GenEO selected no global coarse modes");
    if(n>E->control.ala_geneo_max_global_modes)
        myerror(E,"ALA GenEO global mode count exceeds configured maximum");
    cache->geneo_basis_count=n;
    cache->geneo_local_modes=local_modes;
    cache->geneo_local_offset=mode_offsets[E->parallel.me];
    cache->geneo_eigenvalues=(double *)calloc(n,sizeof(double));
    if(cache->geneo_eigenvalues==NULL)
        myerror(E,"Unable to allocate ALA GenEO eigenvalue diagnostics");
    MPI_Allgatherv(local_selected,local_modes,MPI_DOUBLE,
                   cache->geneo_eigenvalues,mode_counts,mode_offsets,
                   MPI_DOUBLE,E->parallel.world);

    clev=lev-E->control.ala_two_level_offset;
    factor=1 << E->control.ala_two_level_offset;
    celx=E->lmesh.ELX[clev];
    cely=E->lmesh.ELY[clev];
    celz=E->lmesh.ELZ[clev];
    cnpno=E->lmesh.NPNO[clev];
    stride=cnpno+1;
    npno=E->lmesh.NPNO[lev];
    neq=E->lmesh.NEQ[lev];
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        cache->geneo_basis[m]=(double *)calloc(n*stride,sizeof(double));
        if(cache->geneo_basis[m]==NULL)
            myerror(E,"Unable to allocate ALA GenEO basis");
    }
    for(i=0;i<local_modes;i++) {
        col=cache->geneo_local_offset+i;
        norm=0.0;
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(cy=1;cy<=cely;cy++)
                for(cx=1;cx<=celx;cx++)
                    for(cz=1;cz<=celz;cz++) {
                        ce=cz+(cx-1)*celz+(cy-1)*celz*celx;
                        depth=(1.0-E->ECO[clev][m][ce].centre[3])
                            *E->data.radius_km;
                        if(depth>E->control.ala_shallow_patch_depth_km)
                            continue;
                        ex1=min((cx-1)*factor+1,elx);
                        ey1=min((cy-1)*factor+1,ely);
                        ez1=min(cz*factor,elz);
                        bin1=ala_geneo_bin(ex1,ey1,ez1,elx,ely,
                            shallow_min_ez,shallow_layers,hbins,rbins);
                        if(bin1<0 || active_map[bin1]<0)
                            continue;
                        active1=active_map[bin1];
                        value=vectors[active1*nactive+i]
                            /sqrt(aggregate[bin1*nbins+bin1]);
                        cache->geneo_basis[m][col*stride+ce]=value;
                        norm += value*value;
                    }
        if(!isfinite(norm) || norm<=1.0e-30)
            myerror(E,"ALA GenEO selected basis has zero norm");
        norm=sqrt(norm);
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(ce=1;ce<=cnpno;ce++)
                cache->geneo_basis[m][col*stride+ce] /= norm;
    }

    cache->geneo_matrix=(double *)calloc(n*n,sizeof(double));
    cache->geneo_chol=(double *)calloc(n*n,sizeof(double));
    local_column=(double *)calloc(n,sizeof(double));
    global_column=(double *)calloc(n,sizeof(double));
    if(cache->geneo_matrix==NULL || cache->geneo_chol==NULL ||
       local_column==NULL || global_column==NULL)
        myerror(E,"Unable to allocate ALA GenEO global matrix");
    matrix=cache->geneo_matrix;
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        coarse_p[m]=(double *)calloc(cnpno+1,sizeof(double));
        coarse_Ap[m]=(double *)calloc(cnpno+1,sizeof(double));
        fine_p[m]=(double *)calloc(npno+1,sizeof(double));
        fine_Ap[m]=(double *)calloc(npno+1,sizeof(double));
        fine_velocity[m]=(double *)calloc(neq+1,sizeof(double));
        fine_velocity_rhs[m]=(double *)calloc(neq+1,sizeof(double));
        fine_velocity_Ax[m]=(double *)calloc(neq+1,sizeof(double));
        fine_velocity_direction[m]=(double *)calloc(neq+1,sizeof(double));
        if(coarse_p[m]==NULL || coarse_Ap[m]==NULL || fine_p[m]==NULL ||
           fine_Ap[m]==NULL || fine_velocity[m]==NULL ||
           fine_velocity_rhs[m]==NULL || fine_velocity_Ax[m]==NULL ||
           fine_velocity_direction[m]==NULL)
            myerror(E,"Unable to allocate ALA GenEO Galerkin workspace");
    }
    for(col=0;col<n;col++) {
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(ce=1;ce<=cnpno;ce++)
                coarse_p[m][ce]=cache->geneo_basis[m][col*stride+ce];
        apply_ala_galerkin_fixed_schur(
            E,coarse_p,coarse_Ap,fine_p,fine_Ap,fine_velocity,
            fine_velocity_rhs,fine_velocity_Ax,fine_velocity_direction,
            cache,clev,lev,factor,NULL);
        for(row=0;row<n;row++)
            local_column[row]=0.0;
        for(row=cache->geneo_local_offset;
            row<cache->geneo_local_offset+local_modes;row++)
            for(m=1;m<=E->sphere.caps_per_proc;m++)
                for(ce=1;ce<=cnpno;ce++)
                    local_column[row] +=
                        cache->geneo_basis[m][row*stride+ce]
                        *coarse_Ap[m][ce];
        MPI_Allreduce(local_column,global_column,n,MPI_DOUBLE,MPI_SUM,
                      E->parallel.world);
        for(row=0;row<n;row++)
            matrix[row*n+col]=global_column[row];
        if(col==0 || (col+1)%25==0 || col+1==n) {
            setup_seconds=CPU_time0()-setup_start;
            MPI_Allreduce(&setup_seconds,&global_setup_seconds,1,MPI_DOUBLE,
                          MPI_MAX,E->parallel.world);
            if(E->parallel.me==0) {
                fprintf(E->fp,"ALA_GENEO_SETUP_PROGRESS column=%d/%d "
                        "elapsed_seconds_max=%e\n",
                        col+1,n,global_setup_seconds);
                fprintf(stderr,"ALA_GENEO_SETUP_PROGRESS column=%d/%d "
                        "elapsed_seconds_max=%e\n",
                        col+1,n,global_setup_seconds);
                fflush(E->fp);
                fflush(stderr);
            }
        }
    }
    cache->geneo_schur_applications=n;

    anti2=0.0;
    total2=0.0;
    diagonal_min=1.0e300;
    diagonal_max=0.0;
    for(i=0;i<n;i++) {
        diagonal_min=min(diagonal_min,matrix[i*n+i]);
        diagonal_max=max(diagonal_max,matrix[i*n+i]);
        for(j=0;j<n;j++) {
            anti2 += (matrix[i*n+j]-matrix[j*n+i])
                *(matrix[i*n+j]-matrix[j*n+i]);
            total2 += matrix[i*n+j]*matrix[i*n+j];
            cache->geneo_chol[i*n+j]=
                0.5*(matrix[i*n+j]+matrix[j*n+i]);
        }
    }
    symmetry=sqrt(anti2/max(total2,1.0e-300));
    if(!isfinite(symmetry) || symmetry>1.0e-8)
        myerror(E,"ALA GenEO Galerkin matrix is not symmetric");
    if(!isfinite(diagonal_min) || diagonal_min<=0.0)
        myerror(E,"ALA GenEO Galerkin diagonal is not positive");
    shift=E->control.ala_geneo_regularization*diagonal_max;
    for(i=0;i<n;i++)
        cache->geneo_chol[i*n+i] += shift;
    min_pivot=1.0e300;
    for(i=0;i<n;i++) {
        for(j=0;j<=i;j++) {
            sum=cache->geneo_chol[i*n+j];
            for(k=0;k<j;k++)
                sum -= cache->geneo_chol[i*n+k]
                    *cache->geneo_chol[j*n+k];
            if(i==j) {
                pivot=sum;
                if(!isfinite(pivot) ||
                   pivot<=1.0e-14*max(diagonal_max,1.0e-300))
                    myerror(E,"ALA GenEO Galerkin Cholesky failed");
                cache->geneo_chol[i*n+j]=sqrt(pivot);
                min_pivot=min(min_pivot,pivot);
            }
            else
                cache->geneo_chol[i*n+j]=
                    sum/cache->geneo_chol[j*n+j];
        }
        for(j=i+1;j<n;j++)
            cache->geneo_chol[i*n+j]=0.0;
    }
    local_eigen_min=(local_modes>0) ? local_selected[0] : 1.0e300;
    local_eigen_max=(local_modes>0) ? local_selected[local_modes-1] : 0.0;
    local_second_min=(nactive>1) ? eigenvalues[1] : 1.0e300;
    local_second_max=(nactive>1) ? eigenvalues[1] : 0.0;
    MPI_Allreduce(&local_eigen_min,&global_eigen_min,1,MPI_DOUBLE,MPI_MIN,
                  E->parallel.world);
    MPI_Allreduce(&local_eigen_max,&global_eigen_max,1,MPI_DOUBLE,MPI_MAX,
                  E->parallel.world);
    MPI_Allreduce(&local_second_min,&global_second_min,1,MPI_DOUBLE,MPI_MIN,
                  E->parallel.world);
    MPI_Allreduce(&local_second_max,&global_second_max,1,MPI_DOUBLE,MPI_MAX,
                  E->parallel.world);
    local_active_rank=(local_modes>0) ? 1 : 0;
    local_mode_min=(local_modes>0) ? local_modes : 1000000000;
    local_mode_max=local_modes;
    MPI_Allreduce(&local_active_rank,&global_active_ranks,1,MPI_INT,MPI_SUM,
                  E->parallel.world);
    MPI_Allreduce(&local_mode_min,&global_mode_min,1,MPI_INT,MPI_MIN,
                  E->parallel.world);
    MPI_Allreduce(&local_mode_max,&global_mode_max,1,MPI_INT,MPI_MAX,
                  E->parallel.world);
    MPI_Allreduce(&iterations,&global_iterations,1,MPI_INT,MPI_MAX,
                  E->parallel.world);
    setup_seconds=CPU_time0()-setup_start;
    MPI_Allreduce(&setup_seconds,&global_setup_seconds,1,MPI_DOUBLE,MPI_MAX,
                  E->parallel.world);
    if(E->parallel.me==0) {
        fprintf(E->fp,"ALA GenEO coarse space global_modes=%d "
                "active_ranks=%d local_mode_range=(%d,%d) "
                "aggregate_bins=%dx%dx%d "
                "threshold=%e selected_eigen_range=(%e,%e) "
                "second_eigen_range=(%e,%e) "
                "local_jacobi_iterations=%d Galerkin_applications=%d "
                "symmetry_defect=%e diagonal_range=(%e,%e) "
                "regularization=%e shift=%e min_pivot=%e "
                "setup_seconds_max=%e status=pass\n",
                n,global_active_ranks,global_mode_min,global_mode_max,
                hbins,hbins,rbins,threshold,global_eigen_min,global_eigen_max,
                global_second_min,global_second_max,global_iterations,n,
                symmetry,diagonal_min,diagonal_max,
                E->control.ala_geneo_regularization,shift,min_pivot,
                global_setup_seconds);
        fprintf(stderr,"ALA GenEO coarse space global_modes=%d "
                "aggregate_bins=%dx%dx%d threshold=%e "
                "selected_eigen_range=(%e,%e) second_eigen_range=(%e,%e) "
                "Galerkin_applications=%d "
                "symmetry_defect=%e setup_seconds_max=%e status=pass\n",
                n,hbins,hbins,rbins,threshold,global_eigen_min,
                global_eigen_max,global_second_min,global_second_max,n,
                symmetry,global_setup_seconds);
        fflush(E->fp);
        fflush(stderr);
    }
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        free((void *)coarse_p[m]);
        free((void *)coarse_Ap[m]);
        free((void *)fine_p[m]);
        free((void *)fine_Ap[m]);
        free((void *)fine_velocity[m]);
        free((void *)fine_velocity_rhs[m]);
        free((void *)fine_velocity_Ax[m]);
        free((void *)fine_velocity_direction[m]);
    }
    free((void *)aggregate);
    free((void *)active_map);
    free((void *)normalized);
    free((void *)vectors);
    free((void *)eigenvalues);
    free((void *)local_selected);
    free((void *)mode_counts);
    free((void *)mode_offsets);
    free((void *)local_column);
    free((void *)global_column);
}


/* Build a pressure coarse space whose spectral aggregates span rectangular
   groups of horizontal MPI ranks.  Each group owns a disjoint set of pressure
   cells, so its indicator functions form an exact partition of unity.  The
   selected group modes are present on every member rank and the final coarse
   matrix is still assembled by the complete fixed strict-ALA Galerkin map.

   The group spectral operator is the principal restriction of the complete
   diagonal-velocity Schur approximation G diag(K)^-1 G^T.  Its rank-local
   terms use assemble_Ahatp_jacobi_entry(); its cross-rank terms use the exact
   G rows and diagonal velocity inverse weights already exchanged for
   MPI-overlap Schwarz.  Thus mode selection sees the real
   processor-interface coupling instead of a block-diagonal sum of rank-local
   principal matrices. */
static void build_ala_cross_rank_geneo_coarse_cache(
    struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache, int lev)
{
    int m,i,j,k,p,col,row,e1,e2,ex1,ey1,ez1,ex2,ey2,ez2;
    int dx,dy,dz,bin1,bin2,active1,active2,nactive;
    int hbins,rbins,elx,ely,elz,shallow_min_ez,shallow_layers;
    int rank_group_x,rank_group_y,groups_x,groups_y;
    int group_ix,group_iy,group_start_x,group_start_y;
    int group_size_x,group_size_y,rank_offset_x,rank_offset_y;
    int group_hbins_x,group_hbins_y,nbins,color,group_rank,group_size;
    int local_has_support,supported_members,face,ghost_index;
    int tangent_count,ghost_q,ghost_t,ghost_ez;
    int ghost_ex,ghost_ey,ghost_rank_offset_x,ghost_rank_offset_y;
    int face_is_internal,local_boundary;
    int local_interface_entries,global_interface_entries;
    int clev,factor,celx,cely,celz,cnpno,stride,npno,neq;
    int cx,cy,cz,ce,local_modes,n,iterations,global_iterations;
    int desired_modes,candidate,accepted;
    int owned_modes,owned_cursor,local_active,global_active;
    int local_mode_min,local_mode_max,global_mode_min,global_mode_max;
    int local_group_min,local_group_max,global_group_min,global_group_max;
    int *active_map,*mode_counts,*mode_offsets;
    int group_modes[NCS],group_offsets[NCS],group_nbins[NCS];
    int group_hx[NCS],group_rank_offset_x[NCS];
    int group_rank_offset_y[NCS];
    MPI_Comm group_comm[NCS];
    double depth,value,norm,group_norm,threshold;
    double interface_value;
    const double *ghost_record;
    double coefficient,denominator,rayleigh_numerator,rayleigh_denominator;
    double local_eigen_min,local_eigen_max,global_eigen_min,global_eigen_max;
    double local_second_min,local_second_max;
    double global_second_min,global_second_max;
    double selection_anti2,selection_total2,group_selection_symmetry;
    double local_selection_symmetry_max,global_selection_symmetry_max;
    double *local_aggregate,*aggregate,*normalized,*vectors,*eigenvalues;
    double *selected_values[NCS],*selected_shapes[NCS];
    double *local_eigenvalues,*local_column,*global_column,*matrix;
    double *coarse_p[NCS],*coarse_Ap[NCS];
    double *fine_p[NCS],*fine_Ap[NCS],*fine_velocity[NCS];
    double *fine_velocity_rhs[NCS],*fine_velocity_Ax[NCS];
    double *fine_velocity_direction[NCS];
    double anti2,total2,symmetry,diagonal_min,diagonal_max;
    double shift,sum,pivot,min_pivot,setup_start,setup_seconds;
    double global_setup_seconds;
    double CPU_time0();

    setup_start=CPU_time0();
    hbins=E->control.ala_geneo_horizontal_bins;
    rbins=E->control.ala_geneo_radial_bins;
    rank_group_x=E->control.ala_geneo_rank_group_x;
    rank_group_y=E->control.ala_geneo_rank_group_y;
    groups_x=(E->parallel.nprocx+rank_group_x-1)/rank_group_x;
    groups_y=(E->parallel.nprocy+rank_group_y-1)/rank_group_y;
    elx=E->lmesh.ELX[lev];
    ely=E->lmesh.ELY[lev];
    elz=E->lmesh.ELZ[lev];
    shallow_min_ez=elz+1;
    for(ez1=elz;ez1>=1;ez1--) {
        depth=(1.0-0.5*(E->sx[1][3][ez1]+E->sx[1][3][ez1+1]))
            *E->data.radius_km;
        if(depth>E->control.ala_shallow_patch_depth_km)
            break;
        shallow_min_ez=ez1;
    }
    shallow_layers=(shallow_min_ez<=elz) ? elz-shallow_min_ez+1 : 0;
    threshold=E->control.ala_geneo_eigenvalue_threshold;
    owned_modes=0;
    local_active=0;
    local_mode_min=1000000000;
    local_mode_max=0;
    local_group_min=1000000000;
    local_group_max=0;
    local_eigen_min=1.0e300;
    local_eigen_max=0.0;
    local_second_min=1.0e300;
    local_second_max=0.0;
    global_iterations=0;
    local_interface_entries=0;
    local_selection_symmetry_max=0.0;
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        group_comm[m]=MPI_COMM_NULL;
        group_modes[m]=0;
        group_offsets[m]=0;
        selected_values[m]=NULL;
        selected_shapes[m]=NULL;
        group_ix=E->parallel.me_loc[1]/rank_group_x;
        group_iy=E->parallel.me_loc[2]/rank_group_y;
        group_start_x=group_ix*rank_group_x;
        group_start_y=group_iy*rank_group_y;
        group_size_x=min(rank_group_x,E->parallel.nprocx-group_start_x);
        group_size_y=min(rank_group_y,E->parallel.nprocy-group_start_y);
        rank_offset_x=E->parallel.me_loc[1]-group_start_x;
        rank_offset_y=E->parallel.me_loc[2]-group_start_y;
        group_hbins_x=hbins*group_size_x;
        group_hbins_y=hbins*group_size_y;
        nbins=group_hbins_x*group_hbins_y*rbins;
        if(nbins>ALA_GENEO_MAX_BINS)
            myerror(E,"ALA cross-rank GenEO aggregate exceeds compiled maximum");
        color=(((E->sphere.capid[m]-1)*E->parallel.nprocz
               +E->parallel.me_loc[3])*groups_y+group_iy)
               *groups_x+group_ix;
        MPI_Comm_split(E->parallel.world,color,E->parallel.me,
                       &group_comm[m]);
        MPI_Comm_rank(group_comm[m],&group_rank);
        MPI_Comm_size(group_comm[m],&group_size);
        if(group_size!=group_size_x*group_size_y)
            myerror(E,"ALA cross-rank GenEO processor group is incomplete");
        group_nbins[m]=nbins;
        group_hx[m]=group_hbins_x;
        group_rank_offset_x[m]=rank_offset_x;
        group_rank_offset_y[m]=rank_offset_y;
        local_aggregate=(double *)calloc(nbins*nbins,sizeof(double));
        aggregate=(double *)calloc(nbins*nbins,sizeof(double));
        active_map=(int *)calloc(nbins,sizeof(int));
        normalized=(double *)calloc(nbins*nbins,sizeof(double));
        vectors=(double *)calloc(nbins*nbins,sizeof(double));
        eigenvalues=(double *)calloc(nbins,sizeof(double));
        selected_values[m]=(double *)calloc(
            E->control.ala_geneo_max_modes_per_rank,sizeof(double));
        selected_shapes[m]=(double *)calloc(
            E->control.ala_geneo_max_modes_per_rank*nbins,sizeof(double));
        if(local_aggregate==NULL || aggregate==NULL || active_map==NULL ||
           normalized==NULL || vectors==NULL || eigenvalues==NULL ||
           selected_values[m]==NULL || selected_shapes[m]==NULL)
            myerror(E,"Unable to allocate cross-rank GenEO eigenproblem");
        if(shallow_layers>0)
            for(ey1=1;ey1<=ely;ey1++)
                for(ex1=1;ex1<=elx;ex1++)
                    for(ez1=shallow_min_ez;ez1<=elz;ez1++) {
                        e1=ez1+(ex1-1)*elz+(ey1-1)*elz*elx;
                        bin1=(ez1-shallow_min_ez)*rbins/shallow_layers
                            +rbins*((rank_offset_x*hbins
                            +min((ex1-1)*hbins/elx,hbins-1))
                            +group_hbins_x*(rank_offset_y*hbins
                            +min((ey1-1)*hbins/ely,hbins-1)));
                        for(dy=-1;dy<=1;dy++) {
                            ey2=ey1+dy;
                            if(ey2<1 || ey2>ely)
                                continue;
                            for(dx=-1;dx<=1;dx++) {
                                ex2=ex1+dx;
                                if(ex2<1 || ex2>elx)
                                    continue;
                                for(dz=-1;dz<=1;dz++) {
                                    ez2=ez1+dz;
                                    if(ez2<shallow_min_ez || ez2>elz)
                                        continue;
                                    e2=ez2+(ex2-1)*elz
                                        +(ey2-1)*elz*elx;
                                    bin2=(ez2-shallow_min_ez)*rbins
                                            /shallow_layers
                                        +rbins*((rank_offset_x*hbins
                                        +min((ex2-1)*hbins/elx,hbins-1))
                                        +group_hbins_x*(rank_offset_y*hbins
                                        +min((ey2-1)*hbins/ely,hbins-1)));
                                    local_aggregate[bin1*nbins+bin2] +=
                                        assemble_Ahatp_jacobi_entry(
                                            E,e1,e2,lev,m);
                                }
                            }
                        }
                    }

        /* Complete the group principal matrix across MPI interfaces.  Halo
           records contain the full strict-ALA G=D+C element row and the
           diagonal velocity inverse weights at every element node.  A local
           boundary row contributes A(local,ghost); the neighboring rank
           contributes the transpose entry, and the group reduction below
           combines both.  Only the first ghost layer can share velocity nodes
           with a local pressure element, but all tangential and radial
           neighbors are tested geometrically to preserve the element stencil
           exactly. */
        if(shallow_layers>0)
            for(face=0;face<ALA_PATCH_MPI_FACES;face++) {
                face_is_internal=
                    (face==0 && rank_offset_x>0)
                    || (face==1 && rank_offset_x<group_size_x-1)
                    || (face==2 && rank_offset_y>0)
                    || (face==3 && rank_offset_y<group_size_y-1);
                if(!face_is_internal)
                    continue;
                tangent_count=(face<2) ? ely : elx;
                local_boundary=(face==0 || face==2) ? 1
                                                     : ((face==1) ? elx
                                                                  : ely);
                for(ghost_index=0;
                    ghost_index<cache->halo_recv_count[m][face];
                    ghost_index++) {
                    ghost_record=cache->halo_recv_records[m][face]
                        +ghost_index*ALA_HALO_ELEMENT_RECORD_DOUBLES;
                    ghost_q=(int)(ghost_record[0]+0.5);
                    if(ghost_q!=0)
                        continue;
                    ghost_t=(ghost_index/shallow_layers)%tangent_count+1;
                    ghost_ez=shallow_min_ez
                        +ghost_index%shallow_layers;
                    ghost_rank_offset_x=rank_offset_x;
                    ghost_rank_offset_y=rank_offset_y;
                    if(face==0) {
                        ghost_rank_offset_x--;
                        ghost_ex=elx;
                        ghost_ey=ghost_t;
                    }
                    else if(face==1) {
                        ghost_rank_offset_x++;
                        ghost_ex=1;
                        ghost_ey=ghost_t;
                    }
                    else if(face==2) {
                        ghost_rank_offset_y--;
                        ghost_ex=ghost_t;
                        ghost_ey=ely;
                    }
                    else {
                        ghost_rank_offset_y++;
                        ghost_ex=ghost_t;
                        ghost_ey=1;
                    }
                    bin2=(ghost_ez-shallow_min_ez)*rbins/shallow_layers
                        +rbins*((ghost_rank_offset_x*hbins
                        +min((ghost_ex-1)*hbins/elx,hbins-1))
                        +group_hbins_x*(ghost_rank_offset_y*hbins
                        +min((ghost_ey-1)*hbins/ely,hbins-1)));
                    if(face<2) {
                        ex1=local_boundary;
                        for(ey1=max(1,ghost_t-1);
                            ey1<=min(ely,ghost_t+1);ey1++)
                            for(ez1=max(shallow_min_ez,ghost_ez-1);
                                ez1<=min(elz,ghost_ez+1);ez1++) {
                                bin1=(ez1-shallow_min_ez)*rbins
                                        /shallow_layers
                                    +rbins*((rank_offset_x*hbins
                                    +min((ex1-1)*hbins/elx,hbins-1))
                                    +group_hbins_x*(rank_offset_y*hbins
                                    +min((ey1-1)*hbins/ely,hbins-1)));
                                interface_value=ala_halo_record_coupling(
                                    E,
                                    cache->halo_send_records[m][face]
                                      +(ez1-shallow_min_ez
                                        +(ey1-1)*shallow_layers)
                                       *ALA_HALO_ELEMENT_RECORD_DOUBLES,
                                    ghost_record);
                                if(interface_value!=0.0) {
                                    local_aggregate[bin1*nbins+bin2]
                                        +=interface_value;
                                    local_interface_entries++;
                                }
                            }
                    }
                    else {
                        ey1=local_boundary;
                        for(ex1=max(1,ghost_t-1);
                            ex1<=min(elx,ghost_t+1);ex1++)
                            for(ez1=max(shallow_min_ez,ghost_ez-1);
                                ez1<=min(elz,ghost_ez+1);ez1++) {
                                bin1=(ez1-shallow_min_ez)*rbins
                                        /shallow_layers
                                    +rbins*((rank_offset_x*hbins
                                    +min((ex1-1)*hbins/elx,hbins-1))
                                    +group_hbins_x*(rank_offset_y*hbins
                                    +min((ey1-1)*hbins/ely,hbins-1)));
                                interface_value=ala_halo_record_coupling(
                                    E,
                                    cache->halo_send_records[m][face]
                                      +(ez1-shallow_min_ez
                                        +(ex1-1)*shallow_layers)
                                       *ALA_HALO_ELEMENT_RECORD_DOUBLES,
                                    ghost_record);
                                if(interface_value!=0.0) {
                                    local_aggregate[bin1*nbins+bin2]
                                        +=interface_value;
                                    local_interface_entries++;
                                }
                            }
                    }
                }
            }
        MPI_Allreduce(local_aggregate,aggregate,nbins*nbins,MPI_DOUBLE,
                      MPI_SUM,group_comm[m]);
        iterations=0;
        nactive=0;
        if(group_rank==0) {
            selection_anti2=0.0;
            selection_total2=0.0;
            for(i=0;i<nbins;i++)
                for(j=0;j<nbins;j++) {
                    selection_anti2 +=
                        (aggregate[i*nbins+j]-aggregate[j*nbins+i])
                        *(aggregate[i*nbins+j]-aggregate[j*nbins+i]);
                    selection_total2 += aggregate[i*nbins+j]
                        *aggregate[i*nbins+j];
                }
            group_selection_symmetry=sqrt(
                selection_anti2/max(selection_total2,1.0e-300));
            if(!isfinite(group_selection_symmetry))
                myerror(E,
                    "ALA cross-rank GenEO selection symmetry is not finite");
            local_selection_symmetry_max=max(local_selection_symmetry_max,
                                               group_selection_symmetry);
            for(i=0;i<nbins;i++) {
                active_map[i]=-1;
                if(isfinite(aggregate[i*nbins+i]) &&
                   aggregate[i*nbins+i]>0.0)
                    active_map[i]=nactive++;
            }
            if(strcmp(E->control.ala_geneo_basis_type,"spectral")==0) {
                for(i=0;i<nbins;i++)
                    if(active_map[i]>=0)
                        for(j=0;j<nbins;j++)
                            if(active_map[j]>=0) {
                                active1=active_map[i];
                                active2=active_map[j];
                                value=0.5*(aggregate[i*nbins+j]
                                         +aggregate[j*nbins+i]);
                                normalized[active1*nactive+active2]=value
                                    /sqrt(aggregate[i*nbins+i]
                                          *aggregate[j*nbins+j]);
                                if(!isfinite(
                                    normalized[active1*nactive+active2]))
                                    myerror(E,
                                        "ALA cross-rank GenEO operator is not finite");
                            }
                if(nactive>0)
                    iterations=ala_geneo_jacobi_eigensolve(
                        normalized,vectors,eigenvalues,nactive);
                if(nactive>0 &&
                   iterations>=max(50*nactive*nactive,1))
                    myerror(E,
                        "ALA cross-rank GenEO eigensolve did not converge");
                for(i=0;i<nactive;i++)
                    if(!isfinite(eigenvalues[i]))
                        myerror(E,
                            "ALA cross-rank GenEO eigenvalue is not finite");
                for(i=0;i<nactive;i++)
                    for(j=i+1;j<nactive;j++)
                        if(eigenvalues[j]<eigenvalues[i]) {
                            value=eigenvalues[i];
                            eigenvalues[i]=eigenvalues[j];
                            eigenvalues[j]=value;
                            for(k=0;k<nactive;k++) {
                                value=vectors[k*nactive+i];
                                vectors[k*nactive+i]
                                    =vectors[k*nactive+j];
                                vectors[k*nactive+j]=value;
                            }
                        }
            }
            if(strcmp(E->control.ala_geneo_basis_type,
                      "radial_partition")==0) {
                /* Processor groups are also split by the radial rank
                   coordinate.  Only the upper radial groups intersect the
                   configured shallow layer; deep groups must own no shallow
                   coarse modes rather than two identically zero modes. */
                desired_modes=(shallow_layers>0) ? rbins : 0;
                accepted=(desired_modes>0)
                    ? ala_build_radial_partition_shapes(
                        aggregate,nbins,rbins,selected_shapes[m],
                        selected_values[m])
                    : 0;
            }
            else {
                desired_modes=0;
                while(desired_modes<nactive &&
                      desired_modes<E->control.ala_geneo_max_modes_per_rank &&
                      eigenvalues[desired_modes]<=threshold)
                    desired_modes++;
                if(nactive>0)
                    desired_modes=max(desired_modes,
                        E->control.ala_geneo_min_modes_per_rank);
                desired_modes=min(desired_modes,nactive);
                desired_modes=min(desired_modes,
                    E->control.ala_geneo_max_modes_per_rank);

                /* The first mode is the normalized indicator of the complete
                   processor aggregate.  These disjoint indicators are an exact
                   partition of unity and guarantee that the spectral path
                   contains a genuinely cross-rank component. */
                accepted=0;
                if(desired_modes>0) {
                    denominator=0.0;
                    for(j=0;j<nbins;j++)
                        if(active_map[j]>=0) {
                            selected_shapes[m][j]=1.0;
                            denominator += aggregate[j*nbins+j];
                        }
                    if(!isfinite(denominator) || denominator<=1.0e-30)
                        myerror(E,
                            "ALA cross-rank GenEO partition mode has zero mass");
                    denominator=sqrt(denominator);
                    for(j=0;j<nbins;j++)
                        selected_shapes[m][j] /= denominator;
                    rayleigh_numerator=0.0;
                    rayleigh_denominator=0.0;
                    for(i=0;i<nbins;i++)
                        if(active_map[i]>=0) {
                            rayleigh_denominator += aggregate[i*nbins+i]
                                *selected_shapes[m][i]
                                *selected_shapes[m][i];
                            for(j=0;j<nbins;j++)
                                if(active_map[j]>=0)
                                    rayleigh_numerator +=
                                        selected_shapes[m][i]
                                        *aggregate[i*nbins+j]
                                        *selected_shapes[m][j];
                        }
                    selected_values[m][0]=rayleigh_numerator
                        /max(rayleigh_denominator,1.0e-300);
                    accepted=1;
                }

                /* Add low-energy spectral shapes after D-orthogonalizing them
                   against the partition mode and previously accepted shapes.
                   Near-dependent candidates are skipped instead of creating
                   an ill-conditioned global Galerkin matrix. */
                candidate=0;
                while(accepted<desired_modes && candidate<nactive) {
                    for(j=0;j<nbins;j++)
                        selected_shapes[m][accepted*nbins+j]
                            =(active_map[j]>=0)
                            ? vectors[active_map[j]*nactive+candidate]
                              /sqrt(aggregate[j*nbins+j])
                            : 0.0;
                    for(p=0;p<accepted;p++) {
                        coefficient=0.0;
                        denominator=0.0;
                        for(j=0;j<nbins;j++)
                            if(active_map[j]>=0) {
                                coefficient += aggregate[j*nbins+j]
                                    *selected_shapes[m][accepted*nbins+j]
                                    *selected_shapes[m][p*nbins+j];
                                denominator += aggregate[j*nbins+j]
                                    *selected_shapes[m][p*nbins+j]
                                    *selected_shapes[m][p*nbins+j];
                            }
                        coefficient /= max(denominator,1.0e-300);
                        for(j=0;j<nbins;j++)
                            selected_shapes[m][accepted*nbins+j]
                                -=coefficient
                                  *selected_shapes[m][p*nbins+j];
                    }
                    denominator=0.0;
                    for(j=0;j<nbins;j++)
                        if(active_map[j]>=0)
                            denominator += aggregate[j*nbins+j]
                                *selected_shapes[m][accepted*nbins+j]
                                *selected_shapes[m][accepted*nbins+j];
                    candidate++;
                    if(!isfinite(denominator) || denominator<=1.0e-20) {
                        for(j=0;j<nbins;j++)
                            selected_shapes[m][accepted*nbins+j]=0.0;
                        continue;
                    }
                    denominator=sqrt(denominator);
                    for(j=0;j<nbins;j++)
                        selected_shapes[m][accepted*nbins+j] /= denominator;
                    rayleigh_numerator=0.0;
                    rayleigh_denominator=0.0;
                    for(i=0;i<nbins;i++)
                        if(active_map[i]>=0) {
                            rayleigh_denominator += aggregate[i*nbins+i]
                                *selected_shapes[m][accepted*nbins+i]
                                *selected_shapes[m][accepted*nbins+i];
                            for(j=0;j<nbins;j++)
                                if(active_map[j]>=0)
                                    rayleigh_numerator +=
                                      selected_shapes[m][accepted*nbins+i]
                                      *aggregate[i*nbins+j]
                                      *selected_shapes[m][accepted*nbins+j];
                        }
                    selected_values[m][accepted]=rayleigh_numerator
                        /max(rayleigh_denominator,1.0e-300);
                    accepted++;
                }
            }
            if(accepted<desired_modes) {
                group_modes[m]=(accepted<0) ? accepted : -1;
                local_modes=0;
            }
            else {
                local_modes=accepted;
                group_modes[m]=local_modes;
            }
            if(group_modes[m]>0) {
                for(i=0;i<local_modes;i++) {
                    local_eigen_min=min(local_eigen_min,
                                        selected_values[m][i]);
                    local_eigen_max=max(local_eigen_max,
                                        selected_values[m][i]);
                }
                owned_modes += local_modes;
                local_active++;
                local_mode_min=min(local_mode_min,local_modes);
                local_mode_max=max(local_mode_max,local_modes);
                local_group_min=min(local_group_min,group_size);
                local_group_max=max(local_group_max,group_size);
            }
            if(local_modes>1) {
                local_second_min=min(local_second_min,
                                     selected_values[m][1]);
                local_second_max=max(local_second_max,
                                     selected_values[m][1]);
            }
        }
        MPI_Bcast(&group_modes[m],1,MPI_INT,0,group_comm[m]);
        if(group_modes[m]<0) {
            if(group_rank==0) {
                fprintf(stderr,
                    "ALA cross-rank basis construction failed cap=%d "
                    "failure_code=%d basis_type=%s\n",
                    E->sphere.capid[m],group_modes[m],
                    E->control.ala_geneo_basis_type);
                fflush(stderr);
            }
            myerror(E,strcmp(E->control.ala_geneo_basis_type,
                            "radial_partition")==0
                    ? "ALA Stage 6d radial partition basis is incomplete"
                    : "ALA cross-rank GenEO basis is incomplete");
        }
        MPI_Bcast(selected_values[m],
                  E->control.ala_geneo_max_modes_per_rank,
                  MPI_DOUBLE,0,group_comm[m]);
        MPI_Bcast(selected_shapes[m],
                  E->control.ala_geneo_max_modes_per_rank*nbins,
                  MPI_DOUBLE,0,group_comm[m]);
        MPI_Allreduce(&iterations,&i,1,MPI_INT,MPI_MAX,group_comm[m]);
        if(group_rank==0)
            global_iterations=max(global_iterations,i);
        free((void *)local_aggregate);
        free((void *)aggregate);
        free((void *)active_map);
        free((void *)normalized);
        free((void *)vectors);
        free((void *)eigenvalues);
    }

    MPI_Allreduce(&local_interface_entries,&global_interface_entries,1,
                  MPI_INT,MPI_SUM,E->parallel.world);
    MPI_Allreduce(&local_selection_symmetry_max,
                  &global_selection_symmetry_max,1,MPI_DOUBLE,MPI_MAX,
                  E->parallel.world);
    if(global_interface_entries<=0)
        myerror(E,
            "ALA cross-rank GenEO selection found no MPI interface coupling");
    if(!isfinite(global_selection_symmetry_max) ||
       global_selection_symmetry_max>1.0e-8)
        myerror(E,
            "ALA cross-rank GenEO selection operator is not symmetric");

    mode_counts=(int *)calloc(E->parallel.nproc,sizeof(int));
    mode_offsets=(int *)calloc(E->parallel.nproc,sizeof(int));
    if(mode_counts==NULL || mode_offsets==NULL)
        myerror(E,"Unable to allocate cross-rank GenEO mode ownership");
    MPI_Allgather(&owned_modes,1,MPI_INT,mode_counts,1,MPI_INT,
                  E->parallel.world);
    n=0;
    for(i=0;i<E->parallel.nproc;i++) {
        mode_offsets[i]=n;
        n += mode_counts[i];
    }
    if(n<1)
        myerror(E,"ALA cross-rank GenEO selected no global coarse modes");
    if(n>E->control.ala_geneo_max_global_modes)
        myerror(E,
            "ALA cross-rank GenEO global mode count exceeds configured maximum");
    cache->geneo_basis_count=n;
    cache->geneo_local_modes=owned_modes;
    cache->geneo_local_offset=mode_offsets[E->parallel.me];
    owned_cursor=mode_offsets[E->parallel.me];
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        MPI_Comm_rank(group_comm[m],&group_rank);
        if(group_rank==0) {
            group_offsets[m]=owned_cursor;
            owned_cursor += group_modes[m];
        }
        MPI_Bcast(&group_offsets[m],1,MPI_INT,0,group_comm[m]);
    }
    cache->geneo_eigenvalues=(double *)calloc(n,sizeof(double));
    local_eigenvalues=(double *)calloc(n,sizeof(double));
    if(cache->geneo_eigenvalues==NULL || local_eigenvalues==NULL)
        myerror(E,
            "Unable to allocate cross-rank GenEO eigenvalue diagnostics");
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        MPI_Comm_rank(group_comm[m],&group_rank);
        if(group_rank==0)
            for(i=0;i<group_modes[m];i++)
                local_eigenvalues[group_offsets[m]+i]
                    =selected_values[m][i];
    }
    MPI_Allreduce(local_eigenvalues,cache->geneo_eigenvalues,n,MPI_DOUBLE,
                  MPI_SUM,E->parallel.world);

    clev=lev-E->control.ala_two_level_offset;
    factor=1 << E->control.ala_two_level_offset;
    celx=E->lmesh.ELX[clev];
    cely=E->lmesh.ELY[clev];
    celz=E->lmesh.ELZ[clev];
    cnpno=E->lmesh.NPNO[clev];
    stride=cnpno+1;
    npno=E->lmesh.NPNO[lev];
    neq=E->lmesh.NEQ[lev];
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        cache->geneo_basis[m]=(double *)calloc(n*stride,sizeof(double));
        if(cache->geneo_basis[m]==NULL)
            myerror(E,"Unable to allocate cross-rank GenEO basis");
        rank_offset_x=group_rank_offset_x[m];
        rank_offset_y=group_rank_offset_y[m];
        group_hbins_x=group_hx[m];
        nbins=group_nbins[m];
        for(i=0;i<group_modes[m];i++) {
            col=group_offsets[m]+i;
            norm=0.0;
            for(cy=1;cy<=cely;cy++)
                for(cx=1;cx<=celx;cx++)
                    for(cz=1;cz<=celz;cz++) {
                        ce=cz+(cx-1)*celz+(cy-1)*celz*celx;
                        depth=(1.0-E->ECO[clev][m][ce].centre[3])
                            *E->data.radius_km;
                        if(depth>E->control.ala_shallow_patch_depth_km)
                            continue;
                        ex1=min((cx-1)*factor+1,elx);
                        ey1=min((cy-1)*factor+1,ely);
                        ez1=min(cz*factor,elz);
                        if(ez1<shallow_min_ez)
                            continue;
                        bin1=(ez1-shallow_min_ez)*rbins/shallow_layers
                            +rbins*((rank_offset_x*hbins
                            +min((ex1-1)*hbins/elx,hbins-1))
                            +group_hbins_x*(rank_offset_y*hbins
                            +min((ey1-1)*hbins/ely,hbins-1)));
                        value=selected_shapes[m][i*nbins+bin1];
                        cache->geneo_basis[m][col*stride+ce]=value;
                        norm += value*value;
                    }
            MPI_Allreduce(&norm,&group_norm,1,MPI_DOUBLE,MPI_SUM,
                          group_comm[m]);
            local_has_support=(norm>1.0e-30) ? 1 : 0;
            MPI_Allreduce(&local_has_support,&supported_members,1,MPI_INT,
                          MPI_SUM,group_comm[m]);
            if(strcmp(E->control.ala_geneo_basis_type,
                      "radial_partition")==0 &&
               supported_members!=group_size)
                myerror(E,
                    "ALA Stage 6d radial mode misses a member rank");
            if(strcmp(E->control.ala_geneo_basis_type,"spectral")==0 &&
               i==0 && supported_members!=group_size)
                myerror(E,
                    "ALA cross-rank GenEO partition mode misses a member rank");
            if(!isfinite(group_norm) || group_norm<=1.0e-30)
                myerror(E,
                    "ALA cross-rank GenEO selected basis has zero norm");
            group_norm=sqrt(group_norm);
            for(ce=1;ce<=cnpno;ce++)
                cache->geneo_basis[m][col*stride+ce] /= group_norm;
        }
    }

    cache->geneo_matrix=(double *)calloc(n*n,sizeof(double));
    cache->geneo_chol=(double *)calloc(n*n,sizeof(double));
    local_column=(double *)calloc(n,sizeof(double));
    global_column=(double *)calloc(n,sizeof(double));
    if(cache->geneo_matrix==NULL || cache->geneo_chol==NULL ||
       local_column==NULL || global_column==NULL)
        myerror(E,"Unable to allocate cross-rank GenEO global matrix");
    matrix=cache->geneo_matrix;
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        coarse_p[m]=(double *)calloc(cnpno+1,sizeof(double));
        coarse_Ap[m]=(double *)calloc(cnpno+1,sizeof(double));
        fine_p[m]=(double *)calloc(npno+1,sizeof(double));
        fine_Ap[m]=(double *)calloc(npno+1,sizeof(double));
        fine_velocity[m]=(double *)calloc(neq+1,sizeof(double));
        fine_velocity_rhs[m]=(double *)calloc(neq+1,sizeof(double));
        fine_velocity_Ax[m]=(double *)calloc(neq+1,sizeof(double));
        fine_velocity_direction[m]=(double *)calloc(neq+1,sizeof(double));
        if(coarse_p[m]==NULL || coarse_Ap[m]==NULL || fine_p[m]==NULL ||
           fine_Ap[m]==NULL || fine_velocity[m]==NULL ||
           fine_velocity_rhs[m]==NULL || fine_velocity_Ax[m]==NULL ||
           fine_velocity_direction[m]==NULL)
            myerror(E,
                "Unable to allocate cross-rank GenEO Galerkin workspace");
    }
    for(col=0;col<n;col++) {
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(ce=1;ce<=cnpno;ce++)
                coarse_p[m][ce]=cache->geneo_basis[m][col*stride+ce];
        apply_ala_galerkin_fixed_schur(
            E,coarse_p,coarse_Ap,fine_p,fine_Ap,fine_velocity,
            fine_velocity_rhs,fine_velocity_Ax,fine_velocity_direction,
            cache,clev,lev,factor,NULL);
        for(row=0;row<n;row++)
            local_column[row]=0.0;
        for(row=0;row<n;row++)
            for(m=1;m<=E->sphere.caps_per_proc;m++)
                for(ce=1;ce<=cnpno;ce++)
                    local_column[row] +=
                        cache->geneo_basis[m][row*stride+ce]
                        *coarse_Ap[m][ce];
        MPI_Allreduce(local_column,global_column,n,MPI_DOUBLE,MPI_SUM,
                      E->parallel.world);
        for(row=0;row<n;row++)
            matrix[row*n+col]=global_column[row];
        if(col==0 || (col+1)%25==0 || col+1==n) {
            setup_seconds=CPU_time0()-setup_start;
            MPI_Allreduce(&setup_seconds,&global_setup_seconds,1,MPI_DOUBLE,
                          MPI_MAX,E->parallel.world);
            if(E->parallel.me==0) {
                fprintf(E->fp,
                    "ALA_CROSS_RANK_GENEO_SETUP_PROGRESS column=%d/%d "
                    "elapsed_seconds_max=%e\n",
                    col+1,n,global_setup_seconds);
                fprintf(stderr,
                    "ALA_CROSS_RANK_GENEO_SETUP_PROGRESS column=%d/%d "
                    "elapsed_seconds_max=%e\n",
                    col+1,n,global_setup_seconds);
                fflush(E->fp);
                fflush(stderr);
            }
        }
    }
    cache->geneo_schur_applications=n;
    anti2=0.0;
    total2=0.0;
    diagonal_min=1.0e300;
    diagonal_max=0.0;
    for(i=0;i<n;i++) {
        diagonal_min=min(diagonal_min,matrix[i*n+i]);
        diagonal_max=max(diagonal_max,matrix[i*n+i]);
        for(j=0;j<n;j++) {
            anti2 += (matrix[i*n+j]-matrix[j*n+i])
                *(matrix[i*n+j]-matrix[j*n+i]);
            total2 += matrix[i*n+j]*matrix[i*n+j];
            cache->geneo_chol[i*n+j]
                =0.5*(matrix[i*n+j]+matrix[j*n+i]);
        }
    }
    symmetry=sqrt(anti2/max(total2,1.0e-300));
    if(!isfinite(symmetry) || symmetry>1.0e-8)
        myerror(E,"ALA cross-rank GenEO Galerkin matrix is not symmetric");
    if(!isfinite(diagonal_min) || diagonal_min<=0.0)
        myerror(E,"ALA cross-rank GenEO Galerkin diagonal is not positive");
    shift=E->control.ala_geneo_regularization*diagonal_max;
    for(i=0;i<n;i++)
        cache->geneo_chol[i*n+i] += shift;
    min_pivot=1.0e300;
    for(i=0;i<n;i++) {
        for(j=0;j<=i;j++) {
            sum=cache->geneo_chol[i*n+j];
            for(k=0;k<j;k++)
                sum -= cache->geneo_chol[i*n+k]
                    *cache->geneo_chol[j*n+k];
            if(i==j) {
                pivot=sum;
                if(!isfinite(pivot) ||
                   pivot<=1.0e-14*max(diagonal_max,1.0e-300))
                    myerror(E,"ALA cross-rank GenEO Cholesky failed");
                cache->geneo_chol[i*n+j]=sqrt(pivot);
                min_pivot=min(min_pivot,pivot);
            }
            else
                cache->geneo_chol[i*n+j]
                    =sum/cache->geneo_chol[j*n+j];
        }
        for(j=i+1;j<n;j++)
            cache->geneo_chol[i*n+j]=0.0;
    }

    MPI_Allreduce(&local_active,&global_active,1,MPI_INT,MPI_SUM,
                  E->parallel.world);
    MPI_Allreduce(&local_mode_min,&global_mode_min,1,MPI_INT,MPI_MIN,
                  E->parallel.world);
    MPI_Allreduce(&local_mode_max,&global_mode_max,1,MPI_INT,MPI_MAX,
                  E->parallel.world);
    MPI_Allreduce(&local_group_min,&global_group_min,1,MPI_INT,MPI_MIN,
                  E->parallel.world);
    MPI_Allreduce(&local_group_max,&global_group_max,1,MPI_INT,MPI_MAX,
                  E->parallel.world);
    MPI_Allreduce(&local_eigen_min,&global_eigen_min,1,MPI_DOUBLE,MPI_MIN,
                  E->parallel.world);
    MPI_Allreduce(&local_eigen_max,&global_eigen_max,1,MPI_DOUBLE,MPI_MAX,
                  E->parallel.world);
    MPI_Allreduce(&local_second_min,&global_second_min,1,MPI_DOUBLE,MPI_MIN,
                  E->parallel.world);
    MPI_Allreduce(&local_second_max,&global_second_max,1,MPI_DOUBLE,MPI_MAX,
                  E->parallel.world);
    MPI_Allreduce(&global_iterations,&iterations,1,MPI_INT,MPI_MAX,
                  E->parallel.world);
    setup_seconds=CPU_time0()-setup_start;
    MPI_Allreduce(&setup_seconds,&global_setup_seconds,1,MPI_DOUBLE,MPI_MAX,
                  E->parallel.world);
    if(E->parallel.me==0) {
        fprintf(E->fp,
            "ALA cross-rank GenEO coarse space global_modes=%d "
            "basis_type=%s "
            "active_aggregates=%d modes_per_aggregate_range=(%d,%d) "
            "rank_group=%dx%d group_member_range=(%d,%d) "
            "aggregate_bins_max=%dx%dx%d partition_of_unity=disjoint_exact "
            "basis_construction=%s "
            "interface_entries=%d selection_symmetry_defect=%e "
            "threshold=%e selected_rayleigh_range=(%e,%e) "
            "second_mode_rayleigh_range=(%e,%e) "
            "local_jacobi_iterations=%d Galerkin_applications=%d "
            "symmetry_defect=%e diagonal_range=(%e,%e) "
            "regularization=%e shift=%e min_pivot=%e "
            "setup_seconds_max=%e status=pass\n",
            n,E->control.ala_geneo_basis_type,
            global_active,global_mode_min,global_mode_max,
            rank_group_x,rank_group_y,global_group_min,global_group_max,
            hbins*rank_group_x,hbins*rank_group_y,rbins,
            strcmp(E->control.ala_geneo_basis_type,"radial_partition")==0
                ? "geometric_radial_partition"
                : "interface_coupled_jacobi",
            global_interface_entries,global_selection_symmetry_max,threshold,
            global_eigen_min,global_eigen_max,
            global_second_min,global_second_max,iterations,n,
            symmetry,diagonal_min,diagonal_max,
            E->control.ala_geneo_regularization,shift,min_pivot,
            global_setup_seconds);
        fprintf(stderr,
            "ALA cross-rank GenEO coarse space global_modes=%d "
            "basis_type=%s "
            "active_aggregates=%d rank_group=%dx%d "
            "aggregate_bins_max=%dx%dx%d "
            "partition_of_unity=disjoint_exact "
            "basis_construction=%s "
            "interface_entries=%d selection_symmetry_defect=%e "
            "symmetry_defect=%e setup_seconds_max=%e status=pass\n",
            n,E->control.ala_geneo_basis_type,
            global_active,rank_group_x,rank_group_y,
            hbins*rank_group_x,hbins*rank_group_y,rbins,
            strcmp(E->control.ala_geneo_basis_type,"radial_partition")==0
                ? "geometric_radial_partition"
                : "interface_coupled_jacobi",
            global_interface_entries,global_selection_symmetry_max,
            symmetry,global_setup_seconds);
        fflush(E->fp);
        fflush(stderr);
    }
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        free((void *)coarse_p[m]);
        free((void *)coarse_Ap[m]);
        free((void *)fine_p[m]);
        free((void *)fine_Ap[m]);
        free((void *)fine_velocity[m]);
        free((void *)fine_velocity_rhs[m]);
        free((void *)fine_velocity_Ax[m]);
        free((void *)fine_velocity_direction[m]);
        free((void *)selected_values[m]);
        free((void *)selected_shapes[m]);
        MPI_Comm_free(&group_comm[m]);
    }
    free((void *)mode_counts);
    free((void *)mode_offsets);
    free((void *)local_eigenvalues);
    free((void *)local_column);
    free((void *)global_column);
}


/* Calibrate the Jacobi-scaled velocity operator used by the fixed polynomial
   inverse.  The similarity transform sqrt(B) K sqrt(B) is symmetric; a fixed
   Chebyshev polynomial of it therefore yields an SPD approximation to K^-1. */
static void calibrate_ala_velocity_spectrum(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache, int lev)
{
    int m,j,k,node,d,eqn,neq,nno;
    double theta,phi,radius,norm,rayleigh,upper,scale;
    double *q[NCS],*scaled_q[NCS],*Kq[NCS],*transformed_q[NCS];
    double global_vdot();
    void assemble_del2_u();
    void strip_bcs_from_residual();

    neq=E->lmesh.NEQ[lev];
    nno=E->lmesh.NNO[lev];
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        q[m]=(double *)calloc(neq+1,sizeof(double));
        scaled_q[m]=(double *)calloc(neq+1,sizeof(double));
        Kq[m]=(double *)calloc(neq+1,sizeof(double));
        transformed_q[m]=(double *)calloc(neq+1,sizeof(double));
        if(q[m]==NULL || scaled_q[m]==NULL || Kq[m]==NULL ||
           transformed_q[m]==NULL)
            myerror(E,"Unable to allocate ALA velocity spectrum workspace");
        for(node=1;node<=nno;node++) {
            theta=E->SX[lev][m][1][node];
            phi=E->SX[lev][m][2][node];
            radius=E->SX[lev][m][3][node];
            for(d=1;d<=E->mesh.nsd;d++) {
                eqn=E->ID[lev][m][node].doff[d];
                q[m][eqn]=sin(11.324718*theta+7.193147*phi
                              +5.731921*radius+0.417*d
                              +0.113*E->sphere.capid[m]);
            }
        }
    }
    strip_bcs_from_residual(E,q,lev);
    norm=sqrt(max(global_vdot(E,q,q,lev),1.0e-300));
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(j=0;j<neq;j++)
            q[m][j] /= norm;

    rayleigh=0.0;
    for(k=0;k<ALA_VELOCITY_POWER_ITERATIONS;k++) {
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(j=0;j<neq;j++) {
                scale=sqrt(E->ALA_velocity_BI[lev][m][j]);
                scaled_q[m][j]=scale*q[m][j];
            }
        assemble_del2_u(E,scaled_q,Kq,lev,1);
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(j=0;j<neq;j++)
                transformed_q[m][j]
                    =sqrt(E->ALA_velocity_BI[lev][m][j])*Kq[m][j];
        rayleigh=global_vdot(E,q,transformed_q,lev);
        norm=sqrt(max(global_vdot(E,transformed_q,transformed_q,lev),
                      1.0e-300));
        if(!isfinite(rayleigh) || rayleigh<=0.0 || !isfinite(norm) ||
           norm<=1.0e-150)
            myerror(E,"ALA velocity spectral calibration failed");
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(j=0;j<neq;j++)
                q[m][j]=transformed_q[m][j]/norm;
    }

    upper=max(2.0*rayleigh,
              2.0*E->control.ala_two_level_velocity_eigenvalue_min);
    if(upper>E->control.ala_two_level_velocity_eigenvalue_max)
        myerror(E,"ALA velocity spectral safety bound exceeds configured maximum");
    cache->velocity_eigenvalue_max=upper;
    cache->velocity_eigenvalue_min=
        E->control.ala_two_level_velocity_eigenvalue_min;
    if(cache->velocity_eigenvalue_max<=cache->velocity_eigenvalue_min)
        myerror(E,"ALA calibrated velocity interval is empty");
    if(E->parallel.me==0) {
        fprintf(E->fp,"ALA velocity spectrum power_iterations=%d "
                "lambda_estimate=%e safety_factor=2.0 "
                "chebyshev_interval=(%e,%e) configured_upper=%e\n",
                ALA_VELOCITY_POWER_ITERATIONS,rayleigh,
                cache->velocity_eigenvalue_min,
                cache->velocity_eigenvalue_max,
                E->control.ala_two_level_velocity_eigenvalue_max);
        fprintf(stderr,"ALA velocity spectrum power_iterations=%d "
                "lambda_estimate=%e safety_factor=2.0 "
                "chebyshev_interval=(%e,%e) configured_upper=%e\n",
                ALA_VELOCITY_POWER_ITERATIONS,rayleigh,
                cache->velocity_eigenvalue_min,
                cache->velocity_eigenvalue_max,
                E->control.ala_two_level_velocity_eigenvalue_max);
        fflush(E->fp);
        fflush(stderr);
    }
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        free((void *)q[m]);
        free((void *)scaled_q[m]);
        free((void *)Kq[m]);
        free((void *)transformed_q[m]);
    }
}


/* Estimate the largest eigenvalue of M^(1/2) Sc M^(1/2), where M is the
   aggregate Jacobi inverse.  This uses the exact same fixed Galerkin map as
   the coarse polynomial, so the resulting interval remains fixed throughout
   the outer PCG solve. */
static void calibrate_ala_two_level_spectrum(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache, int lev)
{
    int m,ce,k,clev,factor,cnpno,npno,neq;
    double local[2],global[2],norm,rayleigh,upper,scale;
    double *q[NCS],*coarse_p[NCS],*coarse_Ap[NCS];
    double *fine_p[NCS],*fine_Ap[NCS],*fine_velocity[NCS];
    double *fine_velocity_rhs[NCS],*fine_velocity_Ax[NCS];
    double *fine_velocity_direction[NCS];

    clev=lev-E->control.ala_two_level_offset;
    factor=1 << E->control.ala_two_level_offset;
    cnpno=E->lmesh.NPNO[clev];
    npno=E->lmesh.NPNO[lev];
    neq=E->lmesh.NEQ[lev];
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        q[m]=(double *)calloc(cnpno+1,sizeof(double));
        coarse_p[m]=(double *)calloc(cnpno+1,sizeof(double));
        coarse_Ap[m]=(double *)calloc(cnpno+1,sizeof(double));
        fine_p[m]=(double *)calloc(npno+1,sizeof(double));
        fine_Ap[m]=(double *)calloc(npno+1,sizeof(double));
        fine_velocity[m]=(double *)calloc(neq+1,sizeof(double));
        fine_velocity_rhs[m]=(double *)calloc(neq+1,sizeof(double));
        fine_velocity_Ax[m]=(double *)calloc(neq+1,sizeof(double));
        fine_velocity_direction[m]=(double *)calloc(neq+1,sizeof(double));
        if(q[m]==NULL || coarse_p[m]==NULL || coarse_Ap[m]==NULL ||
           fine_p[m]==NULL || fine_Ap[m]==NULL || fine_velocity[m]==NULL ||
           fine_velocity_rhs[m]==NULL || fine_velocity_Ax[m]==NULL ||
           fine_velocity_direction[m]==NULL)
            myerror(E,"Unable to allocate ALA Galerkin spectrum workspace");
        for(ce=1;ce<=cnpno;ce++)
            q[m][ce]=sin(0.7548776662466927*(double)ce
                         +0.5698402909980532*(double)(E->parallel.me+1));
    }

    local[0]=0.0;
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(ce=1;ce<=cnpno;ce++)
            local[0] += q[m][ce]*q[m][ce];
    MPI_Allreduce(local,global,1,MPI_DOUBLE,MPI_SUM,E->parallel.world);
    norm=sqrt(max(global[0],1.0e-300));
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(ce=1;ce<=cnpno;ce++)
            q[m][ce] /= norm;

    rayleigh=0.0;
    for(k=0;k<ALA_COARSE_POWER_ITERATIONS;k++) {
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(ce=1;ce<=cnpno;ce++) {
                scale=sqrt(cache->coarse_bpi[m][ce]);
                coarse_p[m][ce]=scale*q[m][ce];
            }
        apply_ala_galerkin_fixed_schur(
            E,coarse_p,coarse_Ap,fine_p,fine_Ap,fine_velocity,
            fine_velocity_rhs,fine_velocity_Ax,fine_velocity_direction,
            cache,clev,lev,factor,NULL);
        local[0]=0.0;
        local[1]=0.0;
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(ce=1;ce<=cnpno;ce++) {
                coarse_Ap[m][ce] *= sqrt(cache->coarse_bpi[m][ce]);
                local[0] += q[m][ce]*coarse_Ap[m][ce];
                local[1] += coarse_Ap[m][ce]*coarse_Ap[m][ce];
            }
        MPI_Allreduce(local,global,2,MPI_DOUBLE,MPI_SUM,
                      E->parallel.world);
        rayleigh=global[0];
        norm=sqrt(max(global[1],1.0e-300));
        if(!isfinite(rayleigh) || !isfinite(norm) || norm<=1.0e-150)
            myerror(E,"ALA Galerkin spectral calibration failed");
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(ce=1;ce<=cnpno;ce++)
                q[m][ce]=coarse_Ap[m][ce]/norm;
    }

    upper=max(2.0*rayleigh,
              2.0*E->control.ala_two_level_coarse_eigenvalue_min);
    if(upper>E->control.ala_two_level_coarse_eigenvalue_max)
        myerror(E,"ALA Galerkin spectral safety bound exceeds configured maximum");
    cache->coarse_eigenvalue_max=upper;
    if(cache->coarse_eigenvalue_max<=
       E->control.ala_two_level_coarse_eigenvalue_min)
        myerror(E,"ALA Galerkin calibrated interval is empty");
    cache->coarse_eigenvalue_min=
        E->control.ala_two_level_coarse_eigenvalue_min;
    if(E->parallel.me==0) {
        fprintf(E->fp,"ALA Galerkin spectrum power_iterations=%d "
                "lambda_estimate=%e safety_factor=2.0 "
                "chebyshev_interval=(%e,%e) configured_upper=%e\n",
                ALA_COARSE_POWER_ITERATIONS,rayleigh,
                cache->coarse_eigenvalue_min,cache->coarse_eigenvalue_max,
                E->control.ala_two_level_coarse_eigenvalue_max);
        fprintf(stderr,"ALA Galerkin spectrum power_iterations=%d "
                "lambda_estimate=%e safety_factor=2.0 "
                "chebyshev_interval=(%e,%e) configured_upper=%e\n",
                ALA_COARSE_POWER_ITERATIONS,rayleigh,
                cache->coarse_eigenvalue_min,cache->coarse_eigenvalue_max,
                E->control.ala_two_level_coarse_eigenvalue_max);
        fflush(E->fp);
        fflush(stderr);
    }
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        free((void *)q[m]);
        free((void *)coarse_p[m]);
        free((void *)coarse_Ap[m]);
        free((void *)fine_p[m]);
        free((void *)fine_Ap[m]);
        free((void *)fine_velocity[m]);
        free((void *)fine_velocity_rhs[m]);
        free((void *)fine_velocity_Ax[m]);
        free((void *)fine_velocity_direction[m]);
    }
}


static void apply_ala_pressure_preconditioner(struct All_variables *E,
                                              double **r, double **z,
                                              double **work, int lev,
                                              int iteration,
    struct ala_pressure_preconditioner_cache *cache)
{
    int m,j,col,k,e,ex,ey,ez,elz,ncolumns,npno,b,n,i;
    int face,ref,ghost_index;
    int clev,factor,celx,celz,cnpno,ce,cx,cy,cz,neq;
    double damping,theta,delta,sigma,rho_cheb,rho_new;
    double velocity_residual_reduction;
    double rhs[ALA_PATCH_MAX_ELEMENTS],solution[ALA_PATCH_MAX_ELEMENTS];
    double constant_mode[ALA_PATCH_MAX_ELEMENTS];
    double sum,weight,numerator,denominator,*L;
    double local_patch_energy[2],global_patch_energy[2];
    double local_energy[4],global_energy[4];
    double global_local_rhs[ALA_GLOBAL_BASIS_COUNT];
    double global_rhs[ALA_GLOBAL_BASIS_COUNT];
    double global_solution[ALA_GLOBAL_BASIS_COUNT];
    double global_residual2,global_rhs2,global_scale;
    double *coarse_rhs[NCS],*coarse_x[NCS],*coarse_residual[NCS];
    double *coarse_Ax[NCS],*coarse_direction[NCS];
    double *fine_p[NCS],*fine_Ap[NCS],*fine_velocity[NCS];
    double *fine_velocity_rhs[NCS],*fine_velocity_Ax[NCS];
    double *fine_velocity_direction[NCS];
    double *ghost_r[NCS][ALA_PATCH_MPI_FACES];
    double *ghost_work[NCS][ALA_PATCH_MPI_FACES];
    double **effective_r,depth_km,depth_weight;

    npno=E->lmesh.NPNO[lev];
    effective_r=r;
    if(fabs(E->control.ala_pressure_shallow_action_scale-1.0)>1.0e-12) {
        elz=E->lmesh.ELZ[lev];
        depth_weight=sqrt(E->control.ala_pressure_shallow_action_scale);
        for(m=1;m<=E->sphere.caps_per_proc;m++) {
            if(cache->pressure_depth_weighted_r[m]==NULL) {
                cache->pressure_depth_weighted_r[m]=(double *)calloc(
                    npno+1,sizeof(double));
                if(cache->pressure_depth_weighted_r[m]==NULL)
                    myerror(E,"Unable to allocate ALA depth-weighted pressure residual");
            }
            cache->pressure_depth_weighted_r[m][0]=0.0;
            for(e=1;e<=npno;e++) {
                ez=(e-1)%elz+1;
                depth_km=(1.0-0.5*(E->sx[m][3][ez]
                                  +E->sx[m][3][ez+1]))*E->data.radius_km;
                cache->pressure_depth_weighted_r[m][e]=r[m][e]
                    *(depth_km<E->control.ala_pressure_shallow_depth_km
                      ? depth_weight : 1.0);
            }
        }
        effective_r=cache->pressure_depth_weighted_r;
        if(iteration==0 && E->parallel.me==0) {
            fprintf(E->fp,"ALA PRESSURE DEPTH SIMILARITY depth_km=%e "
                    "shallow_action_scale=%e shallow_weight=%e "
                    "operator=W*P^-1*W\n",
                    E->control.ala_pressure_shallow_depth_km,
                    E->control.ala_pressure_shallow_action_scale,
                    depth_weight);
            fflush(E->fp);
        }
    }
    if(E->control.ala_radial_line_preconditioner) {
        elz=E->lmesh.ELZ[lev];
        ncolumns=E->lmesh.ELX[lev]*E->lmesh.ELY[lev];
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(col=0;col<ncolumns;col++) {
                if(!E->ALA_BPI_line_valid[lev][m][col+1]) {
                    for(k=0;k<elz;k++) {
                        e=col*elz+k+1;
                        z[m][e]=E->BPI[lev][m][e]*effective_r[m][e];
                    }
                    continue;
                }

                e=col*elz+1;
                work[m][e]=effective_r[m][e];
                for(k=1;k<elz;k++) {
                    e=col*elz+k+1;
                    work[m][e]=effective_r[m][e]-
                        E->ALA_BPI_line_lower[lev][m][e]*work[m][e-1];
                }
                for(k=0;k<elz;k++) {
                    e=col*elz+k+1;
                    z[m][e]=work[m][e]/E->ALA_BPI_line_diag[lev][m][e];
                }
                for(k=elz-2;k>=0;k--) {
                    e=col*elz+k+1;
                    z[m][e]-=E->ALA_BPI_line_lower[lev][m][e+1]*z[m][e+1];
                }
            }
    }
    else {
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(j=1;j<=npno;j++)
                z[m][j]=E->BPI[lev][m][j]*effective_r[m][j];
    }
    /* The diagonal/radial-line Schur map has its own scale, independent of
     * the velocity solve and optional additive pressure corrections.  This
     * lets the current-rheology block-triangular path correct an over-strong
     * pressure action without altering K, G, or the physical residual. */
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(j=1;j<=npno;j++)
            z[m][j] *= E->control.ala_pressure_bpi_weight;

    if(E->control.ala_shallow_patch_preconditioner) {
        weight=E->control.ala_shallow_patch_weight;
        local_patch_energy[0]=0.0;
        local_patch_energy[1]=0.0;
        for(m=1;m<=E->sphere.caps_per_proc;m++) {
            for(face=0;face<ALA_PATCH_MPI_FACES;face++) {
                ghost_r[m][face]=(double *)calloc(
                    max(cache->halo_recv_count[m][face],1),sizeof(double));
                ghost_work[m][face]=(double *)calloc(
                    max(cache->halo_recv_count[m][face],1),sizeof(double));
                if(ghost_r[m][face]==NULL || ghost_work[m][face]==NULL)
                    myerror(E,"Unable to allocate ALA MPI-overlap work");
            }
            exchange_ala_shallow_halo_values(E,cache,lev,m,effective_r,
                                              ghost_r[m],0);
            for(e=1;e<=npno;e++) {
                work[m][e]=0.0;
                if(cache->multiplicity[m][e]>0)
                    local_patch_energy[0] += effective_r[m][e]
                        *E->BPI[lev][m][e]*effective_r[m][e];
            }
        }
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(b=0;b<cache->blocks[m];b++) {
                n=cache->size[m][b];
                if(n==0)
                    continue;
                L=cache->chol[m]+(size_t)b*cache->patch_capacity
                    *cache->patch_capacity;
                for(i=0;i<n;i++) {
                    e=cache->elements[m][b*cache->patch_capacity+i];
                    constant_mode[i]=1.0
                        /sqrt((double)cache->multiplicity[m][e]);
                    rhs[i]=effective_r[m][e]*constant_mode[i];
                }
                if(strcmp(E->control.ala_shallow_patch_velocity_solver,
                          "element_vanka")==0) {
                    numerator=denominator=0.0;
                    for(i=0;i<n;i++) {
                        numerator += constant_mode[i]*rhs[i];
                        denominator += constant_mode[i]*constant_mode[i];
                    }
                    sum=numerator/max(denominator,1.0e-300);
                    for(i=0;i<n;i++)
                        rhs[i] -= sum*constant_mode[i];
                }
                for(i=0;i<n;i++) {
                    sum=rhs[i];
                    for(j=0;j<i;j++)
                        sum -= L[i*cache->patch_capacity+j]*solution[j];
                    solution[i]=sum/L[i*cache->patch_capacity+i];
                }
                for(i=n-1;i>=0;i--) {
                    sum=solution[i];
                    for(j=i+1;j<n;j++)
                        sum -= L[j*cache->patch_capacity+i]*solution[j];
                    solution[i]=sum/L[i*cache->patch_capacity+i];
                }
                if(strcmp(E->control.ala_shallow_patch_velocity_solver,
                          "element_vanka")==0) {
                    numerator=denominator=0.0;
                    for(i=0;i<n;i++) {
                        numerator += constant_mode[i]*solution[i];
                        denominator += constant_mode[i]*constant_mode[i];
                    }
                    sum=numerator/max(denominator,1.0e-300);
                    for(i=0;i<n;i++)
                        solution[i] -= sum*constant_mode[i];
                }
                for(i=0;i<n;i++) {
                    e=cache->elements[m][b*cache->patch_capacity+i];
                    work[m][e] += solution[i]
                        /sqrt((double)cache->multiplicity[m][e]);
                }
            }
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(b=0;b<cache->interface_blocks[m];b++) {
                n=cache->interface_size[m][b];
                if(n==0)
                    continue;
                face=cache->interface_face[m][b];
                L=cache->interface_chol[m]
                    +(size_t)b*cache->interface_capacity
                        *cache->interface_capacity;
                for(i=0;i<n;i++) {
                    ref=cache->interface_elements[m]
                        [b*cache->interface_capacity+i];
                    if(ref>0) {
                        if(cache->multiplicity[m][ref]==0)
                            myerror(E,"ALA local partition weight is zero");
                        constant_mode[i]=1.0
                            /sqrt((double)cache->multiplicity[m][ref]);
                        rhs[i]=effective_r[m][ref]*constant_mode[i];
                    }
                    else {
                        ghost_index=-ref-1;
                        if(cache->halo_multiplicity[m][face][ghost_index]==0)
                            myerror(E,"ALA ghost partition weight is zero");
                        constant_mode[i]=1.0
                            /sqrt((double)cache->halo_multiplicity[m][face]
                                                [ghost_index]);
                        rhs[i]=ghost_r[m][face][ghost_index]
                            *constant_mode[i];
                    }
                }
                if(strcmp(E->control.ala_shallow_patch_velocity_solver,
                          "element_vanka")==0) {
                    numerator=denominator=0.0;
                    for(i=0;i<n;i++) {
                        numerator += constant_mode[i]*rhs[i];
                        denominator += constant_mode[i]*constant_mode[i];
                    }
                    sum=numerator/max(denominator,1.0e-300);
                    for(i=0;i<n;i++)
                        rhs[i] -= sum*constant_mode[i];
                }
                for(i=0;i<n;i++) {
                    sum=rhs[i];
                    for(j=0;j<i;j++)
                        sum -= L[i*cache->interface_capacity+j]*solution[j];
                    solution[i]=sum/L[i*cache->interface_capacity+i];
                }
                for(i=n-1;i>=0;i--) {
                    sum=solution[i];
                    for(j=i+1;j<n;j++)
                        sum -= L[j*cache->interface_capacity+i]*solution[j];
                    solution[i]=sum/L[i*cache->interface_capacity+i];
                }
                if(strcmp(E->control.ala_shallow_patch_velocity_solver,
                          "element_vanka")==0) {
                    numerator=denominator=0.0;
                    for(i=0;i<n;i++) {
                        numerator += constant_mode[i]*solution[i];
                        denominator += constant_mode[i]*constant_mode[i];
                    }
                    sum=numerator/max(denominator,1.0e-300);
                    for(i=0;i<n;i++)
                        solution[i] -= sum*constant_mode[i];
                }
                for(i=0;i<n;i++) {
                    ref=cache->interface_elements[m]
                        [b*cache->interface_capacity+i];
                    if(ref>0)
                        work[m][ref] += solution[i]
                            /sqrt((double)cache->multiplicity[m][ref]);
                    else {
                        ghost_index=-ref-1;
                        ghost_work[m][face][ghost_index] += solution[i]
                            /sqrt((double)cache->halo_multiplicity[m][face]
                                                [ghost_index]);
                    }
                }
            }
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            exchange_ala_shallow_halo_values(E,cache,lev,m,work,
                                              ghost_work[m],1);
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=1;e<=npno;e++)
                if(cache->multiplicity[m][e]>0) {
                    local_patch_energy[1] += effective_r[m][e]*work[m][e];
                    if(strcmp(E->control.ala_shallow_patch_velocity_solver,
                              "element_vanka")==0)
                        z[m][e] += weight*work[m][e];
                    else
                        z[m][e]=(1.0-weight)*z[m][e]+weight*work[m][e];
                }
        if(iteration==0 ||
           iteration%E->control.ala_coarse_residual_interval==0) {
            MPI_Allreduce(local_patch_energy,global_patch_energy,2,MPI_DOUBLE,
                          MPI_SUM,E->parallel.world);
            if(E->parallel.me==0) {
                fprintf(E->fp,"ALA_SHALLOW_PATCH_ENERGY iteration=%d "
                        "diagonal=%e block=%e block_to_diagonal=%e "
                        "application=%s weight=%e\n",iteration,
                        global_patch_energy[0],global_patch_energy[1],
                        global_patch_energy[1]
                        /max(global_patch_energy[0],1.0e-300),
                        strcmp(E->control.ala_shallow_patch_velocity_solver,
                               "element_vanka")==0
                            ? "additive_weighted_highpass" : "convex_blend",
                        weight);
                fflush(E->fp);
                }
        }
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(face=0;face<ALA_PATCH_MPI_FACES;face++) {
                free((void *)ghost_r[m][face]);
                free((void *)ghost_work[m][face]);
            }
    }

    if(E->control.ala_geneo_preconditioner)
        apply_ala_geneo_correction(E,effective_r,z,lev,iteration,cache);

    if(E->control.ala_pressure_multigrid)
        apply_ala_pressure_multigrid_correction(
            E,effective_r,z,lev,iteration,cache);

    if(!E->control.ala_two_level_preconditioner)
        goto pressure_depth_prolongation;

    clev=lev-E->control.ala_two_level_offset;
    factor=1 << E->control.ala_two_level_offset;
    celx=E->lmesh.ELX[clev];
    celz=E->lmesh.ELZ[clev];
    cnpno=E->lmesh.NPNO[clev];
    neq=E->lmesh.NEQ[lev];
    damping=E->control.ala_two_level_coarse_damping;

    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        coarse_rhs[m]=(double *)calloc(cnpno+1,sizeof(double));
        coarse_x[m]=(double *)calloc(cnpno+1,sizeof(double));
        coarse_residual[m]=(double *)calloc(cnpno+1,sizeof(double));
        coarse_Ax[m]=(double *)calloc(cnpno+1,sizeof(double));
        coarse_direction[m]=(double *)calloc(cnpno+1,sizeof(double));
        fine_p[m]=(double *)calloc(npno+1,sizeof(double));
        fine_Ap[m]=(double *)calloc(npno+1,sizeof(double));
        fine_velocity[m]=(double *)calloc(neq+1,sizeof(double));
        fine_velocity_rhs[m]=(double *)calloc(neq+1,sizeof(double));
        fine_velocity_Ax[m]=(double *)calloc(neq+1,sizeof(double));
        fine_velocity_direction[m]=(double *)calloc(neq+1,sizeof(double));
        if(coarse_rhs[m]==NULL || coarse_x[m]==NULL ||
           coarse_residual[m]==NULL || coarse_Ax[m]==NULL ||
           coarse_direction[m]==NULL || fine_p[m]==NULL ||
           fine_Ap[m]==NULL || fine_velocity[m]==NULL ||
           fine_velocity_rhs[m]==NULL || fine_velocity_Ax[m]==NULL ||
           fine_velocity_direction[m]==NULL)
            myerror(E,"Unable to allocate ALA two-level preconditioner");
    }

    /* P^T restriction: sum every factor^3 fine residuals into the
       geometrically coincident coarse pressure cell. */
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(ey=1;ey<=E->lmesh.ELY[lev];ey++)
            for(ex=1;ex<=E->lmesh.ELX[lev];ex++)
                for(ez=1;ez<=E->lmesh.ELZ[lev];ez++) {
                    e=ez+(ex-1)*E->lmesh.ELZ[lev]
                        +(ey-1)*E->lmesh.ELZ[lev]*E->lmesh.ELX[lev];
                    cx=(ex-1)/factor+1;
                    cy=(ey-1)/factor+1;
                    cz=(ez-1)/factor+1;
                    ce=cz+(cx-1)*celz+(cy-1)*celz*celx;
                    coarse_rhs[m][ce] += effective_r[m][e];
                }

    /* Both coarse solvers are fixed polynomials of the Galerkin ALA operator.
       With zero initial state and the BPI similarity transform, they define
       linear symmetric maps for the outer standard PCG. */
    if(strcmp(E->control.ala_two_level_coarse_solver,
              "scaled_diagonal")==0) {
        /* Balance the diagonal fine map before adding the recurrence-free
           coarse map.  With Q=P*(P^T*P)^-1*P^T and diagonal M, replace M by
           (I-Q)*M*(I-Q); this prevents the factor-two coarse component from
           being counted once by BPI and again by P*D_c^-1*P^T. */
        n=factor*factor*factor;
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(ey=1;ey<=E->lmesh.ELY[lev];ey++)
                for(ex=1;ex<=E->lmesh.ELX[lev];ex++)
                    for(ez=1;ez<=E->lmesh.ELZ[lev];ez++) {
                        e=ez+(ex-1)*E->lmesh.ELZ[lev]
                            +(ey-1)*E->lmesh.ELZ[lev]*E->lmesh.ELX[lev];
                        cx=(ex-1)/factor+1;
                        cy=(ey-1)/factor+1;
                        cz=(ez-1)/factor+1;
                        ce=cz+(cx-1)*celz+(cy-1)*celz*celx;
                        coarse_residual[m][ce] +=
                            E->control.ala_pressure_bpi_weight
                            *E->BPI[lev][m][e]*effective_r[m][e];
                        coarse_Ax[m][ce] +=
                            E->control.ala_pressure_bpi_weight
                            *E->BPI[lev][m][e]
                            *coarse_rhs[m][ce]/(double)n;
                    }
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(ey=1;ey<=E->lmesh.ELY[lev];ey++)
                for(ex=1;ex<=E->lmesh.ELX[lev];ex++)
                    for(ez=1;ez<=E->lmesh.ELZ[lev];ez++) {
                        e=ez+(ex-1)*E->lmesh.ELZ[lev]
                            +(ey-1)*E->lmesh.ELZ[lev]*E->lmesh.ELX[lev];
                        cx=(ex-1)/factor+1;
                        cy=(ey-1)/factor+1;
                        cz=(ez-1)/factor+1;
                        ce=cz+(cx-1)*celz+(cy-1)*celz*celx;
                        z[m][e] += -coarse_residual[m][ce]/(double)n
                            -E->control.ala_pressure_bpi_weight
                             *E->BPI[lev][m][e]
                             *coarse_rhs[m][ce]/(double)n
                            +coarse_Ax[m][ce]/(double)n;
                    }
        /* Repeated diagonal actions are accumulated algebraically so the
           coarse result is exactly alpha*P*D_c^-1*P^T. */
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(ce=1;ce<=cnpno;ce++)
                coarse_x[m][ce]=
                    E->control.ala_two_level_coarse_iterations*damping
                    *cache->coarse_bpi[m][ce]*coarse_rhs[m][ce];
    }
    else if(strcmp(E->control.ala_two_level_coarse_solver,"jacobi")==0) {
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(ce=1;ce<=cnpno;ce++)
                coarse_residual[m][ce]=coarse_rhs[m][ce];
        for(k=0;k<E->control.ala_two_level_coarse_iterations;k++) {
            for(m=1;m<=E->sphere.caps_per_proc;m++)
                for(ce=1;ce<=cnpno;ce++)
                    coarse_x[m][ce] += damping*cache->coarse_bpi[m][ce]
                        *coarse_residual[m][ce];
            if(k+1<E->control.ala_two_level_coarse_iterations) {
                apply_ala_galerkin_fixed_schur(
                    E,coarse_x,coarse_Ax,fine_p,fine_Ap,fine_velocity,
                    fine_velocity_rhs,fine_velocity_Ax,
                    fine_velocity_direction,cache,clev,lev,factor,NULL);
                for(m=1;m<=E->sphere.caps_per_proc;m++)
                    for(ce=1;ce<=cnpno;ce++)
                        coarse_residual[m][ce]=coarse_rhs[m][ce]
                            -coarse_Ax[m][ce];
            }
        }
    }
    else {
        theta=0.5*(cache->coarse_eigenvalue_max+
                   cache->coarse_eigenvalue_min);
        delta=0.5*(cache->coarse_eigenvalue_max-
                   cache->coarse_eigenvalue_min);
        sigma=theta/delta;
        rho_cheb=1.0/sigma;
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(ce=1;ce<=cnpno;ce++) {
                coarse_direction[m][ce]=cache->coarse_bpi[m][ce]
                    *coarse_rhs[m][ce]/theta;
                coarse_x[m][ce]=coarse_direction[m][ce];
            }
        for(k=1;k<E->control.ala_two_level_coarse_iterations;k++) {
            apply_ala_galerkin_fixed_schur(
                E,coarse_x,coarse_Ax,fine_p,fine_Ap,fine_velocity,
                fine_velocity_rhs,fine_velocity_Ax,
                fine_velocity_direction,cache,clev,lev,factor,NULL);
            rho_new=1.0/(2.0*sigma-rho_cheb);
            for(m=1;m<=E->sphere.caps_per_proc;m++)
                for(ce=1;ce<=cnpno;ce++) {
                    coarse_residual[m][ce]=coarse_rhs[m][ce]
                        -coarse_Ax[m][ce];
                    coarse_direction[m][ce]=rho_new*rho_cheb
                        *coarse_direction[m][ce]
                        +(2.0*rho_new/delta)*cache->coarse_bpi[m][ce]
                        *coarse_residual[m][ce];
                    coarse_x[m][ce] += coarse_direction[m][ce];
                }
            rho_cheb=rho_new;
        }
    }

    if(E->control.ala_global_coarse_preconditioner) {
        n=cache->global_basis_count;
        for(i=0;i<n;i++) {
            global_local_rhs[i]=0.0;
            for(m=1;m<=E->sphere.caps_per_proc;m++)
                for(ce=1;ce<=cnpno;ce++)
                    global_local_rhs[i] +=
                        cache->global_basis[m][i*(cnpno+1)+ce]
                        *coarse_rhs[m][ce];
        }
        MPI_Allreduce(global_local_rhs,global_rhs,n,MPI_DOUBLE,MPI_SUM,
                      E->parallel.world);
        for(i=0;i<n;i++) {
            sum=global_rhs[i];
            for(j=0;j<i;j++)
                sum -= cache->global_chol[i*n+j]*global_solution[j];
            global_solution[i]=sum/cache->global_chol[i*n+i];
        }
        for(i=n-1;i>=0;i--) {
            sum=global_solution[i];
            for(j=i+1;j<n;j++)
                sum -= cache->global_chol[j*n+i]*global_solution[j];
            global_solution[i]=sum/cache->global_chol[i*n+i];
        }
        global_scale=E->control.ala_global_coarse_weight
            /E->control.ala_two_level_coarse_weight;
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(ce=1;ce<=cnpno;ce++)
                for(i=0;i<n;i++)
                    coarse_x[m][ce] += global_scale
                        *cache->global_basis[m][i*(cnpno+1)+ce]
                        *global_solution[i];

        if(iteration==0 ||
           iteration%E->control.ala_coarse_residual_interval==0) {
            global_residual2=0.0;
            global_rhs2=0.0;
            for(i=0;i<n;i++) {
                sum=0.0;
                for(j=0;j<n;j++)
                    sum += 0.5*(cache->global_matrix[i*n+j]
                                +cache->global_matrix[j*n+i])
                        *global_solution[j];
                global_residual2 += (global_rhs[i]-sum)
                    *(global_rhs[i]-sum);
                global_rhs2 += global_rhs[i]*global_rhs[i];
            }
            if(E->parallel.me==0) {
                fprintf(E->fp,
                        "ALA_GLOBAL_COARSE_SOLVE iteration=%d basis=%d "
                        "projected_rhs_norm=%e "
                        "projected_residual_reduction=%e weight=%e\n",
                        iteration,n,sqrt(global_rhs2),
                        sqrt(global_residual2/max(global_rhs2,1.0e-300)),
                        E->control.ala_global_coarse_weight);
                fflush(E->fp);
            }
        }
    }

    if(iteration==0 ||
       iteration%E->control.ala_coarse_residual_interval==0) {
        apply_ala_galerkin_fixed_schur(
            E,coarse_x,coarse_Ax,fine_p,fine_Ap,fine_velocity,
            fine_velocity_rhs,fine_velocity_Ax,
            fine_velocity_direction,cache,clev,lev,factor,
            &velocity_residual_reduction);
        local_energy[0]=0.0;
        local_energy[1]=0.0;
        local_energy[2]=0.0;
        local_energy[3]=0.0;
        for(m=1;m<=E->sphere.caps_per_proc;m++) {
            for(e=1;e<=npno;e++)
                local_energy[0] += effective_r[m][e]*z[m][e];
            for(ce=1;ce<=cnpno;ce++) {
                local_energy[1] += E->control.ala_two_level_coarse_weight
                    *coarse_rhs[m][ce]*coarse_x[m][ce];
                local_energy[2] += coarse_rhs[m][ce]*coarse_rhs[m][ce];
                local_energy[3] += (coarse_rhs[m][ce]-coarse_Ax[m][ce])
                    *(coarse_rhs[m][ce]-coarse_Ax[m][ce]);
            }
        }
        MPI_Allreduce(local_energy,global_energy,4,MPI_DOUBLE,MPI_SUM,
                      E->parallel.world);
        if(E->parallel.me==0) {
            fprintf(E->fp,
                    "ALA_TWO_LEVEL_ENERGY iteration=%d fine=%e coarse=%e "
                    "coarse_to_fine=%e coarse_residual_reduction=%e "
                    "velocity_residual_reduction=%e\n",
                    iteration,global_energy[0],
                    global_energy[1],global_energy[1]
                    /max(global_energy[0],1.0e-300),
                    sqrt(global_energy[3]
                         /max(global_energy[2],1.0e-300)),
                    velocity_residual_reduction);
            fflush(E->fp);
        }
    }

    /* Constant P prolongation, added to the fine diagonal smoother. */
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(ey=1;ey<=E->lmesh.ELY[lev];ey++)
            for(ex=1;ex<=E->lmesh.ELX[lev];ex++)
                for(ez=1;ez<=E->lmesh.ELZ[lev];ez++) {
                    e=ez+(ex-1)*E->lmesh.ELZ[lev]
                        +(ey-1)*E->lmesh.ELZ[lev]*E->lmesh.ELX[lev];
                    cx=(ex-1)/factor+1;
                    cy=(ey-1)/factor+1;
                    cz=(ez-1)/factor+1;
                    ce=cz+(cx-1)*celz+(cy-1)*celz*celx;
                    z[m][e] += E->control.ala_two_level_coarse_weight
                        *coarse_x[m][ce];
                }

    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        free((void *)coarse_rhs[m]);
        free((void *)coarse_x[m]);
        free((void *)coarse_residual[m]);
        free((void *)coarse_Ax[m]);
        free((void *)coarse_direction[m]);
        free((void *)fine_p[m]);
        free((void *)fine_Ap[m]);
        free((void *)fine_velocity[m]);
        free((void *)fine_velocity_rhs[m]);
        free((void *)fine_velocity_Ax[m]);
        free((void *)fine_velocity_direction[m]);
    }

pressure_depth_prolongation:
    if(fabs(E->control.ala_pressure_shallow_action_scale-1.0)>1.0e-12) {
        elz=E->lmesh.ELZ[lev];
        depth_weight=sqrt(E->control.ala_pressure_shallow_action_scale);
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=1;e<=npno;e++) {
                ez=(e-1)%elz+1;
                depth_km=(1.0-0.5*(E->sx[m][3][ez]
                                  +E->sx[m][3][ez+1]))*E->data.radius_km;
                if(depth_km<E->control.ala_pressure_shallow_depth_km)
                    z[m][e] *= depth_weight;
            }
    }
}


static void audit_ala_shallow_patch_preconditioner(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache, int lev)
{
    int m,e,npno,pressure_multigrid_enabled;
    double local[4],global[4],left,right,defect;
    double *x[NCS],*y[NCS],*Mx[NCS],*My[NCS],*work[NCS];

    npno=E->lmesh.NPNO[lev];
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        x[m]=(double *)calloc(npno+1,sizeof(double));
        y[m]=(double *)calloc(npno+1,sizeof(double));
        Mx[m]=(double *)calloc(npno+1,sizeof(double));
        My[m]=(double *)calloc(npno+1,sizeof(double));
        work[m]=(double *)calloc(npno+1,sizeof(double));
        if(x[m]==NULL || y[m]==NULL || Mx[m]==NULL || My[m]==NULL ||
           work[m]==NULL)
            myerror(E,"Unable to allocate ALA Schwarz audit workspace");
        for(e=1;e<=npno;e++) {
            x[m][e]=sin(0.371*(double)e
                         +0.113*(double)(E->parallel.me+1));
            y[m][e]=cos(0.529*(double)e
                         +0.197*(double)(E->parallel.me+1));
        }
    }
    /* Certify the symmetric Schwarz baseline in isolation.  Stage 8's
       rediscretized V-cycle is FGMRES-only and uses its own positive-energy
       rollback, so the legacy PCG symmetry threshold does not apply to it. */
    pressure_multigrid_enabled=E->control.ala_pressure_multigrid;
    E->control.ala_pressure_multigrid=0;
    apply_ala_pressure_preconditioner(E,x,Mx,work,lev,-1,cache);
    apply_ala_pressure_preconditioner(E,y,My,work,lev,-1,cache);
    E->control.ala_pressure_multigrid=pressure_multigrid_enabled;
    for(e=0;e<4;e++)
        local[e]=0.0;
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(e=1;e<=npno;e++) {
            local[0] += x[m][e]*My[m][e];
            local[1] += y[m][e]*Mx[m][e];
            local[2] += x[m][e]*Mx[m][e];
            local[3] += y[m][e]*My[m][e];
        }
    MPI_Allreduce(local,global,4,MPI_DOUBLE,MPI_SUM,E->parallel.world);
    left=global[0];
    right=global[1];
    defect=fabs(left-right)/max(fabs(left)+fabs(right),1.0e-300);
    if(E->parallel.me==0) {
        fprintf(E->fp,"ALA MPI-overlap Schwarz audit "
                "scope=base_without_pressure_mg symmetry_defect=%e "
                "energy_x=%e energy_y=%e status=%s\n",defect,
                global[2],global[3],
                (defect<=1.0e-10 && global[2]>0.0 && global[3]>0.0)
                    ? "pass" : "fail");
        fprintf(stderr,"ALA MPI-overlap Schwarz audit "
                "scope=base_without_pressure_mg symmetry_defect=%e "
                "energy_x=%e energy_y=%e status=%s\n",defect,
                global[2],global[3],
                (defect<=1.0e-10 && global[2]>0.0 && global[3]>0.0)
                    ? "pass" : "fail");
        fflush(E->fp);
        fflush(stderr);
    }
    if(!isfinite(defect) || defect>1.0e-10 ||
       !isfinite(global[2]) || !isfinite(global[3]) ||
       global[2]<=0.0 || global[3]<=0.0)
        myerror(E,"ALA MPI-overlap Schwarz failed SPD audit");
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        free((void *)x[m]);
        free((void *)y[m]);
        free((void *)Mx[m]);
        free((void *)My[m]);
        free((void *)work[m]);
    }
}


/* Master loop for pressure and (hence) velocity field */

void solve_constrained_flow_iterative(E)
     struct All_variables *E;

{
    void v_from_vector();
    void p_to_nodes();

    int cycles;

    cycles=E->control.p_iterations;

    /* Solve for velocity and pressure, correct for bc's */

    solve_Ahat_p_fhat(E,E->U,E->P,E->F,E->control.accuracy,&cycles);

    v_from_vector(E);
    p_to_nodes(E,E->P,E->NP,E->mesh.levmax);

    return;
}


static void strict_ala_continuity_metrics(struct All_variables *E,
                                          double **V, double **r,
                                          double **div_u, int lev,
                                          double *mass_norm,
                                          double *cancellation_l2)
{
    int m, e, nel;
    double volume, div_e, c_e, scale_e;
    double local[2], global[2];
    void assemble_div_u();

    assemble_div_u(E, V, div_u, lev);
    local[0] = 0.0;
    local[1] = 0.0;
    nel = E->lmesh.NEL[lev];
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(e=1; e<=nel; e++) {
            volume = E->eco[m][e].area;
            if(volume <= 0.0)
                continue;
            div_e = div_u[m][e];
            c_e = r[m][e] - div_e;
            scale_e = fabs(div_e) + fabs(c_e);
            local[0] += r[m][e] * r[m][e] / volume;
            local[1] += scale_e * scale_e / volume;
        }

    MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_SUM,
                  E->parallel.world);
    *mass_norm = sqrt(max(global[0], 0.0));
    *cancellation_l2 = sqrt(global[0] / max(global[1], 1.0e-300));
}


/* The cancellation denominator is a physical-strength diagnostic, not an
 * algebraic residual norm. Track it explicitly so a shrinking velocity
 * field cannot be mistaken for improved D/C balance. */
static double strict_ala_continuity_term_strength(struct All_variables *E,
                                                  double **r,
                                                  double **div_u, int lev)
{
    int m, e, nel;
    double volume, div_e, c_e, scale_e, local, global;

    nel = E->lmesh.NEL[lev];
    local = 0.0;
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(e=1; e<=nel; e++) {
            volume = E->eco[m][e].area;
            if(volume <= 0.0)
                continue;
            div_e = div_u[m][e];
            c_e = r[m][e] - div_e;
            scale_e = fabs(div_e) + fabs(c_e);
            local += scale_e * scale_e / volume;
        }
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM,
                  E->parallel.world);
    return sqrt(max(global, 0.0));
}


/* Audit the original, unaugmented momentum equation at the current (V,P)
 * iterate. Reuse FGMRES work arrays and do not feed the result into any
 * Krylov recurrence. */
static void strict_ala_momentum_residual_audit(struct All_variables *E,
                                               double **V, double **P,
                                               double **grad_work,
                                               double **residual_work,
                                               int lev, double *rms,
                                               double *relative)
{
    int m, i, neq, gneq;
    double force_norm, residual_norm;
    void assemble_unaugmented_del2_u();
    void assemble_forces();
    void assemble_grad_p();
    void strip_bcs_from_residual();
    double global_vdot();

    neq = E->lmesh.neq;
    gneq = E->mesh.neq;
    /* Rebuild the force so its strict-ALA C^T P contribution matches the
     * current pressure iterate, exactly as in the post-solve audit. */
    assemble_forces(E,0);
    assemble_unaugmented_del2_u(E,V,residual_work,lev,1);
    assemble_grad_p(E,P,grad_work,lev);
    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        for(i=0; i<neq; i++)
            residual_work[m][i] = E->F[m][i] - grad_work[m][i]
                                 - residual_work[m][i];
    }
    strip_bcs_from_residual(E,residual_work,lev);
    force_norm = sqrt(max(global_vdot(E,E->F,E->F,lev),0.0));
    residual_norm = sqrt(max(global_vdot(E,residual_work,residual_work,lev),
                             0.0));
    *rms = residual_norm / sqrt(max((double)gneq,1.0));
    *relative = residual_norm / max(force_norm,1.0e-32);
}


/* Decompose the audited momentum defect without allocating another global
 * velocity field.  In residual-sign convention,
 *
 *   r_raw = f - G^T p - K u,
 *   r_aug = f - G^T p - K_gamma u = r_raw - r_pen,
 *   r_pen = gamma G^T M_p^-1 G u.
 *
 * The split defect independently checks K_gamma*u-K*u-r_pen.  This audit is
 * intentionally restricted to FGMRES startup/restart boundaries; candidate
 * acceptance retains the cheaper raw audit above. */
static void strict_ala_momentum_decomposition_audit(
    struct All_variables *E, double **V, double **P,
    double **work_a, double **work_b, int lev, const char *stage,
    int iteration, double *raw_rms, double *raw_relative)
{
    int m,i,neq,gneq;
    double force_norm,raw_norm,augmented_norm,penalty_norm;
    double predicted_penalty_norm,split_defect_norm,split_defect;
    double raw_penalty_dot,raw_penalty_cosine,predicted_augmented_norm;
    double residual_norm_defect;
    void assemble_forces();
    void assemble_del2_u();
    void assemble_unaugmented_del2_u();
    void assemble_ala_augmented_u();
    void assemble_grad_p();
    double global_vdot();

    neq=E->lmesh.neq;
    gneq=E->mesh.neq;
    assemble_forces(E,0);
    force_norm=sqrt(max(global_vdot(E,E->F,E->F,lev),0.0));

    /* Independently verify K_gamma-K equals the explicit AL action. */
    assemble_del2_u(E,V,work_a,lev,1);
    assemble_unaugmented_del2_u(E,V,work_b,lev,1);
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(i=0;i<neq;i++)
            work_a[m][i]-=work_b[m][i];
    predicted_penalty_norm=sqrt(max(global_vdot(
        E,work_a,work_a,lev),0.0));
    assemble_ala_augmented_u(E,V,work_b,lev,1);
    penalty_norm=sqrt(max(global_vdot(E,work_b,work_b,lev),0.0));
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(i=0;i<neq;i++)
            work_a[m][i]-=work_b[m][i];
    split_defect_norm=sqrt(max(global_vdot(E,work_a,work_a,lev),0.0));
    split_defect=split_defect_norm/max(max(predicted_penalty_norm,
                                           penalty_norm),1.0e-300);

    /* Assemble the augmented residual independently. */
    assemble_del2_u(E,V,work_a,lev,1);
    assemble_grad_p(E,P,work_b,lev);
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(i=0;i<neq;i++)
            work_a[m][i]=E->F[m][i]-work_b[m][i]-work_a[m][i];
    augmented_norm=sqrt(max(global_vdot(E,work_a,work_a,lev),0.0));

    /* Assemble the original unaugmented residual independently. */
    assemble_unaugmented_del2_u(E,V,work_a,lev,1);
    assemble_grad_p(E,P,work_b,lev);
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(i=0;i<neq;i++)
            work_a[m][i]=E->F[m][i]-work_b[m][i]-work_a[m][i];
    raw_norm=sqrt(max(global_vdot(E,work_a,work_a,lev),0.0));

    /* Quantify whether the raw defect is the AL penalty component. */
    assemble_ala_augmented_u(E,V,work_b,lev,1);
    raw_penalty_dot=global_vdot(E,work_a,work_b,lev);
    raw_penalty_cosine=raw_penalty_dot
        /max(raw_norm*penalty_norm,1.0e-300);
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(i=0;i<neq;i++)
            work_a[m][i]-=work_b[m][i];
    predicted_augmented_norm=sqrt(max(global_vdot(
        E,work_a,work_a,lev),0.0));
    residual_norm_defect=fabs(predicted_augmented_norm-augmented_norm)
        /max(max(predicted_augmented_norm,augmented_norm),1.0e-300);

    *raw_rms=raw_norm/sqrt(max((double)gneq,1.0));
    *raw_relative=raw_norm/max(force_norm,1.0e-32);
    if(E->parallel.me==0) {
        fprintf(E->fp,"ALA FGMRES momentum decomposition stage=%s "
                "iteration=%d raw_rms=%e raw_relative=%e "
                "augmented_rms=%e augmented_relative=%e "
                "penalty_rms=%e penalty_relative=%e "
                "raw_penalty_cosine=%e split_defect=%e "
                "residual_norm_defect=%e\n",stage,iteration,
                *raw_rms,*raw_relative,
                augmented_norm/sqrt(max((double)gneq,1.0)),
                augmented_norm/max(force_norm,1.0e-32),
                penalty_norm/sqrt(max((double)gneq,1.0)),
                penalty_norm/max(force_norm,1.0e-32),
                raw_penalty_cosine,split_defect,residual_norm_defect);
        fflush(E->fp);
    }
}


static void strict_ala_depth_diagnostics(struct All_variables *E,
                                         double **r, double **div_u,
                                         int lev, int iteration)
{
    int m, e, b, bins, nel, local_ez, global_ez, global_elz;
    int fields;
    double volume, div_e, c_e, scale_e, depth_km, r2;
    double cancellation, correlation, d_over_c, residual_fraction;
    double *local, *global;
    double local_summary[5], global_summary[5];

    if(!E->control.ala_depth_diagnostics)
        return;
    if(iteration != 0 &&
       (iteration % E->control.ala_depth_diagnostic_interval) != 0)
        return;

    global_elz = E->mesh.ELZ[lev];
    bins = min(E->control.ala_depth_diagnostic_bins, global_elz);
    fields = 7;
    local = (double *)calloc(bins * fields, sizeof(double));
    global = (double *)calloc(bins * fields, sizeof(double));
    if(local == NULL || global == NULL)
        myerror(E, "Unable to allocate ALA depth diagnostics");

    for(b=0; b<5; b++)
        local_summary[b] = 0.0;

    nel = E->lmesh.NEL[lev];
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(e=1; e<=nel; e++) {
            volume = E->eco[m][e].area;
            if(volume <= 0.0)
                continue;
            local_ez = (e - 1) % E->lmesh.ELZ[lev] + 1;
            global_ez = E->lmesh.EZS[lev] + local_ez;
            b = ((global_ez - 1) * bins) / global_elz;
            b = min(max(b, 0), bins - 1);

            depth_km = (1.0 - 0.5 *
                        (E->sx[m][3][local_ez] +
                         E->sx[m][3][local_ez+1]))
                * E->data.radius_km;
            div_e = div_u[m][e];
            c_e = r[m][e] - div_e;
            scale_e = fabs(div_e) + fabs(c_e);
            r2 = r[m][e] * r[m][e] / volume;

            local[b*fields] += div_e * div_e / volume;
            local[b*fields+1] += c_e * c_e / volume;
            local[b*fields+2] += r2;
            local[b*fields+3] += div_e * c_e / volume;
            local[b*fields+4] += scale_e * scale_e / volume;
            local[b*fields+5] += volume;
            local[b*fields+6] += depth_km * volume;

            local_summary[0] += r2;
            local_summary[1] += depth_km * r2;
            if(depth_km <= 200.0)
                local_summary[2] += r2;
            if(depth_km <= 410.0)
                local_summary[3] += r2;
            if(depth_km <= 660.0)
                local_summary[4] += r2;
        }

    MPI_Allreduce(local, global, bins * fields, MPI_DOUBLE, MPI_SUM,
                  E->parallel.world);
    MPI_Allreduce(local_summary, global_summary, 5, MPI_DOUBLE, MPI_SUM,
                  E->parallel.world);

    if(E->parallel.me == 0) {
        fprintf(E->fp,
                "ALA_DEPTH_CONTINUITY iteration=%d solver=%s bins=%d "
                "residual_mean_depth_km=%e top200_fraction=%e "
                "top410_fraction=%e top660_fraction=%e\n",
                iteration, E->control.uzawa, bins,
                global_summary[1] / max(global_summary[0], 1.0e-300),
                global_summary[2] / max(global_summary[0], 1.0e-300),
                global_summary[3] / max(global_summary[0], 1.0e-300),
                global_summary[4] / max(global_summary[0], 1.0e-300));
        fprintf(E->fp,
                "ALA_DEPTH_BIN bin depth_km residual_fraction cancellation "
                "D_over_C D_C_correlation\n");
        for(b=bins-1; b>=0; b--) {
            residual_fraction = global[b*fields+2]
                / max(global_summary[0], 1.0e-300);
            cancellation = sqrt(global[b*fields+2]
                                / max(global[b*fields+4], 1.0e-300));
            d_over_c = sqrt(global[b*fields]
                            / max(global[b*fields+1], 1.0e-300));
            correlation = global[b*fields+3]
                / sqrt(max(global[b*fields] * global[b*fields+1],
                           1.0e-300));
            correlation = min(max(correlation, -1.0), 1.0);
            depth_km = global[b*fields+6]
                / max(global[b*fields+5], 1.0e-300);
            fprintf(E->fp,
                    "ALA_DEPTH_BIN %d %e %e %e %e %e\n",
                    b, depth_km, residual_fraction, cancellation,
                    d_over_c, correlation);
        }
        fflush(E->fp);
    }

    free((void *)local);
    free((void *)global);
}


static void strict_ala_beta_causal_diagnostics(
    struct All_variables *E, double **V, double **active_r,
    double **alternate_r, double **div_u, int lev, int iteration)
{
    int m,e,b,bins,fields,nel,local_ez,global_ez,global_elz;
    double volume,depth_km,supplied,density,delta,div_e;
    double supplied_c,density_c,supplied_scale,density_scale;
    double supplied_cancellation,density_cancellation,correlation;
    double top410_correlation;
    double *local,*global;
    double local_summary[12],global_summary[12];
    double **supplied_r,**density_r;
    void assemble_div_rho_u_with_beta();

    if(!E->control.ala_beta_causal_diagnostics)
        return;
    if(iteration!=0 &&
       iteration%E->control.ala_depth_diagnostic_interval!=0)
        return;

    if(strcmp(E->control.ala_beta_element_source,
              "supplied_average")==0) {
        assemble_div_rho_u_with_beta(E,V,alternate_r,lev,
                                     E->refstate.ala_beta_density);
        supplied_r=active_r;
        density_r=alternate_r;
    }
    else {
        assemble_div_rho_u_with_beta(E,V,alternate_r,lev,
                                     E->refstate.ala_beta_supplied);
        supplied_r=alternate_r;
        density_r=active_r;
    }

    global_elz=E->mesh.ELZ[lev];
    bins=min(E->control.ala_depth_diagnostic_bins,global_elz);
    fields=3;
    local=(double *)calloc(bins*fields,sizeof(double));
    global=(double *)calloc(bins*fields,sizeof(double));
    if(local==NULL || global==NULL)
        myerror(E,"Unable to allocate ALA beta causal diagnostics");
    for(b=0;b<12;b++)
        local_summary[b]=0.0;

    nel=E->lmesh.NEL[lev];
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(e=1;e<=nel;e++) {
            volume=E->eco[m][e].area;
            if(volume<=0.0)
                continue;
            local_ez=(e-1)%E->lmesh.ELZ[lev]+1;
            global_ez=E->lmesh.EZS[lev]+local_ez;
            b=((global_ez-1)*bins)/global_elz;
            b=min(max(b,0),bins-1);
            depth_km=(1.0-0.5*(E->sx[m][3][local_ez]
                              +E->sx[m][3][local_ez+1]))
                *E->data.radius_km;
            supplied=supplied_r[m][e];
            density=density_r[m][e];
            delta=supplied-density;
            div_e=div_u[m][e];
            supplied_c=supplied-div_e;
            density_c=density-div_e;
            supplied_scale=fabs(div_e)+fabs(supplied_c);
            density_scale=fabs(div_e)+fabs(density_c);

            local[b*fields] += supplied*supplied/volume;
            local[b*fields+1] += density*density/volume;
            local[b*fields+2] += delta*delta/volume;
            local_summary[0] += supplied*supplied/volume;
            local_summary[1] += density*density/volume;
            local_summary[2] += delta*delta/volume;
            local_summary[3] += supplied*density/volume;
            local_summary[4] += supplied_scale*supplied_scale/volume;
            local_summary[5] += density_scale*density_scale/volume;
            if(depth_km<=200.0)
                local_summary[6] += delta*delta/volume;
            if(depth_km<=410.0)
                local_summary[7] += delta*delta/volume;
            if(depth_km<=660.0)
                local_summary[8] += delta*delta/volume;
            if(depth_km<=410.0) {
                local_summary[9] += supplied*supplied/volume;
                local_summary[10] += density*density/volume;
                local_summary[11] += supplied*density/volume;
            }
        }

    MPI_Allreduce(local,global,bins*fields,MPI_DOUBLE,MPI_SUM,
                  E->parallel.world);
    MPI_Allreduce(local_summary,global_summary,12,MPI_DOUBLE,MPI_SUM,
                  E->parallel.world);
    supplied_cancellation=sqrt(global_summary[0]
        /max(global_summary[4],1.0e-300));
    density_cancellation=sqrt(global_summary[1]
        /max(global_summary[5],1.0e-300));
    correlation=global_summary[3]
        /sqrt(max(global_summary[0]*global_summary[1],1.0e-300));
    correlation=min(max(correlation,-1.0),1.0);
    top410_correlation=global_summary[11]
        /sqrt(max(global_summary[9]*global_summary[10],1.0e-300));
    top410_correlation=min(max(top410_correlation,-1.0),1.0);

    if(E->parallel.me==0) {
        fprintf(E->fp,
                "ALA_BETA_CAUSAL iteration=%d active=%s "
                "supplied_cancellation=%e density_cancellation=%e "
                "delta_over_supplied=%e delta_over_density=%e "
                "residual_correlation=%e top200_delta_fraction=%e "
                "top410_delta_fraction=%e top660_delta_fraction=%e "
                "top410_supplied_to_density=%e "
                "top410_residual_correlation=%e\n",
                iteration,E->control.ala_beta_element_source,
                supplied_cancellation,density_cancellation,
                sqrt(global_summary[2]/max(global_summary[0],1.0e-300)),
                sqrt(global_summary[2]/max(global_summary[1],1.0e-300)),
                correlation,
                global_summary[6]/max(global_summary[2],1.0e-300),
                global_summary[7]/max(global_summary[2],1.0e-300),
                global_summary[8]/max(global_summary[2],1.0e-300),
                sqrt(global_summary[9]/max(global_summary[10],1.0e-300)),
                top410_correlation);
        fprintf(E->fp,
                "ALA_BETA_DEPTH_BIN bin supplied_fraction density_fraction "
                "delta_fraction delta_over_supplied\n");
        for(b=bins-1;b>=0;b--)
            fprintf(E->fp,"ALA_BETA_DEPTH_BIN %d %e %e %e %e\n",
                    b,
                    global[b*fields]/max(global_summary[0],1.0e-300),
                    global[b*fields+1]/max(global_summary[1],1.0e-300),
                    global[b*fields+2]/max(global_summary[2],1.0e-300),
                    sqrt(global[b*fields+2]
                         /max(global[b*fields],1.0e-300)));
        fflush(E->fp);
    }
    free((void *)local);
    free((void *)global);
}


static void strict_ala_coarse_residual_diagnostics(
    struct All_variables *E, double **r, int lev, int iteration)
{
    int m, e, ex, ey, ez, offset, factor;
    int elx, ely, elz, parent_elx, parent_ely, parent_elz, parent_nel;
    int parent_ex, parent_ey, parent_e;
    double volume, fine_energy, coarse_energy, coarse_fraction;
    double local_fine_energy, global_fine_energy;
    double local_coarse_energy, global_coarse_energy;
    double *parent_r, *parent_volume;

    if(!E->control.ala_coarse_residual_diagnostics)
        return;
    if(iteration != 0 &&
       (iteration % E->control.ala_coarse_residual_interval) != 0)
        return;

    elx = E->lmesh.ELX[lev];
    ely = E->lmesh.ELY[lev];
    elz = E->lmesh.ELZ[lev];
    local_fine_energy = 0.0;
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(e=1; e<=E->lmesh.NEL[lev]; e++) {
            volume = E->eco[m][e].area;
            if(volume > 0.0)
                local_fine_energy += r[m][e] * r[m][e] / volume;
        }
    MPI_Allreduce(&local_fine_energy, &global_fine_energy, 1, MPI_DOUBLE,
                  MPI_SUM, E->parallel.world);

    if(E->parallel.me == 0)
        fprintf(E->fp,
                "ALA_COARSE_RESIDUAL iteration=%d solver=%s fine_energy=%e\n",
                iteration, E->control.uzawa, global_fine_energy);

    factor = 1;
    for(offset=1; offset<=E->control.ala_coarse_residual_levels; offset++) {
        factor *= 2;
        if(elx % factor || ely % factor || elz % factor)
            myerror(E, "ALA coarse residual levels require locally divisible "
                    "element dimensions");

        parent_elx = elx / factor;
        parent_ely = ely / factor;
        parent_elz = elz / factor;
        parent_nel = parent_elx * parent_ely * parent_elz;
        parent_r = (double *)calloc(parent_nel, sizeof(double));
        parent_volume = (double *)calloc(parent_nel, sizeof(double));
        if(parent_r == NULL || parent_volume == NULL)
            myerror(E, "Unable to allocate ALA coarse residual diagnostics");

        local_coarse_energy = 0.0;
        for(m=1; m<=E->sphere.caps_per_proc; m++) {
            memset(parent_r, 0, parent_nel * sizeof(double));
            memset(parent_volume, 0, parent_nel * sizeof(double));
            for(ey=1; ey<=ely; ey++)
                for(ex=1; ex<=elx; ex++)
                    for(ez=1; ez<=elz; ez++) {
                        e = ez + (ex-1)*elz + (ey-1)*elz*elx;
                        volume = E->eco[m][e].area;
                        if(volume <= 0.0)
                            continue;
                        parent_ex = (ex-1) / factor;
                        parent_ey = (ey-1) / factor;
                        parent_e = (ez-1) / factor
                            + parent_ex*parent_elz
                            + parent_ey*parent_elz*parent_elx;
                        parent_r[parent_e] += r[m][e];
                        parent_volume[parent_e] += volume;
                    }
            for(parent_e=0; parent_e<parent_nel; parent_e++)
                if(parent_volume[parent_e] > 0.0)
                    local_coarse_energy += parent_r[parent_e]
                        * parent_r[parent_e] / parent_volume[parent_e];
        }
        MPI_Allreduce(&local_coarse_energy, &global_coarse_energy, 1,
                      MPI_DOUBLE, MPI_SUM, E->parallel.world);

        fine_energy = max(global_fine_energy, 0.0);
        coarse_energy = max(global_coarse_energy, 0.0);
        coarse_fraction = (fine_energy > 0.0)
            ? coarse_energy / fine_energy : 0.0;
        if(coarse_fraction < -1.0e-12 || coarse_fraction > 1.0+1.0e-12)
            myerror(E, "ALA coarse residual fraction is outside [0,1]");
        coarse_fraction = min(max(coarse_fraction, 0.0), 1.0);
        if(E->parallel.me == 0)
            fprintf(E->fp,
                    "ALA_COARSE_LEVEL offset=%d coarse_fraction=%e "
                    "coarse_energy=%e\n",
                    offset, coarse_fraction, coarse_energy);

        free((void *)parent_r);
        free((void *)parent_volume);
    }
    if(E->parallel.me == 0)
        fflush(E->fp);
}

void solve_constrained_flow_iterative_pseudo_surf(E)
     struct All_variables *E;

{
    void v_from_vector_pseudo_surf();
    void p_to_nodes();

    int cycles;

    cycles=E->control.p_iterations;

    /* Solve for velocity and pressure, correct for bc's */

    solve_Ahat_p_fhat(E,E->U,E->P,E->F,E->control.accuracy,&cycles);

    v_from_vector_pseudo_surf(E);
    p_to_nodes(E,E->P,E->NP,E->mesh.levmax);

    return;
}


/* ========================================================================= */

static void print_convergence_progress(struct All_variables *E,
                                       int count, double time0,
                                       double sq_vdotv,
                                       double dvelocity, double dpressure)
{
    double CPU_time0();

    fprintf(E->fp, "AhatP (%03d) after %6.2f s v=%.3e  div/v=%.3e "
            "dv/v=%.3e and dp/p=%.3e for step %d\n",
            count, CPU_time0()-time0, sq_vdotv,E->monitor.incompressibility,
            dvelocity, dpressure, E->monitor.solution_cycles);
    fprintf(stderr, "AhatP (%03d) after %6.2f s v=%.3e div/v=%.3e "
            "dv/v=%.3e and dp/p=%.3e for step %d\n",
            count, CPU_time0()-time0, sq_vdotv, E->monitor.incompressibility,
            dvelocity, dpressure, E->monitor.solution_cycles);

    return;
}


static double strict_ala_inner_relative_accuracy(
    struct All_variables *E, double outer_relative_residual)
{
    double candidate, final_accuracy, relative_accuracy;

    /* The user cap may deliberately request a tighter K solve than the
     * conventional 0.1*outer-tolerance forcing target. */
    final_accuracy = E->control.ala_inner_accuracy_max;
    if(E->control.tole_comp > 0.0)
        final_accuracy = min(final_accuracy,
                             0.1 * E->control.tole_comp);
    candidate = E->control.ala_inner_accuracy_factor
        * max(outer_relative_residual, 0.0);

    if(candidate <= final_accuracy)
        relative_accuracy = final_accuracy;
    else
        /* Quantize by decades so the inexact Schur operator only changes at
         * explicit Krylov restart points. */
        relative_accuracy = pow(10.0, ceil(log10(candidate) - 1.0e-12));

    relative_accuracy = min(relative_accuracy,
                            E->control.ala_inner_accuracy_max);
    return max(relative_accuracy, final_accuracy);
}


static double strict_ala_inner_accuracy(struct All_variables *E,
                                        double **F, int lev,
                                        double relative_accuracy)
{
    double global_vdot();
    double rhs_norm;

    rhs_norm = sqrt(global_vdot(E, F, F, lev) / E->mesh.neq);

    return max(relative_accuracy * rhs_norm, 1.0e-14);
}



static float solve_Ahat_p_fhat(struct All_variables *E,
                               double **V, double **P, double **F,
                               double imp, int *steps_max)
{
    float residual;

    if(E->control.inv_gruneisen == 0)
        residual = solve_Ahat_p_fhat_CG(E, V, P, F, imp, steps_max);
    else {
        if(strcmp(E->control.uzawa, "cg") == 0)
            residual = solve_Ahat_p_fhat_iterCG(E, V, P, F, imp, steps_max);
        else if(strcmp(E->control.uzawa, "ala_cg") == 0)
            residual = solve_Ahat_p_fhat_ALA_PCG(E, V, P, F, imp,
                                                steps_max);
        else if(strcmp(E->control.uzawa, "bicg") == 0)
            residual = solve_Ahat_p_fhat_BiCG(E, V, P, F, imp, steps_max);
        else
            myerror(E, "Error: unknown Uzawa iteration\n");
    }

    return(residual);
}


/* Solve incompressible Stokes flow using
 * conjugate gradient (CG) iterations
 */

static float solve_Ahat_p_fhat_CG(struct All_variables *E,
                                  double **V, double **P, double **FF,
                                  double imp, int *steps_max)
{
    int m, j, count, valid, lev, npno, neq;
    int gnpno, gneq;

    double *r1[NCS], *r2[NCS], *z1[NCS], *s1[NCS], *s2[NCS], *F[NCS];
    double *shuffle[NCS];
    double alpha, delta, r0dotz0, r1dotz1,sq_vdotv;
    double residual, v_res;

    double global_vdot(), global_pdot();

    double time0, CPU_time0();
    float dpressure, dvelocity;

    void assemble_c_u();
    void assemble_div_u();
    void assemble_del2_u();
    void assemble_grad_p();
    void strip_bcs_from_residual();
    int  solve_del2_u();
    void parallel_process_termination();

    gnpno = E->mesh.npno;
    gneq = E->mesh.neq;
    npno = E->lmesh.npno;
    neq = E->lmesh.neq;
    lev = E->mesh.levmax;

    for (m=1; m<=E->sphere.caps_per_proc; m++)   {
        F[m] = (double *)malloc(neq*sizeof(double));
        r1[m] = (double *)malloc((npno+1)*sizeof(double));
        r2[m] = (double *)malloc((npno+1)*sizeof(double));
        z1[m] = (double *)malloc((npno+1)*sizeof(double));
        s1[m] = (double *)malloc((npno+1)*sizeof(double));
        s2[m] = (double *)malloc((npno+1)*sizeof(double));
    }

    time0 = CPU_time0();
    count = 0;


    /* copy the original force vector since we need to keep it intact
       between iterations */
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(j=0;j<neq;j++)
            F[m][j] = FF[m][j];


    /* calculate the contribution of compressibility in the continuity eqn */
    if(E->control.inv_gruneisen != 0) {
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(j=1;j<=npno;j++)
                r2[m][j] = 0.0;

        assemble_c_u(E, V, r2, lev);
    }


    /* calculate the initial velocity residual */
    v_res = initial_vel_residual(E, V, P, F, imp);


    /* initial residual r1 = div(V) */
    assemble_div_u(E, V, r1, lev);


    /* add the contribution of compressibility to the initial residual */
    if(E->control.inv_gruneisen != 0)
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(j=1;j<=npno;j++) {
                r1[m][j] += r2[m][j];
            }

    residual = incompressibility_residual(E, V, r1);


    sq_vdotv = sqrt(E->monitor.vdotv);

    /* pressure and velocity corrections */
    dpressure = 1.0;
    dvelocity = 1.0;


    if (E->control.print_convergence && E->parallel.me==0)  {
        print_convergence_progress(E, count, time0, sq_vdotv,
                                   dvelocity, dpressure);
    }


    r0dotz0 = 0;

    while( (count < *steps_max) &&
           (E->monitor.incompressibility >= E->control.tole_comp) &&
           (dpressure >= imp) && (dvelocity >= imp) )  {


        /* preconditioner BPI ~= inv(K), z1 = BPI*r1 */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=1; j<=npno; j++)
                z1[m][j] = E->BPI[lev][m][j] * r1[m][j];


        /* r1dotz1 = <r1, z1> */
        r1dotz1 = global_pdot(E, r1, z1, lev);
        assert(r1dotz1 != 0.0  /* Division by zero in head of incompressibility iteration */);


        /* update search direction */
        if(count == 0)
            for (m=1; m<=E->sphere.caps_per_proc; m++)
                for(j=1; j<=npno; j++)
                    s2[m][j] = z1[m][j];
        else {
            /* s2 = z1 + s1 * <r1,z1>/<r0,z0> */
            delta = r1dotz1 / r0dotz0;
            for(m=1; m<=E->sphere.caps_per_proc; m++)
                for(j=1; j<=npno; j++)
                    s2[m][j] = z1[m][j] + delta * s1[m][j];
        }


        /* solve K*u1 = grad(s2) for u1 */
        assemble_grad_p(E, s2, F, lev);
        valid = solve_del2_u(E, E->u1, F, imp*v_res, lev);
        if(!valid && (E->parallel.me==0)) {
            fputs("Warning: solver not converging! 1\n", stderr);
            fputs("Warning: solver not converging! 1\n", E->fp);
	    fprintf(stderr,"Solver not converging! \n");
	    fflush(stderr);
	    exit(18);
        }
        strip_bcs_from_residual(E, E->u1, lev);


        /* The outer CG path follows CitcomS' split compressible strategy:
         * C*u is held on the outer-loop right hand side. */
        assemble_div_u(E, E->u1, F, lev);


        /* alpha = <r1, z1> / <s2, F> */
        if(valid)
            /* alpha defined this way is the same as R&W */
            alpha = r1dotz1 / global_pdot(E, s2, F, lev);
        else
            alpha = 0.0;


        /* r2 = r1 - alpha * div(u1) */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=1; j<=npno; j++)
                r2[m][j] = r1[m][j] - alpha * F[m][j];


        /* P = P + alpha * s2 */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=1; j<=npno; j++)
                P[m][j] += alpha * s2[m][j];


        /* V = V - alpha * u1 */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=0; j<neq; j++)
                V[m][j] -= alpha * E->u1[m][j];


        /* compute velocity and incompressibility residual */
        assemble_div_u(E, V, F, lev);
        incompressibility_residual(E, V, F);

        /* compute velocity and pressure corrections */
        dpressure = alpha * sqrt(global_pdot(E, s2, s2, lev)
                                 / (1.0e-32 + global_pdot(E, P, P, lev)));
        dvelocity = alpha * sqrt(global_vdot(E, E->u1, E->u1, lev)
                                 / (1.0e-32 + E->monitor.vdotv));

        count++;

	sq_vdotv = sqrt(E->monitor.vdotv);

        if (E->control.print_convergence && E->parallel.me==0)  {
            print_convergence_progress(E, count, time0, sq_vdotv,
                                       dvelocity, dpressure);
        }


        /* shift array pointers */
        for(m=1; m<=E->sphere.caps_per_proc; m++) {
            shuffle[m] = s1[m];
            s1[m] = s2[m];
            s2[m] = shuffle[m];

            shuffle[m] = r1[m];
            r1[m] = r2[m];
            r2[m] = shuffle[m];
        }

        /* shift <r0, z0> = <r1, z1> */
        r0dotz0 = r1dotz1;

    } /* end loop for conjugate gradient */

    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        free((void *) F[m]);
        free((void *) r1[m]);
        free((void *) r2[m]);
        free((void *) z1[m]);
        free((void *) s1[m]);
        free((void *) s2[m]);
    }

    *steps_max=count;

    return(residual);
}

/* Solve strict ALA pressure correction using preconditioned conjugate
 * gradients on S=(D+C) K^-1 (D+C)^T.  Unlike the legacy compressible CG
 * path, this applies the complete continuity operator and its transpose.
 */

static float solve_Ahat_p_fhat_ALA_PCG(struct All_variables *E,
                                       double **V, double **P, double **FF,
                                       double imp, int *steps_max)
{
    void assemble_div_rho_u();
    void assemble_grad_rho_p();
    void strip_bcs_from_residual();
    void parallel_process_termination();
    int solve_del2_u();

    double global_vdot(), global_pdot();
    double CPU_time0();

    int npno, neq, lev, coarse_lev, coarse_factor, validation_level;
    int m, j, count, valid;
    int restart_search;
    int hybrid_consecutive_count, hybrid_converged;
    int nonpositive_curvature_count;
    int audit_best_iteration, audit_window_iteration;
    int audit_milestone, audit_stagnated;
    int local_invalid_bpi, global_invalid_bpi;
    int local_invalid_velocity_bi, global_invalid_velocity_bi;
    int galerkin_diagnostic_applications, galerkin_applications;
    int coupled_self_contained_preconditioner;
    const char *preconditioner_mode;

    double alpha, beta, rho, rho_old, curvature;
    double min_curvature, sq_vdotv;
    double residual, dpressure, dvelocity;
    double initial_rnorm, rnorm, relative_residual;
    double recursive_rnorm, recursive_relative_residual, drift_ratio;
    double initial_mass_norm, mass_norm, mass_relative_residual;
    double cancellation_l2;
    double initial_term_strength, term_strength;
    double momentum_rms, momentum_relative;
    double inner_accuracy, inner_relative_accuracy;
    double local_bpi_min, local_bpi_max;
    double global_bpi_min, global_bpi_max;
    double audit_best_cancellation, audit_window_cancellation;
    double audit_window_reduction;
    double audit_milestones[6];
    double time0;

    double *F[NCS];
    double *r[NCS], *z[NCS], *p[NCS], *q[NCS];
    double *explicit_r[NCS], *div_u[NCS];
    double *preconditioner_work[NCS];
    struct ala_pressure_preconditioner_cache preconditioner_cache;

    npno = E->lmesh.npno;
    neq = E->lmesh.neq;
    lev = E->mesh.levmax;
    memset(&preconditioner_cache,0,sizeof(preconditioner_cache));
    coupled_self_contained_preconditioner =
        strcmp(E->control.ala_outer_solver,"coupled_fgmres")==0 &&
        (E->control.ala_coupled_multilevel_vcycle ||
         E->control.ala_coupled_element_vanka);

    if(!coupled_self_contained_preconditioner &&
       (E->control.ala_two_level_preconditioner ||
        E->control.ala_geneo_preconditioner)) {
        coarse_lev=lev-E->control.ala_two_level_offset;
        if(coarse_lev<E->mesh.levmin)
            myerror(E,"ALA two-level offset is below the available mesh hierarchy");
        coarse_factor=1 << E->control.ala_two_level_offset;
        if(coarse_factor>4)
            myerror(E,"ALA two-level aggregate factor above four is unsupported");
        if(E->lmesh.ELX[lev]!=coarse_factor*E->lmesh.ELX[coarse_lev] ||
           E->lmesh.ELY[lev]!=coarse_factor*E->lmesh.ELY[coarse_lev] ||
           E->lmesh.ELZ[lev]!=coarse_factor*E->lmesh.ELZ[coarse_lev])
            myerror(E,"ALA two-level pressure transfer does not match mesh hierarchy");
    }
    if(!coupled_self_contained_preconditioner &&
       E->control.ala_pressure_multigrid)
        for(validation_level=E->control.ala_pressure_multigrid_min_level+1;
            validation_level<=lev;validation_level++)
            if(E->lmesh.ELX[validation_level]
                 !=2*E->lmesh.ELX[validation_level-1] ||
               E->lmesh.ELY[validation_level]
                 !=2*E->lmesh.ELY[validation_level-1] ||
               E->lmesh.ELZ[validation_level]
                 !=2*E->lmesh.ELZ[validation_level-1])
                myerror(E,"ALA pressure multigrid hierarchy is not factor two");

    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        F[m] = (double *)malloc(neq*sizeof(double));
        preconditioner_work[m] =
            (double *)malloc((npno+1)*sizeof(double));
        if(coupled_self_contained_preconditioner) {
            r[m]=z[m]=p[m]=q[m]=NULL;
            explicit_r[m]=div_u[m]=NULL;
        }
        else {
            r[m] = (double *)malloc((npno+1)*sizeof(double));
            z[m] = (double *)malloc((npno+1)*sizeof(double));
            p[m] = (double *)malloc((npno+1)*sizeof(double));
            q[m] = (double *)malloc((npno+1)*sizeof(double));
            explicit_r[m] = (double *)malloc((npno+1)*sizeof(double));
            div_u[m] = (double *)malloc((npno+1)*sizeof(double));
        }
    }

    time0 = CPU_time0();
    count = 0;

    global_bpi_min=global_bpi_max=0.0;
    global_invalid_bpi=global_invalid_velocity_bi=0;
    if(!coupled_self_contained_preconditioner) {
        local_bpi_min = 1.0e300;
        local_bpi_max = 0.0;
        local_invalid_bpi = 0;
        for(validation_level=E->control.ala_pressure_multigrid
                ? E->control.ala_pressure_multigrid_min_level : lev;
            validation_level<=lev;validation_level++)
            for(m=1; m<=E->sphere.caps_per_proc; m++)
                for(j=1; j<=E->lmesh.NPNO[validation_level]; j++) {
                    if(!isfinite(E->BPI[validation_level][m][j]) ||
                       E->BPI[validation_level][m][j] <= 0.0)
                        local_invalid_bpi++;
                    else {
                        local_bpi_min=min(local_bpi_min,
                            E->BPI[validation_level][m][j]);
                        local_bpi_max=max(local_bpi_max,
                            E->BPI[validation_level][m][j]);
                    }
                }
        MPI_Allreduce(&local_bpi_min,&global_bpi_min,1,MPI_DOUBLE,MPI_MIN,
                      E->parallel.world);
        MPI_Allreduce(&local_bpi_max,&global_bpi_max,1,MPI_DOUBLE,MPI_MAX,
                      E->parallel.world);
        MPI_Allreduce(&local_invalid_bpi,&global_invalid_bpi,1,MPI_INT,
                      MPI_SUM,E->parallel.world);
        local_invalid_velocity_bi=0;
        if(E->control.ala_shallow_patch_preconditioner ||
           E->control.ala_two_level_preconditioner ||
           E->control.ala_geneo_preconditioner ||
           E->control.ala_pressure_multigrid)
            for(validation_level=E->control.ala_pressure_multigrid
                    ? E->control.ala_pressure_multigrid_min_level : lev;
                validation_level<=lev;validation_level++)
                for(m=1;m<=E->sphere.caps_per_proc;m++)
                    for(j=0;j<E->lmesh.NEQ[validation_level];j++)
                        if(!isfinite(
                              E->ALA_velocity_BI[validation_level][m][j]) ||
                           E->ALA_velocity_BI[validation_level][m][j]<=0.0)
                            local_invalid_velocity_bi++;
        MPI_Allreduce(&local_invalid_velocity_bi,
                      &global_invalid_velocity_bi,1,MPI_INT,MPI_SUM,
                      E->parallel.world);
    }
    if(E->parallel.me==0 &&
       (global_invalid_bpi || global_invalid_velocity_bi)) {
        fprintf(stderr,"ALA preconditioner diagonal validation failed "
                "fine_pressure_invalid=%d fixed_velocity_invalid=%d\n",
                global_invalid_bpi,global_invalid_velocity_bi);
        fflush(stderr);
    }
    if(global_invalid_bpi)
        myerror(E,"ALA fine pressure preconditioner diagonal is not positive");
    if(global_invalid_velocity_bi)
        myerror(E,"ALA fixed velocity inverse diagonal is not positive");
    if(E->parallel.me==0) {
        fprintf(stderr,"ALA preconditioner startup stage=begin\n");
        fflush(stderr);
    }
    if(!coupled_self_contained_preconditioner &&
       E->control.ala_shallow_patch_preconditioner) {
        if(E->parallel.me==0) {
            fprintf(stderr,"ALA preconditioner startup stage=shallow_patch_build_begin\n");
            fflush(stderr);
        }
        build_ala_shallow_patch_cache(E,&preconditioner_cache,lev);
        if(E->parallel.me==0) {
            fprintf(stderr,"ALA preconditioner startup stage=shallow_patch_build_complete\n");
            fflush(stderr);
        }
    }
    if(!coupled_self_contained_preconditioner &&
       E->control.ala_pressure_multigrid) {
        if(E->parallel.me==0) {
            fprintf(stderr,"ALA preconditioner startup stage=pressure_multigrid_build_begin\n");
            fflush(stderr);
        }
        build_ala_pressure_multigrid_cache(E,&preconditioner_cache,lev);
        if(E->parallel.me==0) {
            fprintf(stderr,"ALA preconditioner startup stage=pressure_multigrid_build_complete\n");
            fflush(stderr);
        }
    }
    if(!coupled_self_contained_preconditioner &&
       (E->control.ala_two_level_preconditioner ||
        E->control.ala_geneo_preconditioner)) {
        if(E->parallel.me==0) {
            fprintf(stderr,"ALA preconditioner startup stage=galerkin_diagonal_build_begin\n");
            fflush(stderr);
        }
        build_ala_two_level_cache(E,&preconditioner_cache,lev);
        if(E->parallel.me==0) {
            fprintf(stderr,"ALA preconditioner startup stage=galerkin_diagonal_build_complete\n");
            fflush(stderr);
        }
        if(strcmp(E->control.ala_two_level_velocity_solver,"chebyshev")==0) {
            if(E->parallel.me==0) {
                fprintf(stderr,"ALA preconditioner startup stage=velocity_spectrum_begin\n");
                fflush(stderr);
            }
            calibrate_ala_velocity_spectrum(
                E,&preconditioner_cache,lev);
            if(E->parallel.me==0) {
                fprintf(stderr,"ALA preconditioner startup stage=velocity_spectrum_complete\n");
                fflush(stderr);
            }
        }
        if(E->control.ala_geneo_preconditioner) {
            if(E->parallel.me==0) {
                fprintf(stderr,"ALA preconditioner startup stage=geneo_build_begin\n");
                fflush(stderr);
            }
            if(E->control.ala_geneo_rank_group_x>1 ||
               E->control.ala_geneo_rank_group_y>1)
                build_ala_cross_rank_geneo_coarse_cache(
                    E,&preconditioner_cache,lev);
            else
                build_ala_geneo_coarse_cache(
                    E,&preconditioner_cache,lev);
            if(E->parallel.me==0) {
                fprintf(stderr,"ALA preconditioner startup stage=geneo_build_complete\n");
                fflush(stderr);
            }
        }
        if(E->control.ala_two_level_preconditioner &&
           E->control.ala_global_coarse_preconditioner) {
            if(E->parallel.me==0) {
                fprintf(stderr,"ALA preconditioner startup stage=global_coarse_build_begin\n");
                fflush(stderr);
            }
            build_ala_global_coarse_cache(
                E,&preconditioner_cache,lev);
            if(E->parallel.me==0) {
                fprintf(stderr,"ALA preconditioner startup stage=global_coarse_build_complete\n");
                fflush(stderr);
            }
        }
        if(E->control.ala_two_level_preconditioner &&
           strcmp(E->control.ala_two_level_coarse_solver,"chebyshev")==0) {
            if(E->parallel.me==0) {
                fprintf(stderr,"ALA preconditioner startup stage=galerkin_spectrum_begin\n");
                fflush(stderr);
            }
            calibrate_ala_two_level_spectrum(
                E,&preconditioner_cache,lev);
            if(E->parallel.me==0) {
                fprintf(stderr,"ALA preconditioner startup stage=galerkin_spectrum_complete\n");
                fflush(stderr);
            }
        }
    }
    if(!coupled_self_contained_preconditioner &&
       E->control.ala_shallow_patch_preconditioner)
        audit_ala_shallow_patch_preconditioner(
            E,&preconditioner_cache,lev);
    if(coupled_self_contained_preconditioner)
        preconditioner_mode=E->control.ala_coupled_multilevel_vcycle
            ? "coupled_multilevel_vcycle" : "coupled_element_vanka";
    else if(E->control.ala_pressure_multigrid)
        preconditioner_mode=E->control.ala_shallow_patch_preconditioner
            ? "mpi_overlap_schwarz_plus_pressure_vcycle"
            : "pressure_vcycle";
    else if(E->control.ala_geneo_preconditioner)
        preconditioner_mode=
            (E->control.ala_geneo_rank_group_x>1 ||
             E->control.ala_geneo_rank_group_y>1)
            ? (strcmp(E->control.ala_geneo_basis_type,
                      "radial_partition")==0
               ? "mpi_overlap_schwarz_plus_cross_rank_radial_aggregate"
               : "mpi_overlap_schwarz_plus_cross_rank_geneo")
            : "mpi_overlap_schwarz_plus_geneo";
    else if(E->control.ala_two_level_preconditioner)
        preconditioner_mode=E->control.ala_shallow_patch_preconditioner
            ? (strcmp(E->control.ala_two_level_velocity_solver,"chebyshev")==0
               ? "mpi_overlap_schwarz_plus_galerkin_kpoly"
               : "mpi_overlap_schwarz_plus_galerkin_diagonal")
            : (strcmp(E->control.ala_two_level_velocity_solver,"chebyshev")==0
               ? "true_galerkin_kpoly" : "diagonal_plus_galerkin");
    else if(E->control.ala_shallow_patch_preconditioner)
        preconditioner_mode="mpi_overlap_schwarz";
    else if(E->control.ala_radial_line_preconditioner)
        preconditioner_mode="radial_line";
    else
        preconditioner_mode="diagonal";
    if(E->parallel.me == 0) {
        fprintf(E->fp,
                "ALA pressure preconditioner = %s outer_solver=%s mode=%s "
                "BPI_range=(%e,%e) BPI_weight=%e invalid=%d restart_interval=%d "
                "beta_source=%s beta_causal_diagnostics=%s "
                "gamma=%e global_coarse=%s global_basis=%d "
                "global_weight=%e geneo=%s geneo_basis_type=%s "
                "geneo_basis=%d geneo_weight=%e "
                "geneo_rank_group=%dx%d cache_scope=%s\n",
                E->control.precondition ? "on" : "off",
                E->control.ala_outer_solver,
                preconditioner_mode,
                global_bpi_min, global_bpi_max,
                E->control.ala_pressure_bpi_weight,global_invalid_bpi,
                E->control.ala_pcg_restart_interval,
                E->control.ala_beta_element_source,
                E->control.ala_beta_causal_diagnostics ? "on" : "off",
                E->control.ala_augmented_lagrangian_gamma,
                E->control.ala_global_coarse_preconditioner ? "on" : "off",
                preconditioner_cache.global_basis_count,
                E->control.ala_global_coarse_weight,
                E->control.ala_geneo_preconditioner ? "on" : "off",
                E->control.ala_geneo_basis_type,
                preconditioner_cache.geneo_basis_count,
                E->control.ala_geneo_weight,
                E->control.ala_geneo_rank_group_x,
                E->control.ala_geneo_rank_group_y,
                coupled_self_contained_preconditioner
                    ? "self_contained_skip_legacy_pressure_cache"
                    : "legacy_pressure_cache_active");
        fprintf(stderr,
                "ALA pressure preconditioner = %s outer_solver=%s mode=%s "
                "BPI_range=(%e,%e) BPI_weight=%e invalid=%d restart_interval=%d "
                "beta_source=%s beta_causal_diagnostics=%s "
                "gamma=%e global_coarse=%s global_basis=%d "
                "global_weight=%e geneo=%s geneo_basis_type=%s "
                "geneo_basis=%d geneo_weight=%e "
                "geneo_rank_group=%dx%d cache_scope=%s\n",
                E->control.precondition ? "on" : "off",
                E->control.ala_outer_solver,
                preconditioner_mode,
                global_bpi_min, global_bpi_max,
                E->control.ala_pressure_bpi_weight,global_invalid_bpi,
                E->control.ala_pcg_restart_interval,
                E->control.ala_beta_element_source,
                E->control.ala_beta_causal_diagnostics ? "on" : "off",
                E->control.ala_augmented_lagrangian_gamma,
                E->control.ala_global_coarse_preconditioner ? "on" : "off",
                preconditioner_cache.global_basis_count,
                E->control.ala_global_coarse_weight,
                E->control.ala_geneo_preconditioner ? "on" : "off",
                E->control.ala_geneo_basis_type,
                preconditioner_cache.geneo_basis_count,
                E->control.ala_geneo_weight,
                E->control.ala_geneo_rank_group_x,
                E->control.ala_geneo_rank_group_y,
                coupled_self_contained_preconditioner
                    ? "self_contained_skip_legacy_pressure_cache"
                    : "legacy_pressure_cache_active");
        if(E->control.ala_feasibility_audit) {
            fprintf(E->fp,
                    "ALA_FEASIBILITY_AUDIT enabled inner_rel=%e "
                    "pressure_iterations=%d restart_interval=%d "
                    "window=%d minimum_window_reduction=%e "
                    "preconditioner=%s\n",
                    E->control.ala_inner_accuracy_max,
                    *steps_max,E->control.ala_pcg_restart_interval,
                    E->control.ala_feasibility_window,
                    E->control.ala_feasibility_min_reduction,
                    preconditioner_mode);
            fprintf(stderr,
                    "ALA_FEASIBILITY_AUDIT enabled inner_rel=%e "
                    "pressure_iterations=%d restart_interval=%d "
                    "window=%d minimum_window_reduction=%e\n",
                    E->control.ala_inner_accuracy_max,
                    *steps_max,E->control.ala_pcg_restart_interval,
                    E->control.ala_feasibility_window,
                    E->control.ala_feasibility_min_reduction);
        }
        if(E->control.ala_pressure_multigrid) {
            fprintf(E->fp,"ALA pressure V-cycle levels=%d range=%d:%d "
                    "pre_smooth=%d post_smooth=%d coarse_iterations=%d "
                    "damping=%e blend_weight=%e "
                    "level_operator=G*diag(Kgamma)^-1*Gt\n",
                    lev-E->control.ala_pressure_multigrid_min_level+1,
                    E->control.ala_pressure_multigrid_min_level,lev,
                    E->control.ala_pressure_multigrid_pre_smooth,
                    E->control.ala_pressure_multigrid_post_smooth,
                    E->control.ala_pressure_multigrid_coarse_iterations,
                    E->control.ala_pressure_multigrid_damping,
                    E->control.ala_pressure_multigrid_weight);
            fprintf(stderr,"ALA pressure V-cycle levels=%d range=%d:%d "
                    "pre_smooth=%d post_smooth=%d coarse_iterations=%d "
                    "damping=%e blend_weight=%e "
                    "level_operator=G*diag(Kgamma)^-1*Gt\n",
                    lev-E->control.ala_pressure_multigrid_min_level+1,
                    E->control.ala_pressure_multigrid_min_level,lev,
                    E->control.ala_pressure_multigrid_pre_smooth,
                    E->control.ala_pressure_multigrid_post_smooth,
                    E->control.ala_pressure_multigrid_coarse_iterations,
                    E->control.ala_pressure_multigrid_damping,
                    E->control.ala_pressure_multigrid_weight);
        }
        if(E->control.ala_two_level_preconditioner) {
            fprintf(E->fp,
                    "ALA two-level pressure correction offset=%d level=%d "
                    "coarse_solver=%s coarse_iterations=%d damping=%e "
                    "eigen_interval=(%e,%e) coarse_weight=%e "
                    "velocity_solver=%s velocity_iterations=%d "
                    "velocity_eigen_interval=(%e,%e) "
                    "invalid_diagonals=%d "
                    "operator=Sc=Pt*((D+C)Kapprox^-1(D+C)^T)_fine*P "
                    "velocity_inverse=%s "
                    "transfer=Pt/P matrix_free_galerkin=on\n",
                    E->control.ala_two_level_offset,coarse_lev,
                    E->control.ala_two_level_coarse_solver,
                    E->control.ala_two_level_coarse_iterations,
                    E->control.ala_two_level_coarse_damping,
                    preconditioner_cache.coarse_eigenvalue_min,
                    preconditioner_cache.coarse_eigenvalue_max,
                    E->control.ala_two_level_coarse_weight,
                    E->control.ala_two_level_velocity_solver,
                    E->control.ala_two_level_velocity_iterations,
                    preconditioner_cache.velocity_eigenvalue_min,
                    preconditioner_cache.velocity_eigenvalue_max,
                    global_invalid_velocity_bi,
                    strcmp(E->control.ala_two_level_velocity_solver,
                           "chebyshev")==0
                        ? "fixed_chebyshev_polynomial" : "diagonal");
            fprintf(stderr,
                    "ALA two-level pressure correction offset=%d level=%d "
                    "coarse_solver=%s coarse_iterations=%d damping=%e "
                    "eigen_interval=(%e,%e) coarse_weight=%e "
                    "velocity_solver=%s velocity_iterations=%d "
                    "velocity_eigen_interval=(%e,%e) "
                    "invalid_diagonals=%d "
                    "operator=Sc=Pt*((D+C)Kapprox^-1(D+C)^T)_fine*P "
                    "velocity_inverse=%s "
                    "transfer=Pt/P matrix_free_galerkin=on\n",
                    E->control.ala_two_level_offset,coarse_lev,
                    E->control.ala_two_level_coarse_solver,
                    E->control.ala_two_level_coarse_iterations,
                    E->control.ala_two_level_coarse_damping,
                    preconditioner_cache.coarse_eigenvalue_min,
                    preconditioner_cache.coarse_eigenvalue_max,
                    E->control.ala_two_level_coarse_weight,
                    E->control.ala_two_level_velocity_solver,
                    E->control.ala_two_level_velocity_iterations,
                    preconditioner_cache.velocity_eigenvalue_min,
                    preconditioner_cache.velocity_eigenvalue_max,
                    global_invalid_velocity_bi,
                    strcmp(E->control.ala_two_level_velocity_solver,
                           "chebyshev")==0
                        ? "fixed_chebyshev_polynomial" : "diagonal");
        }
    }
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(j=0; j<neq; j++)
            F[m][j] = FF[m][j];

    if(strcmp(E->control.ala_outer_solver,"fgmres")==0 ||
       strcmp(E->control.ala_outer_solver,"coupled_fgmres")==0) {
        if(strcmp(E->control.ala_outer_solver,"coupled_fgmres")==0) {
            /* A bounded-accuracy K_gamma momentum equilibration is useful on
             * both the cold solve and evolving warm starts.  It preserves the
             * monolithic coupled acceptance gate while avoiding many costly
             * block iterations whose only role is removing momentum error.
             * Zero retains the direct-block Stage 9d rollback. */
            if(E->control.ala_coupled_initial_velocity_relative_tolerance
               >0.0) {
                if(E->parallel.me==0) {
                    fprintf(E->fp,"ALA COUPLED MOMENTUM PREBALANCE "
                            "solution_cycle=%d relative_tolerance=%e\n",
                            E->monitor.solution_cycles,
                            E->control.
                              ala_coupled_initial_velocity_relative_tolerance);
                    fprintf(stderr,"ALA COUPLED MOMENTUM PREBALANCE "
                            "solution_cycle=%d relative_tolerance=%e\n",
                            E->monitor.solution_cycles,
                            E->control.
                              ala_coupled_initial_velocity_relative_tolerance);
                    fflush(E->fp);
                    fflush(stderr);
                }
                initial_vel_residual(
                    E,V,P,F,
                    E->control.
                      ala_coupled_initial_velocity_relative_tolerance);
            }
            else if(E->parallel.me==0) {
                fprintf(E->fp,"ALA COUPLED FGMRES direct block startup "
                        "initial_velocity_solve=skipped\n");
                fprintf(stderr,"ALA COUPLED FGMRES direct block startup "
                        "initial_velocity_solve=skipped\n");
                fflush(E->fp);
                fflush(stderr);
            }
            residual=solve_ala_coupled_fgmres_core(
                E,V,P,steps_max,lev,&preconditioner_cache,
                preconditioner_work);
        }
        else {
            /* The pressure-only Stage 7 path still requires the historical
             * momentum-consistent velocity before forming its Schur RHS. */
            initial_vel_residual(E,V,P,F,imp);
            residual=solve_ala_fgmres_core(
                E,V,P,steps_max,lev,&preconditioner_cache,r,
                explicit_r,div_u,preconditioner_work);
        }
        if(!coupled_self_contained_preconditioner &&
           E->control.ala_pressure_multigrid && E->parallel.me==0) {
            fprintf(E->fp,"ALA_PRESSURE_MG_COST "
                    "operator_applications=%d levels=%d\n",
                    preconditioner_cache.pressure_mg_operator_applications,
                    lev-preconditioner_cache.pressure_mg_min_level+1);
            fprintf(stderr,"ALA_PRESSURE_MG_COST "
                    "operator_applications=%d levels=%d\n",
                    preconditioner_cache.pressure_mg_operator_applications,
                    lev-preconditioner_cache.pressure_mg_min_level+1);
        }
        for(m=1;m<=E->sphere.caps_per_proc;m++) {
            free((void *)F[m]);
            if(!coupled_self_contained_preconditioner) {
                free((void *)r[m]);
                free((void *)z[m]);
                free((void *)p[m]);
                free((void *)q[m]);
                free((void *)explicit_r[m]);
                free((void *)div_u[m]);
            }
            free((void *)preconditioner_work[m]);
        }
        if(!coupled_self_contained_preconditioner &&
           (E->control.ala_shallow_patch_preconditioner ||
            E->control.ala_two_level_preconditioner ||
            E->control.ala_geneo_preconditioner ||
            E->control.ala_pressure_multigrid))
            free_ala_pressure_preconditioner_cache(E,&preconditioner_cache);
        return(residual);
    }

    /* FF contains the current -C^T*P forcing.  Pressure increments below
     * apply the complete strict-ALA transpose explicitly. */
    initial_vel_residual(E, V, P, F, imp);

    assemble_div_rho_u(E, V, r, lev);
    residual = incompressibility_residual(E, V, r);
    initial_rnorm = sqrt(global_pdot(E, r, r, lev));
    relative_residual = (initial_rnorm > 0.0) ? 1.0 : 0.0;
    strict_ala_continuity_metrics(E, V, r, div_u, lev,
                                  &initial_mass_norm, &cancellation_l2);
    initial_term_strength = strict_ala_continuity_term_strength(
        E,r,div_u,lev);
    mass_norm = initial_mass_norm;
    mass_relative_residual = (initial_mass_norm > 0.0) ? 1.0 : 0.0;

    sq_vdotv = sqrt(E->monitor.vdotv);
    dpressure = 1.0;
    dvelocity = 1.0;
    if(E->control.print_convergence && E->parallel.me == 0) {
        print_convergence_progress(E, count, time0, sq_vdotv,
                                   dvelocity, dpressure);
        fprintf(E->fp,
                "ALA PCG continuity residuals: cancellation=%e "
                "mass_norm=%e mass_relative=%e algebraic_relative=%e "
                "Q=%e Q_relative=%e "
                "term_strength=%e term_strength_relative=%e\n",
                cancellation_l2, mass_norm, mass_relative_residual,
                relative_residual, initial_term_strength, 1.0,
                initial_term_strength, 1.0);
    }
    strict_ala_depth_diagnostics(E, r, div_u, lev, count);
    strict_ala_beta_causal_diagnostics(
        E,V,r,preconditioner_work,div_u,lev,count);
    strict_ala_coarse_residual_diagnostics(E, r, lev, count);

    /* A fixed inner target keeps the approximate Schur operator as
     * stationary as the RHS-relative multigrid stopping rule permits. */
    inner_relative_accuracy = E->control.ala_inner_accuracy_max;
    rho_old = 0.0;
    restart_search = 1;
    hybrid_consecutive_count = 0;
    hybrid_converged = 0;
    nonpositive_curvature_count = 0;
    min_curvature = 1.0e300;
    audit_milestones[0] = 0.5;
    audit_milestones[1] = 0.2;
    audit_milestones[2] = 0.1;
    audit_milestones[3] = 0.05;
    audit_milestones[4] = 0.02;
    audit_milestones[5] = 0.01;
    audit_milestone = 0;
    audit_best_cancellation = cancellation_l2;
    audit_best_iteration = 0;
    audit_window_cancellation = cancellation_l2;
    audit_window_iteration = 0;
    audit_window_reduction = 0.0;
    audit_stagnated = 0;

    while(count < *steps_max &&
          (cancellation_l2 >= E->control.tole_comp ||
           (E->control.ala_hybrid_convergence && !hybrid_converged))) {

        apply_ala_pressure_preconditioner(
            E,r,z,preconditioner_work,lev,count,&preconditioner_cache);

        rho = global_pdot(E, r, z, lev);
        if(!isfinite(rho) || rho <= 1.0e-300) {
            if(E->parallel.me == 0) {
                fprintf(E->fp,
                        "ALA PCG breakdown: invalid preconditioned residual "
                        "<r,z>=%e\n", rho);
                fprintf(stderr,
                        "ALA PCG breakdown: invalid preconditioned residual "
                        "<r,z>=%e\n", rho);
            }
            parallel_process_termination();
        }

        if(restart_search) {
            for(m=1; m<=E->sphere.caps_per_proc; m++)
                for(j=1; j<=npno; j++)
                    p[m][j] = z[m][j];
            restart_search = 0;
        }
        else {
            beta = rho / rho_old;
            if(!isfinite(beta)) {
                if(E->parallel.me == 0) {
                    fprintf(E->fp, "ALA PCG breakdown: non-finite beta\n");
                    fprintf(stderr, "ALA PCG breakdown: non-finite beta\n");
                }
                parallel_process_termination();
            }
            for(m=1; m<=E->sphere.caps_per_proc; m++)
                for(j=1; j<=npno; j++)
                    p[m][j] = z[m][j] + beta * p[m][j];
        }

        assemble_grad_rho_p(E, p, F, lev);
        inner_accuracy = strict_ala_inner_accuracy(
            E, F, lev, inner_relative_accuracy);
        valid = solve_del2_u(E, E->u1, F, inner_accuracy, lev);
        if(!valid) {
            if(E->parallel.me == 0) {
                fputs("ALA PCG inner velocity solve failed\n", stderr);
                fputs("ALA PCG inner velocity solve failed\n", E->fp);
            }
            parallel_process_termination();
        }
        strip_bcs_from_residual(E, E->u1, lev);
        assemble_div_rho_u(E, E->u1, q, lev);

        curvature = global_pdot(E, p, q, lev);
        if(!isfinite(curvature) || curvature <= 1.0e-300) {
            nonpositive_curvature_count++;
            if(E->parallel.me == 0) {
                fprintf(E->fp,
                        "ALA PCG breakdown: nonpositive Schur curvature=%e "
                        "at iteration %d\n", curvature, count);
                fprintf(stderr,
                        "ALA PCG breakdown: nonpositive Schur curvature=%e "
                        "at iteration %d\n", curvature, count);
            }
            parallel_process_termination();
        }
        min_curvature = min(min_curvature, curvature);
        alpha = rho / curvature;
        if(!isfinite(alpha)) {
            if(E->parallel.me == 0) {
                fprintf(E->fp, "ALA PCG breakdown: non-finite alpha\n");
                fprintf(stderr, "ALA PCG breakdown: non-finite alpha\n");
            }
            parallel_process_termination();
        }

        for(m=1; m<=E->sphere.caps_per_proc; m++) {
            for(j=1; j<=npno; j++) {
                P[m][j] += alpha * p[m][j];
                r[m][j] -= alpha * q[m][j];
            }
            for(j=0; j<neq; j++)
                V[m][j] -= alpha * E->u1[m][j];
        }

        recursive_rnorm = sqrt(global_pdot(E, r, r, lev));
        recursive_relative_residual = (initial_rnorm > 0.0)
            ? recursive_rnorm / initial_rnorm : 0.0;

        assemble_div_rho_u(E, V, explicit_r, lev);
        residual = incompressibility_residual(E, V, explicit_r);
        rnorm = sqrt(global_pdot(E, explicit_r, explicit_r, lev));
        relative_residual = (initial_rnorm > 0.0)
            ? rnorm / initial_rnorm : 0.0;
        strict_ala_continuity_metrics(E, V, explicit_r, div_u, lev,
                                      &mass_norm, &cancellation_l2);
        term_strength = strict_ala_continuity_term_strength(
            E,explicit_r,div_u,lev);
        mass_relative_residual = (initial_mass_norm > 0.0)
            ? mass_norm / initial_mass_norm : 0.0;
        drift_ratio = (relative_residual > 0.0)
            ? recursive_relative_residual / relative_residual : 1.0;
        residual = cancellation_l2;

        dpressure = fabs(alpha) * sqrt(
            global_pdot(E, p, p, lev) /
            (1.0e-32 + global_pdot(E, P, P, lev)));
        dvelocity = fabs(alpha) * sqrt(
            global_vdot(E, E->u1, E->u1, lev) /
            (1.0e-32 + E->monitor.vdotv));

        count++;
        if(cancellation_l2 < audit_best_cancellation) {
            audit_best_cancellation = cancellation_l2;
            audit_best_iteration = count;
        }
        if(E->control.ala_feasibility_audit) {
            while(audit_milestone < 6 &&
                  cancellation_l2 <= audit_milestones[audit_milestone]) {
                if(E->parallel.me == 0)
                    fprintf(E->fp,
                            "ALA_FEASIBILITY_MILESTONE cancellation=%e "
                            "iteration=%d elapsed_seconds=%e "
                            "mass_relative=%e algebraic_relative=%e\n",
                            audit_milestones[audit_milestone],count,
                            CPU_time0()-time0,mass_relative_residual,
                            relative_residual);
                audit_milestone++;
            }
            if(count-audit_window_iteration >=
               E->control.ala_feasibility_window) {
                audit_window_reduction = 1.0
                    - cancellation_l2
                    / max(audit_window_cancellation,1.0e-300);
                audit_stagnated =
                    audit_window_reduction <
                    E->control.ala_feasibility_min_reduction;
                if(E->parallel.me == 0) {
                    fprintf(E->fp,
                            "ALA_FEASIBILITY_WINDOW from=%d to=%d "
                            "cancellation_start=%e cancellation_end=%e "
                            "relative_reduction=%e best=%e "
                            "best_iteration=%d status=%s\n",
                            audit_window_iteration,count,
                            audit_window_cancellation,cancellation_l2,
                            audit_window_reduction,
                            audit_best_cancellation,audit_best_iteration,
                            audit_stagnated ? "stagnated" : "progressing");
                    fflush(E->fp);
                }
                audit_window_cancellation = cancellation_l2;
                audit_window_iteration = count;
            }
        }
        if(E->control.ala_hybrid_convergence) {
            if(E->monitor.incompressibility <
                   E->control.ala_div_v_tolerance &&
               dvelocity < E->control.ala_update_tolerance &&
               dpressure < E->control.ala_update_tolerance)
                hybrid_consecutive_count++;
            else
                hybrid_consecutive_count = 0;
            if(hybrid_consecutive_count >= E->control.ala_consecutive_steps)
                hybrid_converged = 1;
        }

        sq_vdotv = sqrt(E->monitor.vdotv);
        if(E->control.print_convergence && E->parallel.me == 0) {
            print_convergence_progress(E, count, time0, sq_vdotv,
                                       dvelocity, dpressure);
            fprintf(E->fp,
                    "ALA PCG continuity residuals: cancellation=%e "
                    "mass_norm=%e mass_relative=%e algebraic_relative=%e "
                    "recursive_algebraic=%e drift=%e inner_rel=%e "
                    "curvature=%e Q=%e Q_relative=%e term_strength=%e "
                    "term_strength_relative=%e\n",
                    cancellation_l2, mass_norm, mass_relative_residual,
                    relative_residual, recursive_relative_residual,
                    drift_ratio, inner_relative_accuracy, curvature,
                    term_strength,
                    term_strength/max(initial_term_strength,1.0e-300),
                    term_strength,
                    term_strength/max(initial_term_strength,1.0e-300));
            if(E->control.ala_hybrid_convergence)
                fprintf(E->fp,
                        "ALA hybrid convergence streak = %d/%d "
                        "limits: div/v<%e updates<%e\n",
                        hybrid_consecutive_count,
                        E->control.ala_consecutive_steps,
                        E->control.ala_div_v_tolerance,
                        E->control.ala_update_tolerance);
        }
        strict_ala_depth_diagnostics(E, explicit_r, div_u, lev, count);
        strict_ala_beta_causal_diagnostics(
            E,V,explicit_r,preconditioner_work,div_u,lev,count);
        strict_ala_coarse_residual_diagnostics(E, explicit_r, lev, count);

        rho_old = rho;

        /* Explicit replacement protects a long inexact PCG run from residual
         * drift.  Replacement also restarts the conjugate direction. */
        if((count % E->control.ala_pcg_restart_interval) == 0 ||
           drift_ratio > 10.0 || drift_ratio < 0.1) {
            for(m=1; m<=E->sphere.caps_per_proc; m++)
                for(j=1; j<=npno; j++)
                    r[m][j] = explicit_r[m][j];
            restart_search = 1;
            if(E->parallel.me == 0)
                fprintf(E->fp,
                        "ALA PCG restarted from explicit residual at "
                        "iteration %d\n", count);
        }
    }

    strict_ala_momentum_residual_audit(
        E,V,P,F,E->u1,lev,&momentum_rms,&momentum_relative);

    if(E->parallel.me == 0) {
        fprintf(E->fp,
                "ALA PCG momentum audit final iteration=%d rms=%e "
                "relative=%e equation=unaugmented_Ku_plus_GtP_minus_f\n",
                count,momentum_rms,momentum_relative);
        fprintf(stderr,
                "ALA PCG momentum audit final iteration=%d rms=%e "
                "relative=%e equation=unaugmented_Ku_plus_GtP_minus_f\n",
                count,momentum_rms,momentum_relative);
        fprintf(E->fp,
                "ALA PCG operator audit: minimum_curvature=%e "
                "nonpositive_curvature=%d iterations=%d "
                "Schur_applications=%d pressure_velocity_solves=%d\n",
                min_curvature, nonpositive_curvature_count, count, count,
                count);
        fprintf(stderr,
                "ALA PCG operator audit: minimum_curvature=%e "
                "nonpositive_curvature=%d iterations=%d "
                "Schur_applications=%d pressure_velocity_solves=%d\n",
                min_curvature, nonpositive_curvature_count, count, count,
                count);
        if(E->control.ala_two_level_preconditioner) {
            galerkin_diagnostic_applications=(count>0)
                ? 1+(count-1)/E->control.ala_coarse_residual_interval : 0;
            galerkin_applications=count
                *(E->control.ala_two_level_coarse_iterations-1)
                +galerkin_diagnostic_applications;
            fprintf(E->fp,
                    "ALA_TWO_LEVEL_COST pressure_iterations=%d "
                    "galerkin_applications=%d diagnostic_applications=%d "
                    "coarse_polynomial_steps=%d "
                    "velocity_polynomial_steps=%d\n",
                    count,galerkin_applications,
                    galerkin_diagnostic_applications,
                    E->control.ala_two_level_coarse_iterations,
                    E->control.ala_two_level_velocity_iterations);
            fprintf(stderr,
                    "ALA_TWO_LEVEL_COST pressure_iterations=%d "
                    "galerkin_applications=%d diagnostic_applications=%d "
                    "coarse_polynomial_steps=%d "
                    "velocity_polynomial_steps=%d\n",
                    count,galerkin_applications,
                    galerkin_diagnostic_applications,
                    E->control.ala_two_level_coarse_iterations,
                    E->control.ala_two_level_velocity_iterations);
        }
        if(E->control.ala_geneo_preconditioner) {
            fprintf(E->fp,"ALA_GENEO_COST setup_galerkin_applications=%d "
                    "pressure_iterations=%d dense_coarse_solves=%d "
                    "global_modes=%d\n",
                    preconditioner_cache.geneo_schur_applications,count,count,
                    preconditioner_cache.geneo_basis_count);
            fprintf(stderr,"ALA_GENEO_COST setup_galerkin_applications=%d "
                    "pressure_iterations=%d dense_coarse_solves=%d "
                    "global_modes=%d\n",
                    preconditioner_cache.geneo_schur_applications,count,count,
                    preconditioner_cache.geneo_basis_count);
        }
        if(E->control.ala_feasibility_audit) {
            fprintf(E->fp,
                    "ALA_FEASIBILITY_SUMMARY status=%s final=%e best=%e "
                    "best_iteration=%d target=%e iterations=%d "
                    "last_complete_window_reduction=%e "
                    "mass_norm=%e mass_relative=%e Q=%e Q_relative=%e "
                    "momentum_relative=%e\n",
                    cancellation_l2 < E->control.tole_comp
                        ? "discrete_target_reached"
                        : (audit_stagnated
                           ? "preconditioned_iteration_stagnated"
                           : "iteration_budget_exhausted_while_progressing"),
                    cancellation_l2,audit_best_cancellation,
                    audit_best_iteration,E->control.tole_comp,count,
                    audit_window_reduction,mass_norm,mass_relative_residual,
                    term_strength,
                    term_strength/max(initial_term_strength,1.0e-300),
                    momentum_relative);
            fprintf(stderr,
                    "ALA_FEASIBILITY_SUMMARY status=%s final=%e best=%e "
                    "best_iteration=%d target=%e iterations=%d\n",
                    cancellation_l2 < E->control.tole_comp
                        ? "discrete_target_reached"
                        : (audit_stagnated
                           ? "preconditioned_iteration_stagnated"
                           : "iteration_budget_exhausted_while_progressing"),
                    cancellation_l2,audit_best_cancellation,
                    audit_best_iteration,E->control.tole_comp,count);
        }
    }

    if(E->parallel.me == 0 &&
       cancellation_l2 < E->control.tole_comp &&
       (!E->control.ala_hybrid_convergence || hybrid_converged)) {
        fprintf(E->fp,
                "Strict ALA PCG converged by physical continuity: "
                "cancellation=%e tolerance=%e mass_relative=%e "
                "algebraic_relative=%e iterations=%d\n",
                cancellation_l2, E->control.tole_comp,
                mass_relative_residual, relative_residual, count);
        fflush(E->fp);
    }

    if(cancellation_l2 >= E->control.tole_comp ||
       (E->control.ala_hybrid_convergence && !hybrid_converged)) {
        if(E->parallel.me == 0) {
            fprintf(E->fp,
                    "Strict ALA PCG failed physical continuity: "
                    "cancellation=%e tolerance=%e mass_relative=%e "
                    "algebraic_relative=%e hybrid_streak=%d/%d "
                    "iterations=%d\n",
                    cancellation_l2, E->control.tole_comp,
                    mass_relative_residual, relative_residual,
                    hybrid_consecutive_count,
                    E->control.ala_consecutive_steps, count);
            fprintf(stderr,
                    "Strict ALA PCG failed physical continuity: "
                    "cancellation=%e tolerance=%e mass_relative=%e "
                    "algebraic_relative=%e hybrid_streak=%d/%d "
                    "iterations=%d\n",
                    cancellation_l2, E->control.tole_comp,
                    mass_relative_residual, relative_residual,
                    hybrid_consecutive_count,
                    E->control.ala_consecutive_steps, count);
            fflush(E->fp);
        }
        parallel_process_termination();
    }

    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        free((void *)F[m]);
        free((void *)r[m]);
        free((void *)z[m]);
        free((void *)p[m]);
        free((void *)q[m]);
        free((void *)explicit_r[m]);
        free((void *)div_u[m]);
        free((void *)preconditioner_work[m]);
    }
    if(E->control.ala_shallow_patch_preconditioner ||
       E->control.ala_two_level_preconditioner ||
       E->control.ala_geneo_preconditioner ||
       E->control.ala_pressure_multigrid)
        free_ala_pressure_preconditioner_cache(E,&preconditioner_cache);

    *steps_max = count;
    return(residual);
}


/* Solve compressible Stokes flow using
 * bi-conjugate gradient stablized (BiCG-stab) iterations
 */

static float solve_Ahat_p_fhat_BiCG(struct All_variables *E,
                                    double **V, double **P, double **FF,
                                    double imp, int *steps_max)
{
    void assemble_div_rho_u();
    void assemble_del2_u();
    void assemble_grad_p();
    void assemble_grad_rho_p();
    void strip_bcs_from_residual();
    int  solve_del2_u();
    void parallel_process_termination();

    double global_vdot(), global_pdot();
    double CPU_time0();

    int gnpno, gneq;
    int npno, neq;
    int m, j, count, lev;
    int valid;
    int restart_search;
    int hybrid_consecutive_count, hybrid_converged;

    double alpha, beta, omega,sq_vdotv;
    double r0dotrt, r1dotrt;
    double residual, dpressure, dvelocity;
    double initial_rnorm, relative_residual, rnorm;
    double initial_mass_norm, mass_norm, mass_relative_residual;
    double cancellation_l2;
    double recursive_rnorm, recursive_relative_residual, drift_ratio;
    double denominator, numerator, inner_accuracy;
    double inner_relative_accuracy, previous_inner_relative_accuracy;
    double symmetry_left, symmetry_right, symmetry_defect;
    double max_symmetry_defect, curvature_p, curvature_s;
    int symmetry_sample_count, nonpositive_curvature_count;

    double *F[NCS];
    double *r1[NCS], *r2[NCS], *pt[NCS], *p1[NCS], *p2[NCS];
    double *rt[NCS], *v0[NCS], *s0[NCS], *st[NCS], *t0[NCS];
    double *u0[NCS], *div_u[NCS];
    double *shuffle[NCS];

    double time0, v_res;

    gnpno = E->mesh.npno;
    gneq = E->mesh.neq;
    npno = E->lmesh.npno;
    neq = E->lmesh.neq;
    lev = E->mesh.levmax;

    for (m=1; m<=E->sphere.caps_per_proc; m++)   {
        F[m] = (double *)malloc(neq*sizeof(double));
        r1[m] = (double *)malloc((npno+1)*sizeof(double));
        r2[m] = (double *)malloc((npno+1)*sizeof(double));
        pt[m] = (double *)malloc((npno+1)*sizeof(double));
        p1[m] = (double *)malloc((npno+1)*sizeof(double));
        p2[m] = (double *)malloc((npno+1)*sizeof(double));
        rt[m] = (double *)malloc((npno+1)*sizeof(double));
        v0[m] = (double *)malloc((npno+1)*sizeof(double));
        s0[m] = (double *)malloc((npno+1)*sizeof(double));
        st[m] = (double *)malloc((npno+1)*sizeof(double));
        t0[m] = (double *)malloc((npno+1)*sizeof(double));
        div_u[m] = (double *)malloc((npno+1)*sizeof(double));

        u0[m] = (double *)malloc(neq*sizeof(double));
    }

    time0 = CPU_time0();
    count = 0;

    /* copy the original force vector since we need to keep it intact
       between iterations */
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(j=0;j<neq;j++)
            F[m][j] = FF[m][j];


    /* FF already contains the current -C^T*P force contribution.  The
     * pressure increment operator below must therefore apply both the
     * standard gradient and the strict-ALA C^T response. */
    /* calculate the initial velocity residual */
    v_res = initial_vel_residual(E, V, P, F, imp);


    /* initial residual r1 = div(rho_ref*V) */
    assemble_div_rho_u(E, V, r1, lev);
    residual = incompressibility_residual(E, V, r1);
    initial_rnorm = sqrt(global_pdot(E, r1, r1, lev));
    relative_residual = (initial_rnorm > 0.0) ? 1.0 : 0.0;
    if(E->control.ala_pressure_buoyancy) {
        strict_ala_continuity_metrics(E, V, r1, div_u, lev,
                                      &initial_mass_norm, &cancellation_l2);
        mass_norm = initial_mass_norm;
        mass_relative_residual = (initial_mass_norm > 0.0) ? 1.0 : 0.0;
    }
    else {
        initial_mass_norm = mass_norm = initial_rnorm;
        mass_relative_residual = relative_residual;
        cancellation_l2 = relative_residual;
    }


    /* initial conjugate residual rt = r1 */
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(j=1; j<=npno; j++)
            rt[m][j] = r1[m][j];


    sq_vdotv = sqrt(E->monitor.vdotv);

    /* pressure and velocity corrections */
    dpressure = 1.0;
    dvelocity = 1.0;


    if (E->control.print_convergence && E->parallel.me==0)  {
        print_convergence_progress(E, count, time0, sq_vdotv,
                                   dvelocity, dpressure);
        if(E->control.ala_pressure_buoyancy)
            fprintf(E->fp,
                    "ALA continuity residuals: cancellation=%e "
                    "mass_relative=%e algebraic_relative=%e\n",
                    cancellation_l2, mass_relative_residual,
                    relative_residual);
    }
    if(E->control.ala_pressure_buoyancy)
        strict_ala_depth_diagnostics(E, r1, div_u, lev, count);


    valid = 1;
    r0dotrt = alpha = omega = 0;
    restart_search = 0;
    previous_inner_relative_accuracy = -1.0;
    max_symmetry_defect = 0.0;
    symmetry_sample_count = 0;
    nonpositive_curvature_count = 0;
    hybrid_consecutive_count = 0;
    hybrid_converged = 0;

    while( (count < *steps_max) &&
           ((E->control.ala_pressure_buoyancy &&
             (cancellation_l2 >= E->control.tole_comp ||
              (E->control.ala_hybrid_convergence &&
               !hybrid_converged))) ||
            (!E->control.ala_pressure_buoyancy &&
             E->monitor.incompressibility >= E->control.tole_comp)) &&
           (E->control.ala_pressure_buoyancy ||
            ((dpressure >= imp) && (dvelocity >= imp))) )  {


        inner_relative_accuracy = strict_ala_inner_relative_accuracy(
            E, mass_relative_residual);
        if(previous_inner_relative_accuracy > 0.0 &&
           fabs(inner_relative_accuracy - previous_inner_relative_accuracy)
             > 1.0e-12 * previous_inner_relative_accuracy) {
            /* A changed inner tolerance changes the approximate K inverse.
             * Re-anchor BiCGStab before using the new inexact Schur operator. */
            for(m=1; m<=E->sphere.caps_per_proc; m++)
                for(j=1; j<=npno; j++) {
                    r1[m][j] = t0[m][j];
                    rt[m][j] = t0[m][j];
                }
            restart_search = 1;
            r0dotrt = alpha = omega = 0.0;
            if(E->parallel.me == 0) {
                fprintf(E->fp,
                        "ALA inner relative accuracy changed %.3e -> %.3e; "
                        "restarting BiCGStab at iteration %d\n",
                        previous_inner_relative_accuracy,
                        inner_relative_accuracy, count);
                fprintf(stderr,
                        "ALA inner relative accuracy changed %.3e -> %.3e; "
                        "restarting BiCGStab at iteration %d\n",
                        previous_inner_relative_accuracy,
                        inner_relative_accuracy, count);
            }
        }
        previous_inner_relative_accuracy = inner_relative_accuracy;


        /* r1dotrt = <r1, rt> */
        r1dotrt = global_pdot(E, r1, rt, lev);
        if(!isfinite(r1dotrt) || fabs(r1dotrt) <= 1.0e-300) {
            fprintf(E->fp, "BiCGStab breakdown: invalid <r,rt>\n");
            fprintf(stderr, "BiCGStab breakdown: invalid <r,rt>\n");
            parallel_process_termination();
        }


        /* update search direction */
        if(count == 0 || restart_search) {
            for (m=1; m<=E->sphere.caps_per_proc; m++)
                for(j=1; j<=npno; j++)
                    p2[m][j] = r1[m][j];
            restart_search = 0;
        }
        else {
            /* p2 = r1 + <r1,rt>/<r0,rt> * alpha/omega * (p1 - omega*v0) */
            if(!isfinite(r0dotrt) || fabs(r0dotrt) <= 1.0e-300 ||
               !isfinite(omega) || fabs(omega) <= 1.0e-300) {
                fprintf(E->fp,
                        "BiCGStab breakdown: invalid beta denominator\n");
                fprintf(stderr,
                        "BiCGStab breakdown: invalid beta denominator\n");
                parallel_process_termination();
            }
            beta = (r1dotrt / r0dotrt) * (alpha / omega);
            if(!isfinite(beta)) {
                fprintf(E->fp, "BiCGStab breakdown: non-finite beta\n");
                fprintf(stderr, "BiCGStab breakdown: non-finite beta\n");
                parallel_process_termination();
            }
            for(m=1; m<=E->sphere.caps_per_proc; m++)
                for(j=1; j<=npno; j++)
                    p2[m][j] = r1[m][j] + beta
                        * (p1[m][j] - omega * v0[m][j]);
        }


        /* pressure Schur-diagonal preconditioner */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=1; j<=npno; j++)
                pt[m][j] = E->BPI[lev][m][j] * p2[m][j];


        /* solve K*u0 = (D+C)^T*pt for strict ALA */
        if(E->control.ala_pressure_buoyancy)
            assemble_grad_rho_p(E, pt, F, lev);
        else
            assemble_grad_p(E, pt, F, lev);
        inner_accuracy = E->control.ala_pressure_buoyancy
            ? strict_ala_inner_accuracy(E, F, lev,
                                        inner_relative_accuracy) : imp * v_res;
        valid = solve_del2_u(E, u0, F, inner_accuracy, lev);
        if(!valid && (E->parallel.me==0)) {
            fputs("Warning: solver not converging! 1\n", stderr);
            fputs("Warning: solver not converging! 1\n", E->fp);
        }
        strip_bcs_from_residual(E, u0, lev);


        /* v0 = div(rho_ref*u0) */
        assemble_div_rho_u(E, u0, v0, lev);


        /* alpha = r1dotrt / <rt, v0> */
        denominator = global_pdot(E, rt, v0, lev);
        if(!isfinite(denominator) || fabs(denominator) <= 1.0e-300) {
            fprintf(E->fp, "BiCGStab breakdown: invalid alpha denominator\n");
            fprintf(stderr, "BiCGStab breakdown: invalid alpha denominator\n");
            parallel_process_termination();
        }
        alpha = r1dotrt / denominator;
        if(!isfinite(alpha)) {
            fprintf(E->fp, "BiCGStab breakdown: non-finite alpha\n");
            fprintf(stderr, "BiCGStab breakdown: non-finite alpha\n");
            parallel_process_termination();
        }


        /* s0 = r1 - alpha * v0 */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=1; j<=npno; j++)
                s0[m][j] = r1[m][j] - alpha * v0[m][j];


        /* pressure Schur-diagonal preconditioner */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=1; j<=npno; j++)
                st[m][j] = E->BPI[lev][m][j] * s0[m][j];


        /* solve K*u1 = (D+C)^T*st for strict ALA */
        if(E->control.ala_pressure_buoyancy)
            assemble_grad_rho_p(E, st, F, lev);
        else
            assemble_grad_p(E, st, F, lev);
        inner_accuracy = E->control.ala_pressure_buoyancy
            ? strict_ala_inner_accuracy(E, F, lev,
                                        inner_relative_accuracy) : imp * v_res;
        valid = solve_del2_u(E, E->u1, F, inner_accuracy, lev);
        if(!valid && (E->parallel.me==0)) {
            fputs("Warning: solver not converging! 2\n", stderr);
            fputs("Warning: solver not converging! 2\n", E->fp);
        }
        strip_bcs_from_residual(E, E->u1, lev);


        /* t0 = div(rho_ref * u1) */
        assemble_div_rho_u(E, E->u1, t0, lev);

        if(E->control.ala_pressure_buoyancy &&
           E->control.ala_schur_symmetry_check) {
            /* v0=S*pt and t0=S*st are already available.  This checks the
             * actual inexact Schur applications without another K solve. */
            symmetry_left = global_pdot(E, pt, t0, lev);
            symmetry_right = global_pdot(E, st, v0, lev);
            symmetry_defect = fabs(symmetry_left - symmetry_right)
                / (fabs(symmetry_left) + fabs(symmetry_right) + 1.0e-300);
            curvature_p = global_pdot(E, pt, v0, lev);
            curvature_s = global_pdot(E, st, t0, lev);
            if(!isfinite(symmetry_defect) || !isfinite(curvature_p) ||
               !isfinite(curvature_s)) {
                symmetry_defect = 1.0;
                nonpositive_curvature_count++;
            }
            else if(curvature_p <= 0.0 || curvature_s <= 0.0)
                nonpositive_curvature_count++;
            max_symmetry_defect = max(max_symmetry_defect, symmetry_defect);
            symmetry_sample_count++;
        }


        /* omega = <t0, s0> / <t0, t0> */
        numerator = global_pdot(E, t0, s0, lev);
        denominator = global_pdot(E, t0, t0, lev);
        if(!isfinite(numerator) || !isfinite(denominator) ||
           fabs(denominator) <= 1.0e-300) {
            fprintf(E->fp, "BiCGStab breakdown: invalid omega denominator\n");
            fprintf(stderr, "BiCGStab breakdown: invalid omega denominator\n");
            parallel_process_termination();
        }
        omega = numerator / denominator;
        if(!isfinite(omega)) {
            fprintf(E->fp, "BiCGStab breakdown: non-finite omega\n");
            fprintf(stderr, "BiCGStab breakdown: non-finite omega\n");
            parallel_process_termination();
        }


        /* r2 = s0 - omega * t0 */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=1; j<=npno; j++)
                r2[m][j] = s0[m][j] - omega * t0[m][j];

        recursive_rnorm = sqrt(global_pdot(E, r2, r2, lev));
        recursive_relative_residual = (initial_rnorm > 0.0)
            ? recursive_rnorm / initial_rnorm : 0.0;


        /* P = P + alpha * pt + omega * st */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=1; j<=npno; j++)
                s0[m][j] = alpha * pt[m][j] + omega * st[m][j];

        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=1; j<=npno; j++)
                P[m][j] += s0[m][j];


        /* V = V - alpha * u0 - omega * u1 */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=0; j<neq; j++)
                F[m][j] = alpha * u0[m][j] + omega * E->u1[m][j];

        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=0; j<neq; j++)
                V[m][j] -= F[m][j];


        /* compute velocity and incompressibility residual */
        assemble_div_rho_u(E, V, t0, lev);
        residual = incompressibility_residual(E, V, t0);
        rnorm = sqrt(global_pdot(E, t0, t0, lev));
        relative_residual = (initial_rnorm > 0.0)
            ? rnorm / initial_rnorm : 0.0;
        if(E->control.ala_pressure_buoyancy) {
            strict_ala_continuity_metrics(E, V, t0, div_u, lev,
                                          &mass_norm, &cancellation_l2);
            mass_relative_residual = (initial_mass_norm > 0.0)
                ? mass_norm / initial_mass_norm : 0.0;
        }
        drift_ratio = (relative_residual > 0.0)
            ? recursive_relative_residual / relative_residual : 1.0;
        if(E->control.ala_pressure_buoyancy)
            residual = cancellation_l2;

        /* compute velocity and pressure corrections */
        dpressure = sqrt( global_pdot(E, s0, s0, lev)
                          / (1.0e-32 + global_pdot(E, P, P, lev)) );
        dvelocity = sqrt( global_vdot(E, F, F, lev)
                          / (1.0e-32 + E->monitor.vdotv) );


        count++;

        if(E->control.ala_pressure_buoyancy &&
           E->control.ala_hybrid_convergence) {
            if(E->monitor.incompressibility <
                   E->control.ala_div_v_tolerance &&
               dvelocity < E->control.ala_update_tolerance &&
               dpressure < E->control.ala_update_tolerance)
                hybrid_consecutive_count++;
            else
                hybrid_consecutive_count = 0;
            if(hybrid_consecutive_count >= E->control.ala_consecutive_steps)
                hybrid_converged = 1;
        }

	sq_vdotv = sqrt(E->monitor.vdotv);

        if(E->control.print_convergence && E->parallel.me==0) {
            print_convergence_progress(E, count, time0, sq_vdotv,
                                       dvelocity, dpressure);
            if(E->control.ala_pressure_buoyancy)
                fprintf(E->fp,
                        "ALA continuity residuals: cancellation=%e "
                        "mass_relative=%e algebraic_relative=%e "
                        "recursive_algebraic=%e drift=%e inner_rel=%e\n",
                        cancellation_l2, mass_relative_residual,
                        relative_residual, recursive_relative_residual,
                        drift_ratio, inner_relative_accuracy);
            if(E->control.ala_pressure_buoyancy &&
               E->control.ala_hybrid_convergence)
                fprintf(E->fp,
                        "ALA hybrid convergence streak = %d/%d "
                        "limits: div/v<%e updates<%e\n",
                        hybrid_consecutive_count,
                        E->control.ala_consecutive_steps,
                        E->control.ala_div_v_tolerance,
                        E->control.ala_update_tolerance);
            if(E->control.ala_pressure_buoyancy &&
               E->control.ala_schur_symmetry_check) {
                fprintf(E->fp,
                        "ALA Schur symmetry defect = %e curvature=(%e,%e)\n",
                        symmetry_defect, curvature_p, curvature_s);
                fprintf(stderr,
                        "ALA Schur symmetry defect = %e curvature=(%e,%e)\n",
                        symmetry_defect, curvature_p, curvature_s);
            }
        }
        if(E->control.ala_pressure_buoyancy)
            strict_ala_depth_diagnostics(E, t0, div_u, lev, count);

        /* Re-anchor the recurrence to the explicitly assembled B*u residual.
         * This limits finite-precision drift in long strict-ALA solves. */
        if(E->control.ala_pressure_buoyancy &&
           ((count % 20) == 0 || drift_ratio > 10.0 || drift_ratio < 0.1)) {
            for(m=1; m<=E->sphere.caps_per_proc; m++)
                for(j=1; j<=npno; j++) {
                    r2[m][j] = t0[m][j];
                    rt[m][j] = t0[m][j];
                }
            restart_search = 1;
            if(E->parallel.me == 0)
                fprintf(E->fp,
                        "ALA BiCGStab restarted from explicit residual at "
                        "iteration %d\n", count);
        }


        /* shift array pointers */
        for(m=1; m<=E->sphere.caps_per_proc; m++) {
            shuffle[m] = p1[m];
            p1[m] = p2[m];
            p2[m] = shuffle[m];

            shuffle[m] = r1[m];
            r1[m] = r2[m];
            r2[m] = shuffle[m];
        }

        /* shift <r0, rt> = <r1, rt> */
        r0dotrt = r1dotrt;

    } /* end loop for BiCGStab */

    if(E->control.ala_pressure_buoyancy &&
       E->control.ala_schur_symmetry_check && E->parallel.me == 0) {
        fprintf(E->fp,
                "ALA Schur CG suitability: max_symmetry_defect=%e "
                "tolerance=%e samples=%d nonpositive_curvature=%d result=%s\n",
                max_symmetry_defect,
                E->control.ala_schur_symmetry_tolerance,
                symmetry_sample_count,
                nonpositive_curvature_count,
                (symmetry_sample_count > 0 && max_symmetry_defect <=
                     E->control.ala_schur_symmetry_tolerance &&
                 nonpositive_curvature_count == 0) ? "PASS" : "FAIL");
        fprintf(stderr,
                "ALA Schur CG suitability: max_symmetry_defect=%e "
                "tolerance=%e samples=%d nonpositive_curvature=%d result=%s\n",
                max_symmetry_defect,
                E->control.ala_schur_symmetry_tolerance,
                symmetry_sample_count,
                nonpositive_curvature_count,
                (symmetry_sample_count > 0 && max_symmetry_defect <=
                     E->control.ala_schur_symmetry_tolerance &&
                 nonpositive_curvature_count == 0) ? "PASS" : "FAIL");
    }

    if(E->control.ala_pressure_buoyancy && E->parallel.me == 0 &&
       cancellation_l2 < E->control.tole_comp &&
       (!E->control.ala_hybrid_convergence || hybrid_converged)) {
        fprintf(E->fp,
                "Strict ALA BiCGStab converged by physical continuity: "
                "cancellation=%e tolerance=%e mass_relative=%e "
                "algebraic_relative=%e iterations=%d\n",
                cancellation_l2, E->control.tole_comp,
                mass_relative_residual, relative_residual, count);
        fflush(E->fp);
    }

    if(E->control.ala_pressure_buoyancy &&
       (cancellation_l2 >= E->control.tole_comp ||
        (E->control.ala_hybrid_convergence && !hybrid_converged))) {
        if(E->parallel.me == 0) {
            fprintf(stderr,
                    "Strict ALA BiCGStab failed physical continuity: "
                    "cancellation=%e tolerance=%e mass_relative=%e "
                    "algebraic_relative=%e hybrid_streak=%d/%d iterations=%d\n",
                    cancellation_l2, E->control.tole_comp,
                    mass_relative_residual, relative_residual,
                    hybrid_consecutive_count,
                    E->control.ala_consecutive_steps, count);
            fprintf(E->fp,
                    "Strict ALA BiCGStab failed physical continuity: "
                    "cancellation=%e tolerance=%e mass_relative=%e "
                    "algebraic_relative=%e hybrid_streak=%d/%d iterations=%d\n",
                    cancellation_l2, E->control.tole_comp,
                    mass_relative_residual, relative_residual,
                    hybrid_consecutive_count,
                    E->control.ala_consecutive_steps, count);
            fflush(E->fp);
        }
        parallel_process_termination();
    }


    for(m=1; m<=E->sphere.caps_per_proc; m++) {
    	free((void *) F[m]);
        free((void *) r1[m]);
        free((void *) r2[m]);
        free((void *) pt[m]);
        free((void *) p1[m]);
        free((void *) p2[m]);
        free((void *) rt[m]);
        free((void *) v0[m]);
        free((void *) s0[m]);
        free((void *) st[m]);
        free((void *) t0[m]);
        free((void *) div_u[m]);

        free((void *) u0[m]);
    }

    *steps_max=count;

    return(residual);

}


/* Solve compressible Stokes flow using
 * conjugate gradient (CG) iterations with an outer iteration
 */

static float solve_Ahat_p_fhat_iterCG(struct All_variables *E,
                                      double **V, double **P, double **F,
                                      double imp, int *steps_max)
{
    int m, i;
    int cycles, num_of_loop;
    double residual;
    double relative_err_v, relative_err_p;
    double *old_v[NCS], *old_p[NCS],*diff_v[NCS],*diff_p[NCS];

    const int npno = E->lmesh.npno;
    const int neq = E->lmesh.neq;
    const int lev = E->mesh.levmax;

    double global_vdot(),global_pdot();
    void assemble_forces();

    for (m=1;m<=E->sphere.caps_per_proc;m++)   {
    	old_v[m] = (double *)malloc(neq*sizeof(double));
    	diff_v[m] = (double *)malloc(neq*sizeof(double));
    	old_p[m] = (double *)malloc((npno+1)*sizeof(double));
    	diff_p[m] = (double *)malloc((npno+1)*sizeof(double));
    }

    cycles = E->control.p_iterations;

    residual = 1.0;
    relative_err_v = 1.0;
    relative_err_p = 1.0;
    num_of_loop = 0;

    while((relative_err_v >= E->control.relative_err_accuracy ||
           relative_err_p >= E->control.relative_err_accuracy) &&
          num_of_loop <= E->control.compress_iter_maxstep) {

        for (m=1;m<=E->sphere.caps_per_proc;m++) {
            for(i=0;i<neq;i++) old_v[m][i] = V[m][i];
            for(i=1;i<=npno;i++) old_p[m][i] = P[m][i];
        }

        if(E->control.ala_pressure_buoyancy)
            assemble_forces(E,0);

        residual = solve_Ahat_p_fhat_CG(E,V,P,F,E->control.accuracy,&cycles);

        for (m=1;m<=E->sphere.caps_per_proc;m++)
            for(i=0;i<neq;i++) diff_v[m][i] = V[m][i] - old_v[m][i];

        relative_err_v = sqrt( global_vdot(E,diff_v,diff_v,lev) /
                               (1.0e-32 + global_vdot(E,V,V,lev)) );

        for (m=1;m<=E->sphere.caps_per_proc;m++)
            for(i=1;i<=npno;i++) diff_p[m][i] = P[m][i] - old_p[m][i];

        relative_err_p = sqrt( global_pdot(E,diff_p,diff_p,lev) /
                               (1.0e-32 + global_pdot(E,P,P,lev)) );

        num_of_loop++;

        if(E->parallel.me == 0) {
            fprintf(stderr, "Relative error err_v / v = %e and err_p / p = %e after %d loops\n\n", relative_err_v, relative_err_p, num_of_loop);
            fprintf(E->fp, "Relative error err_v / v = %e and err_p / p = %e after %d loops\n\n", relative_err_v, relative_err_p, num_of_loop);
        }

    } /* end of while */

    for (m=1;m<=E->sphere.caps_per_proc;m++)   {
    	free((void *) old_v[m]);
    	free((void *) old_p[m]);
	free((void *) diff_v[m]);
	free((void *) diff_p[m]);
    }

    return(residual);
}


static double initial_vel_residual(struct All_variables *E,
                                   double **V, double **P, double **F,
                                   double imp)
{
    void assemble_del2_u();
    void assemble_grad_p();
    void strip_bcs_from_residual();
    int  solve_del2_u();
    double global_vdot();

    int neq = E->lmesh.neq;
    int gneq = E->mesh.neq;
    int lev = E->mesh.levmax;
    int i, m, valid;
    double v_res;

    v_res = sqrt(global_vdot(E, F, F, lev) / gneq);

    if (E->parallel.me==0) {
        fprintf(E->fp, "initial residue of momentum equation F %.9e %d\n",
                v_res, gneq);
        fprintf(stderr, "initial residue of momentum equation F %.9e %d\n",
                v_res, gneq);
    }


    /* F = F - grad(P) - K*V */
    assemble_grad_p(E, P, E->u1, lev);
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(i=0; i<neq; i++)
            F[m][i] = F[m][i] - E->u1[m][i];

    assemble_del2_u(E, V, E->u1, lev, 1);
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(i=0; i<neq; i++)
            F[m][i] = F[m][i] - E->u1[m][i];

    strip_bcs_from_residual(E, F, lev);


    /* solve K*u1 = F for u1 */
    valid = solve_del2_u(E, E->u1, F, imp*v_res, lev);
    if(!valid && (E->parallel.me==0)) {
        fputs("Warning: solver not converging! 0\n", stderr);
        fputs("Warning: solver not converging! 0\n", E->fp);
    }
    strip_bcs_from_residual(E, E->u1, lev);


    /* V = V + u1 */
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(i=0; i<neq; i++)
            V[m][i] += E->u1[m][i];

    return(v_res);
}



static double incompressibility_residual(struct All_variables *E,
                                         double **V, double **r)
{
    double global_pdot();
    double global_vdot();

    int gnpno = E->mesh.npno;
    int gneq = E->mesh.neq;
    int lev = E->mesh.levmax;
    double tmp1, tmp2;

    /* incompressiblity residual = norm(r) / norm(V) */

    tmp1 = global_vdot(E, V, V, lev);
    tmp2 = global_pdot(E, r, r, lev);
    E->monitor.incompressibility = sqrt((gneq / gnpno)
                                        *( (1.0e-32 + tmp2)
                                              / (1.0e-32 + tmp1) ));

    E->monitor.vdotv = tmp1;

    return(sqrt(tmp2/gnpno));;
}
