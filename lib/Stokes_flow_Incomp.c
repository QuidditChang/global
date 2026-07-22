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
#include <stdlib.h>
#include <string.h>

void myerror(struct All_variables *,char *);

static float solve_Ahat_p_fhat(struct All_variables *E,
                               double **V, double **P, double **F,
                               double imp, int *steps_max);
static float solve_Ahat_p_fhat_CG(struct All_variables *E,
                                  double **V, double **P, double **F,
                                  double imp, int *steps_max);
static float solve_Ahat_p_fhat_BiCG(struct All_variables *E,
                                    double **V, double **P, double **F,
                                    double imp, int *steps_max);
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
static double strict_ala_inner_relative_accuracy(struct All_variables *E,
                                                 double outer_relative_residual);
static double strict_ala_inner_accuracy(struct All_variables *E,
                                        double **F, int lev,
                                        double relative_accuracy);


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
