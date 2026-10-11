/* Fixed-composition conductive initialization, on an independent thermal clock.
 * Backward Euler with old-state lumped effective capacity and nonlinear new-state
 * conductivity. This is first-order in time (including latent capacity lag), not
 * an exact enthalpy integrator. Validate max_dt by step refinement.
 *
 * Only E->T is written. No particle projection/advection, production clock,
 * reconstruction, boundary flags/TB arrays, heating, assimilation or phase state
 * cache is touched. Geometry and material functions match the production FE.
 */
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "phase_change.h"
#include "pices.h"
#include "thermal_preage.h"

void get_global_shape_fn();

struct Preage {
    struct All_variables *E;
    int nn, ne;
    double *K, *mass, *copies, *diag;
    double *capacity, *conductivity, *oldcapacity, *oldconductivity;
    double *fraction, *oldfraction;
    double *r, *z, *direction, *work;
};

static double *alloc(struct All_variables *E, size_t n)
{
    double *p=calloc(n,sizeof(double));
    if(!p) pices_fail(E,"thermal_preage allocation failed");
    return p;
}
static double global(struct All_variables *E,double x,MPI_Op op)
{
    double y; MPI_Allreduce(&x,&y,1,MPI_DOUBLE,op,E->parallel.world);return y;
}
static void share(struct All_variables *E,double *x)
{
    double *a[NCS];a[1]=x;E->exchange_node_d(E,a,E->mesh.levmax);
}
/* Local radial endpoints are often MPI interfaces, not physical boundaries. */
static int radial_index(struct All_variables *E,int n)
{ return (n-1)%E->lmesh.noz+E->lmesh.nzs; }
static int fixed(struct All_variables *E,int n)
{
    int k=radial_index(E,n);return k==1 || k==E->mesh.noz;
}
static void boundary(struct All_variables *E,double *t)
{
    int n;
    for(n=1;n<=E->lmesh.nno;n++) {
        int k=radial_index(E,n);
        if(k==1)t[n]=1.0;
        else if(k==E->mesh.noz)t[n]=E->refstate.temperature_surface;
    }
}

void thermal_preage_validate(struct All_variables *E)
{
    if(!E->control.thermal_preage)return;
    if(!isfinite(E->control.thermal_preage_Ma) || E->control.thermal_preage_Ma<0 ||
       !isfinite(E->control.thermal_preage_max_dt_Ma) || E->control.thermal_preage_max_dt_Ma<=0)
        pices_fail(E,"thermal_preage requires finite age>=0 and max_dt_Ma>0");
    if(!E->pices.enabled || !E->pices.eba || !E->control.eba_formulation ||
       E->sphere.caps!=12 || E->sphere.caps_per_proc!=1 || !E->control.lith_age ||
       E->control.restart || E->control.post_p || E->convection.tic_method==-1)
        pices_fail(E,"thermal_preage requires fresh full-sphere PICES EBA lith_age initialization");
    if(!E->refstate.has_temperature || !E->refstate.Tref ||
       !isfinite(E->refstate.temperature_surface) || !isfinite(E->refstate.temperature_cmb) ||
       !isfinite(E->data.scalet) || E->data.scalet<=0 ||
       E->mesh.toptbc!=1 || E->mesh.bottbc!=1 ||
       fabs(E->control.TBCtopval)>1e-6 || fabs(E->control.TBCbotval-1)>1e-6)
        pices_fail(E,"thermal_preage requires Tref, positive time scale and physical radial Dirichlet 0/1");
    if(E->data.Ttop<=0 || !isfinite(E->data.Ttop) ||
       !isfinite(E->data.ref_temperature) || E->data.ref_temperature<=0)
        pices_fail(E,"thermal_preage requires positive physical Ttop and DeltaT");
}

/* Assemble production quadrature K. When requested, freeze the row-sum mass
 * at the accepted old field. The phase capacity matches Pices.c exactly;
 * its pressure source is absent because this pre-age has no flow. */
static void assemble(struct Preage *A,const double *t,int freeze_mass)
{
    struct All_variables *E=A->E;
    struct Shape_function N;struct Shape_function_dx dx;struct Shape_function_dA omega;
    double rtf[4][9],tg[9],rho[9],cp[9],k[9],kap[9],grad[9][3];
    int e,g,a,b,d,n;
    memset(A->K,0,(A->ne+1)*64*sizeof(double));
    if(freeze_mass)memset(A->mass,0,(A->nn+1)*sizeof(double));
    for(e=1;e<=A->ne;e++) {
        get_global_shape_fn(E,e,&N,&dx,&omega,0,1,rtf,E->mesh.levmax,1);
        for(g=1;g<=8;g++) {
            tg[g]=0;
            for(a=1;a<=8;a++)tg[g]+=E->N.vpt[GNVINDEX(a,g)]*t[E->ien[1][e].node[a]];
        }
        thermal_transport_at_gp(E,1,e,tg,E->control.reference_conductivity,rho,cp,k,kap);
        for(g=1;g<=8;g++) {
            double w=omega.vpt[g]*g_point[g].weight[2],capacity=rho[g]*cp[g],rg=0;
            int j,idx=e*8+g-1;
            for(a=1;a<=8;a++) {
                int nz=(E->ien[1][e].node[a]-1)%E->lmesh.noz+1;
                rg+=E->refstate.rho[nz]*E->refstate.gravity[nz]*E->N.vpt[GNVINDEX(a,g)];
            }
            for(j=0;j<PHASE_TRANSITIONS;j++) {
                double q,x,dT,dr;
                const struct Phase_transition *phase=&E->control.phase[j];
                phase_change_state(phase,E->sphere.ro-1/rtf[3][g],tg[g],1,rg,0,&q,&x,&dT,&dr);
                capacity+=rho[g]*(tg[g]+E->control.surface_temp)*phase->entropy_jump*dT;
                A->fraction[idx*PHASE_TRANSITIONS+j]=x;
            }
            if(!(w>0) || !(k[g]>0) || !(capacity>0) || !isfinite(k[g]) || !isfinite(capacity))
                pices_fail(E,"thermal_preage nonpositive geometry/conductivity/effective capacity");
            A->capacity[idx]=capacity;A->conductivity[idx]=k[g];
            for(a=1;a<=8;a++) {
                n=E->ien[1][e].node[a];
                if(freeze_mass)A->mass[n]+=w*capacity*E->N.vpt[GNVINDEX(a,g)];
                grad[a][0]=dx.vpt[GNVXINDEX(0,a,g)]*rtf[3][g];
                grad[a][1]=dx.vpt[GNVXINDEX(1,a,g)]*rtf[3][g]/sin(rtf[1][g]);
                grad[a][2]=dx.vpt[GNVXINDEX(2,a,g)];
            }
            for(a=1;a<=8;a++)for(b=1;b<=8;b++)for(d=0;d<3;d++)
                A->K[e*64+(a-1)*8+b-1]+=w*k[g]*grad[a][d]*grad[b][d];
        }
    }
    if(freeze_mass) {
        share(E,A->mass);
        for(n=1;n<=A->nn;n++)if(!(A->mass[n]>0) || !isfinite(A->mass[n]))
            pices_fail(E,"thermal_preage invalid assembled mass");
    }
}

static void apply(struct Preage *A,const double *x,double *y,double dt)
{
    int e,a,b,n;
    memset(y,0,(A->nn+1)*sizeof(double));
    for(e=1;e<=A->ne;e++)for(a=1;a<=8;a++)for(b=1;b<=8;b++)
        y[A->E->ien[1][e].node[a]]+=A->K[e*64+(a-1)*8+b-1]*x[A->E->ien[1][e].node[b]];
    share(A->E,y);
    for(n=1;n<=A->nn;n++)y[n]=dt*y[n]+A->mass[n]*x[n];
}
static double dot(struct Preage *A,const double *x,const double *y)
{
    int n;double s=0;
    for(n=1;n<=A->nn;n++)s+=x[n]*y[n]/A->copies[n];
    return global(A->E,s,MPI_SUM);
}
static double residual(struct Preage *A,const double *old,const double *x,double dt)
{
    int n;double error=0;
    apply(A,x,A->work,dt);
    for(n=1;n<=A->nn;n++) {
        A->r[n]=fixed(A->E,n)?0:A->mass[n]*old[n]-A->work[n];
        if(!isfinite(A->r[n]))error=DBL_MAX;
        else error=fmax(error,fabs(A->r[n])/A->mass[n]);
    }
    return global(A->E,error,MPI_MAX);
}
/* Homogeneous corrections retain symmetry despite inhomogeneous boundary T. */
static int solve(struct Preage *A,const double *old,double *x,double dt,int *iterations)
{
    int n,e,a,it;double rz,next,pap,alpha,beta,error;
    const double tolerance=1e-11;
    memset(A->diag,0,(A->nn+1)*sizeof(double));
    for(e=1;e<=A->ne;e++)for(a=1;a<=8;a++)
        A->diag[A->E->ien[1][e].node[a]]+=A->K[e*64+(a-1)*8+a-1];
    share(A->E,A->diag);
    for(n=1;n<=A->nn;n++)A->diag[n]=A->mass[n]+dt*A->diag[n];
    error=residual(A,old,x,dt);
    if(error<=tolerance)return 1;
    for(n=1;n<=A->nn;n++) {
        A->z[n]=A->r[n]/A->diag[n];A->direction[n]=A->z[n];
    }
    rz=dot(A,A->r,A->z);
    for(it=0;it<2000;it++) {
        apply(A,A->direction,A->work,dt);
        pap=dot(A,A->direction,A->work);
        if(!(pap>0) || !isfinite(pap) || !(rz>0) || !isfinite(rz))return 0;
        alpha=rz/pap;
        for(n=1;n<=A->nn;n++)if(!fixed(A->E,n))x[n]+=alpha*A->direction[n];
        /* True residual on every iteration avoids recursive false convergence. */
        error=residual(A,old,x,dt);
        (*iterations)++;
        if(!isfinite(error))return 0;
        if(error<=tolerance)return 1;
        for(n=1;n<=A->nn;n++)A->z[n]=A->r[n]/A->diag[n];
        next=dot(A,A->r,A->z);beta=next/rz;rz=next;
        for(n=1;n<=A->nn;n++)A->direction[n]=A->z[n]+beta*A->direction[n];
    }
    return 0;
}

void thermal_preage_run(struct All_variables *E)
{
    struct Preage A;
    double age=0,duration,maxdt,dt,last_error=0,ledger_storage=0,ledger_boundary=0;
    double *old,*trial;int n,j,step=0,rejected=0,total_cg=0,total_picard=0;
    if(!E->control.thermal_preage)return;
    thermal_preage_validate(E);
    memset(&A,0,sizeof(A));A.E=E;A.nn=E->lmesh.nno;A.ne=E->lmesh.nel;
    A.K=alloc(E,(A.ne+1)*64);A.mass=alloc(E,A.nn+1);A.copies=alloc(E,A.nn+1);
    A.diag=alloc(E,A.nn+1);A.r=alloc(E,A.nn+1);A.z=alloc(E,A.nn+1);
    A.direction=alloc(E,A.nn+1);A.work=alloc(E,A.nn+1);
    A.capacity=alloc(E,(A.ne+1)*8);A.conductivity=alloc(E,(A.ne+1)*8);
    A.oldcapacity=alloc(E,(A.ne+1)*8);A.oldconductivity=alloc(E,(A.ne+1)*8);
    A.fraction=alloc(E,(A.ne+1)*8*PHASE_TRANSITIONS);
    A.oldfraction=alloc(E,(A.ne+1)*8*PHASE_TRANSITIONS);
    old=alloc(E,A.nn+1);trial=alloc(E,A.nn+1);
    for(n=1;n<=A.nn;n++)A.copies[n]=1;
    share(E,A.copies);
    duration=E->control.thermal_preage_Ma/E->data.scalet;
    maxdt=E->control.thermal_preage_max_dt_Ma/E->data.scalet;
    if(!isfinite(duration) || !isfinite(maxdt) || !(maxdt>0))
        pices_fail(E,"thermal_preage invalid scaled duration or timestep");
    boundary(E,E->T[1]);
    {
        double bad=0;
        for(n=1;n<=A.nn;n++)if(!isfinite(E->T[1][n]) ||
            E->data.Ttop+E->data.ref_temperature*E->T[1][n]<=0)bad=1;
        if(global(E,bad,MPI_MAX))pices_fail(E,"thermal_preage invalid initial absolute temperature");
    }
    dt=fmin(maxdt,duration);
    fprintf(E->fp,"THERMAL_PREAGE_BEGIN age_Ma=%.17g max_dt_Ma=%.17g top=%.17g bottom=1 method=BE_frozen_phase_capacity sources=none composition=fixed production_time=%.17g\n",
        E->control.thermal_preage_Ma,E->control.thermal_preage_max_dt_Ma,E->refstate.temperature_surface,(double)E->monitor.elapsed_time);
    while(age<duration) {
        int attempt,accepted=0;double variation=0;
        if(step>=1000000)pices_fail(E,"thermal_preage exceeded one million accepted steps");
        memcpy(old,E->T[1],(A.nn+1)*sizeof(double));
        assemble(&A,old,1);
        memcpy(A.oldcapacity,A.capacity,(A.ne+1)*8*sizeof(double));
        memcpy(A.oldconductivity,A.conductivity,(A.ne+1)*8*sizeof(double));
        memcpy(A.oldfraction,A.fraction,(A.ne+1)*8*PHASE_TRANSITIONS*sizeof(double));
        dt=fmin(dt,duration-age);
        for(attempt=0;attempt<40;attempt++) {
            int it,converged=0;double invalid=0;
            memcpy(trial,old,(A.nn+1)*sizeof(double));
            if(!(dt>0) || age+dt==age)pices_fail(E,"thermal_preage timestep underflow");
            for(it=0;it<50;it++) {
                assemble(&A,trial,0);
                if(!solve(&A,old,trial,dt,&total_cg))break;
                for(n=1;n<=A.nn;n++)if(!isfinite(trial[n]) ||
                    E->data.Ttop+E->data.ref_temperature*trial[n]<=0)invalid=1;
                if(global(E,invalid,MPI_MAX))break;
                assemble(&A,trial,0);
                last_error=residual(&A,old,trial,dt);total_picard++;
                if(last_error<=1e-10) {converged=1;break;}
            }
            variation=0;
            if(converged)for(j=8;j<(A.ne+1)*8;j++) {
                variation=fmax(variation,fabs(A.capacity[j]/A.oldcapacity[j]-1)/.05);
                variation=fmax(variation,fabs(A.conductivity[j]/A.oldconductivity[j]-1)/.1);
            }
            /* Capacity endpoints can both be small after crossing a narrow phase.
             * Bound the phase fraction as production does, not just capacity. */
            if(converged)for(j=8*PHASE_TRANSITIONS;j<(A.ne+1)*8*PHASE_TRANSITIONS;j++)
                variation=fmax(variation,fabs(A.fraction[j]-A.oldfraction[j])/.05);
            variation=global(E,variation,MPI_MAX);
            if(converged && variation<=1) {accepted=1;break;}
            dt*=.5;rejected++;
        }
        if(!accepted)pices_fail(E,"thermal_preage nonlinear/property adaptation did not converge");
        /* Numerical BE ledger with frozen capacity; shared-node multiplicity
         * is removed. This is NOT an exact nonlinear enthalpy diagnostic. */
        {
            double storage=0,reaction=0;
            apply(&A,trial,A.work,dt);
            for(n=1;n<=A.nn;n++) {
                storage+=A.mass[n]*(trial[n]-old[n])/A.copies[n];
                if(fixed(E,n))reaction+=(A.work[n]-A.mass[n]*trial[n])/A.copies[n];
            }
            ledger_storage+=global(E,storage,MPI_SUM);
            ledger_boundary+=global(E,reaction,MPI_SUM);
        }
        memcpy(E->T[1],trial,(A.nn+1)*sizeof(double));
        age+=dt;step++;
        if(step==1 || step%25==0 || age>=duration) {
            fprintf(E->fp,"THERMAL_PREAGE_PROGRESS step=%d age_Ma=%.17g dt_Ma=%.17g residual=%.17g property_fraction=%.17g\n",step,age*E->data.scalet,dt*E->data.scalet,last_error,variation);
            fflush(E->fp);
        }
        dt=fmin(maxdt,dt*1.5);
    }
    fprintf(E->fp,"THERMAL_PREAGE_END age_Ma=%.17g steps=%d rejected=%d cg=%d picard=%d residual=%.17g frozen_storage=%.17g boundary_input=%.17g balance=%.17g production_time=%.17g\n",age*E->data.scalet,step,rejected,total_cg,total_picard,last_error,ledger_storage,ledger_boundary,ledger_storage-ledger_boundary,(double)E->monitor.elapsed_time);
    fflush(E->fp);
    free(A.K);free(A.mass);free(A.copies);free(A.diag);free(A.r);free(A.z);free(A.direction);free(A.work);
    free(A.capacity);free(A.conductivity);free(A.oldcapacity);free(A.oldconductivity);
    free(A.fraction);free(A.oldfraction);free(old);free(trial);
}
