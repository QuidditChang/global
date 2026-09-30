/* Route B P1: constant-coefficient PIC advection + lumped FE diffusion.
 * No PG formula is replaced. Particle temperatures are persistent extraq.
 * Restart, sources, variable coefficients and legacy CBF are deferred. */
#include <math.h>
#include <float.h>
#include <string.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "parsing.h"
#include "pices.h"
#include "temperature_audit.h"

void temperatures_conform_bcs(struct All_variables *);
void tracer_advection(struct All_variables *);
void get_global_shape_fn();

void pices_fail(struct All_variables *E, const char *why)
{
    fprintf(stderr,"PICES_ERROR rank=%d step=%d: %s\n",
            E->parallel.me,E->monitor.solution_cycles,why);
    fflush(stderr);
    MPI_Abort(E->parallel.world, 72);
    abort();
}
static double *array(struct All_variables *E, size_t n)
{
    double *p=(double *)calloc(n,sizeof(double));
    if(!p) pices_fail(E,"allocation failed");
    return p;
}
static double reduce(struct All_variables *E, double v, MPI_Op op)
{
    double r; MPI_Allreduce(&v,&r,1,MPI_DOUBLE,op,E->parallel.world); return r;
}
static void exchange(struct All_variables *E, double *p)
{
    double *a[NCS]; a[1]=p; E->exchange_node_d(E,a,E->mesh.levmax);
}
static int fixed(struct All_variables *E,int n)
{ return (E->node[1][n] & (TBX|TBY|TBZ))!=0; }

void pices_parameters(struct All_variables *E)
{
    char method[32]; int restart=0,post=0;
    input_string("energy_solver",method,"pg",E->parallel.me);
    if(strcmp(method,"pg") && strcmp(method,"pices"))
        pices_fail(E,"energy_solver must be pg or pices");
    E->pices.enabled=!strcmp(method,"pices");
    E->pices.initialized=E->pices.moving=0;
    E->pices.slot=-1;
    input_boolean("pices_test_no_diffusion",&E->pices.no_diffusion,"off",E->parallel.me);
    input_int("pices_max_substeps",&E->pices.max_substeps,"10000,1,1000000",E->parallel.me);
    input_boolean("pices_checkpoint",&E->pices.checkpoint,"off",E->parallel.me);
    if(E->pices.enabled) {
        input_boolean("restart",&restart,"off",E->parallel.me);
        input_boolean("post_processing",&post,"off",E->parallel.me);
        if((restart && !E->pices.checkpoint) || post) pices_fail(E,"restart requires pices_checkpoint=on; postprocessing unsupported");
    }
}

void pices_validate(struct All_variables *E)
{
    int i;
    if(!E->pices.enabled) return;
    if(E->pices.checkpoint) {
        if(!E->viscosity.update_allowed || E->viscosity.RHEOL!=1 || E->viscosity.SDEPV ||
           E->viscosity.PDEPV || E->viscosity.CDEPV || E->viscosity.FREEZE ||
           E->viscosity.channel || E->viscosity.wedge || E->viscosity.weak_blobs || E->viscosity.weak_zones)
            pices_fail(E,"P2 checkpoint requires constant Newtonian viscosity and rebuilds");
        for(i=0;i<E->viscosity.num_mat;i++)
            if(E->viscosity.N0[i]!=E->viscosity.N0[0] || E->viscosity.E[i]!=0 || E->viscosity.Z[i]!=0)
                pices_fail(E,"P2 checkpoint requires uniform temperature-independent viscosity");
    }
    if(E->pices.max_substeps<1 || E->pices.max_substeps>1000000 ||
       E->trace.itperel<1 || strcmp(E->output.format,"ascii-gz") ||
       E->control.ala_pressure_buoyancy)
        pices_fail(E,"P1 requires positive tracer density, valid substep limit, ascii-gz and BA/EBA");
    if(E->sphere.caps!=12 || E->sphere.caps_per_proc!=1 || !E->control.tracer)
        pices_fail(E,"P1 requires full sphere, one cap per rank, tracer=on");
    if((E->control.restart && !E->pices.checkpoint) || E->control.post_p || E->control.stokes ||
       E->control.pseudo_free_surf || E->control.lith_age || E->composition.on ||
       E->control.mat_control || E->control.vbcs_file || E->trace.ic_method!=0 ||
       E->trace.reclassify_flavors)
        pices_fail(E,"unsupported P1 restart/forcing/composition/tracer initialization");
    if(E->control.Q0!=0 || E->control.tracer_enriched || E->control.disptn_number!=0)
        pices_fail(E,"P1 requires Q0=0, Di=0 and no enriched heating");
    for(i=0;i<PHASE_TRANSITIONS;i++)
        if(E->control.phase[i].entropy_jump!=0 || E->control.phase[i].density_jump!=0)
            pices_fail(E,"P1 phase transitions are unsupported");
    if(E->advection.filter_temperature || !E->advection.ADVECTION ||
       E->advection.fixed_timestep<=0 || E->mesh.toptbc!=1 || E->mesh.bottbc!=1)
        pices_fail(E,"P1 requires fixed positive dt, ADV=on, no filter, radial Dirichlet");
    if(E->output.CBF_frequency!=0 || E->output.output_q_surf_CBF || E->output.output_q_botm_CBF || E->output.write_q_files)
        pices_fail(E,"P1 legacy CBF/heat-flux output is unsupported");
    if(E->control.kT_exponent!=0 || E->control.kC_ratio!=1 ||
       E->control.kd_upper_linear!=0 || E->control.kd_upper_quadratic!=0 ||
       E->control.kd_lower_linear!=0 || E->control.kd_lower_quadratic!=0 ||
       E->control.kd_upper_prefactor!=E->control.kd_lower_prefactor)
        pices_fail(E,"P1 requires spatially constant conductivity");
    for(i=0;i<E->convection.heat_sources.number;i++)
        if(E->convection.heat_sources.Q[i]!=0) pices_fail(E,"P1 radiogenic sources unsupported");
}

/* Matrix element data are cached only under P1's constant-physics guards.
 * The accumulated element absolute row sums bound the assembled |A| rows:
 * shared physical elements occur once; opposite-signed entries may cancel
 * after assembly, so this bound is conservative rather than an equality. */
static void assemble(struct All_variables *E)
{
    struct Shape_function GN;
    struct Shape_function_dx dx;
    struct Shape_function_dA omega;
    double rtf[4][9],tg[9],rho[9],cp[9],k[9],kap[9],grad[9][3];
    double lo[3]={DBL_MAX,DBL_MAX,DBL_MAX},hi[3]={0,0,0};
    double bound=0,edge=DBL_MAX;
    static const int edges[3][4][2]={{{1,2},{4,3},{5,6},{8,7}},
                                   {{1,4},{2,3},{5,8},{6,7}},
                                   {{1,5},{2,6},{3,7},{4,8}}};
    int e,a,b,g,d,q,n;
    struct PICES_STATE *p=&E->pices;
    p->K=array(E,(E->lmesh.nel+1)*64);
    p->emass=array(E,(E->lmesh.nel+1)*8);
    p->mass=array(E,E->lmesh.nno+1); p->rate=array(E,E->lmesh.nno+1);
    p->length=array(E,E->lmesh.nel+1);
    for(e=1;e<=E->lmesh.nel;e++) {
        double *ke=p->K+e*64;
        get_global_shape_fn(E,e,&GN,&dx,&omega,0,1,rtf,E->mesh.levmax,1);
        for(g=1;g<=8;g++) {
            tg[g]=0;
            for(a=1;a<=8;a++) tg[g]+=E->N.vpt[GNVINDEX(a,g)]*E->T[1][E->ien[1][e].node[a]];
        }
        thermal_transport_at_gp(E,1,e,tg,E->control.reference_conductivity,rho,cp,k,kap);
        for(g=1;g<=8;g++) {
            double w=omega.vpt[g]*g_point[g].weight[2];
            double props[3]={rho[g],cp[g],k[g]};
            if(!(w>0) || !isfinite(w)) pices_fail(E,"nonpositive element quadrature");
            for(d=0;d<3;d++) { if(props[d]<lo[d]) lo[d]=props[d]; if(props[d]>hi[d]) hi[d]=props[d]; }
            if(p->no_diffusion) k[g]=0;
            for(a=1;a<=8;a++) {
                double mass=w*rho[g]*cp[g]*E->N.vpt[GNVINDEX(a,g)];
                n=E->ien[1][e].node[a]; p->mass[n]+=mass; p->emass[e*8+a-1]+=mass;
                grad[a][0]=dx.vpt[GNVXINDEX(0,a,g)]*rtf[3][g];
                grad[a][1]=dx.vpt[GNVXINDEX(1,a,g)]*rtf[3][g]/sin(rtf[1][g]);
                grad[a][2]=dx.vpt[GNVXINDEX(2,a,g)];
            }
            for(a=1;a<=8;a++) for(b=1;b<=8;b++)
                for(d=0;d<3;d++) ke[(a-1)*8+b-1]+=w*k[g]*grad[a][d]*grad[b][d];
        }
        p->length[e]=DBL_MAX;
        for(d=0;d<3;d++) {
            double rms=0;
            for(q=0;q<4;q++) {
                double l2=0;
                int n1=E->ien[1][e].node[edges[d][q][0]],n2=E->ien[1][e].node[edges[d][q][1]];
                for(a=1;a<=3;a++) { double v=E->x[1][a][n1]-E->x[1][a][n2]; l2+=v*v; }
                if(!(l2>0) || !isfinite(l2)) pices_fail(E,"invalid physical edge length");
                rms+=l2/4; if(sqrt(l2)<edge) edge=sqrt(l2);
            }
            if(sqrt(rms)<p->length[e]) p->length[e]=sqrt(rms);
        }
        for(a=1;a<=8;a++) for(b=1;b<=8;b++) p->rate[E->ien[1][e].node[a]]+=fabs(ke[(a-1)*8+b-1]);
    }
    for(d=0;d<3;d++) {
        lo[d]=reduce(E,lo[d],MPI_MIN); hi[d]=reduce(E,hi[d],MPI_MAX);
        if(!isfinite(hi[d]) || hi[d]<=0 || fabs(hi[d]-lo[d])>1e-10*hi[d])
            pices_fail(E,"P1 requires constant positive rho, Cp and conductivity");
    }
    p->kappa=p->no_diffusion ? 0 : lo[2]/(lo[0]*lo[1]);
    exchange(E,p->mass); exchange(E,p->rate);
    for(n=1;n<=E->lmesh.nno;n++) {
        if(!(p->mass[n]>0) || !isfinite(p->mass[n])) pices_fail(E,"nonpositive lumped capacity");
        if(!fixed(E,n) && p->rate[n]/p->mass[n]>bound) bound=p->rate[n]/p->mass[n];
    }
    bound=reduce(E,bound,MPI_MAX);
    if(!isfinite(bound)) pices_fail(E,"nonfinite diffusion stability bound");
    p->dt_heat=bound>0 ? .8/bound : DBL_MAX;
    p->min_edge=reduce(E,edge,MPI_MIN);
}

static double interp(const int *nodes,const double *w,const double *field)
{ int a; double t=0; for(a=1;a<=8;a++) t+=w[a]*field[nodes[a]]; return t; }

/* P maps either absolute Tp or a supplied, nonpersistent particle increment. */
static void project(struct All_variables *E,const double *values,double *out,int boundary)
{
    int p,a,n,nodes[9],e; double w[9];
    double *den=array(E,E->lmesh.nno+1);
    int *covered=(int *)calloc(E->lmesh.nel+1,sizeof(int));
    int empty=0,zero=0;
    if(!covered) pices_fail(E,"allocation failed");
    memset(out,0,(E->lmesh.nno+1)*sizeof(double));
    for(p=1;p<=E->trace.ntracers[1];p++) {
        tracer_temperature_weights(E,1,p,nodes,w);
        e=E->trace.ielement[1][p]; covered[e]++;
        for(a=1;a<=8;a++) { n=nodes[a]; den[n]+=w[a]; out[n]+=w[a]*values[p]; }
    }
    for(e=1;e<=E->lmesh.nel;e++) if(!covered[e]) empty++;
    exchange(E,out); exchange(E,den);
    for(n=1;n<=E->lmesh.nno;n++) {
        if(!isfinite(den[n]) || !isfinite(out[n])) pices_fail(E,"nonfinite projection");
        if(den[n]<=0) {
            zero++;
            if(!fixed(E,n)) {
                fprintf(stderr,"PICES_COVERAGE rank=%d node=%d particles=%d empty_elements=%d\n",E->parallel.me,n,E->trace.ntracers[1],empty);
                pices_fail(E,"zero particle coverage at free node");
            }
            out[n]=0;
        } else out[n]/=den[n];
        if(boundary && fixed(E,n)) out[n]=E->T[1][n];
    }
    if(boundary) fprintf(E->fp,"PICES_COVERAGE step=%d empty_elements=%d zero_boundary_nodes=%d\n",E->monitor.solution_cycles,empty,zero);
    free(den); free(covered);
}

static double energy(struct All_variables *E,const double *field)
{
    int e,a;double s=0;
    for(e=1;e<=E->lmesh.nel;e++) for(a=1;a<=8;a++)
        s+=E->pices.emass[e*8+a-1]*field[E->ien[1][e].node[a]];
    return reduce(E,s,MPI_SUM);
}
void pices_initialize(struct All_variables *E)
{
    int p,n,a,nodes[9]; double w[9],v;
    if(E->pices.initialized) pices_fail(E,"duplicate PICES initialization");
    pices_validate(E);
    if(E->pices.slot<0) pices_fail(E,"Tp slot was not registered before tracer allocation");
    for(p=1;p<=E->trace.ntracers[1];p++) {
        tracer_temperature_weights(E,1,p,nodes,w);v=0;
        for(a=1;a<=8;a++) v+=w[a]*E->T[1][nodes[a]];
        E->trace.extraq[1][E->pices.slot][p]=v;
    }
    for(n=1;n<=E->lmesh.nno;n++) E->Tdot[1][n]=0;
    assemble(E); E->pices.initialized=1;
    fprintf(E->fp,"PICES_INIT method=P1_v1 Tp_slot=%d ntracers=%d kappa=%.17g dt_heat=%.17g min_edge=%.17g no_diffusion=%d restart=unsupported\n",E->pices.slot,E->trace.ntracers[1],E->pices.kappa,E->pices.dt_heat,E->pices.min_edge,E->pices.no_diffusion);
    fflush(E->fp);
}

void pices_advance(struct All_variables *E)
{
    int n,e,a,b,p,s,ns,nodes[9];
    int nn=E->lmesh.nno,np;
    double dt=E->advection.timestep,speed=0,cfl,ds,w[9],remaining;
    double *g,*old,*dg,*sub,*mapped,*rhs,*previous;
    double before,advected,heated,mismatch=0,submax=0,tmin=DBL_MAX,tmax=-DBL_MAX;
    FILE *f;
    if(!E->pices.initialized) pices_fail(E,"PICES not initialized");
    f=fopen("rheo.dat","r"); if(f) { fclose(f); pices_fail(E,"rheo.dat temperature edits unsupported"); }
    if(!(dt>0) || !isfinite(dt)) pices_fail(E,"invalid timestep");
    for(n=1;n<=nn;n++) {
        double v=0;for(a=1;a<=3;a++) v+=E->sphere.cap[1].V[a][n]*E->sphere.cap[1].V[a][n];
        if(!isfinite(v)) pices_fail(E,"nonfinite velocity");
        if(sqrt(v)>speed) speed=sqrt(v);
    }
    speed=reduce(E,speed,MPI_MAX);cfl=dt*speed/E->pices.min_edge;
    if(cfl>.25) pices_fail(E,"fixed dt violates particle CFL<=0.25; reduce fixed_timestep");
    if(dt/E->pices.dt_heat>E->pices.max_substeps) pices_fail(E,"heat substep limit exceeded before particle movement");
    ns=(int)ceil(dt/E->pices.dt_heat);if(ns<1) ns=1;
    g=array(E,nn+1);old=array(E,nn+1);dg=array(E,nn+1);mapped=array(E,nn+1);rhs=array(E,nn+1);previous=array(E,nn+1);
    for(n=1;n<=nn;n++) previous[n]=E->T[1][n];
    before=energy(E,previous);
    E->pices.moving=1;tracer_advection(E);E->pices.moving=0;
    np=E->trace.ntracers[1];sub=array(E,np+1);
    temperatures_conform_bcs(E);
    project(E,E->trace.extraq[1][E->pices.slot],g,1);
    advected=energy(E,g);remaining=dt;
    for(s=0;s<ns;s++) {
        ds=(s==ns-1)?remaining:dt/ns;remaining-=ds;
        memcpy(old,g,(nn+1)*sizeof(double));memset(rhs,0,(nn+1)*sizeof(double));
        for(e=1;e<=E->lmesh.nel;e++) for(a=1;a<=8;a++) for(b=1;b<=8;b++)
            rhs[E->ien[1][e].node[a]]-=E->pices.K[e*64+(a-1)*8+b-1]*old[E->ien[1][e].node[b]];
        exchange(E,rhs);
        for(n=1;n<=nn;n++) {
            if(!fixed(E,n)) g[n]=old[n]+ds*rhs[n]/E->pices.mass[n];
            dg[n]=g[n]-old[n];
            if(!isfinite(g[n]) || E->data.Ttop+g[n]*E->data.ref_temperature<0) pices_fail(E,"invalid grid temperature");
        }
        for(p=1;p<=np;p++) {
            double factor,tau;
            tracer_temperature_weights(E,1,p,nodes,w);e=E->trace.ielement[1][p];
            tau=E->pices.kappa>0 ? E->pices.length[e]*E->pices.length[e]/E->pices.kappa : DBL_MAX;
            factor=E->pices.kappa>0 ? -expm1(-ds/tau):0;
            sub[p]=(interp(nodes,w,g)-E->trace.extraq[1][E->pices.slot][p])*factor;
            if(fabs(sub[p])>submax)submax=fabs(sub[p]);
        }
        project(E,sub,mapped,0);
        for(n=1;n<=nn;n++) mapped[n]=dg[n]-mapped[n];
        for(p=1;p<=np;p++) {
            double *tp=&E->trace.extraq[1][E->pices.slot][p];
            tracer_temperature_weights(E,1,p,nodes,w);
            *tp+=sub[p]+interp(nodes,w,mapped);
            if(!isfinite(*tp) || E->data.Ttop+*tp*E->data.ref_temperature<0) pices_fail(E,"invalid particle temperature");
        }
    }
    heated=energy(E,g);
    project(E,E->trace.extraq[1][E->pices.slot],mapped,0);
    for(n=1;n<=nn;n++) {
        double error=fabs(mapped[n]-g[n]);
        if(!fixed(E,n) && error>mismatch)mismatch=error;
        E->T[1][n]=g[n]; E->Tdot[1][n]=0; /* legacy Eulerian derivative intentionally unavailable */
        if(g[n]<tmin)tmin=g[n];if(g[n]>tmax)tmax=g[n];
    }
    E->monitor.T_interior=reduce(E,tmax,MPI_MAX);
    E->advection.timesteps++;E->advection.total_timesteps++;
    E->advection.last_sub_iterations=ns;E->monitor.elapsed_time+=dt;
    audit_temperature(E,"thermal_exit",ns,-1);
    fprintf(E->fp,"PICES_STEP step=%d rank=%d time=%.17g dt=%.17g substeps=%d particles=%d cfl=%.17g Tmin=%.17g Tmax=%.17g mismatch=%.17g subgrid_max=%.17g remap_energy=%.17g heat_energy=%.17g derivative=heat_stage_only\n",E->monitor.solution_cycles,E->parallel.me,(double)E->monitor.elapsed_time,dt,ns,np,cfl,reduce(E,tmin,MPI_MIN),reduce(E,tmax,MPI_MAX),reduce(E,mismatch,MPI_MAX),reduce(E,submax,MPI_MAX),advected-before,heated-advected);
    fflush(E->fp);
    free(g);free(old);free(dg);free(mapped);free(rhs);free(previous);free(sub);
}

/* Full checkpoint already restored Tp: rebuild only nonpersistent geometry. */
void pices_restore(struct All_variables *E)
{
    int p,n,nodes[9]; double w[9];
    pices_validate(E);
    if(E->pices.initialized || E->pices.slot<0) pices_fail(E,"invalid restore lifecycle");
    for(p=1;p<=E->trace.ntracers[1];p++) {
        double tp=E->trace.extraq[1][E->pices.slot][p];
        E->trace.ielement[1][p]=-99;
        tracer_temperature_weights(E,1,p,nodes,w);
        if(!isfinite(tp) || E->data.Ttop+tp*E->data.ref_temperature<0) pices_fail(E,"invalid restored Tp");
    }
    for(n=1;n<=E->lmesh.nno;n++)
        if(!isfinite(E->T[1][n]) || E->data.Ttop+E->T[1][n]*E->data.ref_temperature<0) pices_fail(E,"invalid restored T");
    assemble(E); E->pices.initialized=1;
    E->advection.total_timesteps=E->monitor.solution_cycles+1;
    E->monitor.T_interior=0;
    for(n=1;n<=E->lmesh.nno;n++) if(E->T[1][n]>E->monitor.T_interior) E->monitor.T_interior=E->T[1][n];
    E->monitor.T_interior=reduce(E,E->monitor.T_interior,MPI_MAX);
    fprintf(E->fp,"PICES_RESTORE step=%d particles=%d Tp_slot=%d time=%.17g\n",E->monitor.solution_cycles,E->trace.ntracers[1],E->pices.slot,(double)E->monitor.elapsed_time);
}
