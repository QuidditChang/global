/* Route B PIC advection + lumped FE heat solve. P3 EBA is opt-in.
 * Particle temperatures persist through migration and same-partition restart.
 * P4 adds one accepted-step TA and a cached material heat residual for CBF. */
#include <math.h>
#include <float.h>
#include <string.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "parsing.h"
#include "pices.h"
#include "pices_benchmark.h"
#include "temperature_audit.h"
#include "phase_change.h"
#include "lith_age.h"
#include "pices_coverage_stats.h"
void thermal_heat_sources(struct All_variables *, int, double *, double *, double *);

void temperatures_conform_bcs(struct All_variables *);
void get_global_shape_fn();

/* One accepted-step snapshot feeds both text and NPZ diagnostics. Q arrays
 * are element-volume averages, time averaged over accepted heat substeps.
 * adi and phase follow the existing profile sink-positive convention. */
static double energy(struct All_variables *,const double *);
static struct Pices_diagnostics diagnostics;
const char *pices_diagnostic_names[PD_NVALUES]={
 "dt_nd","time_nd","energy_before_nd","energy_after_nd","remap_nd",
 "capacity_rebase_nd","heat_storage_nd","boundary_reaction_nd",
 "ta_storage_nd","reconciliation_nd","closure_residual_nd","heat_residual_nd",
 "particle_count_global","particle_Tmin_K","particle_Tmax_K","particle_Tmean_count_K",
 "grid_Tmin_K","grid_Tmax_K","energy_scale_J","source_scale_W_m3","internal_input_nd","adiabatic_sink_nd",
 "viscous_input_nd","phase_sink_nd","viscous_raw_nd"};
const struct Pices_diagnostics *pices_diagnostics(struct All_variables *E)
{ return E->pices.enabled ? &diagnostics : NULL; }
double pices_profile_source(struct All_variables *E,int element,int term)
{ return diagnostics.valid && diagnostics.step==E->monitor.solution_cycles ?
    diagnostics.q[element*PD_NQ+term] : NAN; }

void pices_fail(struct All_variables *E, const char *why)
{
    fprintf(stderr,"PICES_ERROR rank=%d step=%d: %s\n",
            E->parallel.me,E->monitor.solution_cycles,why);
    fflush(stderr);
    MPI_Abort(E->parallel.world, 72);
    abort();
}
int pices_temperature_valid(struct All_variables *E,double t)
{
    double kelvin=E->data.Ttop+E->data.ref_temperature*t;
    return isfinite(t) && isfinite(kelvin) && kelvin>=-1e-9;
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
    char method[32],projection[32]; int restart=0,post=0;
    input_string("energy_solver",method,"pg",E->parallel.me);
    if(strcmp(method,"pg") && strcmp(method,"pices"))
        pices_fail(E,"energy_solver must be pg or pices");
    E->pices.enabled=!strcmp(method,"pices");
    input_string("pices_projection",projection,"lumped",E->parallel.me);
    if(strcmp(projection,"lumped") && strcmp(projection,"bounded_consistent"))pices_fail(E,"unknown pices_projection");
    E->pices.consistent_projection=!strcmp(projection,"bounded_consistent");
    if(E->pices.consistent_projection && !E->pices.enabled)pices_fail(E,"pices_projection requires PICES");
    E->pices.initialized=0;
    E->pices.slot=-1;
    input_boolean("pices_test_no_diffusion",&E->pices.no_diffusion,"off",E->parallel.me);
    input_double("pices_max_timestep_Ma",&E->pices.max_timestep_Ma,"0.1",E->parallel.me);
    input_int("pices_max_substeps",&E->pices.max_substeps,"10000,1,1000000",E->parallel.me);
    input_boolean("pices_ta_before_heat",&E->pices.ta_before_heat,"off",E->parallel.me);
    input_boolean("pices_p4",&E->pices.p4,"off",E->parallel.me);
    input_boolean("pices_eba",&E->pices.eba,"off",E->parallel.me);
    input_boolean("pices_checkpoint",&E->pices.checkpoint,"off",E->parallel.me);
    if(E->pices.ta_before_heat && (!E->pices.enabled || !E->pices.p4))
        pices_fail(E,"pre-heat TA requires PICES and pices_p4");
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
    if(E->pices.p4 && (!E->pices.eba || (E->control.lith_age &&
       (!E->control.lith_age_asml || !E->control.lith_age_time || E->control.temperature_bound_adj))))
        pices_fail(E,"P4 requires P3 EBA and relaxed TA with time-dependent ages and no boundary adjustment");
    if(!isfinite(E->control.Q0) || !isfinite(E->control.disptn_number))
        pices_fail(E,"nonfinite EBA heating parameters");
    if(!isfinite(E->pices.max_timestep_Ma) || E->pices.max_timestep_Ma<=0)
        pices_fail(E,"pices_max_timestep_Ma must be finite and positive");
    if(E->pices.checkpoint) {
        if(pices_checkpoint_coupled(E) && !E->pices.p4)
            pices_fail(E,"coupled checkpoint requires P4 EBA heat state");
        if(E->composition.on && !E->pices.p4)
            pices_fail(E,"composition checkpoint requires P4 heat-state schema");
        if(!E->viscosity.update_allowed || !E->viscosity.FROM_SYSTEM || (E->viscosity.RHEOL!=1 && E->viscosity.RHEOL!=7) || E->viscosity.SDEPV ||
           E->viscosity.PDEPV || E->viscosity.CDEPV || E->viscosity.FREEZE ||
           E->viscosity.channel || E->viscosity.wedge || E->viscosity.weak_blobs || E->viscosity.weak_zones)
            pices_fail(E,"checkpoint requires supported Newtonian rheology and viscosity rebuilds");
        if(E->viscosity.RHEOL==1)for(i=0;i<E->viscosity.num_mat;i++)
            if(E->viscosity.N0[i]!=E->viscosity.N0[0] || E->viscosity.E[i]!=0 || E->viscosity.Z[i]!=0)
                pices_fail(E,"P2 checkpoint requires uniform temperature-independent viscosity");
    }
    if(E->composition.on && (!E->pices.eba || E->composition.ibuoy_type!=1 ||
       E->viscosity.CDEPV || E->composition.ncomp!=E->trace.nflavors-1))
        pices_fail(E,"P8a supports EBA ratio composition without compositional rheology");
    if(E->trace.reclassify_flavors && (E->trace.nflavors!=25 ||
       E->control.kC_primordial_flavor!=24 || !E->composition.on ||
       !E->pices.p4 || !E->control.lith_age))
        pices_fail(E,"P8a reclassification requires 25 flavors, primordial 24 and P4 age forcing");
    if(E->control.vbcs_file && (!E->pices.eba || E->mesh.topvbc!=1 || E->data.timedir!=1 || E->control.remove_rigid_rotation))
        pices_fail(E,"P8b plate forcing requires forward EBA, prescribed top velocities and remove_rigid_rotation=off");
    if(E->pices.max_substeps<1 || E->pices.max_substeps>1000000 ||
       E->trace.itperel<1 ||
       (strcmp(E->output.format,"ascii") && strcmp(E->output.format,"ascii-gz")) ||
       E->control.ala_pressure_buoyancy)
        pices_fail(E,"P1 requires positive tracer density, valid substep limit, ascii/ascii-gz and BA/EBA");
    if(E->sphere.caps!=12 || E->sphere.caps_per_proc!=1 || !E->control.tracer)
        pices_fail(E,"P1 requires full sphere, one cap per rank, tracer=on");
    if((E->control.restart && !E->pices.checkpoint) || E->control.post_p || E->control.stokes ||
       E->control.pseudo_free_surf || (E->control.lith_age && !E->pices.p4) ||
       E->control.mat_control || E->trace.ic_method!=0)
        pices_fail(E,"unsupported P1 restart/forcing/composition/tracer initialization");
    if(E->control.tracer_enriched || (E->pices.eba && !E->control.eba_formulation))
        pices_fail(E,"PICES requires EBA and no enriched heating");
    if(!E->pices.eba && (E->control.Q0!=0 || E->control.disptn_number!=0))
        pices_fail(E,"P1 requires Q0=0, Di=0 and no enriched heating");
    for(i=0;!E->pices.eba && i<PHASE_TRANSITIONS;i++)
        if(E->control.phase[i].entropy_jump!=0 || E->control.phase[i].density_jump!=0)
            pices_fail(E,"P1 phase transitions are unsupported");
    if(E->advection.filter_temperature || !E->advection.ADVECTION ||
       !isfinite(E->advection.fixed_timestep) || E->advection.fixed_timestep<0 || E->mesh.toptbc!=1 || E->mesh.bottbc!=1)
        pices_fail(E,"PICES requires nonnegative fixed dt, ADV=on, no filter, radial Dirichlet");
    if((!E->pices.p4 && (E->output.CBF_frequency!=0 || E->output.output_q_surf_CBF || E->output.output_q_botm_CBF)) || E->output.write_q_files)
        pices_fail(E,"P1 legacy CBF/heat-flux output is unsupported");
    if(!E->pices.eba && (E->control.kT_exponent!=0 || E->control.kC_ratio!=1 ||
       E->control.kd_upper_linear!=0 || E->control.kd_upper_quadratic!=0 ||
       E->control.kd_lower_linear!=0 || E->control.kd_lower_quadratic!=0 ||
       E->control.kd_upper_prefactor!=E->control.kd_lower_prefactor))
        pices_fail(E,"P1 requires spatially constant conductivity");
    for(i=0;i<E->convection.heat_sources.number;i++)
        if(E->convection.heat_sources.Q[i]!=0) pices_fail(E,"P1 radiogenic sources unsupported");
}

/* Split the existing PG latent term, retaining physical radial velocity. */
static void phase_coefficients(struct All_variables *E,int e,int g,
    const struct Shape_function_dx *dx,double inverse_r,double t,double rho,
    double cp,double *capacity,double *pressure,double *fractions)
{
    int a,j,n,nz; double rg=0,drg=0,ur=0;
    *capacity=rho*cp; if(pressure)*pressure=0;
    for(a=1;a<=8;a++) {
        double N=E->N.vpt[GNVINDEX(a,g)];
        n=E->ien[1][e].node[a];nz=(n-1)%E->lmesh.noz+1;
        rg+=E->refstate.rho[nz]*E->refstate.gravity[nz]*N;
        drg+=E->refstate.rho[nz]*E->refstate.gravity[nz]*dx->vpt[GNVXINDEX(2,a,g)];
        if(pressure)ur+=E->sphere.cap[1].V[3][n]*N;
    }
    for(j=0;j<PHASE_TRANSITIONS;j++) {
        double q,x,dT,dr;
        const struct Phase_transition *phase=&E->control.phase[j];
        phase_change_state(phase,E->sphere.ro-1/inverse_r,t,1,rg,drg,&q,&x,&dT,&dr);
        fractions[j]=x;
        *capacity+=rho*(t+E->control.surface_temp)*phase->entropy_jump*dT;
        if(pressure && phase->entropy_jump!=0)*pressure-=rho*(t+E->control.surface_temp)*phase->entropy_jump*dr*ur;
    }
    if(!isfinite(*capacity) || !(*capacity>0) || (pressure && !isfinite(*pressure)))
        pices_fail(E,"nonpositive/nonfinite effective phase capacity or pressure source");
}

/* P1 caches matrices; P3 rebuilds them from each frozen heat-stage field.
 * The accumulated element absolute row sums bound the assembled |A| rows:
 * shared physical elements occur once; opposite-signed entries may cancel
 * after assembly, so this bound is conservative rather than an equality. */
static void assemble(struct All_variables *E,const double *field)
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
    if(!p->K) {
        p->K=array(E,(E->lmesh.nel+1)*64);
        p->emass=array(E,(E->lmesh.nel+1)*8);
        p->mass=array(E,E->lmesh.nno+1); p->rate=array(E,E->lmesh.nno+1);
        p->length=array(E,E->lmesh.nel+1);p->ekappa=array(E,E->lmesh.nel+1);
        p->gp_capacity=array(E,(E->lmesh.nel+1)*8);
        p->gp_k=array(E,(E->lmesh.nel+1)*8);
        p->gp_fraction=array(E,(E->lmesh.nel+1)*8*PHASE_TRANSITIONS);
    }
    if(p->p4) {
        if(!p->heat_residual)p->heat_residual=array(E,(E->lmesh.nel+1)*8);
        if(!p->element_source)p->element_source=array(E,(E->lmesh.nel+1)*8);
    }
    memset(p->K,0,(E->lmesh.nel+1)*64*sizeof(double));
    memset(p->emass,0,(E->lmesh.nel+1)*8*sizeof(double));
    memset(p->mass,0,(E->lmesh.nno+1)*sizeof(double));
    memset(p->rate,0,(E->lmesh.nno+1)*sizeof(double));
    for(e=1;e<=E->lmesh.nel;e++) {
        double *ke=p->K+e*64;
        p->ekappa[e]=0;
        get_global_shape_fn(E,e,&GN,&dx,&omega,0,1,rtf,E->mesh.levmax,1);
        for(g=1;g<=8;g++) {
            tg[g]=0;
            for(a=1;a<=8;a++) tg[g]+=E->N.vpt[GNVINDEX(a,g)]*field[E->ien[1][e].node[a]];
        }
        thermal_transport_at_gp(E,1,e,tg,E->control.reference_conductivity,rho,cp,k,kap);
        for(g=1;g<=8;g++) {
            double w=omega.vpt[g]*g_point[g].weight[2];
            double props[3]={rho[g],cp[g],k[g]};
            double capacity=rho[g]*cp[g],pressure; int idx=e*8+g-1;
            if(p->eba) phase_coefficients(E,e,g,&dx,rtf[3][g],tg[g],rho[g],cp[g],&capacity,NULL,p->gp_fraction+idx*PHASE_TRANSITIONS);
            if(!(k[g]>0) || !isfinite(capacity) || !(capacity>0))
                pices_fail(E,"nonpositive conductivity/effective capacity");
            p->gp_capacity[idx]=capacity;p->gp_k[idx]=k[g];
            if(!p->no_diffusion && k[g]/capacity>p->ekappa[e]) p->ekappa[e]=k[g]/capacity;
            if(!(w>0) || !isfinite(w)) pices_fail(E,"nonpositive element quadrature");
            for(d=0;d<3;d++) { if(props[d]<lo[d]) lo[d]=props[d]; if(props[d]>hi[d]) hi[d]=props[d]; }
            if(p->no_diffusion) k[g]=0;
            for(a=1;a<=8;a++) {
                double mass=p->eba ? w*capacity*E->N.vpt[GNVINDEX(a,g)] : w*rho[g]*cp[g]*E->N.vpt[GNVINDEX(a,g)];
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
        if(p5_length_scale()!=1)p->length[e]*=p5_length_scale();
        for(a=1;a<=8;a++) for(b=1;b<=8;b++) p->rate[E->ien[1][e].node[a]]+=fabs(ke[(a-1)*8+b-1]);
    }
    for(d=0;d<3;d++) {
        lo[d]=reduce(E,lo[d],MPI_MIN); hi[d]=reduce(E,hi[d],MPI_MAX);
        if(!isfinite(hi[d]) || hi[d]<=0 || (!p->eba && fabs(hi[d]-lo[d])>1e-10*hi[d]))
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

/* Element heating follows the established EBA quadrature convention. The
 * shared source routines use actual velocity; only their temperature view is
 * switched to the frozen heat-stage field. No legacy CBF residual is called. */
static void heat_load(struct All_variables *E,const double *field,double *load,double totals[5])
{
    int e,a,g,z,n; struct Shape_function GN;struct Shape_function_dx dx;
    struct Shape_function_dA omega; double rtf[4][9],tg[9],rho[9],cp[9],k[9],kap[9];
    double *adi=array(E,E->lmesh.nel+1),*visc=array(E,E->lmesh.nel+1),*raw=array(E,E->lmesh.nel+1),*save=E->T[1];
    E->T[1]=(double *)field;thermal_heat_sources(E,1,adi,visc,raw);E->T[1]=save;
    if(E->pices.p4)memset(E->pices.element_source,0,(E->lmesh.nel+1)*8*sizeof(double));
    memset(load,0,(E->lmesh.nno+1)*sizeof(double));memset(totals,0,5*sizeof(double));
    memset(diagnostics.trial,0,(E->lmesh.nel+1)*PD_NQ*sizeof(double));
    memset(diagnostics.volume,0,(E->lmesh.nel+1)*sizeof(double));
    for(e=1;e<=E->lmesh.nel;e++) {
        get_global_shape_fn(E,e,&GN,&dx,&omega,0,1,rtf,E->mesh.levmax,1);
        z=(e-1)%E->lmesh.elz+1;
        for(g=1;g<=8;g++) { tg[g]=0;for(a=1;a<=8;a++) tg[g]+=E->N.vpt[GNVINDEX(a,g)]*field[E->ien[1][e].node[a]]; }
        thermal_transport_at_gp(E,1,e,tg,E->control.reference_conductivity,rho,cp,k,kap);
        for(g=1;g<=8;g++) {
            double capacity,pressure,fractions[PHASE_TRANSITIONS],w=omega.vpt[g]*g_point[g].weight[2];
            double internal=.5*(E->refstate.rho[z]+E->refstate.rho[z+1])*E->control.Q0;
            phase_coefficients(E,e,g,&dx,rtf[3][g],tg[g],rho[g],cp[g],&capacity,&pressure,fractions);
            diagnostics.volume[e]+=w;
            diagnostics.trial[e*PD_NQ+PD_INTERNAL]+=w*internal;
            diagnostics.trial[e*PD_NQ+PD_ADI]+=w*adi[e];
            diagnostics.trial[e*PD_NQ+PD_VISC]+=w*visc[e];
            diagnostics.trial[e*PD_NQ+PD_PHASE]-=w*pressure;
            diagnostics.trial[e*PD_NQ+PD_RAW]+=w*raw[e];
            totals[0]+=w*internal;totals[1]-=w*adi[e];totals[2]+=w*visc[e];totals[3]+=w*pressure;totals[4]+=w*raw[e];
            for(a=1;a<=8;a++) {n=E->ien[1][e].node[a];load[n]+=w*E->N.vpt[GNVINDEX(a,g)]*(internal-adi[e]+visc[e]+pressure);
                if(E->pices.p4)E->pices.element_source[e*8+a-1]+=w*E->N.vpt[GNVINDEX(a,g)]*(internal-adi[e]+visc[e]+pressure);}
        }
    }
    free(adi);free(visc);free(raw);
}

/* Trial updates never modify Tp. Each retry starts from the same frozen
 * coefficients; coefficients are rebuilt for the next accepted substep. */
static double eba_heat_stage(struct All_variables *E,double *g,const double *old,
    double *rhs,double remaining,double ledger[7])
{
    int e,a,b,n,j,attempt,nn=E->lmesh.nno,ng=(E->lmesh.nel+1)*8;
    double ds,maxrate=0,change,totals[5],storage=0,boundary=0;
    double *capacity=array(E,ng),*conductivity=array(E,ng),*fraction=array(E,ng*PHASE_TRANSITIONS);
    double energy_old_mass=energy(E,old);
    assemble(E,old);
    diagnostics.values[PD_REBASE]+=energy(E,old)-energy_old_mass;
    memcpy(capacity,E->pices.gp_capacity,ng*sizeof(double));memcpy(conductivity,E->pices.gp_k,ng*sizeof(double));
    memcpy(fraction,E->pices.gp_fraction,ng*PHASE_TRANSITIONS*sizeof(double));
    heat_load(E,old,rhs,totals);
    for(e=1;e<=E->lmesh.nel;e++) for(a=1;a<=8;a++) for(b=1;b<=8;b++)
        rhs[E->ien[1][e].node[a]]-=E->pices.K[e*64+(a-1)*8+b-1]*old[E->ien[1][e].node[b]];
    exchange(E,rhs);
    for(n=1;n<=nn;n++) if(!fixed(E,n)) {
        double rate=fabs(rhs[n]/E->pices.mass[n]);
        if(!isfinite(rate)) pices_fail(E,"nonfinite heat/source rate");
        if(rate>maxrate)maxrate=rate;
    }
    maxrate=reduce(E,maxrate,MPI_MAX);ds=fmin(remaining,E->pices.dt_heat);
    if(maxrate>0)ds=fmin(ds,.01/maxrate);
    for(attempt=0;attempt<50;attempt++) {
        change=0;
        for(n=1;n<=nn;n++) {
            g[n]=fixed(E,n)?old[n]:old[n]+ds*rhs[n]/E->pices.mass[n];
            if(!pices_temperature_valid(E,g[n]))change=DBL_MAX;
        }
        change=reduce(E,change,MPI_MAX);
        if(change==0) {
            assemble(E,g);
            for(j=8;j<ng;j++) {
                double c=fabs(E->pices.gp_capacity[j]/capacity[j]-1)/.1;
                double k=fabs(E->pices.gp_k[j]/conductivity[j]-1)/.1;
                change=fmax(change,fmax(c,k));
            }
            for(j=8*PHASE_TRANSITIONS;j<ng*PHASE_TRANSITIONS;j++)
                change=fmax(change,fabs(E->pices.gp_fraction[j]-fraction[j])/.05);
            change=reduce(E,change,MPI_MAX);
        }
        /* Restore frozen M and kappa for update accounting and Tp relaxation. */
        assemble(E,old);
        if(change<=1)break;
        ds*=.5;
    }
    if(attempt==50 || !(ds>0) || remaining-ds==remaining) pices_fail(E,"heat/source adaptation stalled");
    for(e=1;e<=E->lmesh.nel;e++) for(a=1;a<=8;a++) {
        n=E->ien[1][e].node[a];storage+=E->pices.emass[e*8+a-1]*(g[n]-old[n]);
        if(fixed(E,n))boundary-=ds*E->pices.emass[e*8+a-1]/E->pices.mass[n]*rhs[n];
        if(E->pices.p4) {
            double residual=E->pices.element_source[e*8+a-1];
            for(b=1;b<=8;b++)residual-=E->pices.K[e*64+(a-1)*8+b-1]*old[E->ien[1][e].node[b]];
            E->pices.heat_residual[e*8+a-1]+=(ds*residual-E->pices.emass[e*8+a-1]*(g[n]-old[n]))/E->advection.timestep;
        }
    }
    for(e=1;e<=E->lmesh.nel;e++)for(j=0;j<PD_TA;j++)
        diagnostics.q[e*PD_NQ+j]+=ds*diagnostics.trial[e*PD_NQ+j];
    ledger[6]+=ds*totals[4];ledger[0]+=storage;
    for(j=0;j<4;j++)ledger[j+1]+=ds*totals[j];
    /* Boundary reaction closes the frozen-capacity FE ledger. It is not CBF. */
    ledger[5]+=boundary;
    free(capacity);free(conductivity);free(fraction);return ds;
}

static double interp(const int *nodes,const double *w,const double *field)
{ int a; double t=0; for(a=1;a<=8;a++) t+=w[a]*field[nodes[a]]; return t; }

/* P maps either absolute Tp or a supplied, nonpersistent particle increment. */
static void project(struct All_variables *E,const double *values,double *out,int boundary)
{
    int p,a,n,nodes[9],e; double w[9];
    if(E->pices.consistent_projection){pices_consistent_project(E,values,out,boundary);return;}
    double *den=array(E,E->lmesh.nno+1);
    /* TA may project in-place onto E->T: preserve prescribed values before
     * clearing the accumulation buffer. */
    double *boundary_values=E->T[1],*saved_boundary=NULL;
    if(boundary && out==boundary_values) {
        saved_boundary=array(E,E->lmesh.nno+1);
        memcpy(saved_boundary,boundary_values,(E->lmesh.nno+1)*sizeof(double));
        boundary_values=saved_boundary;
    }
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
        if(boundary && fixed(E,n)) out[n]=boundary_values[n];
    }
    if(boundary) fprintf(E->fp,"PICES_COVERAGE step=%d empty_elements=%d zero_boundary_nodes=%d\n",E->monitor.solution_cycles,empty,zero);
    free(saved_boundary); free(den); free(covered);
}

static double energy(struct All_variables *E,const double *field)
{
    int e,a;double s=0;
    for(e=1;e<=E->lmesh.nel;e++) for(a=1;a<=8;a++)
        s+=E->pices.emass[e*8+a-1]*field[E->ien[1][e].node[a]];
    return reduce(E,s,MPI_SUM);
}
/* Convex interpolation of the nodal reaction operator, not of a grid
 * temperature increment. For constant alpha this is exactly the same
 * exponential relaxation applied to Tp. No particle temperature clipping. */
static double ta_particle_relax(const int *nodes,const double *w,
                                const double *alpha,const double *target,
                                double old)
{
    int a;double fraction=0,source=0;
    for(a=1;a<=8;a++) {
        double f=w[a]*alpha[nodes[a]];
        fraction+=f;source+=f*target[nodes[a]];
    }
    return (1.0-fraction)*old+source;
}

static void pices_assimilate(struct All_variables *E,double *g,double dt)
{
    int n,p,nodes[9],nn=E->lmesh.nno;double w[9],maximum=0,active=0;
    double *before=array(E,nn+1),*projected=array(E,nn+1);
    double *alpha=array(E,nn+1),*target=array(E,nn+1);
    double *save=E->T[1],time=E->monitor.elapsed_time,stored=0,reconciliation=0;
    memcpy(before,g,(nn+1)*sizeof(double));
    double before_rebase=energy(E,g);
    assemble(E,g); /* same frozen heat capacity for both diagnostic ledgers */
    diagnostics.values[PD_REBASE]+=energy(E,g)-before_rebase;
    E->T[1]=g;E->monitor.elapsed_time=time+dt;
    lith_age_conform_tbc(E);
    lith_age_temperature_bound_adj(E,E->mesh.levmax);
    temperatures_conform_bcs(E);
    memset(E->assim_delta_T[1],0,(nn+1)*sizeof(double));
    for(n=1;n<=nn;n++) {
        lith_age_relaxation_node(E,1,n,&alpha[n],&target[n]);
        if(!isfinite(alpha[n]) || alpha[n]<0 || alpha[n]>1 ||
           !pices_temperature_valid(E,target[n]))
            pices_fail(E,"invalid particle TA target/fraction");
        active=fmax(active,alpha[n]);
    }
    active=reduce(E,active,MPI_MAX);
    if(active>0) {
        /* Separate the pre-existing particle/grid reconciliation from TA.
         * TA storage below includes only the projection change caused by TA.
         * Use the configured absolute projection on both sides of TA. */
        project(E,E->trace.extraq[1][E->pices.slot],projected,1);
        for(n=1;n<=nn;n++)E->assim_delta_T[1][n]=projected[n]-before[n];
        reconciliation=energy(E,E->assim_delta_T[1]);
        for(p=1;p<=E->trace.ntracers[1];p++) {
            double old=E->trace.extraq[1][E->pices.slot][p],after;
            tracer_temperature_weights(E,1,p,nodes,w);
            after=ta_particle_relax(nodes,w,alpha,target,old);
            if(!pices_temperature_valid(E,after)) {
                fprintf(stderr,"PICES_TA_INVALID rank=%d step=%d particle=%d before_nd=%.17g after_nd=%.17g after_K=%.17g operator=particle_relaxation_v1\n",
                    E->parallel.me,E->monitor.solution_cycles,p,old,after,E->data.Ttop+after*E->data.ref_temperature);
                pices_fail(E,"invalid TA particle temperature");
            }
            E->trace.extraq[1][E->pices.slot][p]=after;
        }
        /* Use one identical projection boundary/initial state for both maps. */
        project(E,E->trace.extraq[1][E->pices.slot],g,1);
        for(n=1;n<=nn;n++) {
            E->assim_delta_T[1][n]=g[n]-projected[n];
            maximum=fmax(maximum,fabs(E->assim_delta_T[1][n]));
            if(!pices_temperature_valid(E,g[n]))
                pices_fail(E,"invalid TA grid temperature");
        }
        stored=energy(E,E->assim_delta_T[1]);
    }
    fprintf(E->fp,"PICES_TA step=%d calls=1 time=%.17g dt=%.17g operator=particle_relaxation_v1 stage=%s delta_max=%.17g storage=%.17g reconciliation_storage=%.17g total_grid_storage=%.17g\n",
        E->monitor.solution_cycles,(double)E->monitor.elapsed_time,dt,
        E->pices.ta_before_heat?"before_heat":"after_heat",
        reduce(E,maximum,MPI_MAX),stored,reconciliation,stored+reconciliation);
    diagnostics.values[PD_TASTORE]=stored;
    diagnostics.values[PD_RECON]=reconciliation;
    for(int e=1;e<=E->lmesh.nel;e++)for(int a=1;a<=8;a++)
        diagnostics.q[e*PD_NQ+PD_TA]+=E->pices.emass[e*8+a-1]*
            E->assim_delta_T[1][E->ien[1][e].node[a]];
    E->monitor.elapsed_time=time;E->T[1]=save;
    free(before);free(projected);free(alpha);free(target);
}

/* Opt-in sparse audit, at the post-advection composition stage. Sharing a
 * diagnostic node mask includes neighbors across MPI/cap seams without
 * modifying any physical field. Geometry is written once, compositions every
 * observed step, including the first occupied step after an empty episode. */
static void empty_audit(struct All_variables *E,const int *count)
{
    static FILE *out;
    static unsigned char *seen,*was_empty;
    static int last=-1;
    int e,a,d,step=E->monitor.solution_cycles;
    double *mask,active=0;
    if(!getenv("PICES_EMPTY_AUDIT") || step==last)return;
    if(!out) {
        char path[128];
        snprintf(path,sizeof(path),"PICES_Audit/rank_%06d.jsonl",E->parallel.me);
        out=fopen(path,"w");
        seen=calloc(E->lmesh.nel+1,1);was_empty=calloc(E->lmesh.nel+1,1);
        if(!out || !seen || !was_empty)pices_fail(E,"cannot initialize sparse empty audit");
        fprintf(out,"{\"type\":\"meta\",\"schema\":1,\"rank\":%d,\"stage\":\"post_advection_composition\",\"coordinates\":\"cartesian_radius_nd\",\"composition\":\"flavors_0_through_ncomp\",\"ri\":%.17g,\"ro\":%.17g,\"radius_km\":%.17g}\n",E->parallel.me,E->sphere.ri,E->sphere.ro,E->data.radius_km);
    }
    for(e=1;e<=E->lmesh.nel;e++)if(!count[e] || was_empty[e])active++;
    active=reduce(E,active,MPI_SUM);
    fprintf(out,"{\"type\":\"frame\",\"step\":%d,\"time_Ma\":%.17g}\n",step,(double)E->monitor.elapsed_time*E->data.scalet);
    if(active>0) {
        mask=calloc(E->lmesh.nno+1,sizeof(double));
        if(!mask)pices_fail(E,"sparse audit mask allocation failed");
        for(e=1;e<=E->lmesh.nel;e++)if(!count[e] || was_empty[e])
            for(a=1;a<=8;a++)mask[E->ien[1][e].node[a]]=1;
        exchange(E,mask);
        for(e=1;e<=E->lmesh.nel;e++) {
            int selected=0;double sum=0;
            for(a=1;a<=8;a++)if(mask[E->ien[1][e].node[a]]>0)selected=1;
            if(!selected)continue;
            if(!seen[e]) {
                fprintf(out,"{\"type\":\"mesh\",\"element\":%d,\"volume_nd\":%.17g,\"xyz\":[",e,E->eco[1][e].area);
                for(a=1;a<=8;a++)for(d=1;d<=3;d++)
                    fprintf(out,"%s%.17g",a==1 && d==1?"":",",E->x[1][d][E->ien[1][e].node[a]]);
                fprintf(out,"]}\n");seen[e]=1;
            }
            for(a=0;a<E->composition.ncomp;a++)sum+=E->composition.comp_el[1][a][e];
            fprintf(out,"{\"type\":\"cell\",\"step\":%d,\"element\":%d,\"count\":%d,\"c\":[%.17g",step,e,count[e],1-sum);
            for(a=0;a<E->composition.ncomp;a++)fprintf(out,",%.17g",E->composition.comp_el[1][a][e]);
            fprintf(out,"]}\n");
        }
        free(mask);
    }
    for(e=1;e<=E->lmesh.nel;e++)was_empty[e]=!count[e];
    if(fflush(out) || ferror(out))pices_fail(E,"sparse empty audit write failed");
    last=step;
}

/* Read-only coverage snapshots around movement. Never relocate, replenish,
 * or change composition while diagnosing an empty element. */
static int *coverage_before;
static int coverage_before_valid, coverage_before_particles;
static void capture_particle_coverage(struct All_variables *E)
{
    int p,e;
    if(!E->composition.on)return;
    if(!coverage_before)coverage_before=calloc(E->lmesh.nel+1,sizeof(int));
    if(!coverage_before)pices_fail(E,"coverage diagnostic allocation failed");
    memset(coverage_before,0,(E->lmesh.nel+1)*sizeof(int));
    for(p=1;p<=E->trace.ntracers[1];p++) {
        e=E->trace.ielement[1][p];
        if(e>0 && e<=E->lmesh.nel)coverage_before[e]++;
    }
    coverage_before_particles=E->trace.ntracers[1];coverage_before_valid=1;
}

/* Counters describe this process lifetime, not checkpoint-persisted history. */
/* Optional global histogram and one-step transition census, not a sampler.
 * Restart seeds the previous counts from the restored state; no state edits. */
static void coverage_statistics(struct All_variables *E,const int *count)
{
    static int *previous,last=-1;
    int e,i,step=E->monitor.solution_cycles,valid;
    long long local[PICES_COVERAGE_STATS_SIZE]={0},global[PICES_COVERAGE_STATS_SIZE];
    if(!getenv("PICES_COVERAGE_STATS") || step==last)return;
    valid=previous && step==last+1;
    if(!previous)previous=calloc(E->lmesh.nel+1,sizeof(int));
    if(!previous)pices_fail(E,"coverage statistics allocation failed");
    for(e=1;e<=E->lmesh.nel;e++) {
        pices_coverage_count(local,count[e],previous[e],valid);
        previous[e]=count[e];
    }
    MPI_Allreduce(local,global,PICES_COVERAGE_STATS_SIZE,MPI_LONG_LONG_INT,MPI_SUM,E->parallel.world);
    if(E->parallel.me==0) {
        fprintf(E->fp,"PICES_COVERAGE_STATS step=%d stage=post_advection_composition transition_valid=%d bins=0,1,2,3,4,5,6,7,8,9plus histogram=",step,valid);
        for(i=0;i<10;i++)fprintf(E->fp,"%s%lld",i?",":"",global[i]);
        fprintf(E->fp," transitions=");
        for(i=0;i<100;i++)fprintf(E->fp,"%s%lld",i?",":"",global[10+i]);
        for(i=0;i<3;i++)fprintf(E->fp," low%d=%lld add%d=%lld",2<<i,global[110+2*i],2<<i,global[111+2*i]);
        fprintf(E->fp,"\n");
    }
    last=step;
}

static void report_empty_elements(struct All_variables *E)
{
    static int *streak,last_step=-1,peak;
    static long long fallback_total;
    int e,p,i,step=E->monitor.solution_cycles,longest=0;
    long long local[5]={0},global[5];
    double volume=0,domain=0,global_volume,global_domain,global_longest;
    int *direct=calloc(E->lmesh.nel+1,sizeof(int));
    if(!streak)streak=calloc(E->lmesh.nel+1,sizeof(int));
    if(!direct || !streak)pices_fail(E,"coverage monitor allocation failed");
    for(p=1;p<=E->trace.ntracers[1];p++) {
        e=E->trace.ielement[1][p];
        if(e<1 || e>E->lmesh.nel)local[2]++;else direct[e]++;
    }
    for(e=1;e<=E->lmesh.nel;e++) {
        int count=0;double v=E->eco[1][e].area;
        for(i=0;i<E->trace.nflavors;i++)count+=E->trace.ntracer_flavor[1][i][e];
        if(count!=direct[e])local[3]++;
        if(!isfinite(v) || v<=0)pices_fail(E,"invalid coverage element volume");
        domain+=v;
        if(!direct[e]) {local[4]++;volume+=v;}
        if(step!=last_step) {
            if(!direct[e]) {
                if(streak[e]==0) {
                    double x=0,y=0,z=0,rad=0;
                    for(int a=1;a<=8;a++){int n=E->ien[1][e].node[a];x+=E->x[1][1][n];y+=E->x[1][2][n];z+=E->x[1][3][n];rad+=E->sx[1][3][n];}
                    fprintf(E->trace.fpt,"PICES_EMPTY_EVENT step=%d element=%d event=begin lon_deg=%.17g lat_deg=%.17g depth_km=%.17g\n",step,e,atan2(y,x)*180/M_PI,atan2(z,hypot(x,y))*180/M_PI,(E->sphere.ro-rad/8)*E->data.radius_km);
                }
                streak[e]++;
            } else if(streak[e]) {
                fprintf(E->trace.fpt,"PICES_EMPTY_EVENT step=%d element=%d event=recovered empty_steps=%d count=%d\n",step,e,streak[e],direct[e]);streak[e]=0;
            }
        }
        if(streak[e]>longest)longest=streak[e];
    }
    local[0]=coverage_before_valid?coverage_before_particles:E->trace.ntracers[1];
    local[1]=E->trace.ntracers[1];
    MPI_Allreduce(local,global,5,MPI_LONG_LONG_INT,MPI_SUM,E->parallel.world);
    if(global[0]!=global[1] || global[2] || global[3])pices_fail(E,"particle conservation/host/count failure");
    global_volume=reduce(E,volume,MPI_SUM);global_domain=reduce(E,domain,MPI_SUM);
    global_longest=reduce(E,longest,MPI_MAX);if(global_longest>peak)peak=(int)global_longest;
    if(step!=last_step)fallback_total+=global[4];
    if(E->parallel.me==0) {
        fprintf(E->fp,"PICES_PARTICLE_BALANCE step=%d before_valid=%d before=%lld after=%lld delta=%lld invalid_hosts=%lld count_mismatch_elements=%lld empty_elements=%lld\n",step,coverage_before_valid,global[0],global[1],global[1]-global[0],global[2],global[3],global[4]);
        fprintf(E->fp,"PICES_EMPTY_MONITOR step=%d empty=%lld volume_nd=%.17g volume_fraction=%.17g longest_current=%d longest_since_start=%d fallback_element_steps=%lld policy=retain_previous_composition\n",step,global[4],global_volume,global_volume/global_domain,(int)global_longest,peak,fallback_total);
        fflush(E->fp);
    }
    if(E->parallel.me==161 && E->lmesh.nel>=27151)
        fprintf(E->trace.fpt,"PICES_EMPTY_TARGET step=%d element=27151 count=%d streak=%d\n",step,direct[27151],streak[27151]);
    if(getenv("PICES_COVERAGE_PROBE")) {
        extern void pices_coverage_probe(struct All_variables *,int *,int,int);
        pices_coverage_probe(E,direct,coverage_before_valid,global[4]>0);
    }
    coverage_statistics(E,direct);
    empty_audit(E,direct);
    fflush(E->trace.fpt);last_step=step;free(direct);coverage_before_valid=0;
}

static void composition_check(struct All_variables *E)
{
    int e,i,p,total,primordial=0,global_primordial;
    double error=0,global_error,cmin=1,cmax=0,global_min,global_max;
    if(!E->composition.on)return;
    report_empty_elements(E);
    for(p=1;p<=E->trace.ntracers[1];p++)
        if(E->trace.extraq[1][0][p]==E->control.kC_primordial_flavor)primordial++;
    for(e=1;e<=E->lmesh.nel;e++) {
        total=0;
        for(i=0;i<E->trace.nflavors;i++)total+=E->trace.ntracer_flavor[1][i][e];
        /* Ratio method retains the last valid composition when empty. */
        for(i=0;i<E->composition.ncomp;i++) {
            double c=E->composition.comp_el[1][i][e];
            if(!isfinite(c) || c<0 || c>1)pices_fail(E,"invalid composition fraction");
            if(total>0)error=fmax(error,fabs(c-E->trace.ntracer_flavor[1][i+1][e]/(double)total));
        }
        if(E->control.kC_primordial_flavor>0 && E->control.kC_primordial_flavor<=E->composition.ncomp) {
            double c=E->composition.comp_el[1][E->control.kC_primordial_flavor-1][e];
            cmin=fmin(cmin,c);cmax=fmax(cmax,c);
        }
    }
    MPI_Allreduce(&primordial,&global_primordial,1,MPI_INT,MPI_SUM,E->parallel.world);
    global_error=reduce(E,error,MPI_MAX);
    global_min=reduce(E,cmin,MPI_MIN);global_max=reduce(E,cmax,MPI_MAX);
    if(global_error>1e-12)pices_fail(E,"composition/count mismatch");
    fprintf(E->fp,"PICES_COMPOSITION step=%d primordial=%d error=%.17g Cmin=%.17g Cmax=%.17g kC_ratio=%.17g\n",
        E->monitor.solution_cycles,global_primordial,global_error,global_min,global_max,E->control.kC_ratio);
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
    p5_particle_initial(E);
    composition_check(E);
    assemble(E,E->T[1]); E->pices.initialized=1;
    fprintf(E->fp,"PICES_INIT method=P1_v1 Tp_slot=%d ntracers=%d kappa=%.17g dt_heat=%.17g min_edge=%.17g no_diffusion=%d checkpoint=%d eba=%d\n",E->pices.slot,E->trace.ntracers[1],E->pices.kappa,E->pices.dt_heat,E->pices.min_edge,E->pices.no_diffusion,E->pices.checkpoint,E->pices.eba);
    fflush(E->fp);
}

/* Bound the outer step; heat/source subcycling has its own stability limit.
 * Existing clocks are float: use an explicit resolution bound at forcing knots. */
int pices_time_finished(struct All_variables *E)
{
    double age=E->control.start_age-(double)E->monitor.elapsed_time*E->data.scalet;
    return E->pices.enabled && E->advection.fixed_timestep==0 &&
        (E->control.lith_age || E->control.vbcs_file) &&
        age<=4*FLT_EPSILON*fmax(1.0,fabs(E->control.start_age));
}

void pices_limit_timestep(struct All_variables *E)
{
    int n,d;
    double speed=0,limit=E->advection.timestep,candidate=limit,particle,age=-1,knot=0;
    if(!E->pices.initialized || !(E->data.scalet>0))pices_fail(E,"automatic step before PICES setup");
    for(n=1;n<=E->lmesh.nno;n++) {
        double v=0;
        for(d=1;d<=3;d++)v+=E->sphere.cap[1].V[d][n]*E->sphere.cap[1].V[d][n];
        if(!isfinite(v))pices_fail(E,"nonfinite automatic-step velocity");
        speed=fmax(speed,sqrt(v));
    }
    speed=reduce(E,speed,MPI_MAX);
    particle=speed>0?.24*E->pices.min_edge/speed:DBL_MAX;
    limit=fmin(limit,fmin(particle,E->pices.max_timestep_Ma/E->data.scalet));
    if(E->control.lith_age || E->control.vbcs_file) {
        double resolution=4*FLT_EPSILON*fmax(1.0,fabs(E->control.start_age));
        age=E->control.start_age-(double)E->monitor.elapsed_time*E->data.scalet;
        if(pices_time_finished(E))pices_fail(E,"automatic step requested past present day; stop the time controller");
        knot=age-floor(age);
        if(knot<=resolution)knot=fmin(age,1.0);
        limit=fmin(limit,fmin(age,knot)/E->data.scalet);
    }
    if(!(limit>0) || !isfinite(limit))pices_fail(E,"invalid automatic timestep");
    E->advection.timestep=(float)limit;
    if(E->advection.timestep>limit)E->advection.timestep=nextafterf(E->advection.timestep,0.f);
    if(E->monitor.elapsed_time+E->advection.timestep==E->monitor.elapsed_time)
        pices_fail(E,"automatic timestep below clock resolution");
    fprintf(E->fp,"PICES_TIMESTEP step=%d candidate=%.17g dt=%.17g particle_limit=%.17g age_Ma=%.17g knot_Ma=%.17g\n",
        E->monitor.solution_cycles,candidate,(double)E->advection.timestep,particle,age,knot);
}

/* Read-only operator-stage minima. Locations use the owning rank and local
 * index; no persistent particle ID is implied. K units match profile output. */
static void stage_minima(struct All_variables *E,const double *g,const char *stage,int substep)
{
    struct {double value;int rank;} local,global;
    int kind,i,index=0,nn=E->lmesh.nno;
    double loc[4];
    if(!E->pices.ta_before_heat)return;
    for(kind=0;kind<2;kind++) {
        local.value=DBL_MAX;local.rank=E->parallel.me;
        for(i=1;i<=(kind?E->trace.ntracers[1]:nn);i++) {
            double v=kind?E->trace.extraq[1][E->pices.slot][i]:g[i];
            if(v<local.value){local.value=v;index=i;}
        }
        MPI_Allreduce(&local,&global,1,MPI_DOUBLE_INT,MPI_MINLOC,E->parallel.world);
        if(E->parallel.me==global.rank) {
            loc[0]=index;
            for(i=0;i<3;i++)loc[i+1]=kind?E->trace.basicq[1][i][index]:E->sx[1][i+1][index];
        }
        MPI_Bcast(loc,4,MPI_DOUBLE,global.rank,E->parallel.world);
        if(E->parallel.me==0)fprintf(E->fp,
            "PICES_STAGE step=%d stage=%s substep=%d field=%s Tmin_K=%.17g rank=%d index=%.0f lon_deg=%.17g lat_deg=%.17g depth_km=%.17g\n",
            E->monitor.solution_cycles,stage,substep,kind?"particle":"grid",
            E->data.Ttop+E->data.ref_temperature*global.value,global.rank,loc[0],
            loc[2]*180/M_PI,90-loc[1]*180/M_PI,(E->sphere.ro-loc[3])*E->data.radius_km);
    }
}

void pices_advance(struct All_variables *E)
{
    int n,e,a,b,p,s,ns,nodes[9];
    int nn=E->lmesh.nno,np;
    double dt=E->advection.timestep,speed=0,cfl,ds,w[9],remaining;
    double *g,*old,*dg,*sub,*mapped,*rhs,*previous;
    double ledger[7]={0,0,0,0,0,0,0};
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
    if(!E->pices.eba && dt/E->pices.dt_heat>E->pices.max_substeps) pices_fail(E,"heat substep limit exceeded before particle movement");
    ns=E->pices.eba?1:(int)ceil(dt/E->pices.dt_heat);if(ns<1) ns=1;
    g=array(E,nn+1);old=array(E,nn+1);dg=array(E,nn+1);mapped=array(E,nn+1);rhs=array(E,nn+1);previous=array(E,nn+1);
    for(n=1;n<=nn;n++) previous[n]=E->T[1][n];
    if(E->pices.eba)assemble(E,previous);
    before=energy(E,previous);
    diagnostics.valid=0;diagnostics.step=E->monitor.solution_cycles;
    memset(diagnostics.values,0,sizeof(diagnostics.values));
    if(!diagnostics.q) {
        diagnostics.q=array(E,(E->lmesh.nel+1)*PD_NQ);
        diagnostics.trial=array(E,(E->lmesh.nel+1)*PD_NQ);
        diagnostics.volume=array(E,E->lmesh.nel+1);
    }
    memset(diagnostics.q,0,(E->lmesh.nel+1)*PD_NQ*sizeof(double));

    capture_particle_coverage(E);
    tracer_move_particles(E);
    /* Reclassify once at the arrival age, before conductivity is evaluated.
     * Like TA, this temporary clock view does not commit the accepted step. */
    {
        double time=E->monitor.elapsed_time;
        E->monitor.elapsed_time=time+dt;
        if(E->trace.reclassify_flavors) lith_age_conform_tbc(E);
        tracer_update_composition(E);
        composition_check(E);
        E->monitor.elapsed_time=time;
    }
    np=E->trace.ntracers[1];sub=array(E,np+1);
    temperatures_conform_bcs(E);
    project(E,E->trace.extraq[1][E->pices.slot],g,1);
    advected=energy(E,g);remaining=dt;
    stage_minima(E,g,"after_remap",-1);
    if(E->pices.p4 && E->pices.ta_before_heat && E->control.lith_age) {
        pices_assimilate(E,g,dt);
        stage_minima(E,g,"after_ta_before_heat",-1);
    }
    if(E->pices.p4) {E->pices.cbf_valid=0;memset(E->pices.heat_residual,0,(E->lmesh.nel+1)*8*sizeof(double));}
    for(s=0;E->pices.eba ? remaining>0 : s<ns;s++) {
        if(s>=E->pices.max_substeps)pices_fail(E,"heat substep limit exceeded");
        memcpy(old,g,(nn+1)*sizeof(double));memset(rhs,0,(nn+1)*sizeof(double));
        if(E->pices.eba) {ds=eba_heat_stage(E,g,old,rhs,remaining,ledger);remaining-=ds;}
        else {
        ds=(s==ns-1)?remaining:dt/ns;remaining-=ds;
        for(e=1;e<=E->lmesh.nel;e++) for(a=1;a<=8;a++) for(b=1;b<=8;b++)
            rhs[E->ien[1][e].node[a]]-=E->pices.K[e*64+(a-1)*8+b-1]*old[E->ien[1][e].node[b]];
        exchange(E,rhs);
        for(n=1;n<=nn;n++) {
            if(!fixed(E,n)) g[n]=old[n]+ds*rhs[n]/E->pices.mass[n];
            dg[n]=g[n]-old[n];
            if(!pices_temperature_valid(E,g[n])) pices_fail(E,"invalid grid temperature");
        }
        }
        stage_minima(E,g,"heat_before_particle",s);
        for(n=1;n<=nn;n++)dg[n]=g[n]-old[n];
        for(p=1;p<=np;p++) {
            double factor,tau,kappa;
            tracer_temperature_weights(E,1,p,nodes,w);e=E->trace.ielement[1][p];
            kappa=E->pices.eba?E->pices.ekappa[e]:E->pices.kappa;
            tau=kappa>0 ? E->pices.length[e]*E->pices.length[e]/kappa : DBL_MAX;
            factor=kappa>0 ? -expm1(-ds/tau):0;
            sub[p]=(interp(nodes,w,g)-E->trace.extraq[1][E->pices.slot][p])*factor;
            if(fabs(sub[p])>submax)submax=fabs(sub[p]);
        }
        project(E,sub,mapped,0);
        for(n=1;n<=nn;n++) mapped[n]=dg[n]-mapped[n];
        for(p=1;p<=np;p++) {
            double *tp=&E->trace.extraq[1][E->pices.slot][p];
            tracer_temperature_weights(E,1,p,nodes,w);
            *tp+=sub[p]+interp(nodes,w,mapped);
            if(!pices_temperature_valid(E,*tp)) pices_fail(E,"invalid particle temperature");
        }
        stage_minima(E,g,"after_particle_heat",s);
    }
    if(E->pices.eba) {
        ns=s;for(a=0;a<7;a++)ledger[a]=reduce(E,ledger[a],MPI_SUM);
        fprintf(E->fp,"PICES_EBA step=%d storage=%.17g internal=%.17g adiabatic=%.17g viscous=%.17g phase_pressure=%.17g boundary_reaction=%.17g\n",E->monitor.solution_cycles,ledger[0],ledger[1],ledger[2],ledger[3],ledger[4],ledger[5]);
    }
    if(E->pices.eba)fprintf(E->fp,"PICES_VISCOUS step=%d mode=%d raw=%.17g applied=%.17g state=accepted_stokes\n",
        E->monitor.solution_cycles,E->control.qvis_mode,ledger[6],ledger[3]);
    heated=E->pices.eba?advected+ledger[0]:energy(E,g);
    if(E->pices.p4) {
        E->pices.cbf_valid=1;
        if(E->control.lith_age && !E->pices.ta_before_heat)pices_assimilate(E,g,dt);
    }
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
    if(E->pices.eba) {
        double *v=diagnostics.values,lo=DBL_MAX,hi=-DBL_MAX,sum=0;
        const double length=E->data.radius_km*1000.0;
        v[PD_DT]=dt;v[PD_TIME]=E->monitor.elapsed_time;
        v[PD_BEFORE]=before;v[PD_AFTER]=energy(E,g);v[PD_REMAP]=advected-before;
        v[PD_HEAT]=ledger[0];v[PD_BOUNDARY]=ledger[5];
        v[PD_QINTERNAL]=ledger[1];v[PD_QADI]=-ledger[2];v[PD_QVISC]=ledger[3];
        v[PD_QPHASE]=-ledger[4];v[PD_QRAW]=ledger[6];
        v[PD_HEAT_RESIDUAL]=ledger[0]-ledger[1]-ledger[2]-ledger[3]-ledger[4]-ledger[5];
        v[PD_RESIDUAL]=v[PD_AFTER]-before-v[PD_REMAP]-v[PD_REBASE]
            -ledger[1]-ledger[2]-ledger[3]-ledger[4]-ledger[5]-v[PD_TASTORE]-v[PD_RECON];
        for(p=1;p<=np;p++) {
            double t=E->data.Ttop+E->data.ref_temperature*E->trace.extraq[1][E->pices.slot][p];
            lo=fmin(lo,t);hi=fmax(hi,t);sum+=t;
        }
        v[PD_PARTICLES]=reduce(E,(double)np,MPI_SUM);
        v[PD_TPMIN]=reduce(E,lo,MPI_MIN);v[PD_TPMAX]=reduce(E,hi,MPI_MAX);
        v[PD_TPMEAN]=reduce(E,sum,MPI_SUM)/v[PD_PARTICLES];
        v[PD_GRIDMIN]=E->data.Ttop+E->data.ref_temperature*reduce(E,tmin,MPI_MIN);
        v[PD_GRIDMAX]=E->data.Ttop+E->data.ref_temperature*reduce(E,tmax,MPI_MAX);
        v[PD_ENERGY_SCALE]=E->data.k0*E->data.ref_temperature*length*E->data.scalet*1e6*365.25*86400;
        v[PD_POWER_DENSITY_SCALE]=E->data.k0*E->data.ref_temperature/(length*length);
        for(e=1;e<=E->lmesh.nel;e++)for(a=0;a<PD_NQ;a++)
            diagnostics.q[e*PD_NQ+a]/=dt*diagnostics.volume[e];
        diagnostics.valid=1;
        if(E->parallel.me==0) {
            fprintf(E->fp,"PICES_DIAGNOSTICS step=%d stage=accepted_step scope=global energy=frozen_capacity_T source=accepted_step_average",diagnostics.step);
            for(a=0;a<PD_NVALUES;a++)fprintf(E->fp," %s=%.17g",pices_diagnostic_names[a],v[a]);
            fprintf(E->fp,"\n");
        }
    }

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
        if(!pices_temperature_valid(E,tp)) pices_fail(E,"invalid restored Tp");
    }
    for(n=1;n<=E->lmesh.nno;n++)
        if(!pices_temperature_valid(E,E->T[1][n])) pices_fail(E,"invalid restored T");
    if(E->pices.p4 && E->control.lith_age) {
        lith_age_conform_tbc(E);
        lith_age_temperature_bound_adj(E,E->mesh.levmax);
        temperatures_conform_bcs(E);
    }
    composition_check(E);
    assemble(E,E->T[1]); E->pices.initialized=1;
    E->advection.total_timesteps=E->monitor.solution_cycles+1;
    E->monitor.T_interior=0;
    for(n=1;n<=E->lmesh.nno;n++) if(E->T[1][n]>E->monitor.T_interior) E->monitor.T_interior=E->T[1][n];
    E->monitor.T_interior=reduce(E,E->monitor.T_interior,MPI_MAX);
    fprintf(E->fp,"PICES_RESTORE step=%d particles=%d Tp_slot=%d time=%.17g\n",E->monitor.solution_cycles,E->trace.ntracers[1],E->pices.slot,(double)E->monitor.elapsed_time);
}
