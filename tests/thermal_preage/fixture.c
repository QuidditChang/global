/* Test-only linker wrapper. Production initialization, operators and parser
 * are linked unchanged. Manufactured composition is explicitly opt-in. */
#include <math.h>
#include <stdint.h>
#include <string.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "pices.h"
#include "material_properties.h"

void __real_thermal_preage_run(struct All_variables *E);
void get_global_shape_fn();

static uint64_t digest(uint64_t h, const void *p, size_t n)
{
    const unsigned char *s=p;
    if(!p) return h;
    while(n--) { h^=*s++; h*=UINT64_C(1099511628211); }
    return h;
}
#define HASH(h,x) ((h)=digest((h),&(x),sizeof(x)))
#define ARRAY(h,p,n) ((h)=digest((h),(p),(n)*sizeof(*(p))))

/* Keep groups separate so mutation failures identify the affected contract. */
static void state_hashes(struct All_variables *E, uint64_t out[8])
{
    int m,n,d,k;
    for(k=0;k<8;k++)out[k]=UINT64_C(14695981039346656037);
    HASH(out[0],E->monitor); HASH(out[0],E->advection);
    HASH(out[1],E->control);
    HASH(out[2],E->trace);
    HASH(out[3],E->composition);
    HASH(out[4],E->pices);
    HASH(out[5],E->refstate);
    ARRAY(out[5],E->refstate.Tref,E->lmesh.noz+1);
    ARRAY(out[5],E->refstate.rho,E->lmesh.noz+1);
    ARRAY(out[5],E->refstate.heat_capacity,E->lmesh.noz+1);
    ARRAY(out[5],E->refstate.gravity,E->lmesh.noz+1);
    ARRAY(out[5],E->refstate.thermal_expansivity,E->lmesh.noz+1);
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        ARRAY(out[6],E->node[m],E->lmesh.nno+1);
        for(d=1;d<=3;d++) {
            ARRAY(out[6],E->sphere.cap[m].TB[d],E->lmesh.nno+1);
            ARRAY(out[7],E->sphere.cap[m].V[d],E->lmesh.nno+1);
        }
        ARRAY(out[7],E->Tdot[m],E->lmesh.nno+1);
        ARRAY(out[7],E->DataT[m],E->lmesh.nno+1);
        ARRAY(out[7],E->P[m],E->lmesh.npno+1);
        ARRAY(out[7],E->assim_delta_T[m],E->lmesh.nno+1);
        ARRAY(out[7],E->heating_adi[m],E->lmesh.nel+1);
        ARRAY(out[7],E->heating_adi_base[m],E->lmesh.nel+1);
        ARRAY(out[7],E->heating_phase[m],E->lmesh.nel+1);
        ARRAY(out[7],E->heating_visc[m],E->lmesh.nel+1);
        ARRAY(out[7],E->heating_visc_raw[m],E->lmesh.nel+1);
        ARRAY(out[7],E->heating_visc_capped[m],E->lmesh.nel+1);
        ARRAY(out[7],E->heating_latent[m],E->lmesh.nel+1);
        ARRAY(out[7],E->heating_internal[m],E->lmesh.nel+1);
        ARRAY(out[7],E->heating_assim[m],E->lmesh.nel+1);
        for(d=0;d<PHASE_TRANSITIONS;d++) {
            ARRAY(out[7],E->phase_B[d][m],E->lmesh.nno+1);
            ARRAY(out[7],E->phase_boundary[d][m],E->lmesh.nsf+1);
        }
        ARRAY(out[2],E->trace.ielement[m],E->trace.ntracers[m]+1);
        for(d=0;d<E->trace.nflavors;d++)ARRAY(out[2],E->trace.ntracer_flavor[m][d],E->lmesh.nel+1);
        for(d=0;d<E->trace.number_of_basic_quantities;d++)
            ARRAY(out[2],E->trace.basicq[m][d],E->trace.ntracers[m]+1);
        for(d=0;d<E->trace.number_of_extra_quantities;d++)
            ARRAY(out[2],E->trace.extraq[m][d],E->trace.ntracers[m]+1);
        for(d=0;d<E->composition.ncomp;d++) {
            ARRAY(out[3],E->composition.comp_el[m][d],E->lmesh.nel+1);
            ARRAY(out[3],E->composition.comp_node[m][d],E->lmesh.nno+1);
        }
    }
    ARRAY(out[7],E->age_t,E->mesh.nox*E->mesh.noy+1);
    ARRAY(out[7],E->flag_depth2,E->mesh.nox*E->mesh.noy+1);
}

static void dump_nodes(struct All_variables *E,const char *stage)
{
    char path[128];FILE *f;int m,n,d;
    snprintf(path,sizeof(path),"nodes.%s.%d.txt",stage,E->parallel.me);
    f=fopen(path,"w");assert(f);
    fprintf(f,"# cap node x y z radial_index radius T DataT Tref Tdot TB1 TB2 TB3 flags C\n");
    for(m=1;m<=E->sphere.caps_per_proc;m++)for(n=1;n<=E->lmesh.nno;n++) {
        int nz=(n-1)%E->lmesh.noz+1;
        fprintf(f,"%d %d %.17g %.17g %.17g %d %.17g %.17g %.17g %.17g %.17g",
            E->sphere.capid[m],n,E->x[m][1][n],E->x[m][2][n],E->x[m][3][n],
            nz+E->lmesh.nzs-1,E->sx[m][3][n],E->T[m][n],E->DataT[m][n],
            E->refstate.Tref[nz],E->Tdot[m][n]);
        for(d=1;d<=3;d++)fprintf(f," %.17g",(double)E->sphere.cap[m].TB[d][n]);
        fprintf(f," %u %.17g\n",E->node[m][n],E->composition.on ? E->composition.comp_node[m][E->control.kC_primordial_flavor-1][n] : 0.);
    }
    fclose(f);
}

static void manufactured_composition(struct All_variables *E)
{
    const char *mode=getenv("PREAGE_TEST_COMPOSITION");int m,p,e,a;
    void init_composition(struct All_variables *);
    if(!mode)return;
    assert(E->composition.on && E->composition.ncomp>=E->control.kC_primordial_flavor);
    /* Elementwise binary hemispheres make actual tracer flavor counts and
     * reconstructed composition consistent and partition-independent. */
    for(m=1;m<=E->sphere.caps_per_proc;m++)for(p=1;p<=E->trace.ntracers[m];p++) {
        double x=0;e=E->trace.ielement[m][p];
        for(a=1;a<=8;a++)x+=E->x[m][1][E->ien[m][e].node[a]]/8.;
        E->trace.extraq[m][0][p]=!strcmp(mode,"uniform") || x>1e-8 ? E->control.kC_primordial_flavor : 0.;
    }
    recount_tracers_of_flavors(E);
    init_composition(E);
}

/* Export geometry, not initializer matrices: Python independently assembles
 * the entire global FE mass/stiffness and solves the backward-Euler system. */
static void dump_quadrature(struct All_variables *E)
{
    struct Shape_function GN;struct Shape_function_dx dx;struct Shape_function_dA omega;
    double rtf[4][9];int e,g,a;char path[128];FILE *f;
    if(!getenv("PREAGE_TEST_QUADRATURE"))return;
    snprintf(path,sizeof(path),"quadrature.%d.txt",E->parallel.me);f=fopen(path,"w");assert(f);
    fprintf(f,"# element gp weight inverse_radius theta then 8 * (node N dtheta dphi dr rho cp) then k_prefactor\n");
    for(e=1;e<=E->lmesh.nel;e++) {
        get_global_shape_fn(E,e,&GN,&dx,&omega,0,1,rtf,E->mesh.levmax,1);
        for(g=1;g<=8;g++) {
            fprintf(f,"%d %d %.17g %.17g %.17g",e,g,omega.vpt[g]*g_point[g].weight[2],rtf[3][g],rtf[1][g]);
            for(a=1;a<=8;a++) {
                int n=E->ien[1][e].node[a],nz=(n-1)%E->lmesh.noz+1;
                fprintf(f," %d %.17g %.17g %.17g %.17g %.17g %.17g",n,E->N.vpt[GNVINDEX(a,g)],
                    dx.vpt[GNVXINDEX(0,a,g)],dx.vpt[GNVXINDEX(1,a,g)],dx.vpt[GNVXINDEX(2,a,g)],
                    E->refstate.rho[nz],E->refstate.heat_capacity[nz]);
            }
            fprintf(f," %.17g\n",conductivity_element_prefactor(E,1,e,E->control.reference_conductivity));
        }
    }
    fclose(f);
}

void __wrap_thermal_preage_run(struct All_variables *E)
{
    static const char *names[8]={"clock","controls","particles","composition","pices","reference","boundary","other_fields"};
    uint64_t before[8],after[8];int k,bad=0;char path[128];FILE *f;
    manufactured_composition(E);
    if(getenv("PREAGE_TEST_BAD_T"))E->T[1][2]=!strcmp(getenv("PREAGE_TEST_BAD_T"),"nan") ? NAN : -1.;
    dump_nodes(E,"before");dump_quadrature(E);state_hashes(E,before);
    __real_thermal_preage_run(E);
    state_hashes(E,after);dump_nodes(E,"preaged");
    snprintf(path,sizeof(path),"isolation.%d.txt",E->parallel.me);f=fopen(path,"w");assert(f);
    for(k=0;k<8;k++) {
        fprintf(f,"%s %016llx %016llx\n",names[k],(unsigned long long)before[k],(unsigned long long)after[k]);
        if(before[k]!=after[k])bad=1;
    }
    fclose(f);
    if(bad) { fprintf(stderr,"PREAGE_TEST_ERROR: production state mutated on rank %d\n",E->parallel.me);MPI_Abort(E->parallel.world,86); }
}

void preage_test_final(struct All_variables *E)
{
    char path[128];FILE *f;int m,p,a,nodes[9];double w[9];
    dump_nodes(E,"final");
    snprintf(path,sizeof(path),"particles.%d.txt",E->parallel.me);f=fopen(path,"w");assert(f);
    for(m=1;m<=E->sphere.caps_per_proc;m++)for(p=1;p<=E->trace.ntracers[m];p++) {
        double expected=0;tracer_temperature_weights(E,m,p,nodes,w);
        for(a=1;a<=8;a++)expected+=w[a]*E->T[m][nodes[a]];
        fprintf(f,"%d %d %.17g %.17g\n",m,p,E->trace.extraq[m][E->pices.slot][p],expected);
    }
    fclose(f);
    snprintf(path,sizeof(path),"metadata.%d.txt",E->parallel.me);f=fopen(path,"w");assert(f);
    fprintf(f,"enabled %d\nduration_Ma %.17g\nmax_dt_Ma %.17g\nscalet %.17g\nclock %.17g\nstep %d\nheat_steps %d\ninitialized %d\nsurface %.17g\nnoz %d\nsurface_temp %.17g\nouter_radius %.17g\n",
        E->control.thermal_preage,E->control.thermal_preage_Ma,E->control.thermal_preage_max_dt_Ma,
        (double)E->data.scalet,(double)E->monitor.elapsed_time,E->monitor.solution_cycles,E->advection.total_timesteps,
        E->pices.initialized,E->refstate.temperature_surface,E->mesh.noz,
        (double)E->control.surface_temp,E->sphere.ro);
    fclose(f);
    snprintf(path,sizeof(path),"phase.%d.txt",E->parallel.me);f=fopen(path,"w");assert(f);
    for(a=0;a<PHASE_TRANSITIONS;a++) {
        const struct Phase_transition *phase=&E->control.phase[a];
        fprintf(f,"%.17g %.17g %.17g %.17g %.17g\n",(double)phase->depth,(double)phase->entropy_jump,
            (double)phase->clapeyron,(double)phase->transT,(double)phase->inv_width);
    }
    fclose(f);fflush(E->fp);
    /* Exercise the production field writer without Stokes/time marching,
       as the separate Pyre initialization-only app will do on the cluster. */
    if(getenv("PREAGE_TEST_WRITE_OUTPUT"))
        (E->problem_output)(E, E->monitor.solution_cycles);
    if(E->parallel.me==0)fprintf(stderr,"PREAGE_TEST_DONE\n");
}
