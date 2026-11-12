/* Explicit transport controls for the EBA ratio-composition volume ledger.
 * Production initialization is used once. No Stokes, heat update, relabeling,
 * volume correction, or production checkpoint is performed in these cases. */
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <stdint.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "composition_related.h"
#include "pices.h"
#include "pices_sha256.h"
#include "pices_volume_audit.h"
#include "pices_volume_control.h"
#include "pices_volume_control_math.h"

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif
#define VC_NF 25
#define VC_VERSION "PICES_VOLUME_CONTROL_V1"
static int reference_slot=-1, running;
static const char *case_names[]={"A_static","B_exact","B_dt","C_half","C_quarter"};
static const int case_factor[]={1,1,1,2,4};
void expand_tracer_arrays(struct All_variables *,int);
void cart_to_sphere(struct All_variables *,double,double,double,double *,double *,double *);
void get_global_shape_fn();

struct Control_metrics {
    long long count[VC_NF],active,active24,particles;
    unsigned long long hash_sum,hash_xor;
    double rms[3],max_error[3],radial_max;
};
struct Control {
    char root[PATH_MAX],snapshot[PATH_MAX],local_hash[65],signature[65],thermal_hash[65];
    struct Pices_control_flow flow;
    int requested_steps,steps;
    double requested_end,end,dt,min_edge,speed_bound,shell,cfl;
    double *volume,baseline_volume[VC_NF];
    struct Control_metrics reference;
};

int pices_volume_control_enabled(struct All_variables *E)
{
    const char *s=getenv("PICES_VOLUME_CONTROL");
    if(!s || !strcmp(s,"0") || !strcmp(s,"off"))return 0;
    if(strcmp(s,"1"))pices_fail(E,"PICES_VOLUME_CONTROL must be 0 or 1");
    return 1;
}
void pices_volume_control_register(struct All_variables *E)
{
    if(!pices_volume_control_enabled(E))return;
    if(E->control.restart || !E->pices.enabled || !E->pices.eba || pices_motion_enabled(E))
        pices_fail(E,"control registration requires fresh EBA PICES and no motion audit");
    if(reference_slot>=0)pices_fail(E,"duplicate control reference registration");
    reference_slot=E->trace.number_of_extra_quantities;
    E->trace.number_of_extra_quantities+=4; /* x0,y0,z0,immutable flavor */
}
void pices_volume_control_check_velocity(struct All_variables *E,const double *v)
{
    double speed;
    if(!running)return;
    speed=sqrt(v[1]*v[1]+v[2]*v[2]+v[3]*v[3]);
    if(!isfinite(speed) || speed>100*E->data.scalev)
        pices_fail(E,"control velocity would trigger production speed clipping");
}
void pices_volume_control_check_bounds(struct All_variables *E,const double raw[3],double x,double y,double z)
{
    if(running && (fabs(raw[0]-x)>1e-12 || fabs(raw[1]-y)>1e-12 || fabs(raw[2]-z)>1e-12))
        pices_fail(E,"control trajectory would trigger shell boundary projection");
}
static FILE *create_file(struct All_variables *E,const char *path)
{
    FILE *f=fopen(path,"wbx");
    if(!f)pices_fail(E,"control output already exists or cannot be created");
    return f;
}
static void directory(struct All_variables *E,const char *path)
{
    if(!E->parallel.me && mkdir(path,0775))pices_fail(E,"control directory exists or cannot be created; use a fresh control run");
    MPI_Barrier(E->parallel.world);
}
static double parameter(struct All_variables *E,const char *name,double fallback)
{
    const char *s=getenv(name);char *end;double x;
    if(!s)return fallback;
    errno=0;x=strtod(s,&end);
    if(errno || *end || !isfinite(x) || x<=0)pices_fail(E,"invalid positive control parameter");
    return x;
}
static void transfer(struct All_variables *E,FILE *f,PicesSHA *sha,void *p,size_t size,size_t n,int writing)
{
    size_t done=writing?fwrite(p,size,n,f):fread(p,size,n,f);
    if(done!=n)pices_fail(E,"control baseline snapshot I/O failed");
    pices_sha_add(sha,p,size*n);
}
/* Native snapshots are local to this executable/run, not restart checkpoints.
 * All populated particle state and the empty-cell composition fallback are
 * restored exactly, even if previous migration resized the particle arrays. */
static void snapshot(struct All_variables *E,struct Control *c,int writing)
{
    FILE *f=writing?create_file(E,c->snapshot):fopen(c->snapshot,"rb");
    PicesSHA sha;char digest[65];int m,q,n,header[8],expected[8];
    if(!f)pices_fail(E,"cannot open common control baseline");
    expected[0]=1;expected[1]=E->parallel.me;expected[2]=E->parallel.nproc;
    expected[3]=E->sphere.caps_per_proc;expected[4]=E->lmesh.nel;expected[5]=E->lmesh.nno;
    expected[6]=E->trace.number_of_basic_quantities;expected[7]=E->trace.number_of_extra_quantities;
    memcpy(header,expected,sizeof(header));pices_sha_init(&sha);
    transfer(E,f,&sha,header,sizeof(int),8,writing);
    if(memcmp(header,expected,sizeof(header)))pices_fail(E,"control snapshot topology mismatch");
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        n=E->trace.ntracers[m];transfer(E,f,&sha,&n,sizeof(int),1,writing);
        if(n<0)pices_fail(E,"invalid control snapshot particle count");
        if(!writing) {
            while(n>=E->trace.max_ntracers[m]-5)expand_tracer_arrays(E,m);
            E->trace.ntracers[m]=n;E->trace.ilater[m]=0;
        }
        for(q=0;q<header[6];q++)transfer(E,f,&sha,E->trace.basicq[m][q]+1,sizeof(double),n,writing);
        for(q=0;q<header[7];q++)transfer(E,f,&sha,E->trace.extraq[m][q]+1,sizeof(double),n,writing);
        transfer(E,f,&sha,E->trace.ielement[m]+1,sizeof(int),n,writing);
        for(q=0;q<E->composition.ncomp;q++) {
            transfer(E,f,&sha,E->composition.comp_el[m][q]+1,sizeof(double),E->lmesh.nel,writing);
            transfer(E,f,&sha,E->composition.comp_node[m][q]+1,sizeof(double),E->lmesh.nno,writing);
        }
        for(q=0;q<VC_NF;q++)transfer(E,f,&sha,E->trace.ntracer_flavor[m][q]+1,sizeof(int),E->lmesh.nel,writing);
        transfer(E,f,&sha,E->T[m]+1,sizeof(double),E->lmesh.nno,writing);
    }
    if(!writing && fgetc(f)!=EOF)pices_fail(E,"trailing data in control baseline");
    if(fclose(f))pices_fail(E,"cannot close control baseline");
    pices_sha_end(&sha,digest);
    if(writing)strcpy(c->local_hash,digest);
    else if(strcmp(c->local_hash,digest))pices_fail(E,"control baseline SHA256 changed");
}
static void common_signature(struct All_variables *E,struct Control *c)
{
    PicesSHA sha;char *all=NULL;
    if(!E->parallel.me){all=malloc((size_t)65*E->parallel.nproc);if(!all)pices_fail(E,"control signature allocation failed");}
    MPI_Gather(c->local_hash,65,MPI_CHAR,all,65,MPI_CHAR,0,E->parallel.world);
    if(!E->parallel.me){pices_sha_init(&sha);pices_sha_add(&sha,all,(size_t)65*E->parallel.nproc);pices_sha_end(&sha,c->signature);free(all);}
    MPI_Bcast(c->signature,65,MPI_CHAR,0,E->parallel.world);
}
static void thermal_signature(struct All_variables *E,char out[65])
{
    PicesSHA sha;int m;pices_sha_init(&sha);
    for(m=1;m<=E->sphere.caps_per_proc;m++)pices_sha_add(&sha,E->T[m]+1,E->lmesh.nno*sizeof(double));
    pices_sha_end(&sha,out);
}
static unsigned long long mix64(unsigned long long x)
{
    x^=x>>30;x*=0xbf58476d1ce4e5b9ULL;x^=x>>27;x*=0x94d049bb133111ebULL;return x^(x>>31);
}
static unsigned long long add_word(unsigned long long hash,double value)
{
    uint64_t word;memcpy(&word,&value,sizeof(word));return mix64(hash^word);
}
/* A permutation-independent double fingerprint follows reference coordinates,
 * original/current flavor and Tp through all MPI transfers. Counts are checked
 * independently for every flavor, not only C24. */
static void metrics(struct All_variables *E,struct Control *c,double time,int static_case,struct Control_metrics *out)
{
    struct Control_metrics local;double sums[3]={0},maxima[4]={0};
    int m,p,d,flavor;unsigned long long hash;double ref[3],pos[3],expected[3];
    memset(&local,0,sizeof(local));memset(out,0,sizeof(*out));
    for(m=1;m<=E->sphere.caps_per_proc;m++)for(p=1;p<=E->trace.ntracers[m];p++) {
        double error=0,r0=0,r=0,w,label=E->trace.extraq[m][0][p];
        if(!isfinite(label) || label<0 || label>=VC_NF || floor(label)!=label || label!=E->trace.extraq[m][reference_slot+3][p])
            pices_fail(E,"control particle flavor changed or is invalid");
        flavor=(int)label;local.count[flavor]++;local.particles++;
        hash=0x9e3779b97f4a7c15ULL;
        for(d=0;d<3;d++) {
            ref[d]=E->trace.extraq[m][reference_slot+d][p];pos[d]=E->trace.basicq[m][d+3][p];
            if(!isfinite(ref[d]) || !isfinite(pos[d]))pices_fail(E,"nonfinite control position");
            hash=add_word(hash,ref[d]);r0+=ref[d]*ref[d];r+=pos[d]*pos[d];
        }
        if(sqrt(r)<E->sphere.ri-1e-12 || sqrt(r)>E->sphere.ro+1e-12)pices_fail(E,"control particle outside spherical shell");
        if(!isfinite(E->trace.extraq[m][E->pices.slot][p]))pices_fail(E,"nonfinite control particle temperature");
        hash=add_word(hash,label);hash=add_word(hash,E->trace.extraq[m][E->pices.slot][p]);
        local.hash_sum+=hash;local.hash_xor^=hash;
        w=pices_control_omega(&c->flow,ref);
        pices_control_exact(&c->flow,ref,static_case?0:time,expected);
        for(d=0;d<3;d++){double delta=pos[d]-expected[d];error+=delta*delta;}
        sums[0]+=error;maxima[0]=fmax(maxima[0],error);
        maxima[3]=fmax(maxima[3],fabs(sqrt(r)-sqrt(r0)));
        if(w>0) {
            local.active++;sums[1]+=error;maxima[1]=fmax(maxima[1],error);
            if(flavor==24){local.active24++;sums[2]+=error;maxima[2]=fmax(maxima[2],error);}
        }
    }
    MPI_Allreduce(local.count,out->count,VC_NF,MPI_LONG_LONG_INT,MPI_SUM,E->parallel.world);
    MPI_Allreduce(&local.active,&out->active,1,MPI_LONG_LONG_INT,MPI_SUM,E->parallel.world);
    MPI_Allreduce(&local.active24,&out->active24,1,MPI_LONG_LONG_INT,MPI_SUM,E->parallel.world);
    MPI_Allreduce(&local.particles,&out->particles,1,MPI_LONG_LONG_INT,MPI_SUM,E->parallel.world);
    MPI_Allreduce(&local.hash_sum,&out->hash_sum,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,E->parallel.world);
    MPI_Allreduce(&local.hash_xor,&out->hash_xor,1,MPI_UNSIGNED_LONG_LONG,MPI_BXOR,E->parallel.world);
    MPI_Allreduce(MPI_IN_PLACE,sums,3,MPI_DOUBLE,MPI_SUM,E->parallel.world);
    MPI_Allreduce(MPI_IN_PLACE,maxima,4,MPI_DOUBLE,MPI_MAX,E->parallel.world);
    out->rms[0]=sqrt(sums[0]/fmax(1,out->particles));out->rms[1]=sqrt(sums[1]/fmax(1,out->active));out->rms[2]=sqrt(sums[2]/fmax(1,out->active24));
    for(d=0;d<3;d++)out->max_error[d]=sqrt(maxima[d]);out->radial_max=maxima[3];
}
static void identity(struct All_variables *E,const struct Control_metrics *a,const struct Control_metrics *b)
{
    if(memcmp(a->count,b->count,sizeof(a->count)) || a->particles!=b->particles ||
       a->active!=b->active || a->active24!=b->active24 || a->hash_sum!=b->hash_sum || a->hash_xor!=b->hash_xor)
        pices_fail(E,"control reference identity/label/temperature/count changed during migration");
}
static void geometry(struct All_variables *E,struct Control *c)
{
    int m,e,g;struct Shape_function GN;struct Shape_function_dx dx;
    struct Shape_function_dA omega;double rtf[4][9],edge=DBL_MAX,shell=0;
    c->volume=calloc((size_t)E->sphere.caps_per_proc*E->lmesh.nel,sizeof(double));
    if(!c->volume)pices_fail(E,"control quadrature allocation failed");
    for(m=1;m<=E->sphere.caps_per_proc;m++)for(e=1;e<=E->lmesh.nel;e++) {
        double v=0;get_global_shape_fn(E,e,&GN,&dx,&omega,0,1,rtf,E->mesh.levmax,m);
        for(g=1;g<=8;g++)v+=omega.vpt[g]*g_point[g].weight[2];
        if(!isfinite(v) || v<=0)pices_fail(E,"invalid control cell volume");
        c->volume[(size_t)(m-1)*E->lmesh.nel+e-1]=v;shell+=v;
        edge=fmin(edge,E->pices.cfl_edge[e]);
    }
    MPI_Allreduce(&shell,&c->shell,1,MPI_DOUBLE,MPI_SUM,E->parallel.world);
    MPI_Allreduce(&edge,&c->min_edge,1,MPI_DOUBLE,MPI_MIN,E->parallel.world);
}
static void integrate(struct All_variables *E,const struct Control *c,double out[VC_NF])
{
    int m,e,f;memset(out,0,VC_NF*sizeof(double));
    for(m=1;m<=E->sphere.caps_per_proc;m++)for(e=1;e<=E->lmesh.nel;e++) {
        double v=c->volume[(size_t)(m-1)*E->lmesh.nel+e-1],background=1;
        for(f=1;f<VC_NF;f++){double x=E->composition.comp_el[m][f-1][e];out[f]+=v*x;background-=x;}
        out[0]+=v*background;
    }
    MPI_Allreduce(MPI_IN_PLACE,out,VC_NF,MPI_DOUBLE,MPI_SUM,E->parallel.world);
}
static void node_velocity(struct All_variables *E,const struct Control *c,int zero)
{
    int m,n,d;
    for(m=1;m<=E->sphere.caps_per_proc;m++)for(n=1;n<=E->lmesh.nno;n++) {
        double x[3],u[3],theta=E->sx[m][1][n],phi=E->sx[m][2][n];
        double st=sin(theta),ct=cos(theta),sp=sin(phi),cp=cos(phi);
        for(d=0;d<3;d++)x[d]=E->x[m][d+1][n];pices_control_velocity(&c->flow,x,u);
        if(zero)u[0]=u[1]=u[2]=0;
        E->sphere.cap[m].V[1][n]=u[0]*ct*cp+u[1]*ct*sp-u[2]*st;
        E->sphere.cap[m].V[2][n]=-u[0]*sp+u[1]*cp;
        E->sphere.cap[m].V[3][n]=u[0]*st*cp+u[1]*st*sp+u[2]*ct;
    }
}
static void exact_move(struct All_variables *E,const struct Control *c,double time)
{
    int m,p,d;double ref[3],x[3];
    for(m=1;m<=E->sphere.caps_per_proc;m++)for(p=1;p<=E->trace.ntracers[m];p++) {
        for(d=0;d<3;d++)ref[d]=E->trace.extraq[m][reference_slot+d][p];
        pices_control_exact(&c->flow,ref,time,x);
        for(d=0;d<3;d++)E->trace.basicq[m][d+3][p]=x[d];
        cart_to_sphere(E,x[0],x[1],x[2],&E->trace.basicq[m][0][p],&E->trace.basicq[m][1][p],&E->trace.basicq[m][2][p]);
    }
    tracer_relocate_particles(E);
}
static void row(FILE *f,const char *name,int step,double time,double end,double dt,const struct Control_metrics *m,double volume,double wall)
{
    fprintf(f,"%s,%d,%.17g,%.17g,%.17g,%lld,%lld,%lld,%lld,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,1,%.17g,%.17g\n",name,step,time,time/end,dt,m->particles,m->count[24],m->active,m->active24,m->rms[0],m->max_error[0],m->rms[1],m->max_error[1],m->rms[2],m->max_error[2],m->radial_max,volume,wall);
    fflush(f);
}
static void metadata(struct All_variables *E,const struct Control *c)
{
    FILE *f;if(E->parallel.me)return;
    f=create_file(E,"volume_control.json");
    fprintf(f,"{\"format_version\":\"%s\",\"baseline_step\":0,\"ranks\":%d,\"initial_signature\":\"%s\",\"requested_end_time_nd\":%.17g,\"end_time_nd\":%.17g,\"requested_base_steps\":%d,\"base_steps\":%d,\"center_nd\":[%.17g,%.17g,%.17g],\"support_radius_nd\":%.17g,\"omega_nd\":%.17g,\"max_angle_rad\":%.17g,\"cfl_limit\":%.17g,\"minimum_edge_nd\":%.17g,\"max_speed_bound_nd\":%.17g,\"shell_volume_nd\":%.17g,\"static_tolerance_nd\":%.17g,\"trajectory_exact_tolerance_nd\":1e-11,\"accept_semantics\":\"completed_control_transport_no_stokes\",\"time_note\":\"Common end time is N times float base dt; half and quarter dt are exact binary subdivisions. Ledger time is legacy float; trajectory time is double.\",\"cases\":[\"A_static\",\"B_exact\",\"B_dt\",\"C_half\",\"C_quarter\"]}\n",VC_VERSION,E->parallel.nproc,c->signature,c->requested_end,c->end,c->requested_steps,c->steps,c->flow.center[0],c->flow.center[1],c->flow.center[2],c->flow.radius,c->flow.omega,c->flow.omega*c->end,c->cfl,c->min_edge,c->speed_bound,c->shell,1e-10*c->shell);
    fclose(f);
}
void pices_volume_control_run(struct All_variables *E)
{
    struct Control c;struct Control_metrics m;int cap,p,d,k,step,n,f;double requested,vol[VC_NF];
    char path[PATH_MAX],thermal[65];FILE *summary=NULL,*trajectory=NULL;double started,time,dt;
    memset(&c,0,sizeof(c));
    if(!pices_volume_control_enabled(E))return;
    if(reference_slot<0 || !E->pices.enabled || !E->pices.eba || E->control.restart ||
       pices_motion_enabled(E) || E->trace.nflavors!=VC_NF || E->composition.ncomp!=24 ||
       !E->composition.on || E->composition.ibuoy_type!=1 || E->control.pseudo_free_surf || E->monitor.solution_cycles!=0)
        pices_fail(E,"volume control requires fresh step-zero EBA PICES ratio composition, 25 flavors, fixed surface and no motion audit");
    if(!getcwd(c.root,sizeof(c.root)))pices_fail(E,"cannot resolve control root");
    directory(E,"Initial");directory(E,"Cases");
    if(snprintf(c.snapshot,sizeof(c.snapshot),"%s/Initial/baseline.%06d.bin",c.root,E->parallel.me)>=(int)sizeof(c.snapshot))pices_fail(E,"control path too long");
    requested=parameter(E,"PICES_CONTROL_BASE_STEPS",100);
    if(requested!=floor(requested) || requested>1000000)pices_fail(E,"invalid control base steps");
    c.requested_steps=(int)requested;c.requested_end=parameter(E,"PICES_CONTROL_END_TIME_ND",4e-6);
    c.flow.center[0]=c.flow.center[1]=0.64/sqrt(2.04);c.flow.center[2]=0.2*c.flow.center[0];
    c.flow.radius=0.08;c.cfl=0.25;
    if(0.64-c.flow.radius<=E->sphere.ri+2*E->trace.box_cushion || 0.64+c.flow.radius>=E->sphere.ro-2*E->trace.box_cushion)
        pices_fail(E,"control vortex support is not strictly inside shell");
    geometry(E,&c);
    if(!isfinite(c.min_edge) || c.min_edge<=0)pices_fail(E,"invalid control CFL geometry");
    /* max |u| occurs in the equatorial plane at s=R/sqrt(5). */
    requested=ceil((16/(25*sqrt(5.0)))*M_PI*c.flow.radius/(c.cfl*c.min_edge));
    if(requested>1000000)pices_fail(E,"control CFL would require more than one million base steps");
    c.steps=max(c.requested_steps,(int)requested);
    c.dt=(double)(float)(c.requested_end/c.steps);c.end=c.dt*c.steps;
    if(!isfinite(c.end) || c.dt<=0)pices_fail(E,"control timestep cannot be represented");
    c.flow.omega=M_PI/c.end;c.speed_bound=(16/(25*sqrt(5.0)))*c.flow.omega*c.flow.radius;
    if(c.speed_bound>=100*E->data.scalev)pices_fail(E,"control flow exceeds production velocity limiter; increase control end time");
    for(cap=1;cap<=E->sphere.caps_per_proc;cap++)for(p=1;p<=E->trace.ntracers[cap];p++) {
        for(d=0;d<3;d++)E->trace.extraq[cap][reference_slot+d][p]=E->trace.basicq[cap][d+3][p];
        E->trace.extraq[cap][reference_slot+3][p]=E->trace.extraq[cap][0][p];
        for(d=6;d<E->trace.number_of_basic_quantities;d++)E->trace.basicq[cap][d][p]=0;
    }
    E->trace.reclassify_flavors=0;thermal_signature(E,c.thermal_hash);
    metrics(E,&c,0,1,&c.reference);
    if(!c.reference.active || !c.reference.active24)pices_fail(E,"control vortex has no active C24 particles; experiment would be uninformative");
    if(c.reference.particles>INT_MAX)pices_fail(E,"production tracer checksum cannot represent total particle count");
    integrate(E,&c,c.baseline_volume);snapshot(E,&c,1);common_signature(E,&c);metadata(E,&c);
    if(setenv("PICES_VOLUME_AUDIT","1",1))pices_fail(E,"cannot enable control ledger");
    if(!E->parallel.me) {
        summary=create_file(E,"cases.csv");
        fprintf(summary,"case,status,steps,dt_nd,end_time_nd,initial_signature,particles,C24_particles,active_particles,active_C24_particles,C24_initial_nd,C24_final_nd,trajectory_rms_nd,trajectory_max_nd,active_rms_nd,active_max_nd,C24_active_rms_nd,C24_active_max_nd,radial_displacement_max_nd,identity_ok,wall_seconds\n");fflush(summary);
    }
    running=1;
    for(k=0;k<5;k++) {
        pices_volume_audit_reset(E);
        if(chdir(c.root))pices_fail(E,"cannot return to control root");
        snprintf(path,sizeof(path),"Cases/%s",case_names[k]);directory(E,path);
        if(chdir(path))pices_fail(E,"cannot enter control case");
        snapshot(E,&c,0);E->trace.ilast_tracer_count=(int)c.reference.particles;
        E->monitor.solution_cycles=0;E->monitor.elapsed_time=0;E->advection.timesteps=0;E->advection.total_timesteps=0;
        dt=c.dt/case_factor[k];n=c.steps*case_factor[k];E->advection.timestep=(float)dt;
        node_velocity(E,&c,k==0);metrics(E,&c,0,k==0,&m);identity(E,&c.reference,&m);
        integrate(E,&c,vol);
        for(f=0;f<VC_NF;f++)if(vol[f]!=c.baseline_volume[f])pices_fail(E,"restored composition baseline differs");
        started=MPI_Wtime();pices_volume_audit_initialize(E);pices_volume_audit_commit(E);
        if(!E->parallel.me) {
            FILE *case_meta=create_file(E,"control_case.json");
            fprintf(case_meta,"{\"case\":\"%s\",\"initial_signature\":\"%s\",\"steps\":%d,\"end_time_nd\":%.17g,\"accept_semantics\":\"completed_control_transport_no_stokes\"}\n",case_names[k],c.signature,n,c.end);fclose(case_meta);
            trajectory=create_file(E,"trajectory.csv");
            fprintf(trajectory,"case,step,time_nd,phase_fraction,dt_nd,particles,C24_particles,active_particles,active_C24_particles,trajectory_rms_nd,trajectory_max_nd,active_rms_nd,active_max_nd,C24_active_rms_nd,C24_active_max_nd,radial_displacement_max_nd,identity_ok,C24_volume_nd,wall_seconds\n");
            row(trajectory,case_names[k],0,0,c.end,dt,&m,vol[24],MPI_Wtime()-started);
        }
        for(step=1;step<=n;step++) {
            time=step*dt;E->monitor.solution_cycles=step;E->monitor.elapsed_time=(float)time;
            E->advection.timesteps=E->advection.total_timesteps=step;
            pices_volume_audit_begin(E);
            if(k==1)exact_move(E,&c,time);else if(k>1)tracer_move_particles(E);
            pices_volume_audit_moved(E);recount_tracers_of_flavors(E);fill_composition(E);
            pices_volume_audit_rebuilt(E);
            metrics(E,&c,time,k==0,&m);identity(E,&c.reference,&m);integrate(E,&c,vol);
            if(k==0) {
                if(m.max_error[0]!=0 || m.radial_max!=0)pices_fail(E,"static control moved particles");
                for(f=0;f<VC_NF;f++)if(!isfinite(vol[f]) || fabs(vol[f]-c.baseline_volume[f])>1e-10*c.shell)pices_fail(E,"static reconstruction control changed species volume");
            }
            if(k==1 && m.max_error[0]>1e-11)pices_fail(E,"exact control trajectory differs from analytic reference");
            pices_volume_audit_commit(E);
            if(!E->parallel.me)row(trajectory,case_names[k],step,time,c.end,dt,&m,vol[24],MPI_Wtime()-started);
        }
        thermal_signature(E,thermal);if(strcmp(thermal,c.thermal_hash))pices_fail(E,"control unexpectedly changed nodal temperature");
        if(!E->parallel.me) {
            fclose(trajectory);trajectory=NULL;
            fprintf(summary,"%s,complete,%d,%.17g,%.17g,%s,%lld,%lld,%lld,%lld,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,1,%.17g\n",case_names[k],n,dt,c.end,c.signature,m.particles,m.count[24],m.active,m.active24,c.baseline_volume[24],vol[24],m.rms[0],m.max_error[0],m.rms[1],m.max_error[1],m.rms[2],m.max_error[2],m.radial_max,MPI_Wtime()-started);fflush(summary);
            fprintf(E->fp,"PICES_VOLUME_CONTROL case=%s complete steps=%d C24=%.17g trajectory_max=%.17g\n",case_names[k],n,vol[24],m.max_error[0]);fflush(E->fp);
        }
        MPI_Barrier(E->parallel.world);
    }
    running=0;pices_volume_audit_reset(E);if(chdir(c.root))pices_fail(E,"cannot return to completed control root");
    if(!E->parallel.me) {
        FILE *fcomplete=create_file(E,"control_complete.json");
        fprintf(fcomplete,"{\"format_version\":\"%s\",\"completed\":true,\"cases\":5,\"initial_signature\":\"%s\",\"A_gate_passed\":true}\n",VC_VERSION,c.signature);fclose(fcomplete);fclose(summary);
    }
    free(c.volume);MPI_Barrier(E->parallel.world);
}
