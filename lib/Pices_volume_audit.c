/* EBA species-volume ledger. This diagnoses the count-ratio reconstruction;
 * Vcell/Ncell is an accounting weight, NOT a persistent particle volume.
 * Owned-element quadrature defines the ledger. Local nodal weights contain
 * only owned-element contributions, so shared MPI nodes are not double counted.
 */
#include <math.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "pices.h"
#include "pices_volume_audit.h"

#define VA_TOP 10
#define VA_BIN_KINDS 4
#define VA_BINS 8
#define VA_METRICS 8
#define VA_RADIAL 8
#define VA_SAMPLE 13
/* comp_el and ledger arithmetic are double; legacy TWW/MASS projection
 * weights are float. Keep their diagnostic roundoff tolerances separate. */
#define VA_ELEMENT_TOL 1e-10
#define VA_NODE_TOL 1e-7

struct Volume_audit {
    int initialized,enabled,active,moved,pending,initial,step,nf,nc,nel,nno,nlayers,interval,detail;
    int baseline_step;
    double time,age,baseline_time,shell,mean_volume,mean_density,started,seconds,bounds_error[2];
    double *volume,*node_weight,*depth,*before24,*moved24,*baseline;
    long long *counts,*total,*old_total;
    long long *nbefore,*nmoved,*nafter,*transitions;
    double *before,*arrival,*after,*node,*transfer,*source,*residual;
    double *radial,*bins,*samples,*gathered;
    FILE *global_file,*transition_file,*radial_file,*bin_file,*sample_file;
};
static struct Volume_audit audit;
void get_global_shape_fn();

static size_t cell(int cap,int e)
{ return (size_t)(cap-1)*audit.nel+(size_t)(e-1); }
static void *allocate(struct All_variables *E,size_t n,size_t size)
{
    void *p=calloc(n,size);
    if(!p)pices_fail(E,"volume audit allocation failed");
    return p;
}
static void sum_double(struct All_variables *E,double *v,int n)
{ MPI_Allreduce(MPI_IN_PLACE,v,n,MPI_DOUBLE,MPI_SUM,E->parallel.world); }
static void sum_integer(struct All_variables *E,long long *v,int n)
{ MPI_Allreduce(MPI_IN_PLACE,v,n,MPI_LONG_LONG_INT,MPI_SUM,E->parallel.world); }
static double element_fraction(struct All_variables *E,int cap,int e,int f)
{
    int k;double c=1;
    if(f)return E->composition.comp_el[cap][f-1][e];
    for(k=0;k<audit.nf-1;k++)c-=E->composition.comp_el[cap][k][e];
    return c;
}
static double nodal_fraction(struct All_variables *E,int cap,int n,int f)
{
    int k;double c=1;
    if(f)return E->composition.comp_node[cap][f-1][n];
    for(k=0;k<audit.nf-1;k++)c-=E->composition.comp_node[cap][k][n];
    return c;
}
static void integrate(struct All_variables *E,double *out,double *nodes,long long *counts)
{
    int m,e,f,n;size_t c;double bounds_error[2]={0,0};
    memset(out,0,audit.nf*sizeof(double));
    if(nodes)memset(nodes,0,audit.nf*sizeof(double));
    if(counts)memset(counts,0,audit.nf*sizeof(long long));
    for(m=1;m<=audit.nc;m++) {
        for(e=1;e<=audit.nel;e++) {
            c=cell(m,e);
            for(f=0;f<audit.nf;f++) {
                double fraction=element_fraction(E,m,e,f);
                bounds_error[0]=fmax(bounds_error[0],isfinite(fraction)?fmax(-fraction,fraction-1):HUGE_VAL);
                out[f]+=audit.volume[c]*fraction;
                if(counts)counts[f]+=E->trace.ntracer_flavor[m][f][e];
            }
        }
        if(nodes)for(n=1;n<=audit.nno;n++)for(f=0;f<audit.nf;f++) {
            double fraction=nodal_fraction(E,m,n,f);
            bounds_error[1]=fmax(bounds_error[1],isfinite(fraction)?fmax(-fraction,fraction-1):HUGE_VAL);
            nodes[f]+=audit.node_weight[(size_t)(m-1)*audit.nno+n-1]*fraction;
        }
    }
    MPI_Allreduce(bounds_error,audit.bounds_error,2,MPI_DOUBLE,MPI_MAX,E->parallel.world);
    sum_double(E,out,audit.nf);
    if(nodes)sum_double(E,nodes,audit.nf);
    if(counts)sum_integer(E,counts,audit.nf);
}
static FILE *open_file(struct All_variables *E,const char *name,const char *header)
{
    char path[256];FILE *f;
    snprintf(path,sizeof(path),"Diagnostic/VolumeAudit/%s",name);
    /* Never append another invocation to an existing baseline silently. */
    f=fopen(path,"wx");
    if(!f)pices_fail(E,"volume audit output exists or cannot be created; choose a fresh run directory");
    fprintf(f,"%s\n",header);fflush(f);return f;
}
static void snapshot_before(struct All_variables *E)
{
    int m,e,f;size_t c;
    integrate(E,audit.before,NULL,audit.nbefore);
    for(m=1;m<=audit.nc;m++)for(e=1;e<=audit.nel;e++) {
        c=cell(m,e);audit.old_total[c]=0;
        for(f=0;f<audit.nf;f++)audit.old_total[c]+=E->trace.ntracer_flavor[m][f][e];
        audit.before24[c]=element_fraction(E,m,e,24);
    }
}
static void reset_step(struct All_variables *E)
{
    size_t matrix=(size_t)audit.nf*audit.nf;
    audit.step=E->monitor.solution_cycles;
    audit.time=E->monitor.elapsed_time;
    audit.age=E->control.start_age-audit.time*E->data.scalet;
    audit.detail=audit.initial || audit.step%audit.interval==0;
    memset(audit.transitions,0,matrix*sizeof(long long));
    memset(audit.transfer,0,matrix*sizeof(double));
    audit.moved=0;audit.pending=0;audit.active=1;
}
void pices_volume_audit_initialize(struct All_variables *E)
{
    int m,e,a,g,f;size_t c,ncell;char *flag,*end;long interval;
    struct Shape_function GN;struct Shape_function_dx dx;
    struct Shape_function_dA omega;double rtf[4][9],local[2]={0,0};
    FILE *meta;
    if(audit.initialized)return;
    audit.initialized=1;flag=getenv("PICES_VOLUME_AUDIT");
    if(!flag || !strcmp(flag,"0") || !strcmp(flag,"off"))return;
    if(strcmp(flag,"1"))pices_fail(E,"PICES_VOLUME_AUDIT must be 0 or 1");
    if(!E->pices.enabled || !E->pices.eba || !E->composition.on || E->composition.ibuoy_type!=1 ||
       E->trace.nflavors!=25 || E->composition.ncomp!=24 || E->control.pseudo_free_surf)
        pices_fail(E,"volume audit requires EBA PICES ratio composition with flavors 0..24 and fixed surface");
    audit.enabled=1;audit.initial=1;audit.started=MPI_Wtime();
    audit.nf=E->trace.nflavors;audit.nc=E->sphere.caps_per_proc;
    audit.nel=E->lmesh.nel;audit.nno=E->lmesh.nno;audit.nlayers=E->mesh.elz;
    audit.interval=50;flag=getenv("PICES_VOLUME_AUDIT_INTERVAL");
    if(flag) {
        errno=0;interval=strtol(flag,&end,10);
        if(errno || *end || interval<1 || interval>1000000)pices_fail(E,"invalid PICES_VOLUME_AUDIT_INTERVAL");
        audit.interval=(int)interval;
    }
    ncell=(size_t)audit.nc*audit.nel;
#define ALLOC(name,n,type) audit.name=allocate(E,(n),sizeof(type))
    ALLOC(volume,ncell,double);ALLOC(depth,ncell,double);
    ALLOC(node_weight,(size_t)audit.nc*audit.nno,double);
    ALLOC(before24,ncell,double);ALLOC(moved24,ncell,double);
    ALLOC(counts,ncell*audit.nf,long long);ALLOC(total,ncell,long long);ALLOC(old_total,ncell,long long);
    ALLOC(nbefore,audit.nf,long long);ALLOC(nmoved,audit.nf,long long);ALLOC(nafter,audit.nf,long long);
    ALLOC(transitions,audit.nf*audit.nf,long long);ALLOC(transfer,audit.nf*audit.nf,double);
    ALLOC(before,audit.nf,double);ALLOC(arrival,audit.nf,double);ALLOC(after,audit.nf,double);
    ALLOC(node,audit.nf,double);ALLOC(source,audit.nf,double);ALLOC(residual,audit.nf,double);ALLOC(baseline,audit.nf,double);
    ALLOC(radial,audit.nlayers*VA_RADIAL,double);ALLOC(bins,VA_BIN_KINDS*VA_BINS*VA_METRICS,double);
    ALLOC(samples,VA_TOP*VA_SAMPLE,double);
    if(E->parallel.me==0)ALLOC(gathered,(size_t)E->parallel.nproc*VA_TOP*VA_SAMPLE,double);
#undef ALLOC
    for(m=1;m<=audit.nc;m++)for(e=1;e<=audit.nel;e++) {
        double rad=0;c=cell(m,e);
        get_global_shape_fn(E,e,&GN,&dx,&omega,0,1,rtf,E->mesh.levmax,m);
        for(g=1;g<=8;g++) {
            double w=omega.vpt[g]*g_point[g].weight[2];audit.volume[c]+=w;
            for(a=1;a<=8;a++)audit.node_weight[(size_t)(m-1)*audit.nno+E->ien[m][e].node[a]-1]+=E->N.vpt[GNVINDEX(a,g)]*w;
        }
        if(!isfinite(audit.volume[c]) || audit.volume[c]<=0)pices_fail(E,"invalid volume audit quadrature");
        for(a=1;a<=8;a++)rad+=E->sx[m][3][E->ien[m][e].node[a]];
        audit.depth[c]=(E->sphere.ro-rad/8)*E->data.radius_km;
        local[0]+=audit.volume[c];local[1]+=1;
    }
    sum_double(E,local,2);audit.shell=local[0];audit.mean_volume=local[0]/local[1];
    reset_step(E);snapshot_before(E);
    audit.baseline_step=audit.step;audit.baseline_time=audit.time;
    memcpy(audit.baseline,audit.before,audit.nf*sizeof(double));
    local[0]=0;for(f=0;f<audit.nf;f++)local[0]+=audit.nbefore[f];
    audit.mean_density=local[0]/audit.shell;
    if(E->parallel.me==0) {
        if(mkdir("Diagnostic",0775) && errno!=EEXIST)pices_fail(E,"cannot create Diagnostic");
        if(mkdir("Diagnostic/VolumeAudit",0775) && errno!=EEXIST)pices_fail(E,"cannot create VolumeAudit");
        audit.global_file=open_file(E,"globals.csv","status,step,time_nd,age_Ma,flavor,count_before,count_moved,count_after,volume_before_nd,volume_moved_nd,volume_after_nd,node_volume_nd,source_nd,residual_nd,baseline_volume_nd,shell_volume_nd");
        audit.transition_file=open_file(E,"transitions.csv","status,step,old_flavor,new_flavor,count,volume_nd");
        audit.radial_file=open_file(E,"radial.csv","status,step,layer,depth_km,shell_volume_nd,C24_before_nd,C24_moved_nd,C24_after_nd,particles,C24_particles,empty_volume_nd");
        audit.bin_file=open_file(E,"bins.csv","status,step,kind,bin,elements,volume_nd,C24_before_nd,C24_moved_nd,C24_after_nd,particles,C24_particles");
        audit.sample_file=open_file(E,"samples.csv","status,step,rank,cap,element,depth_km,volume_nd,count_before,count_moved,C24_before,C24_moved,C24_after,delta_C24_volume_nd");
        meta=open_file(E,"metadata.json", "");
        fprintf(meta,"{\"schema_version\":1,\"format_version\":\"PICES_VOLUME_AUDIT_V1\",\"radial_layers\":%d,\"baseline_kind\":\"diagnostic_start\",\"baseline_step\":%d,\"baseline_time_nd\":%.17g,\"flavors\":%d,\"ranks\":%d,\"interval\":%d,\"radius_km\":%.17g,\"initial_particle_density_nd\":%.17g,\"mean_element_volume_nd\":%.17g,\"element_budget_relative_tolerance\":1e-10,\"element_fraction_tolerance\":1e-10,\"node_relative_tolerance\":1e-7,\"node_fraction_tolerance\":1e-7,\"node_precision_note\":\"Legacy TWW and MASS projection weights are float; raw node differences are retained.\",\"note\":\"No correction. Existing drift before diagnostic_start is not removed. Pending rows require an accepted row to enter evolution curves. Volume transfers use current-cell V/N, not persistent particle weights.\",\"bin_bounds\":{\"volume_ratio\":[0.25,0.5,1,2,4,8,16],\"particle_count\":[1,4,8,16,32,64,128],\"sampling_ratio\":[0.125,0.25,0.5,1,2,4,8]}}\n",audit.nlayers,audit.baseline_step,audit.baseline_time,audit.nf,E->parallel.nproc,audit.interval,(double)E->data.radius_km,audit.mean_density,audit.mean_volume);
        fclose(meta);
    }
    pices_volume_audit_moved(E);
    pices_volume_audit_rebuilt(E);
    audit.seconds=MPI_Wtime()-audit.started;
}
void pices_volume_audit_begin(struct All_variables *E)
{
    if(!audit.enabled)return;
    audit.started=MPI_Wtime();
    if(audit.pending)pices_fail(E,"volume audit previous step lacks successful solve commit");
    audit.initial=0;reset_step(E);snapshot_before(E);
    audit.seconds+=MPI_Wtime()-audit.started;
}
void pices_volume_audit_moved(struct All_variables *E)
{
    int m,e,p,f;size_t c,ncell=(size_t)audit.nc*audit.nel;double started;
    if(!audit.enabled || !audit.active)return;
    started=MPI_Wtime();
    memset(audit.counts,0,ncell*audit.nf*sizeof(long long));memset(audit.total,0,ncell*sizeof(long long));
    memset(audit.nmoved,0,audit.nf*sizeof(long long));memset(audit.arrival,0,audit.nf*sizeof(double));
    for(m=1;m<=audit.nc;m++)for(p=1;p<=E->trace.ntracers[m];p++) {
        double label=E->trace.extraq[m][0][p];e=E->trace.ielement[m][p];
        if(e<1 || e>audit.nel || !isfinite(label) || label<0 || label>=audit.nf || floor(label)!=label)
            pices_fail(E,"invalid host/flavor in volume audit");
        f=(int)label;c=cell(m,e);audit.counts[c*audit.nf+f]++;audit.total[c]++;audit.nmoved[f]++;
    }
    for(m=1;m<=audit.nc;m++)for(e=1;e<=audit.nel;e++) {
        c=cell(m,e);
        for(f=0;f<audit.nf;f++) {
            double fraction=audit.total[c]?(double)audit.counts[c*audit.nf+f]/audit.total[c]:element_fraction(E,m,e,f);
            audit.arrival[f]+=audit.volume[c]*fraction;
            if(f==24)audit.moved24[c]=fraction;
        }
    }
    sum_integer(E,audit.nmoved,audit.nf);sum_double(E,audit.arrival,audit.nf);
    audit.moved=1;audit.seconds+=MPI_Wtime()-started;
}
void pices_volume_audit_transition(struct All_variables *E,int cap,int e,int old,int new)
{
    size_t c,k;
    if(!audit.enabled || !audit.active || !audit.moved || old==new)return;
    if(cap<1 || cap>audit.nc || e<1 || e>audit.nel || old<0 || old>=audit.nf || new<0 || new>=audit.nf)
        pices_fail(E,"invalid volume audit transition");
    c=cell(cap,e);if(!audit.total[c])pices_fail(E,"volume audit transition in empty element");
    k=(size_t)old*audit.nf+new;audit.transitions[k]++;audit.transfer[k]+=audit.volume[c]/audit.total[c];
}
static int choose_bin(double x,const double limits[7])
{ int b=0;while(b<7 && x>=limits[b])b++;return b; }
static void add_bin(int kind,int b,double v,double before,double moved,double after,long long n,long long n24)
{
    double *r=audit.bins+(kind*VA_BINS+b)*VA_METRICS;
    r[0]++;r[1]+=v;r[2]+=v*before;r[3]+=v*moved;r[4]+=v*after;r[5]+=n;r[6]+=n24;
}
static void detail_statistics(struct All_variables *E)
{
    const double vb[7]={.25,.5,1,2,4,8,16},nb[7]={1,4,8,16,32,64,128},db[7]={.125,.25,.5,1,2,4,8};
    int m,e,k,b,z,slot;size_t c;double *r;
    memset(audit.radial,0,audit.nlayers*VA_RADIAL*sizeof(double));
    memset(audit.bins,0,VA_BIN_KINDS*VA_BINS*VA_METRICS*sizeof(double));
    memset(audit.samples,0,VA_TOP*VA_SAMPLE*sizeof(double));
    for(k=0;k<VA_TOP;k++)audit.samples[k*VA_SAMPLE]=-1;
    for(m=1;m<=audit.nc;m++)for(e=1;e<=audit.nel;e++) {
        double v,after,delta,score;long long n,n24;
        c=cell(m,e);v=audit.volume[c];n=audit.total[c];
        n24=E->trace.ntracer_flavor[m][24][e];after=element_fraction(E,m,e,24);
        z=E->lmesh.nzs-1+(e-1)%E->lmesh.elz;
        if(z<0 || z>=audit.nlayers)pices_fail(E,"invalid volume audit radial index");
        r=audit.radial+z*VA_RADIAL;r[0]+=v;r[1]+=v*audit.depth[c];
        r[2]+=v*audit.before24[c];r[3]+=v*audit.moved24[c];r[4]+=v*after;
        r[5]+=n;r[6]+=n24;if(!n)r[7]+=v;
        b=n?(audit.old_total[c]?3:2):(audit.old_total[c]?0:1);
        add_bin(0,b,v,audit.before24[c],audit.moved24[c],after,n,n24);
        add_bin(1,choose_bin(v/audit.mean_volume,vb),v,audit.before24[c],audit.moved24[c],after,n,n24);
        add_bin(2,choose_bin((double)n,nb),v,audit.before24[c],audit.moved24[c],after,n,n24);
        add_bin(3,choose_bin(n/v/audit.mean_density,db),v,audit.before24[c],audit.moved24[c],after,n,n24);
        delta=v*(after-audit.before24[c]);score=fabs(delta);slot=0;
        for(k=1;k<VA_TOP;k++)if(audit.samples[k*VA_SAMPLE]<audit.samples[slot*VA_SAMPLE])slot=k;
        if(score>audit.samples[slot*VA_SAMPLE]) {
            r=audit.samples+slot*VA_SAMPLE;r[0]=score;r[1]=E->parallel.me;r[2]=m;r[3]=e;r[4]=audit.depth[c];r[5]=v;
            r[6]=(double)audit.old_total[c];r[7]=(double)n;r[8]=audit.before24[c];r[9]=audit.moved24[c];r[10]=after;r[11]=delta;r[12]=0;
        }
    }
    sum_double(E,audit.radial,audit.nlayers*VA_RADIAL);
    sum_double(E,audit.bins,VA_BIN_KINDS*VA_BINS*VA_METRICS);
    MPI_Gather(audit.samples,VA_TOP*VA_SAMPLE,MPI_DOUBLE,audit.gathered,VA_TOP*VA_SAMPLE,MPI_DOUBLE,0,E->parallel.world);
    if(E->parallel.me==0) {
        memset(audit.samples,0,VA_TOP*VA_SAMPLE*sizeof(double));
        for(k=0;k<VA_TOP;k++)audit.samples[k*VA_SAMPLE]=-1;
        for(b=0;b<E->parallel.nproc*VA_TOP;b++) {
            r=audit.gathered+b*VA_SAMPLE;slot=0;
            for(k=1;k<VA_TOP;k++)if(audit.samples[k*VA_SAMPLE]<audit.samples[slot*VA_SAMPLE])slot=k;
            if(r[0]>audit.samples[slot*VA_SAMPLE])memcpy(audit.samples+slot*VA_SAMPLE,r,VA_SAMPLE*sizeof(double));
        }
    }
}
static void publish(struct All_variables *E,const char *status)
{
    int f,g,k,b,z;double *r;
    const char *kinds[4]={"occupancy","volume_ratio","particle_count","sampling_ratio"};
    const char *occupancy[4]={"new_empty","persistent_empty","recovered","occupied"};
    if(E->parallel.me)return;
    for(f=0;f<audit.nf;f++)for(g=0;g<audit.nf;g++) {
        k=f*audit.nf+g;
        if(audit.transitions[k])fprintf(audit.transition_file,"%s,%d,%d,%d,%lld,%.17g\n",status,audit.step,f,g,audit.transitions[k],audit.transfer[k]);
    }
    if(audit.detail) {
        for(z=0;z<audit.nlayers;z++) {
            r=audit.radial+z*VA_RADIAL;
            fprintf(audit.radial_file,"%s,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.0f,%.0f,%.17g\n",status,audit.step,z+1,r[0]?r[1]/r[0]:0,r[0],r[2],r[3],r[4],r[5],r[6],r[7]);
        }
        for(k=0;k<VA_BIN_KINDS;k++)for(b=0;b<(k?VA_BINS:4);b++) {
            char label[32];r=audit.bins+(k*VA_BINS+b)*VA_METRICS;
            if(k)snprintf(label,sizeof(label),"%d",b);else snprintf(label,sizeof(label),"%s",occupancy[b]);
            fprintf(audit.bin_file,"%s,%d,%s,%s,%.0f,%.17g,%.17g,%.17g,%.17g,%.0f,%.0f\n",status,audit.step,kinds[k],label,r[0],r[1],r[2],r[3],r[4],r[5],r[6]);
        }
        for(k=0;k<VA_TOP;k++) {
            r=audit.samples+k*VA_SAMPLE;if(r[0]<0)continue;
            fprintf(audit.sample_file,"%s,%d,%.0f,%.0f,%.0f,%.17g,%.17g,%.0f,%.0f,%.17g,%.17g,%.17g,%.17g\n",status,audit.step,r[1],r[2],r[3],r[4],r[5],r[6],r[7],r[8],r[9],r[10],r[11]);
        }
    }
    /* Publish auxiliary accepted tables before the global commit rows. */
    fflush(audit.transition_file);fflush(audit.radial_file);fflush(audit.bin_file);fflush(audit.sample_file);
    for(f=0;f<audit.nf;f++)fprintf(audit.global_file,"%s,%d,%.17g,%.17g,%d,%lld,%lld,%lld,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",status,audit.step,audit.time,audit.age,f,audit.nbefore[f],audit.nmoved[f],audit.nafter[f],audit.before[f],audit.arrival[f],audit.after[f],audit.node[f],audit.source[f],audit.residual[f],audit.baseline[f],audit.shell);
    fflush(audit.global_file);
}
void pices_volume_audit_rebuilt(struct All_variables *E)
{
    int f,g;double started,sum_elements=0,sum_nodes=0;
    if(!audit.enabled || !audit.active || !audit.moved)return;
    started=MPI_Wtime();
    integrate(E,audit.after,audit.node,audit.nafter);
    for(f=0;f<audit.nf;f++){sum_elements+=audit.after[f];sum_nodes+=audit.node[f];}
    if(E->parallel.me==0 && (!isfinite(sum_elements) || !isfinite(sum_nodes) ||
       fabs(sum_elements-audit.shell)>VA_ELEMENT_TOL*audit.shell || fabs(sum_nodes-audit.shell)>VA_NODE_TOL*audit.shell ||
       audit.bounds_error[0]>VA_ELEMENT_TOL || audit.bounds_error[1]>VA_NODE_TOL))
        fprintf(stderr,"VOLUME_AUDIT_ANOMALY step=%d element_sum_error=%.17g node_sum_error=%.17g element_fraction_bound_violation=%.17g node_fraction_bound_violation=%.17g\n",audit.step,sum_elements-audit.shell,sum_nodes-audit.shell,audit.bounds_error[0],audit.bounds_error[1]);
    sum_integer(E,audit.transitions,audit.nf*audit.nf);sum_double(E,audit.transfer,audit.nf*audit.nf);
    memset(audit.source,0,audit.nf*sizeof(double));
    for(f=0;f<audit.nf;f++)for(g=0;g<audit.nf;g++) {
        audit.source[g]+=audit.transfer[f*audit.nf+g];audit.source[f]-=audit.transfer[f*audit.nf+g];
    }
    for(f=0;f<audit.nf;f++) {
        long long source_count=0;
        for(g=0;g<audit.nf;g++)source_count+=audit.transitions[g*audit.nf+f]-audit.transitions[f*audit.nf+g];
        audit.residual[f]=audit.after[f]-audit.arrival[f]-audit.source[f];
        if(E->parallel.me==0 && (audit.nbefore[f]!=audit.nmoved[f] || audit.nafter[f]-audit.nmoved[f]!=source_count ||
           !isfinite(audit.residual[f]) || fabs(audit.residual[f])>VA_ELEMENT_TOL*audit.shell))
            fprintf(stderr,"VOLUME_AUDIT_ANOMALY step=%d flavor=%d particle_delta=%lld conversion_count_residual=%lld budget_residual=%.17g\n",audit.step,f,audit.nmoved[f]-audit.nbefore[f],audit.nafter[f]-audit.nmoved[f]-source_count,audit.residual[f]);
    }
    if(E->parallel.me==0)for(f=0;f<audit.nf;f++)if(audit.transitions[24*audit.nf+f] || audit.transitions[f*audit.nf+24])
        fprintf(stderr,"VOLUME_AUDIT_ANOMALY step=%d primordial_label_changed other_flavor=%d\n",audit.step,f);
    audit.time=E->monitor.elapsed_time;audit.age=E->control.start_age-audit.time*E->data.scalet;
    if(audit.detail)detail_statistics(E);
    publish(E,"pending");audit.pending=1;
    audit.seconds+=MPI_Wtime()-started;
}
void pices_volume_audit_commit(struct All_variables *E)
{
    if(!audit.enabled || !audit.pending)return;
    publish(E,"accepted");audit.pending=0;audit.active=0;
    if(E->parallel.me==0) {
        fprintf(E->fp,"PICES_VOLUME_AUDIT step=%d status=accepted C24_volume_nd=%.17g C24_motion_delta_nd=%.17g C24_source_nd=%.17g C24_closure_nd=%.17g node_minus_element_nd=%.17g baseline_step=%d cumulative_local_seconds=%.9g element_fraction_bound_violation=%.17g node_fraction_bound_violation=%.17g\n",audit.step,audit.after[24],audit.arrival[24]-audit.before[24],audit.source[24],audit.residual[24],audit.node[24]-audit.after[24],audit.baseline_step,audit.seconds,audit.bounds_error[0],audit.bounds_error[1]);
        fflush(E->fp);
    }
}
