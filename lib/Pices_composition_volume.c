/* Experimental EBA fixed-volume composition reconstruction.
 *
 * Particles carry immutable positive volume weights. Raw cell mixtures retain
 * as much volume as their capacity permits. A nonnegative, finite-support
 * transport plan moves ONLY surplus into deficits. Its rows/columns are scaled
 * to capacities by iterative proportional fitting; all species use the same
 * plan. No remote/global species scaling is used. Infeasible support is an
 * explicit scientific failure, not permission to enlarge the support silently.
 *
 * K = diag(V) H^hops is symmetric: H maps cell values through shared nodes
 * using positive FE integrals. Thus K is applied without storing a global graph
 * and communicates only via the existing double-precision node exchange.
 * This is a bounded global-volume remap prototype, not a proof of consistency
 * with the local physical transport equation.
 */
#include <errno.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "pices.h"
#include "pices_composition_volume.h"

#define VF 25
#define VOLUME_BUDGET_TOL 1e-10
#define WEIGHT_ROUND_TOL (64*DBL_EPSILON)
struct Volume_state {
    int geometry, nel, nno, pending, files_open;
    double shell;
    double *volume,*weight,*owned_node,*node_volume,*center;
    double *raw,*corrected,*nodal,*q,*surplus,*deficit,*left,*right,*ka,*kb;
    double *work0,*work1,*input,*output,*node_work;
    double local_source[VF],step_source[VF];
    struct Pices_volume_summary last;
    FILE *steps,*flavors;
};
static struct Volume_state volume_state;
void get_global_shape_fn(struct All_variables *,int,struct Shape_function *,
    struct Shape_function_dx *,struct Shape_function_dA *,int,int,double [4][9],int,int);

static void *storage(struct All_variables *E,size_t count,size_t size)
{
    void *p=calloc(count,size);
    if(!p)pices_fail(E,"composition volume allocation failed");
    return p;
}
static void sums(struct All_variables *E,double *values,int count)
{ MPI_Allreduce(MPI_IN_PLACE,values,count,MPI_DOUBLE,MPI_SUM,E->parallel.world); }
static double largest(struct All_variables *E,double value)
{ double out;MPI_Allreduce(&value,&out,1,MPI_DOUBLE,MPI_MAX,E->parallel.world);return out; }
static double smallest(struct All_variables *E,double value)
{ double out;MPI_Allreduce(&value,&out,1,MPI_DOUBLE,MPI_MIN,E->parallel.world);return out; }
static void exchange(struct All_variables *E,double *values)
{
    double *p[NCS]={0};p[1]=values;
    (E->exchange_node_d)(E,p,E->mesh.levmax);
}
static double old_fraction(struct All_variables *E,int e,int flavor)
{
    int f;double x=1;
    if(flavor)return E->composition.comp_el[1][flavor-1][e];
    for(f=0;f<24;f++)x-=E->composition.comp_el[1][f][e];
    return x;
}
int pices_composition_volume_enabled(struct All_variables *E)
{ return E->pices.volume_mode!=0; }
void pices_composition_volume_register(struct All_variables *E)
{
    if(!pices_composition_volume_enabled(E))return;
    if(!E->pices.enabled || !E->pices.eba || E->trace.nflavors!=VF)
        pices_fail(E,"composition volume mode requires EBA PICES and 25 flavors");
    E->pices.volume_slot=E->trace.number_of_extra_quantities++;
}
double pices_composition_volume_weight(struct All_variables *E,int cap,int particle)
{
    double value;
    if(!pices_composition_volume_enabled(E) || E->pices.volume_slot<0 ||
       E->pices.volume_slot>=E->trace.number_of_extra_quantities)
        pices_fail(E,"missing composition volume particle slot");
    value=E->trace.extraq[cap][E->pices.volume_slot][particle];
    if(!isfinite(value) || value<=0)pices_fail(E,"nonpositive or nonfinite particle volume weight");
    return value;
}
static void geometry(struct All_variables *E)
{
    struct Volume_state *v=&volume_state;
    struct Shape_function GN;struct Shape_function_dx dx;
    struct Shape_function_dA omega;double rtf[4][9];int e,n,a,g,d;
    if(v->geometry)return;
    if(E->sphere.caps_per_proc!=1 || !E->composition.on || E->composition.ncomp!=24 ||
       E->composition.ibuoy_type!=1 || E->control.pseudo_free_surf)
        pices_fail(E,"composition volume requires one-cap EBA ratio composition and fixed boundaries");
    v->nel=E->lmesh.nel;v->nno=E->lmesh.nno;
#define ALLOC(name,count) v->name=storage(E,(count),sizeof(double))
    ALLOC(volume,v->nel);ALLOC(weight,(size_t)v->nel*8);ALLOC(center,(size_t)v->nel*3);
    ALLOC(owned_node,v->nno+1);ALLOC(node_volume,v->nno+1);ALLOC(node_work,v->nno+1);
    ALLOC(raw,(size_t)v->nel*VF);ALLOC(corrected,(size_t)v->nel*VF);ALLOC(nodal,(size_t)v->nno*VF);
    ALLOC(q,v->nel);ALLOC(surplus,v->nel);ALLOC(deficit,v->nel);ALLOC(left,v->nel);ALLOC(right,v->nel);
    ALLOC(ka,v->nel);ALLOC(kb,v->nel);ALLOC(work0,v->nel);ALLOC(work1,v->nel);ALLOC(input,v->nel);ALLOC(output,v->nel);
#undef ALLOC
    for(e=1;e<=v->nel;e++) {
        get_global_shape_fn(E,e,&GN,&dx,&omega,0,1,rtf,E->mesh.levmax,1);
        for(g=1;g<=8;g++) {
            double measure=omega.vpt[g]*g_point[g].weight[2];
            if(!isfinite(measure) || measure<=0)pices_fail(E,"invalid volume quadrature measure");
            v->volume[e-1]+=measure;
            for(a=1;a<=8;a++)v->weight[(size_t)(e-1)*8+a-1]+=E->N.vpt[GNVINDEX(a,g)]*measure;
        }
        for(a=1;a<=8;a++) {
            double w=v->weight[(size_t)(e-1)*8+a-1];n=E->ien[1][e].node[a];
            if(!isfinite(w) || w<=0)pices_fail(E,"nonpositive composition node weight");
            v->owned_node[n]+=w;
            for(d=0;d<3;d++)v->center[(size_t)(e-1)*3+d]+=E->x[1][d+1][n]/8;
        }
        v->shell+=v->volume[e-1];
    }
    sums(E,&v->shell,1);
    memcpy(v->node_volume,v->owned_node,(v->nno+1)*sizeof(double));exchange(E,v->node_volume);
    for(n=1;n<=v->nno;n++)if(!isfinite(v->node_volume[n]) || v->node_volume[n]<=0)
        pices_fail(E,"invalid shared-node volume");
    v->geometry=1;
}
/* Self-adjoint under the V inner product; no composition is modified. */
static void average(struct All_variables *E,const double *input,double *output)
{
    struct Volume_state *v=&volume_state;int e,a,n;
    memset(v->node_work,0,(v->nno+1)*sizeof(double));
    for(e=1;e<=v->nel;e++)for(a=1;a<=8;a++) {
        n=E->ien[1][e].node[a];
        v->node_work[n]+=v->weight[(size_t)(e-1)*8+a-1]*input[e-1];
    }
    exchange(E,v->node_work);
    for(n=1;n<=v->nno;n++)v->node_work[n]/=v->node_volume[n];
    for(e=1;e<=v->nel;e++) {
        double value=0;
        for(a=1;a<=8;a++)value+=v->weight[(size_t)(e-1)*8+a-1]*v->node_work[E->ien[1][e].node[a]];
        output[e-1]=value/v->volume[e-1];
    }
}
static void kernel(struct All_variables *E,const double *input,double *output)
{
    struct Volume_state *v=&volume_state;int e,k;double *from=v->work0,*to=v->work1,*swap;
    memcpy(from,input,v->nel*sizeof(double));
    for(k=0;k<E->pices.volume_hops;k++) {
        average(E,from,to);swap=from;from=to;to=swap;
    }
    for(e=0;e<v->nel;e++)output[e]=v->volume[e]*from[e];
}
static void deposit(struct All_variables *E,double totals[VF])
{
    struct Volume_state *v=&volume_state;int p,e,f;
    memset(v->raw,0,(size_t)v->nel*VF*sizeof(double));memset(totals,0,VF*sizeof(double));
    for(p=1;p<=E->trace.ntracers[1];p++) {
        double label=E->trace.extraq[1][0][p],w=pices_composition_volume_weight(E,1,p);
        e=E->trace.ielement[1][p];
        if(e<1 || e>v->nel || !isfinite(label) || label<0 || label>=VF || floor(label)!=label)
            pices_fail(E,"invalid volume particle label or host");
        f=(int)label;v->raw[(size_t)(e-1)*VF+f]+=w;
    }
    for(e=0;e<v->nel;e++)for(f=0;f<VF;f++)totals[f]+=v->raw[(size_t)e*VF+f];
    sums(E,totals,VF);
}
void pices_composition_volume_initialize(struct All_variables *E,int restored)
{
    struct Volume_state *v=&volume_state;int p,e,f;double totals[VF];
    int *counts;
    if(!pices_composition_volume_enabled(E))return;
    geometry(E);
    if(!restored) {
        if(E->pices.volume_initialized)pices_fail(E,"duplicate composition weight initialization");
        counts=storage(E,v->nel,sizeof(int));
        for(p=1;p<=E->trace.ntracers[1];p++) {
            e=E->trace.ielement[1][p];if(e<1 || e>v->nel)pices_fail(E,"invalid initial volume host");counts[e-1]++;
        }
        for(e=0;e<v->nel;e++)if(!counts[e])pices_fail(E,"fixed-volume initialization requires no empty initial cells");
        for(p=1;p<=E->trace.ntracers[1];p++) {
            e=E->trace.ielement[1][p]-1;E->trace.extraq[1][E->pices.volume_slot][p]=v->volume[e]/counts[e];
        }
        free(counts);deposit(E,totals);
        memcpy(E->pices.volume_baseline,totals,sizeof(totals));memset(E->pices.volume_source,0,sizeof(totals));
    }
    else {
        if(!E->pices.volume_initialized)pices_fail(E,"checkpoint lacks initialized particle volume history");
        deposit(E,totals);
        for(f=0;f<VF;f++)if(!isfinite(E->pices.volume_baseline[f]) || !isfinite(E->pices.volume_source[f]) ||
           fabs(totals[f]-E->pices.volume_baseline[f]-E->pices.volume_source[f])>VOLUME_BUDGET_TOL*v->shell)
            pices_fail(E,"restored species weight budget differs from checkpoint");
    }
    {
        double total=0;for(f=0;f<VF;f++)total+=totals[f];
        if(fabs(total-v->shell)>VOLUME_BUDGET_TOL*v->shell)pices_fail(E,"particle weights do not fill the shell");
    }
    E->pices.volume_initialized=1;
}
void pices_composition_volume_begin(struct All_variables *E)
{
    if(!pices_composition_volume_enabled(E))return;
    if(volume_state.pending)pices_fail(E,"previous composition volume step was not accepted");
    memset(volume_state.local_source,0,sizeof(volume_state.local_source));
}
void pices_composition_volume_transition(struct All_variables *E,int cap,int particle,int old_flavor,int new_flavor)
{
    double w;
    if(!pices_composition_volume_enabled(E) || !E->pices.volume_initialized || old_flavor==new_flavor)return;
    if(old_flavor<0 || old_flavor>=VF || new_flavor<0 || new_flavor>=VF)
        pices_fail(E,"invalid weighted flavor transition");
    if(old_flavor==24 || new_flavor==24)pices_fail(E,"primordial flavor changed in volume experiment");
    w=pices_composition_volume_weight(E,cap,particle);
    volume_state.local_source[old_flavor]-=w;volume_state.local_source[new_flavor]+=w;
}
static FILE *new_output(struct All_variables *E,const char *path,const char *header)
{
    FILE *f=fopen(path,"wx");
    if(!f)pices_fail(E,"volume-conserve diagnostic exists or cannot be created; use a fresh segment directory");
    if(header)fprintf(f,"%s\n",header);return f;
}
static void outputs(struct All_variables *E)
{
    struct Volume_state *v=&volume_state;FILE *meta;
    if(v->files_open)return;v->files_open=1;
    if(E->parallel.me)return;
    if(mkdir("Diagnostic",0775) && errno!=EEXIST)pices_fail(E,"cannot create Diagnostic");
    if(mkdir("Diagnostic/VolumeConserve",0775) && errno!=EEXIST)pices_fail(E,"cannot create VolumeConserve");
    v->steps=new_output(E,"Diagnostic/VolumeConserve/steps.csv","status,step,time_nd,mode,feasible,iterations,hops,shell_volume_nd,capacity_ratio_min,capacity_ratio_max,surplus_cells,deficit_cells,empty_cells,minimum_transfer_nd,row_error,column_error,species_error,node_error,fraction_error,partition_error,rounding_adjustment_nd,transfer_distance_rms_nd,wall_seconds,reason");
    v->flavors=new_output(E,"Diagnostic/VolumeConserve/flavors.csv","status,step,time_nd,flavor,target_volume_nd,raw_volume_nd,corrected_volume_nd,node_volume_nd,original_volume_nd,step_source_nd,cumulative_source_nd,baseline_volume_nd,shell_volume_nd,feasible");
    meta=new_output(E,"Diagnostic/VolumeConserve/metadata.json",NULL);
    fprintf(meta,"{\"format_version\":\"PICES_VOLUME_CONSERVE_V1\",\"flavors\":25,\"baseline_step\":%d,\"mode\":%d,\"ranks\":%d,\"support_hops\":%d,\"iterations_limit\":%d,\"row_tolerance\":%.17g,\"species_budget_tolerance_shell\":1e-10,\"support\":\"shared-node element neighbors, at most support_hops\",\"method\":\"surplus-only finite-support positive transport, common to all species\",\"node_mapping\":\"double positive finite-element lumped-volume convex projection\",\"rounding\":\"Rows are closed only after convergence within row_tolerance; the amount and species-budget change are recorded.\",\"note\":\"Prototype global volume conservation and boundedness do not certify local physical transport accuracy. Pending rows require accepted records.\"}\n",E->monitor.solution_cycles,E->pices.volume_mode,E->parallel.nproc,E->pices.volume_hops,E->pices.volume_iterations,E->pices.volume_tolerance);
    fclose(meta);
}
static void write_rows(struct All_variables *E,const char *status)
{
    struct Volume_state *v=&volume_state;struct Pices_volume_summary *s=&v->last;int f;
    if(E->parallel.me)return;
    fprintf(v->steps,"%s,%d,%.17g,%d,%d,%d,%d,%.17g,%.17g,%.17g,%lld,%lld,%lld,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%s\n",status,s->step,(double)E->monitor.elapsed_time,E->pices.volume_mode,s->feasible,s->iterations,s->hops,s->shell,s->capacity_min,s->capacity_max,s->surplus_cells,s->deficit_cells,s->empty_cells,s->minimum_transfer_nd,s->row_error,s->column_error,s->species_error,s->node_error,s->fraction_error,s->partition_error,s->rounding_adjustment_nd,s->transfer_distance_rms_nd,s->wall_seconds,s->reason);
    for(f=0;f<VF;f++)fprintf(v->flavors,"%s,%d,%.17g,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d\n",status,s->step,(double)E->monitor.elapsed_time,f,s->target[f],s->raw[f],s->corrected[f],s->node[f],s->original[f],s->source[f],E->pices.volume_source[f],E->pices.volume_baseline[f],s->shell,s->feasible);
    fflush(v->steps);fflush(v->flavors);
}
static int failed(struct All_variables *E,const char *reason)
{
    snprintf(volume_state.last.reason,sizeof(volume_state.last.reason),"%s",reason);
    volume_state.last.feasible=0;return 0;
}
static int transport(struct All_variables *E)
{
    struct Volume_state *v=&volume_state;struct Pices_volume_summary *s=&v->last;
    int e,f,it;double maxerror=0,tol=E->pices.volume_tolerance;
    double total_surplus=0,total_deficit=0;long long counts[3]={0};
    s->capacity_min=DBL_MAX;
    for(e=0;e<v->nel;e++) {
        double q=0,delta,ratio;
        for(f=0;f<VF;f++)q+=v->raw[(size_t)e*VF+f];v->q[e]=q;
        ratio=q/v->volume[e];delta=q-v->volume[e];
        maxerror=fmax(maxerror,fabs(delta)/v->volume[e]);
        s->capacity_min=fmin(s->capacity_min,ratio);s->capacity_max=fmax(s->capacity_max,ratio);
        v->surplus[e]=delta>WEIGHT_ROUND_TOL*v->volume[e]?delta:0;
        v->deficit[e]=delta< -WEIGHT_ROUND_TOL*v->volume[e]?-delta:0;
        total_surplus+=v->surplus[e];total_deficit+=v->deficit[e];
        counts[0]+=v->surplus[e]>0;counts[1]+=v->deficit[e]>0;counts[2]+=q==0;
        v->left[e]=0;v->right[e]=v->surplus[e]>0?1:0;
    }
    sums(E,&total_surplus,1);sums(E,&total_deficit,1);
    MPI_Allreduce(MPI_IN_PLACE,counts,3,MPI_LONG_LONG_INT,MPI_SUM,E->parallel.world);
    s->surplus_cells=counts[0];s->deficit_cells=counts[1];s->empty_cells=counts[2];
    s->capacity_min=smallest(E,s->capacity_min);s->capacity_max=largest(E,s->capacity_max);
    s->minimum_transfer_nd=.5*(total_surplus+total_deficit);
    if(fabs(total_surplus-total_deficit)>VOLUME_BUDGET_TOL*v->shell)
        return failed(E,"particle total volume differs from shell capacity");
    maxerror=largest(E,maxerror);
    if(maxerror<=tol) {
        memcpy(v->corrected,v->raw,(size_t)v->nel*VF*sizeof(double));
        s->row_error=maxerror;return 1;
    }
    for(it=1;it<=E->pices.volume_iterations;it++) {
        int bad=0;double scale=0,row_error=0,column_error=0;
        kernel(E,v->right,v->kb);
        for(e=0;e<v->nel;e++) {
            if(v->deficit[e]>0 && (!isfinite(v->kb[e]) || v->kb[e]<=0))bad=1;
            v->left[e]=v->deficit[e]>0 && v->kb[e]>0?v->deficit[e]/v->kb[e]:0;
            if(!isfinite(v->left[e]))bad=1;
        }
        if(largest(E,bad))return failed(E,"deficit unreachable within configured support");
        kernel(E,v->left,v->ka);
        for(e=0;e<v->nel;e++) {
            if(v->surplus[e]>0 && (!isfinite(v->ka[e]) || v->ka[e]<=0))bad=1;
            v->right[e]=v->surplus[e]>0 && v->ka[e]>0?v->surplus[e]/v->ka[e]:0;
            if(!isfinite(v->right[e]))bad=1;
            scale=fmax(scale,v->right[e]);
        }
        if(largest(E,bad))return failed(E,"surplus unreachable within configured support");
        scale=largest(E,scale);
        if(!(scale>0) || !isfinite(scale))return failed(E,"invalid capacity scaling");
        /* Remove the arbitrary common scaling before it becomes ill-conditioned. */
        if(scale>1e20 || scale<1e-20)for(e=0;e<v->nel;e++) {
            v->right[e]/=scale;v->left[e]*=scale;v->ka[e]*=scale;
        }
        kernel(E,v->right,v->kb);
        for(e=0;e<v->nel;e++) {
            double row=v->left[e]*v->kb[e],column=v->right[e]*v->ka[e];
            if(!isfinite(row) || !isfinite(column))bad=1;
            row_error=fmax(row_error,fabs(row-v->deficit[e])/v->volume[e]);
            column_error=fmax(column_error,fabs(column-v->surplus[e])/v->volume[e]);
        }
        if(largest(E,bad))return failed(E,"nonfinite capacity transport");
        s->iterations=it;s->row_error=largest(E,row_error);s->column_error=largest(E,column_error);
        if(fmax(s->row_error,s->column_error)<=tol)break;
    }
    if(it>E->pices.volume_iterations)return failed(E,"capacity constraints did not converge within iteration limit");
    for(f=0;f<VF;f++) {
        for(e=0;e<v->nel;e++)v->input[e]=v->q[e]>0?v->right[e]*v->raw[(size_t)e*VF+f]/v->q[e]:0;
        kernel(E,v->input,v->output);
        for(e=0;e<v->nel;e++) {
            double retained=v->q[e]>0?fmin(1,v->volume[e]/v->q[e]):0;
            v->corrected[(size_t)e*VF+f]=v->raw[(size_t)e*VF+f]*retained+v->left[e]*v->output[e];
        }
    }
    /* Second moment of W measures actual donor/destination centre distance. */
    {
        int d;double moment=0,total=0;
        for(e=0;e<v->nel;e++) {
            double r2=0;for(d=0;d<3;d++)r2+=v->center[(size_t)e*3+d]*v->center[(size_t)e*3+d];
            v->input[e]=v->right[e]*r2;
        }
        kernel(E,v->input,v->output);
        for(e=0;e<v->nel;e++) {
            double r2=0;for(d=0;d<3;d++)r2+=v->center[(size_t)e*3+d]*v->center[(size_t)e*3+d];
            moment+=v->left[e]*(v->output[e]+r2*v->kb[e]);total+=v->left[e]*v->kb[e];
        }
        for(d=0;d<3;d++) {
            for(e=0;e<v->nel;e++)v->input[e]=v->right[e]*v->center[(size_t)e*3+d];
            kernel(E,v->input,v->output);
            for(e=0;e<v->nel;e++)moment-=2*v->left[e]*v->center[(size_t)e*3+d]*v->output[e];
        }
        sums(E,&moment,1);sums(E,&total,1);
        if(moment< -1e-12*fmax(total,DBL_MIN))return failed(E,"negative transport distance moment");
        s->transfer_distance_rms_nd=total>0?sqrt(fmax(0,moment)/total):0;
    }
    return 1;
}
/* Only roundoff-sized row closure is allowed after the transport converged.
 * The adjustment and its actual species-budget impact remain visible. */
static int close_rows(struct All_variables *E)
{
    struct Volume_state *v=&volume_state;struct Pices_volume_summary *s=&v->last;
    int e,f,bad=0;double amount=0;
    for(e=0;e<v->nel;e++) {
        double total=0,partition=0,maximum=-1,before[VF];int dominant=0;
        for(f=0;f<VF;f++) {
            double x=v->corrected[(size_t)e*VF+f];
            before[f]=x;if(!isfinite(x) || x<0)bad=1;total+=x;
        }
        if(!isfinite(total) || total<=0 || fabs(total/v->volume[e]-1)>2*E->pices.volume_tolerance) {bad=1;continue;}
        for(f=0;f<VF;f++) {
            double old=v->corrected[(size_t)e*VF+f],fraction=old/total;
            partition+=fraction;
            v->corrected[(size_t)e*VF+f]=fraction;
            if(fraction>maximum){maximum=fraction;dominant=f;}
        }
        v->corrected[(size_t)e*VF+dominant]+=1-partition;
        /* Background is represented as the complement in the production code. */
        partition=0;for(f=1;f<VF;f++)partition+=v->corrected[(size_t)e*VF+f];
        v->corrected[(size_t)e*VF]=1-partition;
        for(f=0;f<VF;f++) {
            double x=v->corrected[(size_t)e*VF+f];
            if(x< -1e-12 || x>1+1e-12)bad=1;
            s->fraction_error=fmax(s->fraction_error,fmax(-x,x-1));
            v->corrected[(size_t)e*VF+f]=x*v->volume[e];
            amount+=fabs(v->corrected[(size_t)e*VF+f]-before[f]);
            s->corrected[f]+=x*v->volume[e];
        }
    }
    sums(E,&amount,1);s->rounding_adjustment_nd=amount;
    s->fraction_error=largest(E,s->fraction_error);sums(E,s->corrected,VF);
    if(largest(E,bad))return failed(E,"capacity row closure exceeds rounding tolerance");
    for(f=0;f<VF;f++)s->species_error=fmax(s->species_error,fabs(s->corrected[f]-s->target[f])/v->shell);
    if(s->species_error>VOLUME_BUDGET_TOL)return failed(E,"corrected species budget does not match particle history");
    return 1;
}
static int node_projection(struct All_variables *E)
{
    struct Volume_state *v=&volume_state;struct Pices_volume_summary *s=&v->last;int e,n,a,f,bad=0;
    for(f=0;f<VF;f++) {
        memset(v->node_work,0,(v->nno+1)*sizeof(double));
        for(e=1;e<=v->nel;e++)for(a=1;a<=8;a++) {
            n=E->ien[1][e].node[a];
            v->node_work[n]+=v->weight[(size_t)(e-1)*8+a-1]*(v->corrected[(size_t)(e-1)*VF+f]/v->volume[e-1]);
        }
        exchange(E,v->node_work);
        for(n=1;n<=v->nno;n++)v->nodal[(size_t)(n-1)*VF+f]=v->node_work[n]/v->node_volume[n];
    }
    for(n=1;n<=v->nno;n++) {
        double partition=0;for(f=1;f<VF;f++)partition+=v->nodal[(size_t)(n-1)*VF+f];
        v->nodal[(size_t)(n-1)*VF]=1-partition;partition=0;
        for(f=0;f<VF;f++) {
            double x=v->nodal[(size_t)(n-1)*VF+f];
            if(!isfinite(x) || x< -1e-12 || x>1+1e-12)bad=1;
            s->fraction_error=fmax(s->fraction_error,fmax(-x,x-1));partition+=x;
            s->node[f]+=v->owned_node[n]*x;
        }
        s->partition_error=fmax(s->partition_error,fabs(partition-1));
    }
    sums(E,s->node,VF);s->fraction_error=largest(E,s->fraction_error);s->partition_error=largest(E,s->partition_error);
    for(f=0;f<VF;f++)s->node_error=fmax(s->node_error,fabs(s->node[f]-s->corrected[f])/v->shell);
    if(largest(E,bad) || s->partition_error>1e-12 || s->node_error>VOLUME_BUDGET_TOL)
        return failed(E,"double nodal projection violates bounds or integral closure");
    return 1;
}
int pices_composition_volume_reconstruct(struct All_variables *E)
{
    struct Volume_state *v=&volume_state;struct Pices_volume_summary *s=&v->last;
    int e,n,f;double started=MPI_Wtime();
    if(!pices_composition_volume_enabled(E))return 1;
    if(!E->pices.volume_initialized || !v->geometry)pices_fail(E,"composition weights not initialized before reconstruction");
    if(v->pending)pices_fail(E,"duplicate unaccepted conservative composition reconstruction");
    memset(s,0,sizeof(*s));s->initialized=1;s->feasible=1;s->shell=v->shell;
    s->step=E->monitor.solution_cycles;s->hops=E->pices.volume_hops;strcpy(s->reason,"ok");
    memcpy(v->step_source,v->local_source,sizeof(v->step_source));sums(E,v->step_source,VF);
    memset(v->local_source,0,sizeof(v->local_source));
    for(f=0;f<VF;f++) {
        E->pices.volume_source[f]+=v->step_source[f];s->source[f]=v->step_source[f];
        s->target[f]=E->pices.volume_baseline[f]+E->pices.volume_source[f];
    }
    deposit(E,s->raw);
    for(e=1;e<=v->nel;e++)for(f=0;f<VF;f++)s->original[f]+=v->volume[e-1]*old_fraction(E,e,f);
    sums(E,s->original,VF);
    for(f=0;f<VF;f++)if(fabs(s->raw[f]-s->target[f])>VOLUME_BUDGET_TOL*v->shell)
        pices_fail(E,"particle volume budget changed without a weighted source");
    if(transport(E) && close_rows(E) && node_projection(E)) {
        if(E->pices.volume_mode==2) {
            for(e=1;e<=v->nel;e++)for(f=1;f<VF;f++)E->composition.comp_el[1][f-1][e]=v->corrected[(size_t)(e-1)*VF+f]/v->volume[e-1];
            for(n=1;n<=v->nno;n++)for(f=1;f<VF;f++)E->composition.comp_node[1][f-1][n]=v->nodal[(size_t)(n-1)*VF+f];
        }
    }
    s->wall_seconds=MPI_Wtime()-started;outputs(E);v->pending=1;write_rows(E,"pending");
    if(!s->feasible && E->pices.volume_mode==2)pices_fail(E,s->reason);
    return s->feasible;
}
void pices_composition_volume_check(struct All_variables *E)
{
    struct Volume_state *v=&volume_state;double global[VF]={0},bounds=0,partition=0;int e,f;
    if(!pices_composition_volume_enabled(E) || !E->pices.volume_initialized)return;
    /* Shadow deliberately retains ratio physics; only active changes its contract. */
    if(E->pices.volume_mode!=2)return;
    /* The controlled driver takes its common initial snapshot before its first
     * reconstruction. Initial ratio fields already match V/N-assigned weights. */
    if(v->last.initialized && !v->last.feasible)pices_fail(E,"active composition has no feasible reconstruction");
    for(e=1;e<=v->nel;e++) {
        double total=0;
        for(f=0;f<VF;f++) {
            double x=old_fraction(E,e,f);if(!isfinite(x))bounds=HUGE_VAL;
            bounds=fmax(bounds,fmax(-x,x-1));total+=x;global[f]+=v->volume[e-1]*x;
        }
        partition=fmax(partition,fabs(total-1));
    }
    sums(E,global,VF);bounds=largest(E,bounds);partition=largest(E,partition);
    if(bounds>1e-12 || partition>1e-12)pices_fail(E,"active species fractions violate their bounds");
    for(f=0;f<VF;f++)if(fabs(global[f]-E->pices.volume_baseline[f]-E->pices.volume_source[f])>VOLUME_BUDGET_TOL*v->shell)
        pices_fail(E,"active species integral violates weighted history");
}
void pices_composition_volume_commit(struct All_variables *E)
{
    if(!pices_composition_volume_enabled(E) || !volume_state.pending)return;
    write_rows(E,"accepted");volume_state.pending=0;
    if(!E->parallel.me) {
        struct Pices_volume_summary *s=&volume_state.last;
        fprintf(E->fp,"PICES_VOLUME_CONSERVE step=%d mode=%d feasible=%d iterations=%d hops=%d C24_target=%.17g C24_corrected=%.17g minimum_transfer=%.17g row_error=%.17g species_error=%.17g reason=%s\n",s->step,E->pices.volume_mode,s->feasible,s->iterations,s->hops,s->target[24],s->corrected[24],s->minimum_transfer_nd,s->row_error,s->species_error,s->reason);fflush(E->fp);
    }
}
void pices_composition_volume_reset(struct All_variables *E,int rebase)
{
    struct Volume_state *v=&volume_state;
    if(!pices_composition_volume_enabled(E))return;
    if(v->steps)fclose(v->steps);if(v->flavors)fclose(v->flavors);
    v->steps=v->flavors=NULL;v->files_open=v->pending=0;memset(&v->last,0,sizeof(v->last));
    memset(v->local_source,0,sizeof(v->local_source));memset(v->step_source,0,sizeof(v->step_source));
    if(rebase) {
        double totals[VF];
        if(!v->geometry || !E->pices.volume_initialized)pices_fail(E,"cannot rebase missing volume state");
        deposit(E,totals);memcpy(E->pices.volume_baseline,totals,sizeof(totals));memset(E->pices.volume_source,0,sizeof(totals));
    }
}
const struct Pices_volume_summary *pices_composition_volume_latest(struct All_variables *E)
{ (void)E;return &volume_state.last; }
double pices_composition_volume_fraction(struct All_variables *E,int element,int flavor)
{
    struct Volume_state *v=&volume_state;(void)E;
    if(!v->last.initialized || !v->last.feasible || element<1 || element>v->nel || flavor<0 || flavor>=VF)return NAN;
    return v->corrected[(size_t)(element-1)*VF+flavor]/v->volume[element-1];
}
