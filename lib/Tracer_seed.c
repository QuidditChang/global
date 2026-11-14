/* Opt-in initialization only. The quota is fixed per cell (B), or proportional
 * to FE volume (C). Both use the same sampler and private per-cell RNG. Positions
 * are uniform in the spherical wedge used by the actual PICES host predicate.
 * These geometries differ; their difference is measured, never hidden. */
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include "global_defs.h"
#include "element_definitions.h"
#include "pices.h"
#include "tracer_seed.h"
#include "tracer_seed_math.h"

void allocate_tracer_arrays(struct All_variables *,int,int);
void get_global_shape_fn(struct All_variables *,int,struct Shape_function *,
    struct Shape_function_dx *,struct Shape_function_dA *,int,int,double [4][9],int,int);

static struct {
    double *volume,*host_volume;
    int *quota;
    long long target;
    double shell;
    int pending;
} seeding;

static void *seed_alloc(struct All_variables *E,size_t n,size_t width) {
    void *p=calloc(n?n:1,width);
    if(!p)pices_fail(E,"tracer seed allocation failed");
    return p;
}

void tracer_seed_configure(struct All_variables *E,const char *mode,int seed) {
    if(!strcmp(mode,"legacy"))E->trace.seed_distribution=0;
    else if(!strcmp(mode,"volume_uniform"))E->trace.seed_distribution=1;
    else if(!strcmp(mode,"cell_uniform"))E->trace.seed_distribution=2;
    else pices_fail(E,"tracer_seed_distribution must be legacy, cell_uniform or volume_uniform");
    if(seed<0)pices_fail(E,"tracer_seed must be a nonnegative integer");
    E->trace.seed=seed;
}

void tracer_seed_check_ic(struct All_variables *E) {
    if(E->trace.seed_distribution && E->trace.ic_method!=0)
        pices_fail(E,"cell_uniform/volume_uniform require tracer_ic_method=0; paired-checkpoint restart keeps this setting and restores saved particles");
}

static const char *seed_mode(struct All_variables *E) {
    return E->trace.seed_distribution==2?"cell_uniform":"volume_uniform";
}
static const char *seed_expectation(struct All_variables *E) {
    return E->trace.seed_distribution==2?"tracers_per_element":"total_times_FE_fraction";
}

static void geometry(struct All_variables *E,int e,double uv[4][2],
                     double area[2],double *lo,double *hi,double *solid) {
    double unit[4][3];int a,d;
    *lo=E->sx[1][3][E->ien[1][e].node[1]];
    *hi=E->sx[1][3][E->ien[1][e].node[5]];
    if(!isfinite(*lo)||!isfinite(*hi)||!(*hi>*lo)||!(*lo>0))
        pices_fail(E,"invalid radial shell in tracer seeder");
    for(a=0;a<4;a++) {
        int n=E->ien[1][e].node[a+1],top=E->ien[1][e].node[a+5];
        double theta=E->sx[1][1][n],phi=E->sx[1][2][n];
        if(fabs(E->sx[1][3][n]-*lo)>64*DBL_EPSILON ||
           fabs(E->sx[1][3][top]-*hi)>64*DBL_EPSILON)
            pices_fail(E,"tracer seeder requires concentric radial faces");
        uv[a][0]=E->trace.UV[1][1][n];uv[a][1]=E->trace.UV[1][2][n];
        unit[a][0]=sin(theta)*cos(phi);unit[a][1]=sin(theta)*sin(phi);unit[a][2]=cos(theta);
        for(d=0;d<2;d++)if(!isfinite(uv[a][d]))pices_fail(E,"invalid seed gnomonic geometry");
    }
    area[0]=tracer_seed_area(uv[0],uv[1],uv[2]);
    area[1]=tracer_seed_area(uv[0],uv[2],uv[3]);
    *solid=tracer_seed_solid_angle(unit[0],unit[1],unit[2])+
           tracer_seed_solid_angle(unit[0],unit[2],unit[3]);
    if(!(area[0]>0)||!(area[1]>0)||!isfinite(*solid)||!(*solid>0))
        pices_fail(E,"degenerate tracer seed cell");
}

void tracer_seed_generate(struct All_variables *E) {
    int e,g,p,a,d,rank=E->parallel.me,nproc=E->parallel.nproc,local_count=0;
    long long cells=E->lmesh.nel,global_cells,generated,low[4]={0,0,0,0};
    double low_volume[3]={0,0,0};
    long double local_volume=0,total_volume=0,prefix=0,cumulative=0,offset;
    long double *rank_volumes=NULL,*starts=NULL,bounds[2];
    uint64_t offset_state=(uint64_t)E->trace.seed;
    long long previous,next;
    double c[3],east[3],north[3],theta_f=E->trace.UV[1][1][0],phi_f=E->trace.UV[1][2][0];
    struct Shape_function GN;struct Shape_function_dx dx;struct Shape_function_dA omega;double rtf[4][9];
    if(!E->pices.enabled||!E->pices.eba||E->sphere.caps!=12||E->sphere.caps_per_proc!=1||E->control.pseudo_free_surf)
        pices_fail(E,"cell_uniform/volume_uniform seeding requires full-sphere EBA PICES with one cap/rank and fixed boundaries");
    if(E->trace.ic_method!=0||E->trace.itperel<1||E->pices.volume_mode)
        pices_fail(E,"cell_uniform/volume_uniform seeding requires fresh random IC, positive tracers_per_element and pices_volume_mode=off");
    if(E->trace.nflavors!=25||!E->composition.on||E->composition.ibuoy_type!=1)
        pices_fail(E,"cell_uniform/volume_uniform seeding currently requires the 25-flavor ratio-composition experiment");
    if(seeding.pending)pices_fail(E,"duplicate tracer seed initialization");
    MPI_Allreduce(&cells,&global_cells,1,MPI_LONG_LONG,MPI_SUM,E->parallel.world);
    seeding.target=global_cells*(long long)E->trace.itperel;
    /* Legacy tracer/checkpoint counters are signed int. Do not silently wrap. */
    if(global_cells<=0||seeding.target<=0||seeding.target>INT_MAX)
        pices_fail(E,"global tracer target exceeds the existing signed-int counter range");
    seeding.volume=seed_alloc(E,E->lmesh.nel,sizeof(double));
    seeding.host_volume=seed_alloc(E,E->lmesh.nel,sizeof(double));
    seeding.quota=seed_alloc(E,E->lmesh.nel,sizeof(int));
    for(e=1;e<=E->lmesh.nel;e++) {
        double uv[4][2],area[2],lo,hi,solid;
        get_global_shape_fn(E,e,&GN,&dx,&omega,0,1,rtf,E->mesh.levmax,1);
        for(g=1;g<=8;g++) {
            double w=omega.vpt[g]*g_point[g].weight[2];
            if(!isfinite(w)||!(w>0))pices_fail(E,"nonpositive FE seed quadrature");
            seeding.volume[e-1]+=w;
        }
        local_volume+=(long double)seeding.volume[e-1];
        geometry(E,e,uv,area,&lo,&hi,&solid);
        seeding.host_volume[e-1]=solid*(hi*hi*hi-lo*lo*lo)/3;
    }
    if(rank==0) {
        rank_volumes=seed_alloc(E,nproc,sizeof(long double));
        starts=seed_alloc(E,2*nproc,sizeof(long double));
    }
    MPI_Gather(&local_volume,1,MPI_LONG_DOUBLE,rank_volumes,1,MPI_LONG_DOUBLE,0,E->parallel.world);
    if(rank==0)for(p=0;p<nproc;p++){starts[2*p]=total_volume;total_volume+=rank_volumes[p];starts[2*p+1]=total_volume;}
    MPI_Bcast(&total_volume,1,MPI_LONG_DOUBLE,0,E->parallel.world);
    MPI_Scatter(starts,2,MPI_LONG_DOUBLE,bounds,2,MPI_LONG_DOUBLE,0,E->parallel.world);
    prefix=bounds[0];
    free(starts);free(rank_volumes);
    if(!isfinite(total_volume)||!(total_volume>0))pices_fail(E,"invalid global seed FE volume");
    seeding.shell=(double)total_volume;
    offset=tracer_seed_unit(&offset_state);
    previous=(long long)floorl((prefix/total_volume)*seeding.target+offset);
    for(e=1;e<=E->lmesh.nel;e++) {
        long long n;
        double lambda;
        if(E->trace.seed_distribution==2) {
            n=E->trace.itperel;lambda=E->trace.itperel;
        } else {
            /* Preserve C's arithmetic and stream exactly. */
            lambda=(double)(seeding.target*(long double)seeding.volume[e-1]/total_volume);
            cumulative+=(long double)seeding.volume[e-1];
            /* The two ranks sharing a boundary use its identical root-owned
             * bit pattern, not differently associated floating-point sums. */
            next=(long long)floorl(((e==E->lmesh.nel?bounds[1]:prefix+cumulative)/total_volume)*seeding.target+offset);
            n=next-previous;previous=next;
            if(fabs(n-lambda)>1+1e-6)pices_fail(E,"seed quota is not a floor/ceil of its FE-volume expectation");
        }
        if(n<0||n>INT_MAX-local_count)pices_fail(E,"invalid integer tracer seed quota");
        if(lambda<1){low[0]++;low_volume[0]+=seeding.volume[e-1];}
        if(lambda<4){low[1]++;low_volume[1]+=seeding.volume[e-1];}
        if(lambda<8){low[2]++;low_volume[2]+=seeding.volume[e-1];}
        if(!n)low[3]++;
        seeding.quota[e-1]=(int)n;local_count+=(int)n;
    }
    generated=local_count;
    MPI_Allreduce(MPI_IN_PLACE,&generated,1,MPI_LONG_LONG,MPI_SUM,E->parallel.world);
    if(generated!=seeding.target)pices_fail(E,"global integer seed quotas do not sum to target");
    MPI_Allreduce(MPI_IN_PLACE,low,4,MPI_LONG_LONG,MPI_SUM,E->parallel.world);
    MPI_Allreduce(MPI_IN_PLACE,low_volume,3,MPI_DOUBLE,MPI_SUM,E->parallel.world);
    if(rank==0) {
        fprintf(stderr,"TRACER_SEED_QUOTA mode=%s expectation=%s target=%lld lambda_lt1_cells=%lld lambda_lt4_cells=%lld lambda_lt8_cells=%lld lambda_lt1_volume_fraction=%.17g lambda_lt4_volume_fraction=%.17g lambda_lt8_volume_fraction=%.17g empty_cells=%lld\n",
            seed_mode(E),seed_expectation(E),seeding.target,low[0],low[1],low[2],low_volume[0]/seeding.shell,low_volume[1]/seeding.shell,low_volume[2]/seeding.shell,low[3]);
        fflush(stderr);
    }
    if(low[3])pices_fail(E,"volume-proportional seed quotas leave empty initial cells; increase total particles explicitly, no minimum quota is added");
    /* The legacy allocator is one-based and adds floor(N/4) spare entries. */
    if(local_count<4||local_count>INT_MAX/5*4-4)
        pices_fail(E,"tracer seed rank count incompatible with legacy tracer allocation");
    allocate_tracer_arrays(E,1,local_count);
    c[0]=sin(theta_f)*cos(phi_f);c[1]=sin(theta_f)*sin(phi_f);c[2]=cos(theta_f);
    east[0]=-sin(phi_f);east[1]=cos(phi_f);east[2]=0;
    north[0]=-cos(theta_f)*cos(phi_f);north[1]=-cos(theta_f)*sin(phi_f);north[2]=sin(theta_f);
    for(e=1;e<=E->lmesh.nel;e++) {
        double uv[4][2],area[2],lo,hi,solid;
        uint64_t state=(uint64_t)E->trace.seed^((uint64_t)(rank+1)<<32)^(uint64_t)e;
        geometry(E,e,uv,area,&lo,&hi,&solid);
        for(a=0;a<seeding.quota[e-1];a++) {
            double u,v,r,theta,phi,xyz[3],scale;
            if(!tracer_seed_direction(&state,uv,area,&u,&v))pices_fail(E,"seed spherical-area rejection did not converge");
            r=tracer_seed_radius(&state,lo,hi);scale=r/sqrt(1+u*u+v*v);
            for(d=0;d<3;d++)xyz[d]=scale*(c[d]+u*east[d]+v*north[d]);
            theta=atan2(hypot(xyz[0],xyz[1]),xyz[2]);phi=atan2(xyz[1],xyz[0]);if(phi<0)phi+=2*M_PI;
            if(!pices_seed_host_valid(E,1,e,theta,phi,r))pices_fail(E,"tracer seed point failed its actual PICES host predicate");
            p=++E->trace.ntracers[1];E->trace.ielement[1][p]=e;
            E->trace.basicq[1][0][p]=theta;E->trace.basicq[1][1][p]=phi;E->trace.basicq[1][2][p]=r;
            for(d=0;d<3;d++)E->trace.basicq[1][d+3][p]=xyz[d];
        }
    }
    if(E->trace.ntracers[1]!=local_count)pices_fail(E,"tracer seed local count mismatch");
    seeding.pending=1;
}

void tracer_seed_report(struct All_variables *E) {
    int e,p,f,rank=E->parallel.me,*count,*flavor24;
    long long generated,local_min=LLONG_MAX,local_max=0,empty=0,rank_min,rank_max,flavors[25]={0};
    double lambda_min=DBL_MAX,lambda_max=0,geometry_error=0,quota_error=0,host_shell=0,volume24=0;
    if(!seeding.pending)return; /* Restarts never regenerate or rewrite provenance. */
    count=seed_alloc(E,E->lmesh.nel,sizeof(int));flavor24=seed_alloc(E,E->lmesh.nel,sizeof(int));
    for(p=1;p<=E->trace.ntracers[1];p++) {
        e=E->trace.ielement[1][p];
        if(e<1||e>E->lmesh.nel)pices_fail(E,"invalid host after seed relocation");
        count[e-1]++;
        if(E->trace.nflavors) {
            double label=E->trace.extraq[1][0][p];f=(int)label;
            if(label!=f||f<0||f>=25)pices_fail(E,"invalid initial seed flavor");
            flavors[f]++;if(f==24)flavor24[e-1]++;
        }
    }
    for(e=0;e<E->lmesh.nel;e++) {
        double lambda=E->trace.seed_distribution==2?E->trace.itperel:seeding.target*seeding.volume[e]/seeding.shell;
        if(count[e]!=seeding.quota[e])pices_fail(E,"initial host relocation changed a tracer seed quota");
        if(count[e]<local_min)local_min=count[e];if(count[e]>local_max)local_max=count[e];
        if(!count[e])empty++;
        lambda_min=fmin(lambda_min,lambda);lambda_max=fmax(lambda_max,lambda);
        geometry_error=fmax(geometry_error,fabs(seeding.volume[e]/seeding.host_volume[e]-1));
        quota_error=fmax(quota_error,fabs(count[e]-lambda));
        host_shell+=seeding.host_volume[e];
        if(count[e])volume24+=seeding.volume[e]*flavor24[e]/count[e];
    }
    fprintf(E->trace.fpt,"TRACER_SEED_RANK mode=%s expectation=%s seed=%d rank=%d particles=%d min_cell=%lld max_cell=%lld empty_cells=%lld lambda_min=%.17g lambda_max=%.17g\n",
        seed_mode(E),seed_expectation(E),E->trace.seed,rank,E->trace.ntracers[1],local_min,local_max,empty,lambda_min,lambda_max);
    fflush(E->trace.fpt);
    generated=rank_min=rank_max=E->trace.ntracers[1];
#define REDUCE(value,type,op) MPI_Allreduce(MPI_IN_PLACE,&value,1,type,op,E->parallel.world)
    REDUCE(generated,MPI_LONG_LONG,MPI_SUM);REDUCE(rank_min,MPI_LONG_LONG,MPI_MIN);REDUCE(rank_max,MPI_LONG_LONG,MPI_MAX);
    REDUCE(local_min,MPI_LONG_LONG,MPI_MIN);REDUCE(local_max,MPI_LONG_LONG,MPI_MAX);REDUCE(empty,MPI_LONG_LONG,MPI_SUM);
    REDUCE(lambda_min,MPI_DOUBLE,MPI_MIN);REDUCE(lambda_max,MPI_DOUBLE,MPI_MAX);
    REDUCE(geometry_error,MPI_DOUBLE,MPI_MAX);REDUCE(quota_error,MPI_DOUBLE,MPI_MAX);
    REDUCE(host_shell,MPI_DOUBLE,MPI_SUM);REDUCE(volume24,MPI_DOUBLE,MPI_SUM);
#undef REDUCE
    MPI_Allreduce(MPI_IN_PLACE,flavors,25,MPI_LONG_LONG,MPI_SUM,E->parallel.world);
    if(generated!=seeding.target||quota_error>(E->trace.seed_distribution==2?0:1+1e-6))
        pices_fail(E,"tracer seed final quota audit failed");
    if(rank==0) {
        FILE *files[2]={stderr,E->fp};int k;
        for(k=0;k<2;k++) {
            fprintf(files[k],"TRACER_SEED_GLOBAL mode=%s seed=%d target=%lld generated=%lld elements=%lld fe_volume=%.17g host_volume=%.17g min_cell=%lld max_cell=%lld empty_cells=%lld rank_min=%lld rank_max=%lld lambda_min=%.17g lambda_max=%.17g quota_error_max=%.17g fe_host_relative_error_max=%.17g quota=%s expectation=%s sampler=gnomonic_area_r3_v1 ordering=rank_element\n",
                seed_mode(E),E->trace.seed,seeding.target,generated,seeding.target/E->trace.itperel,seeding.shell,host_shell,local_min,local_max,empty,rank_min,rank_max,lambda_min,lambda_max,quota_error,geometry_error,
                E->trace.seed_distribution==2?"constant_per_cell_v1":"systematic_prefix_v1",seed_expectation(E));
            fprintf(files[k],"TRACER_SEED_FLAVOR flavor=24 count=%lld count_fraction=%.17g volume_fraction=%.17g scope=initial_geology\n",flavors[24],(double)flavors[24]/generated,volume24/seeding.shell);
            fflush(files[k]);
        }
    }
    free(count);free(flavor24);free(seeding.volume);free(seeding.host_volume);free(seeding.quota);
    memset(&seeding,0,sizeof(seeding));
}
