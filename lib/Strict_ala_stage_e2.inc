/* Strict-ALA Stage E2 is deliberately compiled into Stokes_flow_Incomp.c.
 * It observes selected independently assembled pressure residuals and runs
 * all collective analysis after the production FGMRES solve has stopped. */

#define ALA_E2_SNAPSHOT_COUNT 19
#define ALA_E2_LATE_FIRST 6
#define ALA_E2_LATE_COUNT 13
#define ALA_E2_MAX_MODES 6
#define ALA_E2_MAGIC 0x414c413245325331ULL

struct ala_e2_file_header {
    unsigned long long magic;
    unsigned long long checksum;
    int version;
    int rank;
    int ranks;
    int caps;
    int npno;
    int count;
    int iterations[ALA_E2_SNAPSHOT_COUNT];
};

struct ala_e2_observer {
    int enabled;
    int lev;
    int count;
    int iterations[ALA_E2_SNAPSHOT_COUNT];
    double **snapshot[ALA_E2_SNAPSHOT_COUNT];
    double capture_seconds;
    char case_name[8];
    char shared_dir[512];
    char output_dir[512];
};

static const int ala_e2_schedule[ALA_E2_SNAPSHOT_COUNT]={
    0,5,10,20,30,35,40,42,44,46,48,49,50,51,52,54,56,58,60
};

static unsigned long long ala_e2_checksum_bytes(
    unsigned long long hash,const unsigned char *bytes,size_t count)
{
    size_t i;
    for(i=0;i<count;i++) {
        hash^=(unsigned long long)bytes[i];
        hash*=1099511628211ULL;
    }
    return hash;
}

static unsigned long long ala_e2_local_checksum(
    struct All_variables *E,struct ala_e2_observer *o)
{
    int i,m;
    unsigned long long hash=1469598103934665603ULL;
    for(i=0;i<o->count;i++) {
        hash=ala_e2_checksum_bytes(hash,
            (const unsigned char *)&o->iterations[i],sizeof(int));
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            hash=ala_e2_checksum_bytes(hash,
                (const unsigned char *)&o->snapshot[i][m][1],
                E->lmesh.NPNO[o->lev]*sizeof(double));
    }
    return hash;
}

static int ala_e2_schedule_index(int iteration)
{
    int i;
    for(i=0;i<ALA_E2_SNAPSHOT_COUNT;i++)
        if(ala_e2_schedule[i]==iteration) return i;
    return -1;
}

static FILE *ala_e2_open_root(struct All_variables *E,
    struct ala_e2_observer *o,const char *name,const char *mode)
{
    char path[1024];
    FILE *fp;
    snprintf(path,sizeof(path),"%s/%s",o->output_dir,name);
    fp=fopen(path,mode);
    if(fp==NULL) myerror(E,"Unable to open strict-ALA Stage-E2 output");
    return fp;
}

static void *strict_ala_stage_e2_begin(struct All_variables *E,int lev)
{
    struct ala_e2_observer *o;
    const char *required=getenv("STRICT_ALA_STAGE_E2_REQUIRED");
    const char *shared=getenv("STRICT_ALA_STAGE_E2_SHARED_DIR");
    const char *output=getenv("STRICT_ALA_STAGE_E2_OUTPUT_DIR");
    const char *name=getenv("STRICT_ALA_CASE");
    if(required==NULL) return NULL;
    if(shared==NULL || output==NULL || name==NULL ||
       (strcmp(name,"E0")!=0 && strcmp(name,"E1")!=0))
        myerror(E,"Incomplete strict-ALA Stage-E2 environment contract");
    o=(struct ala_e2_observer *)calloc(1,sizeof(*o));
    if(o==NULL) myerror(E,"Unable to allocate Stage-E2 observer");
    o->enabled=1;
    o->lev=lev;
    snprintf(o->case_name,sizeof(o->case_name),"%s",name);
    snprintf(o->shared_dir,sizeof(o->shared_dir),"%s",shared);
    snprintf(o->output_dir,sizeof(o->output_dir),"%s",output);
    if(E->parallel.me==0) {
        fprintf(E->fp,"STRICT_ALA_STAGE_E2_BEGIN case=%s schedule=predeclared19 "
                "online_collectives=0 residual=true_continuity\n",name);
        fflush(E->fp);
    }
    return o;
}

static void strict_ala_stage_e2_capture(struct All_variables *E,void *state,
    double **residual,int iteration)
{
    struct ala_e2_observer *o=(struct ala_e2_observer *)state;
    int slot,m,e;
    double start;
    if(o==NULL || !o->enabled) return;
    slot=ala_e2_schedule_index(iteration);
    if(slot<0) return;
    if(slot!=o->count)
        myerror(E,"Stage-E2 snapshot schedule is incomplete or out of order");
    start=CPU_time0();
    o->snapshot[slot]=ala_schur_alloc_pressure(E,o->lev);
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(e=1;e<=E->lmesh.NPNO[o->lev];e++)
            o->snapshot[slot][m][e]=residual[m][e];
    o->iterations[slot]=iteration;
    o->count++;
    o->capture_seconds+=CPU_time0()-start;
}

static void ala_e2_snapshot_path(struct ala_e2_observer *o,int rank,
    char *path,size_t size)
{
    snprintf(path,size,"%s/%s.rank%06d.bin",o->shared_dir,o->case_name,rank);
}

static unsigned long long ala_e2_write_snapshots(
    struct All_variables *E,struct ala_e2_observer *o)
{
    struct ala_e2_file_header h;
    char path[1024];
    FILE *fp;
    int i,m;
    unsigned long long local,global;
    if(o->count!=ALA_E2_SNAPSHOT_COUNT)
        myerror(E,"Strict-ALA Stage-E2 did not capture all snapshots");
    memset(&h,0,sizeof(h));
    h.magic=ALA_E2_MAGIC;
    h.version=1;
    h.rank=E->parallel.me;
    h.ranks=E->parallel.nproc;
    h.caps=E->sphere.caps_per_proc;
    h.npno=E->lmesh.NPNO[o->lev];
    h.count=o->count;
    for(i=0;i<o->count;i++) h.iterations[i]=o->iterations[i];
    local=ala_e2_local_checksum(E,o);
    h.checksum=local;
    ala_e2_snapshot_path(o,E->parallel.me,path,sizeof(path));
    fp=fopen(path,"wb");
    if(fp==NULL) myerror(E,"Unable to create Stage-E2 distributed snapshot");
    if(fwrite(&h,sizeof(h),1,fp)!=1)
        myerror(E,"Unable to write Stage-E2 snapshot header");
    for(i=0;i<o->count;i++)
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            if(fwrite(&o->snapshot[i][m][1],sizeof(double),h.npno,fp)
               !=(size_t)h.npno)
                myerror(E,"Unable to write Stage-E2 snapshot payload");
    if(fclose(fp)!=0) myerror(E,"Unable to close Stage-E2 snapshot");
    fp=fopen(path,"rb");
    if(fp==NULL || fread(&h,sizeof(h),1,fp)!=1 || h.checksum!=local)
        myerror(E,"Unable to verify Stage-E2 snapshot header");
    {
        unsigned long long verify=1469598103934665603ULL;
        double *buffer=(double *)calloc(h.npno,sizeof(double));
        if(buffer==NULL) myerror(E,"Unable to allocate Stage-E2 verify buffer");
        for(i=0;i<h.count;i++) {
            verify=ala_e2_checksum_bytes(verify,
                (const unsigned char *)&h.iterations[i],sizeof(int));
            for(m=1;m<=h.caps;m++) {
                if(fread(buffer,sizeof(double),h.npno,fp)!=(size_t)h.npno)
                    myerror(E,"Truncated Stage-E2 snapshot during verification");
                verify=ala_e2_checksum_bytes(verify,
                    (const unsigned char *)buffer,h.npno*sizeof(double));
            }
        }
        free(buffer);
        if(verify!=local || fgetc(fp)!=EOF)
            myerror(E,"Stage-E2 snapshot payload verification failed");
    }
    fclose(fp);
    MPI_Allreduce(&local,&global,1,MPI_UNSIGNED_LONG_LONG,MPI_BXOR,
                  E->parallel.world);
    return global;
}

static struct ala_e2_observer *ala_e2_read_case(
    struct All_variables *E,struct ala_e2_observer *current,const char *name)
{
    struct ala_e2_observer *o;
    struct ala_e2_file_header h;
    char path[1024];
    FILE *fp;
    int i,m;
    unsigned long long checksum;
    o=(struct ala_e2_observer *)calloc(1,sizeof(*o));
    if(o==NULL) myerror(E,"Unable to allocate Stage-E2 archived observer");
    o->enabled=1;
    o->lev=current->lev;
    snprintf(o->case_name,sizeof(o->case_name),"%s",name);
    snprintf(o->shared_dir,sizeof(o->shared_dir),"%s",current->shared_dir);
    snprintf(o->output_dir,sizeof(o->output_dir),"%s",current->output_dir);
    ala_e2_snapshot_path(o,E->parallel.me,path,sizeof(path));
    fp=fopen(path,"rb");
    if(fp==NULL) myerror(E,"Missing Stage-E2 archived E0 snapshot");
    if(fread(&h,sizeof(h),1,fp)!=1 || h.magic!=ALA_E2_MAGIC || h.version!=1 ||
       h.rank!=E->parallel.me || h.ranks!=E->parallel.nproc ||
       h.caps!=E->sphere.caps_per_proc ||
       h.npno!=E->lmesh.NPNO[o->lev] || h.count!=ALA_E2_SNAPSHOT_COUNT)
        myerror(E,"Invalid Stage-E2 archived snapshot header");
    o->count=h.count;
    for(i=0;i<o->count;i++) {
        if(h.iterations[i]!=ala_e2_schedule[i])
            myerror(E,"Stage-E2 archived schedule mismatch");
        o->iterations[i]=h.iterations[i];
        o->snapshot[i]=ala_schur_alloc_pressure(E,o->lev);
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            if(fread(&o->snapshot[i][m][1],sizeof(double),h.npno,fp)
               !=(size_t)h.npno)
                myerror(E,"Truncated Stage-E2 archived snapshot");
    }
    if(fgetc(fp)!=EOF) myerror(E,"Trailing bytes in Stage-E2 snapshot");
    fclose(fp);
    checksum=ala_e2_local_checksum(E,o);
    if(checksum!=h.checksum) myerror(E,"Stage-E2 snapshot checksum mismatch");
    return o;
}

static void ala_e2_free_observer(struct All_variables *E,
    struct ala_e2_observer *o)
{
    int i;
    if(o==NULL) return;
    for(i=0;i<o->count;i++) ala_schur_free_field(E,o->snapshot[i]);
    free(o);
}

static void ala_e2_covariance(struct All_variables *E,
    struct ala_e2_observer **cases,int *indices,int n,double *cov)
{
    int i,j,m,e;
    double *local=(double *)calloc(n*n,sizeof(double));
    if(local==NULL) myerror(E,"Unable to allocate Stage-E2 covariance");
    for(i=0;i<n;i++)
        for(j=i;j<n;j++) {
            struct ala_e2_observer *oi=cases[i];
            struct ala_e2_observer *oj=cases[j];
            double value=0.0;
            for(m=1;m<=E->sphere.caps_per_proc;m++)
                for(e=1;e<=E->lmesh.NPNO[oi->lev];e++)
                    value+=oi->snapshot[indices[i]][m][e]
                          *oj->snapshot[indices[j]][m][e];
            local[i*n+j]=local[j*n+i]=value;
        }
    MPI_Allreduce(local,cov,n*n,MPI_DOUBLE,MPI_SUM,E->parallel.world);
    free(local);
}

static void ala_e2_sort_eigensystem(double *values,double *vectors,int n)
{
    int i,j,k;
    for(i=0;i<n-1;i++) {
        k=i;
        for(j=i+1;j<n;j++) if(values[j]>values[k]) k=j;
        if(k!=i) {
            double t=values[i]; values[i]=values[k]; values[k]=t;
            for(j=0;j<n;j++) {
                t=vectors[j*n+i]; vectors[j*n+i]=vectors[j*n+k];
                vectors[j*n+k]=t;
            }
        }
    }
}

static int ala_e2_pod(double *cov,int n,int normalize,double *values,
    double *vectors,double *norms)
{
    double *a;
    double trace=0.0,scale,total=0.0;
    int i,j;
    a=(double *)calloc(n*n,sizeof(double));
    if(a==NULL) return 0;
    for(i=0;i<n;i++) {
        norms[i]=sqrt(max(cov[i*n+i],0.0));
        trace+=normalize ? (norms[i]>0.0 ? 1.0 : 0.0) : cov[i*n+i];
    }
    scale=max(trace,1.0e-300);
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            a[i*n+j]=(normalize
                ? cov[i*n+j]/max(norms[i]*norms[j],1.0e-300)
                : cov[i*n+j])/scale;
    ala_geneo_jacobi_eigensolve(a,vectors,values,n);
    ala_e2_sort_eigensystem(values,vectors,n);
    for(i=0;i<n;i++) {
        values[i]=max(values[i]*scale,0.0);
        total+=values[i];
    }
    free(a);
    return total>0.0;
}

static void ala_e2_write_manifest(struct All_variables *E,
    struct ala_e2_observer *o,unsigned long long checksum)
{
    FILE *fp;
    int i,m;
    double cov[ALA_E2_SNAPSHOT_COUNT*ALA_E2_SNAPSHOT_COUNT];
    unsigned long long local_hash,global_hash;
    struct ala_e2_observer *cases[ALA_E2_SNAPSHOT_COUNT];
    int indices[ALA_E2_SNAPSHOT_COUNT];
    for(i=0;i<o->count;i++) { cases[i]=o; indices[i]=i; }
    ala_e2_covariance(E,cases,indices,o->count,cov);
    fp=E->parallel.me==0 ? ala_e2_open_root(E,o,
        "strict_ala_stage_E2_snapshots_manifest.csv","a") : NULL;
    for(i=0;i<o->count;i++) {
        local_hash=1469598103934665603ULL;
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            local_hash=ala_e2_checksum_bytes(local_hash,
                (const unsigned char *)&o->snapshot[i][m][1],
                E->lmesh.NPNO[o->lev]*sizeof(double));
        MPI_Allreduce(&local_hash,&global_hash,1,MPI_UNSIGNED_LONG_LONG,
                      MPI_BXOR,E->parallel.world);
        if(fp!=NULL)
            fprintf(fp,"%s,%d,%.17e,%016llx,%s/%s.rank%%06d.bin,%d,%d\n",
                o->case_name,o->iterations[i],
                sqrt(max(cov[i*o->count+i],0.0)),global_hash,o->shared_dir,
                o->case_name,E->parallel.nproc,E->lmesh.NPNO[o->lev]);
    }
    if(fp!=NULL) fclose(fp);
    (void)checksum;
}

static void ala_e2_write_gram(struct All_variables *E,
    struct ala_e2_observer *e0,struct ala_e2_observer *e1)
{
    struct ala_e2_observer *cases[2*ALA_E2_SNAPSHOT_COUNT];
    int indices[2*ALA_E2_SNAPSHOT_COUNT];
    double cov[4*ALA_E2_SNAPSHOT_COUNT*ALA_E2_SNAPSHOT_COUNT];
    FILE *fp;
    int i,j,n=2*ALA_E2_SNAPSHOT_COUNT;
    for(i=0;i<ALA_E2_SNAPSHOT_COUNT;i++) {
        cases[i]=e0; indices[i]=i;
        cases[i+ALA_E2_SNAPSHOT_COUNT]=e1;
        indices[i+ALA_E2_SNAPSHOT_COUNT]=i;
    }
    ala_e2_covariance(E,cases,indices,n,cov);
    if(E->parallel.me!=0) return;
    fp=ala_e2_open_root(E,e1,"strict_ala_stage_E2_snapshot_gram.csv","w");
    fprintf(fp,"case_i,iteration_i,case_j,iteration_j,inner_product\n");
    for(i=0;i<n;i++) for(j=0;j<n;j++)
        fprintf(fp,"%s,%d,%s,%d,%.17e\n",cases[i]->case_name,
            cases[i]->iterations[indices[i]],cases[j]->case_name,
            cases[j]->iterations[indices[j]],cov[i*n+j]);
    fclose(fp);
}

static int ala_e2_write_pod_set(struct All_variables *E,
    struct ala_e2_observer *out,struct ala_e2_observer **cases,int *indices,
    int n,const char *case_label,const char *basis,int normalize,
    double *ret_values,double *ret_vectors,double *ret_norms)
{
    double *cov,*values,*vectors,*norms,total=0.0,cumulative=0.0;
    FILE *fp;
    int i,ok,n50=0,n80=0,n90=0,n95=0;
    cov=(double *)calloc(n*n,sizeof(double));
    values=ret_values ? ret_values : (double *)calloc(n,sizeof(double));
    vectors=ret_vectors ? ret_vectors : (double *)calloc(n*n,sizeof(double));
    norms=ret_norms ? ret_norms : (double *)calloc(n,sizeof(double));
    if(cov==NULL || values==NULL || vectors==NULL || norms==NULL)
        myerror(E,"Unable to allocate Stage-E2 POD workspace");
    ala_e2_covariance(E,cases,indices,n,cov);
    ok=ala_e2_pod(cov,n,normalize,values,vectors,norms);
    for(i=0;i<n;i++) total+=values[i];
    for(i=0;i<n;i++) {
        cumulative+=values[i]/max(total,1.0e-300);
        if(!n50 && cumulative>=0.50) n50=i+1;
        if(!n80 && cumulative>=0.80) n80=i+1;
        if(!n90 && cumulative>=0.90) n90=i+1;
        if(!n95 && cumulative>=0.95) n95=i+1;
    }
    cumulative=0.0;
    if(E->parallel.me==0) {
        fp=ala_e2_open_root(E,out,"strict_ala_stage_E2_pod.csv","a");
        for(i=0;i<n;i++) {
            double fraction=values[i]/max(total,1.0e-300);
            cumulative+=fraction;
            fprintf(fp,"%s,%s,%d,%.17e,%.17e,%.17e,%d,%d,%d,%d,%s\n",
                case_label,basis,i+1,sqrt(max(values[i],0.0)),fraction,
                cumulative,n50,n80,n90,n95,
                normalize ? "direction_case_balanced" : "physical_energy");
        }
        fclose(fp);
    }
    free(cov);
    if(!ret_values) free(values);
    if(!ret_vectors) free(vectors);
    if(!ret_norms) free(norms);
    return ok;
}

static double **ala_e2_build_mode(struct All_variables *E,
    struct ala_e2_observer **cases,int *indices,int n,double *vectors,
    double *norms,int mode)
{
    double **q=ala_schur_alloc_pressure(E,cases[0]->lev);
    double norm;
    int i,m,e;
    for(i=0;i<n;i++) {
        double coefficient=vectors[i*n+mode]/max(norms[i],1.0e-300);
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=1;e<=E->lmesh.NPNO[cases[0]->lev];e++)
                q[m][e]+=coefficient*cases[i]->snapshot[indices[i]][m][e];
    }
    norm=sqrt(max(global_pdot(E,q,q,cases[0]->lev),0.0));
    if(!isfinite(norm) || norm<=1.0e-30)
        myerror(E,"Degenerate Stage-E2 POD mode");
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(e=1;e<=E->lmesh.NPNO[cases[0]->lev];e++) q[m][e]/=norm;
    return q;
}

static int ala_e2_selected_count(double *values,int n)
{
    double total=0.0,cumulative=0.0,tolerance;
    int i,selected=0;
    for(i=0;i<n;i++) total+=values[i];
    tolerance=n*DBL_EPSILON*max(values[0],1.0);
    for(i=0;i<n && selected<ALA_E2_MAX_MODES;i++) {
        if(values[i]<=tolerance) break;
        cumulative+=values[i]/max(total,1.0e-300);
        selected++;
        if(cumulative>=0.90) break;
    }
    return selected;
}

static void ala_e2_write_coefficients(struct All_variables *E,
    struct ala_e2_observer *out,struct ala_e2_observer **late_cases,
    int *late_indices,int n,double *vectors,double *norms,int selected)
{
    FILE *fp=NULL;
    double **q;
    int mode,c,i;
    if(E->parallel.me==0) {
        fp=ala_e2_open_root(E,out,
            "strict_ala_stage_E2_joint_pod_coefficients.csv","w");
        fprintf(fp,"case,iteration,pod_mode,coefficient,contraction,sign_change\n");
    }
    for(mode=0;mode<selected;mode++) {
        q=ala_e2_build_mode(E,late_cases,late_indices,n,vectors,norms,mode);
        for(c=0;c<2;c++) {
            double previous=0.0;
            struct ala_e2_observer *o=(c==0) ? late_cases[0]
                : late_cases[ALA_E2_SNAPSHOT_COUNT-ALA_E2_LATE_FIRST];
            for(i=0;i<o->count;i++) {
                double coefficient=global_pdot(E,q,o->snapshot[i],o->lev);
                if(E->parallel.me==0)
                    fprintf(fp,"%s,%d,%d,%.17e,%.17e,%d\n",o->case_name,
                        o->iterations[i],mode+1,coefficient,
                        i ? fabs(coefficient)/max(fabs(previous),1.0e-300) : 1.0,
                        i && coefficient*previous<0.0);
                previous=coefficient;
            }
        }
        ala_schur_free_field(E,q);
    }
    if(fp!=NULL) fclose(fp);
}

static void ala_e2_write_angles(struct All_variables *E,
    struct ala_e2_observer *e0,struct ala_e2_observer *e1)
{
    const int n=ALA_E2_LATE_COUNT;
    struct ala_e2_observer *c0[n],*c1[n];
    int index[n],i,j,k,dims[2];
    double C0[n*n],C1[n*n],C01[n*n],v0[n*n],v1[n*n],l0[n],l1[n];
    double norms0[n],norms1[n],cross[n*n],ata[n*n],sv[n],vv[n*n];
    double total0=0.0,total1=0.0,cum0=0.0,cum1=0.0;
    FILE *fp=NULL;
    for(i=0;i<n;i++) { c0[i]=e0; c1[i]=e1; index[i]=i+ALA_E2_LATE_FIRST; }
    ala_e2_covariance(E,c0,index,n,C0);
    ala_e2_covariance(E,c1,index,n,C1);
    for(i=0;i<n;i++) for(j=0;j<n;j++) {
        double local=0.0,global;
        int m,e;
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=1;e<=E->lmesh.NPNO[e0->lev];e++)
                local+=e0->snapshot[index[i]][m][e]
                      *e1->snapshot[index[j]][m][e];
        MPI_Allreduce(&local,&global,1,MPI_DOUBLE,MPI_SUM,E->parallel.world);
        C01[i*n+j]=global;
    }
    ala_e2_pod(C0,n,1,l0,v0,norms0);
    ala_e2_pod(C1,n,1,l1,v1,norms1);
    for(i=0;i<n;i++) { total0+=l0[i]; total1+=l1[i]; }
    dims[0]=dims[1]=n;
    for(i=0;i<n;i++) {
        cum0+=l0[i]/max(total0,1.0e-300);
        cum1+=l1[i]/max(total1,1.0e-300);
        if(dims[0]==n && cum0>=0.80 && cum1>=0.80) dims[0]=i+1;
        if(dims[1]==n && cum0>=0.90 && cum1>=0.90) dims[1]=i+1;
    }
    if(E->parallel.me==0) {
        fp=ala_e2_open_root(E,e1,"strict_ala_stage_E2_subspace_angles.csv","w");
        fprintf(fp,"energy_level,subspace_dimension,angle_index,angle_degrees,cosine,threshold_status\n");
    }
    for(k=0;k<2;k++) {
        int d=dims[k];
        memset(cross,0,sizeof(cross)); memset(ata,0,sizeof(ata));
        for(i=0;i<d;i++) for(j=0;j<d;j++) {
            int a,b;
            for(a=0;a<n;a++) for(b=0;b<n;b++)
                cross[i*d+j]+=v0[a*n+i]
                    *(C01[a*n+b]/max(norms0[a]*norms1[b],1.0e-300))
                    *v1[b*n+j]/sqrt(max(l0[i]*l1[j],1.0e-300));
        }
        for(i=0;i<d;i++) for(j=0;j<d;j++) {
            int a;
            for(a=0;a<d;a++) ata[i*d+j]+=cross[a*d+i]*cross[a*d+j];
        }
        ala_geneo_jacobi_eigensolve(ata,vv,sv,d);
        for(i=0;i<d;i++) sv[i]=sqrt(max(sv[i],0.0));
        for(i=0;i<d-1;i++) for(j=i+1;j<d;j++) if(sv[j]>sv[i]) {
            double t=sv[i]; sv[i]=sv[j]; sv[j]=t;
        }
        if(E->parallel.me==0) for(i=0;i<d;i++) {
            double cosine=min(1.0,max(0.0,sv[i]));
            fprintf(fp,"%d,%d,%d,%.17e,%.17e,PROVISIONAL_DIAGNOSTIC\n",
                k ? 90 : 80,d,i+1,acos(cosine)*180.0/M_PI,cosine);
        }
    }
    if(fp!=NULL) fclose(fp);
}

static void ala_e2_projection_energy(struct All_variables *E,double **q,
    int lev,int factor,double energy[ALA_SCHUR_DEPTH_BANDS])
{
    int nx=E->mesh.ELX[lev],ny=E->mesh.ELY[lev];
    int bx=(nx+factor-1)/factor,by=(ny+factor-1)/factor;
    int size=E->sphere.caps*bx*by*ALA_SCHUR_DEPTH_BANDS;
    double *sum=(double *)calloc(size,sizeof(double));
    double *count=(double *)calloc(size,sizeof(double));
    double *gsum=(double *)calloc(size,sizeof(double));
    double *gcount=(double *)calloc(size,sizeof(double));
    int m,ex,ey,ez,e,b,index;
    if(sum==NULL || count==NULL || gsum==NULL || gcount==NULL)
        myerror(E,"Unable to allocate Stage-E2 multiresolution workspace");
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(ey=1;ey<=E->lmesh.ELY[lev];ey++)
            for(ex=1;ex<=E->lmesh.ELX[lev];ex++)
                for(ez=1;ez<=E->lmesh.ELZ[lev];ez++) {
                    int gx=E->lmesh.EXS[lev]+ex-1;
                    int gy=E->lmesh.EYS[lev]+ey-1;
                    e=ez+(ex-1)*E->lmesh.ELZ[lev]
                        +(ey-1)*E->lmesh.ELZ[lev]*E->lmesh.ELX[lev];
                    b=ala_schur_depth_band((1.0-E->ECO[lev][m][e].centre[3])
                                           *E->data.radius_km);
                    index=((((E->sphere.capid[m]-1)*by+gy/factor)*bx
                           +gx/factor)*ALA_SCHUR_DEPTH_BANDS+b);
                    sum[index]+=q[m][e]; count[index]+=1.0;
                }
    MPI_Allreduce(sum,gsum,size,MPI_DOUBLE,MPI_SUM,E->parallel.world);
    MPI_Allreduce(count,gcount,size,MPI_DOUBLE,MPI_SUM,E->parallel.world);
    for(b=0;b<ALA_SCHUR_DEPTH_BANDS;b++) energy[b]=0.0;
    for(index=0;index<size;index++) {
        b=index%ALA_SCHUR_DEPTH_BANDS;
        if(gcount[index]>0.0) energy[b]+=gsum[index]*gsum[index]/gcount[index];
    }
    free(sum); free(count); free(gsum); free(gcount);
}

static void ala_e2_spatial_row(struct All_variables *E,
    struct ala_e2_observer *out,FILE *scale,FILE *mpi,const char *case_name,
    const char *kind,int id,double **q)
{
    int m,ex,ey,ez,e,b,d;
    int lev=out->lev,elx=E->lmesh.ELX[lev],ely=E->lmesh.ELY[lev];
    int elz=E->lmesh.ELZ[lev];
    double local_total_depth[ALA_SCHUR_DEPTH_BANDS],global_total_depth[ALA_SCHUR_DEPTH_BANDS];
    double p128[ALA_SCHUR_DEPTH_BANDS],p16[ALA_SCHUR_DEPTH_BANDS];
    double p4[ALA_SCHUR_DEPTH_BANDS],global_depth[ALA_SCHUR_DEPTH_BANDS][4];
    double local_mpi[4][2],global_mpi[4][2],total=0.0;
    for(b=0;b<ALA_SCHUR_DEPTH_BANDS;b++) local_total_depth[b]=0.0;
    for(d=0;d<4;d++) { local_mpi[d][0]=0.0; local_mpi[d][1]=0.0; }
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(ey=1;ey<=ely;ey++) for(ex=1;ex<=elx;ex++) for(ez=1;ez<=elz;ez++) {
            double energy;
            e=ez+(ex-1)*elz+(ey-1)*elz*elx;
            energy=q[m][e]*q[m][e];
            b=ala_schur_depth_band((1.0-E->ECO[lev][m][e].centre[3])
                                   *E->data.radius_km);
            local_total_depth[b]+=energy;
            d=min(min(ex-1,elx-ex),min(ey-1,ely-ey));
            d=min(d,3);
            local_mpi[d][0]+=1.0;
            local_mpi[d][1]+=energy;
            total+=energy;
        }
    MPI_Allreduce(local_total_depth,global_total_depth,ALA_SCHUR_DEPTH_BANDS,
                  MPI_DOUBLE,MPI_SUM,E->parallel.world);
    ala_e2_projection_energy(E,q,lev,128,p128);
    ala_e2_projection_energy(E,q,lev,16,p16);
    ala_e2_projection_energy(E,q,lev,4,p4);
    for(b=0;b<ALA_SCHUR_DEPTH_BANDS;b++) {
        global_depth[b][0]=max(p128[b],0.0);
        global_depth[b][1]=max(p16[b]-p128[b],0.0);
        global_depth[b][2]=max(p4[b]-p16[b],0.0);
        global_depth[b][3]=max(global_total_depth[b]-p4[b],0.0);
    }
    MPI_Allreduce(local_mpi,global_mpi,8,MPI_DOUBLE,MPI_SUM,E->parallel.world);
    MPI_Allreduce(MPI_IN_PLACE,&total,1,MPI_DOUBLE,MPI_SUM,E->parallel.world);
    if(E->parallel.me==0) {
        static const char *depths[5]={"0_200","200_410","410_660","660_1000","1000_cmb"};
        static const char *bands[4]={"GLOBAL_VERY_LONG","LONG_INTERMEDIATE","PATCH_LOCAL","HIGH_FREQUENCY"};
        for(b=0;b<5;b++) for(d=0;d<4;d++)
            fprintf(scale,"%s,%s,%d,%s,%s,%.17e,%.17e,geometry_native_nested_v1\n",
                case_name,kind,id,depths[b],bands[d],global_depth[b][d],
                global_depth[b][d]/max(total,1.0e-300));
        for(d=0;d<4;d++)
            fprintf(mpi,"%s,%s,%d,%s,%.0f,%.17e,%.17e\n",case_name,kind,id,
                d==3 ? "3+" : (d==2 ? "2" : (d==1 ? "1" : "0")),
                global_mpi[d][0],global_mpi[d][1],
                global_mpi[d][1]/max(global_mpi[d][0],1.0));
    }
}

static void ala_e2_true_mode_audit(struct All_variables *E,
    struct ala_e2_observer *out,struct ala_e2_observer **cases,int *indices,
    int n,double *values,double *vectors,double *norms,int selected,
    struct ala_pressure_preconditioner_cache *cache)
{
    double **q,**z=ala_schur_alloc_pressure(E,out->lev);
    double **action=ala_schur_alloc_pressure(E,out->lev);
    double **work=ala_schur_alloc_pressure(E,out->lev);
    double **aux=ala_schur_alloc_pressure(E,out->lev);
    double **rhs=ala_schur_alloc_velocity(E,out->lev);
    double **u=ala_schur_alloc_velocity(E,out->lev);
    double **vwork=ala_schur_alloc_velocity(E,out->lev);
    FILE *fp=NULL;
    int mode,map,m,e;
    long long saved_inner_calls=E->control.ala_stage_abc_inner_call_count;
    long long saved_inner_cycles=E->control.ala_stage_abc_inner_cycle_count;
    long long saved_schur=E->control.ala_stage_abc_schur_action_count;
    long long saved_preconditioners=E->control.ala_stage_abc_preconditioner_count;
    long long saved_k=E->control.ala_stage_abc_k_application_count;
    long long saved_e_k=E->control.ala_stage_e_k_operator_applications;
    long long saved_sweeps=E->control.ala_stage_e_velocity_smoother_sweeps;
    double saved_inner_seconds=E->control.ala_stage_abc_inner_seconds;
    const char *maps[3]={"BPI","CONFIGURED","PURE_SCHWARZ"};
    double total=0.0,cumulative=0.0;
    for(mode=0;mode<n;mode++) total+=values[mode];
    if(E->parallel.me==0) {
        fp=ala_e2_open_root(E,out,"strict_ala_stage_E2_true_mode_preconditioner.csv","w");
        fprintf(fp,"mode,map,mode_energy_fraction,cumulative_energy,E_P,cosine,alpha_opt,Pq_norm,SPq_norm,qTPq,qTSPq,positive_qTPq,positive_qTSPq,tight_solve_achieved,map_semantics\n");
    }
    for(mode=0;mode<selected;mode++) {
        q=ala_e2_build_mode(E,cases,indices,n,vectors,norms,mode);
        cumulative+=values[mode]/max(total,1.0e-300);
        for(map=0;map<3;map++) {
            double achieved,qz,qa,zz,aa,error2,cosine,alpha;
            ala_schur_preconditioner_action(E,q,z,work,aux,out->lev,
                                            map==0 ? 0 : (map==1 ? 3 : 1),cache);
            ala_schur_apply_true_action(E,z,action,rhs,u,vwork,out->lev,&achieved);
            qz=global_pdot(E,q,z,out->lev);
            qa=global_pdot(E,q,action,out->lev);
            zz=global_pdot(E,z,z,out->lev);
            aa=global_pdot(E,action,action,out->lev);
            error2=1.0+aa-2.0*qa;
            cosine=qa/max(sqrt(max(aa,0.0)),1.0e-300);
            alpha=qa/max(aa,1.0e-300);
            if(E->parallel.me==0)
                fprintf(fp,"%d,%s,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%d,%d,%.17e,%s\n",
                    mode+1,maps[map],values[mode]/max(total,1.0e-300),cumulative,
                    sqrt(max(error2,0.0)),cosine,alpha,sqrt(max(zz,0.0)),
                    sqrt(max(aa,0.0)),qz,qa,qz>0.0,qa>0.0,achieved,
                    map==2 ? "configured_minus_BPI_exact_component" : "complete_map");
            for(m=1;m<=E->sphere.caps_per_proc;m++)
                for(e=1;e<=E->lmesh.NPNO[out->lev];e++) {
                    z[m][e]=0.0; action[m][e]=0.0; work[m][e]=0.0; aux[m][e]=0.0;
                }
        }
        ala_schur_free_field(E,q);
    }
    if(fp!=NULL) fclose(fp);
    E->control.ala_stage_abc_inner_call_count=saved_inner_calls;
    E->control.ala_stage_abc_inner_cycle_count=saved_inner_cycles;
    E->control.ala_stage_abc_schur_action_count=saved_schur;
    E->control.ala_stage_abc_preconditioner_count=saved_preconditioners;
    E->control.ala_stage_abc_k_application_count=saved_k;
    E->control.ala_stage_e_k_operator_applications=saved_e_k;
    E->control.ala_stage_e_velocity_smoother_sweeps=saved_sweeps;
    E->control.ala_stage_abc_inner_seconds=saved_inner_seconds;
    ala_schur_free_field(E,z); ala_schur_free_field(E,action);
    ala_schur_free_field(E,work); ala_schur_free_field(E,aux);
    ala_schur_free_field(E,rhs); ala_schur_free_field(E,u);
    ala_schur_free_field(E,vwork);
}

static void strict_ala_stage_e2_finalize(struct All_variables *E,void *state,
    struct ala_pressure_preconditioner_cache *cache)
{
    struct ala_e2_observer *o=(struct ala_e2_observer *)state,*e0;
    struct ala_e2_observer *late_cases[2*ALA_E2_LATE_COUNT];
    int late_indices[2*ALA_E2_LATE_COUNT];
    double values[2*ALA_E2_LATE_COUNT];
    double vectors[4*ALA_E2_LATE_COUNT*ALA_E2_LATE_COUNT];
    double norms[2*ALA_E2_LATE_COUNT];
    double e0_values[ALA_E2_LATE_COUNT],e1_values[ALA_E2_LATE_COUNT];
    double e0_vectors[ALA_E2_LATE_COUNT*ALA_E2_LATE_COUNT];
    double e1_vectors[ALA_E2_LATE_COUNT*ALA_E2_LATE_COUNT];
    double e0_norms[ALA_E2_LATE_COUNT],e1_norms[ALA_E2_LATE_COUNT];
    unsigned long long checksum;
    FILE *scale=NULL,*mpi=NULL;
    int i,n,selected=0;
    double start;
    if(o==NULL || !o->enabled) return;
    start=CPU_time0();
    checksum=ala_e2_write_snapshots(E,o);
    if(E->parallel.me==0 && strcmp(o->case_name,"E0")==0) {
        FILE *fp=ala_e2_open_root(E,o,"strict_ala_stage_E2_snapshots_manifest.csv","w");
        fprintf(fp,"case,iteration,residual_norm,global_checksum,storage_identifier,ranks,local_npno\n");
        fclose(fp);
        fp=ala_e2_open_root(E,o,"strict_ala_stage_E2_pod.csv","w");
        fprintf(fp,"case,basis_type,mode,singular_value,energy_fraction,cumulative_fraction,N50,N80,N90,N95,weighting\n");
        fclose(fp);
    }
    MPI_Barrier(E->parallel.world);
    ala_e2_write_manifest(E,o,checksum);
    if(strcmp(o->case_name,"E0")==0) {
        struct ala_e2_observer *cases[ALA_E2_SNAPSHOT_COUNT];
        int indices[ALA_E2_SNAPSHOT_COUNT];
        for(i=0;i<o->count;i++) { cases[i]=o; indices[i]=i; }
        ala_e2_write_pod_set(E,o,cases,indices,o->count,"E0","full_physical",0,NULL,NULL,NULL);
        for(i=ALA_E2_LATE_FIRST;i<o->count;i++) {
            cases[i-ALA_E2_LATE_FIRST]=o; indices[i-ALA_E2_LATE_FIRST]=i;
        }
        ala_e2_write_pod_set(E,o,cases,indices,o->count-ALA_E2_LATE_FIRST,
                             "E0","late_physical",0,NULL,NULL,NULL);
        ala_e2_write_pod_set(E,o,cases,indices,o->count-ALA_E2_LATE_FIRST,
                             "E0","late_direction",1,NULL,NULL,NULL);
        if(E->parallel.me==0) fprintf(E->fp,"STRICT_ALA_STAGE_E2_ARCHIVED case=E0 snapshots=19\n");
        ala_e2_free_observer(E,o);
        return;
    }
    e0=ala_e2_read_case(E,o,"E0");
    ala_e2_write_gram(E,e0,o);
    {
        struct ala_e2_observer *cases[ALA_E2_SNAPSHOT_COUNT];
        int indices[ALA_E2_SNAPSHOT_COUNT];
        for(i=0;i<o->count;i++) { cases[i]=o; indices[i]=i; }
        ala_e2_write_pod_set(E,o,cases,indices,o->count,"E1","full_physical",0,NULL,NULL,NULL);
        for(i=ALA_E2_LATE_FIRST;i<o->count;i++) {
            cases[i-ALA_E2_LATE_FIRST]=o; indices[i-ALA_E2_LATE_FIRST]=i;
        }
        ala_e2_write_pod_set(E,o,cases,indices,o->count-ALA_E2_LATE_FIRST,
                             "E1","late_physical",0,NULL,NULL,NULL);
        ala_e2_write_pod_set(E,o,cases,indices,o->count-ALA_E2_LATE_FIRST,
                             "E1","late_direction",1,NULL,NULL,NULL);
    }
    n=2*ALA_E2_LATE_COUNT;
    for(i=0;i<ALA_E2_LATE_COUNT;i++) {
        late_cases[i]=e0; late_indices[i]=i+ALA_E2_LATE_FIRST;
        late_cases[i+ALA_E2_LATE_COUNT]=o;
        late_indices[i+ALA_E2_LATE_COUNT]=i+ALA_E2_LATE_FIRST;
    }
    ala_e2_write_pod_set(E,o,late_cases,late_indices,n,"JOINT",
        "late_physical",0,NULL,NULL,NULL);
    ala_e2_write_pod_set(E,o,late_cases,late_indices,n,"JOINT",
        "late_direction_case_balanced",1,values,vectors,norms);
    {
        struct ala_e2_observer *separate[ALA_E2_LATE_COUNT];
        int separate_indices[ALA_E2_LATE_COUNT];
        double covariance[ALA_E2_LATE_COUNT*ALA_E2_LATE_COUNT];
        for(i=0;i<ALA_E2_LATE_COUNT;i++) {
            separate[i]=e0; separate_indices[i]=i+ALA_E2_LATE_FIRST;
        }
        ala_e2_covariance(E,separate,separate_indices,ALA_E2_LATE_COUNT,covariance);
        ala_e2_pod(covariance,ALA_E2_LATE_COUNT,1,e0_values,e0_vectors,e0_norms);
        for(i=0;i<ALA_E2_LATE_COUNT;i++) separate[i]=o;
        ala_e2_covariance(E,separate,separate_indices,ALA_E2_LATE_COUNT,covariance);
        ala_e2_pod(covariance,ALA_E2_LATE_COUNT,1,e1_values,e1_vectors,e1_norms);
    }
    selected=ala_e2_selected_count(values,n);
    ala_e2_write_coefficients(E,o,late_cases,late_indices,n,vectors,norms,selected);
    ala_e2_write_angles(E,e0,o);
    if(E->parallel.me==0) {
        scale=ala_e2_open_root(E,o,"strict_ala_stage_E2_depth_scale_spectrum.csv","w");
        fprintf(scale,"case,object_type,object_id,depth_band,horizontal_scale_band,energy,normalized_fraction,method\n");
        mpi=ala_e2_open_root(E,o,"strict_ala_stage_E2_mpi_distance.csv","w");
        fprintf(mpi,"case,object_type,object_id,mpi_distance_band,element_count,energy,energy_per_element\n");
    }
    for(i=0;i<ALA_E2_SNAPSHOT_COUNT;i++) {
        ala_e2_spatial_row(E,o,scale,mpi,"E0","snapshot",e0->iterations[i],
                           e0->snapshot[i]);
        ala_e2_spatial_row(E,o,scale,mpi,"E1","snapshot",o->iterations[i],
                           o->snapshot[i]);
    }
    for(i=0;i<selected;i++) {
        double **q=ala_e2_build_mode(E,late_cases,late_indices,n,vectors,norms,i);
        ala_e2_spatial_row(E,o,scale,mpi,"JOINT","pod_mode",i+1,q);
        ala_schur_free_field(E,q);
    }
    {
    int e0_selected=ala_e2_selected_count(e0_values,ALA_E2_LATE_COUNT);
    int e1_selected=ala_e2_selected_count(e1_values,ALA_E2_LATE_COUNT);
    int separate_max=max(e0_selected,e1_selected);
    for(i=0;i<separate_max;i++) {
        struct ala_e2_observer *separate[ALA_E2_LATE_COUNT];
        int separate_indices[ALA_E2_LATE_COUNT],j;
        double **q;
        for(j=0;j<ALA_E2_LATE_COUNT;j++) {
            separate[j]=e0; separate_indices[j]=j+ALA_E2_LATE_FIRST;
        }
        if(i<e0_selected) {
            q=ala_e2_build_mode(E,separate,separate_indices,ALA_E2_LATE_COUNT,
                                e0_vectors,e0_norms,i);
            ala_e2_spatial_row(E,o,scale,mpi,"E0","pod_mode",i+1,q);
            ala_schur_free_field(E,q);
        }
        for(j=0;j<ALA_E2_LATE_COUNT;j++) separate[j]=o;
        if(i<e1_selected) {
            q=ala_e2_build_mode(E,separate,separate_indices,ALA_E2_LATE_COUNT,
                                e1_vectors,e1_norms,i);
            ala_e2_spatial_row(E,o,scale,mpi,"E1","pod_mode",i+1,q);
            ala_schur_free_field(E,q);
        }
    }
    }
    if(scale!=NULL) fclose(scale);
    if(mpi!=NULL) fclose(mpi);
    ala_e2_true_mode_audit(E,o,late_cases,late_indices,n,values,vectors,norms,
                           selected,cache);
    if(E->parallel.me==0) {
        fprintf(E->fp,"STRICT_ALA_STAGE_E2_COMPLETE case=E1 snapshots=38 "
                "selected_modes=%d capture_seconds=%e postsolve_seconds=%e\n",
                selected,o->capture_seconds,CPU_time0()-start);
        fflush(E->fp);
    }
    ala_e2_free_observer(E,e0);
    ala_e2_free_observer(E,o);
}
