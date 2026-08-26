/* Stage F1b selector-only principal-submatrix construction.
 *
 * The legacy shallow-patch cache is immutable.  This file is reached only
 * when ala_shallow_patch_local_operator=operator_consistent.  It builds the
 * exact element-sum principal restriction of K_gamma on the velocity-node
 * union of a pressure patch, including direct-face and two-hop corner element
 * records.  The temporary velocity matrix/factor is discarded after the
 * pressure Schur factor has been formed.
 */

#define ALA_F1B_RECORD_META 4
#define ALA_F1B_RECORD_COORD (ALA_F1B_RECORD_META)
#define ALA_F1B_RECORD_B (ALA_F1B_RECORD_COORD+ALA_VANKA_DOF)
#define ALA_F1B_RECORD_FIXED (ALA_F1B_RECORD_B+ALA_VANKA_DOF)
#define ALA_F1B_RECORD_K (ALA_F1B_RECORD_FIXED+ALA_VANKA_DOF)
#define ALA_F1B_RECORD_DOUBLES (ALA_F1B_RECORD_K+ALA_VANKA_DOF*ALA_VANKA_DOF)
#define ALA_F1B_MAX_VELOCITY_DOF 441

struct ala_f1b_element_pool {
    int local_count;
    int remote_count;
    int support_min_ez;
    int support_layers;
    double *local_records;
    double *remote_records;
    double temporary_bytes;
    double resident_bytes;
    double mpi_payload_bytes;
};

struct ala_f1b_patch_stats {
    int count;
    int rhs_solve_count;
    int max_pdim;
    int max_vdim;
    double max_symmetry_defect;
    double min_k_pivot_ratio;
    double max_k_solve_residual;
    double max_temporary_bytes;
    double max_factor_bytes;
    double estimated_dense_work;
};

struct ala_f1b_replay_snapshot {
    int valid;
    int vdim;
    double *matrix;
    double coordinates[3*ALA_F1B_MAX_VELOCITY_DOF];
    unsigned char components[ALA_F1B_MAX_VELOCITY_DOF];
};

static struct ala_f1b_replay_snapshot ala_f1b_replay_snapshot;

static int ala_f1b_same_point(const double *a,const double *b)
{
    double dx=a[0]-b[0],dy=a[1]-b[1],dz=a[2]-b[2];
    return(dx*dx+dy*dy+dz*dz<=1.0e-18);
}


static void ala_f1b_record_center(const double *record,double center[3])
{
    int a,d;
    for(d=0;d<3;d++) {
        center[d]=0.0;
        for(a=0;a<8;a++)
            center[d] += record[ALA_F1B_RECORD_COORD+3*a+d];
        center[d] *= 0.125;
    }
}


static void ala_f1b_fill_element_record(struct All_variables *E,int e,
    int lev,int m,int ex,int ey,int ez,double *record)
{
    int a,b,da,db,i,j,node;
    double matrix[ALA_VANKA_DOF*ALA_VANKA_DOF];
    void get_elt_k();
    void get_ala_aug_k();
    const int dims=E->mesh.nsd;
    const int ends=enodes[dims];

    memset(record,0,ALA_F1B_RECORD_DOUBLES*sizeof(double));
    record[0]=(double)E->parallel.me;
    record[1]=(double)ex;
    record[2]=(double)ey;
    record[3]=(double)ez;
    get_elt_k(E,e,matrix,lev,m,0);
    if(E->control.augmented_Lagr)
        myerror(E,"F1b operator-consistent candidate does not permit legacy augmentation");
    if(E->control.ala_augmented_lagrangian_gamma>0.0)
        get_ala_aug_k(E,e,matrix,lev,m);
    for(a=1;a<=ends;a++) {
        node=E->IEN[lev][m][e].node[a];
        for(da=0;da<dims;da++) {
            i=(a-1)*dims+da;
            record[ALA_F1B_RECORD_COORD+i]=E->X[lev][m][da+1][node];
            record[ALA_F1B_RECORD_B+i]=ALA_COMBINED_PRESSURE_COEFFICIENT(
                E->elt_del[lev][m][e].g[i][0],
                E->elt_c[lev][m][e].c[i][0]);
            record[ALA_F1B_RECORD_FIXED+i]=
                (da==0 && (E->NODE[lev][m][node]&VBX)) ||
                (da==1 && (E->NODE[lev][m][node]&VBY)) ||
                (da==2 && (E->NODE[lev][m][node]&VBZ));
        }
    }
    for(a=1;a<=ends;a++)
        for(da=0;da<dims;da++) {
            i=(a-1)*dims+da;
            for(b=1;b<=ends;b++)
                for(db=0;db<dims;db++) {
                    j=(b-1)*dims+db;
                    record[ALA_F1B_RECORD_K+i*ALA_VANKA_DOF+j]=
                        (record[ALA_F1B_RECORD_FIXED+i]!=0.0 ||
                         record[ALA_F1B_RECORD_FIXED+j]!=0.0)
                        ? 0.0 : matrix[i*ALA_VANKA_DOF+j];
                }
            if(record[ALA_F1B_RECORD_FIXED+i]!=0.0)
                record[ALA_F1B_RECORD_B+i]=0.0;
        }
}


static int ala_f1b_send_record_to(struct All_variables *E,
    const double *record,int target,
    const int *neighbors,int nproc,int overlap,int elx,int ely)
{
    int f,g,direct=0,two_hop=0;
    int ex=(int)(record[1]+0.5),ey=(int)(record[2]+0.5);
    int face_hit[4];
    face_hit[0]=(ex<=overlap+1);
    face_hit[1]=(ex>elx-overlap-1);
    face_hit[2]=(ey<=overlap+1);
    face_hit[3]=(ey>ely-overlap-1);
    for(f=0;f<4;f++)
        if(neighbors[(E->parallel.me)*4+f]==target && face_hit[f])
            direct=1;
    if(direct)
        return(1);
    for(f=0;f<4;f++)
        for(g=0;g<4;g++)
            if(neighbors[neighbors[E->parallel.me*4+f]*4+g]==target)
                two_hop=1;
    if(!two_hop || target==E->parallel.me)
        return(0);
    /* Only corner blocks are needed from a two-hop horizontal neighbour. */
    return((face_hit[0]||face_hit[1]) && (face_hit[2]||face_hit[3]));
}


static void ala_f1b_build_element_pool(struct All_variables *E,int lev,int m,
    int shallow_min_ez,struct ala_f1b_element_pool *pool)
{
    int ex,ey,ez,e,i,t,f,nproc,elx,ely,elz,overlap;
    int *neighbors,*send_records,*recv_records,*send_counts,*recv_counts;
    int *send_displs,*recv_displs,total_send,total_recv,cursor;
    double *send_buffer;
    const double *record;

    memset(pool,0,sizeof(*pool));
    nproc=E->parallel.nproc;
    elx=E->lmesh.ELX[lev];
    ely=E->lmesh.ELY[lev];
    elz=E->lmesh.ELZ[lev];
    overlap=E->control.ala_shallow_patch_mpi_overlap;
    pool->support_min_ez=max(1,shallow_min_ez-1);
    pool->support_layers=elz-pool->support_min_ez+1;
    pool->local_count=elx*ely*pool->support_layers;
    pool->local_records=(double *)calloc(
        (size_t)max(pool->local_count,1)*ALA_F1B_RECORD_DOUBLES,
        sizeof(double));
    if(pool->local_records==NULL)
        myerror(E,"Unable to allocate F1b local element records");
    i=0;
    for(ey=1;ey<=ely;ey++)
        for(ex=1;ex<=elx;ex++)
            for(ez=pool->support_min_ez;ez<=elz;ez++) {
                e=ez+(ex-1)*elz+(ey-1)*elz*elx;
                ala_f1b_fill_element_record(E,e,lev,m,ex,ey,ez,
                    pool->local_records+(size_t)i*ALA_F1B_RECORD_DOUBLES);
                i++;
            }

    neighbors=(int *)calloc((size_t)nproc*4,sizeof(int));
    send_records=(int *)calloc(nproc,sizeof(int));
    recv_records=(int *)calloc(nproc,sizeof(int));
    send_counts=(int *)calloc(nproc,sizeof(int));
    recv_counts=(int *)calloc(nproc,sizeof(int));
    send_displs=(int *)calloc(nproc,sizeof(int));
    recv_displs=(int *)calloc(nproc,sizeof(int));
    if(neighbors==NULL || send_records==NULL || recv_records==NULL ||
       send_counts==NULL || recv_counts==NULL || send_displs==NULL ||
       recv_displs==NULL)
        myerror(E,"Unable to allocate F1b halo metadata");
    {
        int local_neighbors[4];
        for(f=0;f<4;f++)
            local_neighbors[f]=E->parallel.PROCESSOR[lev][m].pass[f+1];
        MPI_Allgather(local_neighbors,4,MPI_INT,neighbors,4,MPI_INT,
                      E->parallel.world);
    }
    for(t=0;t<nproc;t++)
        if(t!=E->parallel.me)
            for(i=0;i<pool->local_count;i++) {
                record=pool->local_records
                    +(size_t)i*ALA_F1B_RECORD_DOUBLES;
                if(ala_f1b_send_record_to(E,record,t,neighbors,nproc,
                                           overlap,elx,ely))
                    send_records[t]++;
            }
    MPI_Alltoall(send_records,1,MPI_INT,recv_records,1,MPI_INT,
                 E->parallel.world);
    total_send=total_recv=0;
    for(t=0;t<nproc;t++) {
        send_displs[t]=total_send;
        recv_displs[t]=total_recv;
        send_counts[t]=send_records[t]*ALA_F1B_RECORD_DOUBLES;
        recv_counts[t]=recv_records[t]*ALA_F1B_RECORD_DOUBLES;
        total_send += send_counts[t];
        total_recv += recv_counts[t];
    }
    send_buffer=(double *)calloc(max(total_send,1),sizeof(double));
    pool->remote_records=(double *)calloc(max(total_recv,1),sizeof(double));
    if(send_buffer==NULL || pool->remote_records==NULL)
        myerror(E,"Unable to allocate F1b halo records");
    for(t=0;t<nproc;t++) {
        cursor=send_displs[t];
        for(i=0;i<pool->local_count;i++) {
            record=pool->local_records
                +(size_t)i*ALA_F1B_RECORD_DOUBLES;
            if(ala_f1b_send_record_to(E,record,t,neighbors,nproc,
                                       overlap,elx,ely)) {
                memcpy(send_buffer+cursor,record,
                       ALA_F1B_RECORD_DOUBLES*sizeof(double));
                cursor += ALA_F1B_RECORD_DOUBLES;
            }
        }
    }
    MPI_Alltoallv(send_buffer,send_counts,send_displs,MPI_DOUBLE,
                  pool->remote_records,recv_counts,recv_displs,MPI_DOUBLE,
                  E->parallel.world);
    pool->remote_count=total_recv/ALA_F1B_RECORD_DOUBLES;
    pool->temporary_bytes=
        (double)(pool->local_count+pool->remote_count)
        *ALA_F1B_RECORD_DOUBLES*sizeof(double)
        +(double)total_send*sizeof(double);
    pool->resident_bytes=(double)(pool->local_count+pool->remote_count)
        *ALA_F1B_RECORD_DOUBLES*sizeof(double);
    pool->mpi_payload_bytes=(double)total_send*sizeof(double);
    free(send_buffer);
    free(neighbors);
    free(send_records);
    free(recv_records);
    free(send_counts);
    free(recv_counts);
    free(send_displs);
    free(recv_displs);
}


static const double *ala_f1b_local_record(struct All_variables *E,int lev,
    const struct ala_f1b_element_pool *pool,int e)
{
    int ex,ey,ez,index,elx,elz;
    elx=E->lmesh.ELX[lev];
    elz=E->lmesh.ELZ[lev];
    ez=(e-1)%elz+1;
    ex=((e-1)/elz)%elx+1;
    ey=(e-1)/(elz*elx)+1;
    if(ez<pool->support_min_ez)
        return(NULL);
    index=((ey-1)*elx+(ex-1))*pool->support_layers
        +(ez-pool->support_min_ez);
    return(pool->local_records+(size_t)index*ALA_F1B_RECORD_DOUBLES);
}


static const double *ala_f1b_find_remote_record(
    const struct ala_f1b_element_pool *pool,const double center[3])
{
    int i;
    double candidate[3];
    const double *record;
    for(i=0;i<pool->remote_count;i++) {
        record=pool->remote_records+(size_t)i*ALA_F1B_RECORD_DOUBLES;
        ala_f1b_record_center(record,candidate);
        if(ala_f1b_same_point(center,candidate))
            return(record);
    }
    return(NULL);
}


static int ala_f1b_factor(double *matrix,int n,double *min_ratio)
{
    int i,j,k;
    double maxdiag=0.0,sum,pivot;
    for(i=0;i<n;i++)
        maxdiag=max(maxdiag,matrix[i*n+i]);
    if(!isfinite(maxdiag) || maxdiag<=0.0)
        return(0);
    for(i=0;i<n;i++)
        for(j=0;j<=i;j++) {
            sum=matrix[i*n+j];
            for(k=0;k<j;k++)
                sum -= matrix[i*n+k]*matrix[j*n+k];
            if(i==j) {
                pivot=sum;
                if(!isfinite(pivot) || pivot<=1.0e-14*maxdiag)
                    return(0);
                *min_ratio=min(*min_ratio,pivot/maxdiag);
                matrix[i*n+j]=sqrt(pivot);
            }
            else
                matrix[i*n+j]=sum/matrix[j*n+j];
            matrix[j*n+i]=0.0;
        }
    return(1);
}


static void ala_f1b_solve(const double *factor,const double *rhs,double *x,
    int n)
{
    int i,j;
    double sum;
    for(i=0;i<n;i++) {
        sum=rhs[i];
        for(j=0;j<i;j++)
            sum -= factor[i*n+j]*x[j];
        x[i]=sum/factor[i*n+i];
    }
    for(i=n-1;i>=0;i--) {
        sum=x[i];
        for(j=i+1;j<n;j++)
            sum -= factor[j*n+i]*x[j];
        x[i]=sum/factor[i*n+i];
    }
}


static int ala_f1b_record_key_seen(const double *record,int used[][4],int n)
{
    int i,d;
    for(i=0;i<n;i++) {
        for(d=0;d<4;d++)
            if(used[i][d]!=(int)(record[d]+0.5))
                break;
        if(d==4)
            return(1);
    }
    return(0);
}


static int ala_f1b_velocity_index(const double *coord,int component,
    const double *vcoord,const unsigned char *vcomp,int n)
{
    int i;
    for(i=0;i<n;i++)
        if(vcomp[i]==component &&
           ala_f1b_same_point(coord,vcoord+3*i))
            return(i);
    return(-1);
}


static int ala_f1b_build_operator_consistent_patch(
    struct All_variables *E,const struct ala_f1b_element_pool *pool,
    const double **pressure_records,int np,
    double matrix[ALA_PATCH_MAX_ELEMENTS][ALA_PATCH_MAX_ELEMENTS],
    struct ala_f1b_patch_stats *stats)
{
    int p,a,d,i,j,k,r,nv,nused,map[ALA_VANKA_DOF];
    int used[1024][4];
    unsigned char vcomp[ALA_F1B_MAX_VELOCITY_DOF];
    double vcoord[3*ALA_F1B_MAX_VELOCITY_DOF];
    double *K,*factor,*B,*rhs,*solution;
    double bnorm2,rnorm2,sum,numerator,denominator,symmetry,min_pivot;
    const double *record,*coord;

    nv=0;
    for(p=0;p<np;p++) {
        record=pressure_records[p];
        if(record==NULL)
            return(0);
        for(a=0;a<8;a++)
            for(d=0;d<3;d++) {
                i=3*a+d;
                if(record[ALA_F1B_RECORD_FIXED+i]!=0.0)
                    continue;
                coord=record+ALA_F1B_RECORD_COORD+3*a;
                if(ala_f1b_velocity_index(coord,d,vcoord,vcomp,nv)<0) {
                    if(nv>=ALA_F1B_MAX_VELOCITY_DOF)
                        return(0);
                    memcpy(vcoord+3*nv,coord,3*sizeof(double));
                    vcomp[nv]=(unsigned char)d;
                    nv++;
                }
            }
    }
    if(nv<=0 || np<=0)
        return(0);
    K=(double *)calloc((size_t)nv*nv,sizeof(double));
    factor=(double *)calloc((size_t)nv*nv,sizeof(double));
    B=(double *)calloc((size_t)np*nv,sizeof(double));
    rhs=(double *)calloc(nv,sizeof(double));
    solution=(double *)calloc(nv,sizeof(double));
    if(K==NULL || factor==NULL || B==NULL || rhs==NULL || solution==NULL)
        myerror(E,"Unable to allocate F1b dense patch workspace");
    for(p=0;p<np;p++) {
        record=pressure_records[p];
        for(a=0;a<8;a++)
            for(d=0;d<3;d++) {
                i=3*a+d;
                k=ala_f1b_velocity_index(
                    record+ALA_F1B_RECORD_COORD+3*a,d,vcoord,vcomp,nv);
                if(k>=0)
                    B[p*nv+k] += record[ALA_F1B_RECORD_B+i];
            }
    }

    nused=0;
    for(r=0;r<pool->local_count+pool->remote_count;r++) {
        record=(r<pool->local_count)
            ? pool->local_records+(size_t)r*ALA_F1B_RECORD_DOUBLES
            : pool->remote_records
                +(size_t)(r-pool->local_count)*ALA_F1B_RECORD_DOUBLES;
        for(a=0;a<8;a++)
            for(d=0;d<3;d++)
                map[3*a+d]=ala_f1b_velocity_index(
                    record+ALA_F1B_RECORD_COORD+3*a,d,vcoord,vcomp,nv);
        for(i=0;i<ALA_VANKA_DOF;i++)
            if(map[i]>=0)
                break;
        if(i==ALA_VANKA_DOF)
            continue;
        if(ala_f1b_record_key_seen(record,used,nused))
            continue;
        if(nused>=1024)
            myerror(E,"F1b patch incident-element bound exceeded");
        for(d=0;d<4;d++)
            used[nused][d]=(int)(record[d]+0.5);
        nused++;
        for(i=0;i<ALA_VANKA_DOF;i++)
            if(map[i]>=0)
                for(j=0;j<ALA_VANKA_DOF;j++)
                    if(map[j]>=0)
                        K[map[i]*nv+map[j]] +=
                            record[ALA_F1B_RECORD_K+i*ALA_VANKA_DOF+j];
    }

    numerator=denominator=0.0;
    for(i=0;i<nv;i++)
        for(j=0;j<nv;j++) {
            numerator=max(numerator,fabs(K[i*nv+j]-K[j*nv+i]));
            denominator=max(denominator,
                            max(fabs(K[i*nv+j]),fabs(K[j*nv+i])));
        }
    symmetry=numerator/max(denominator,1.0e-300);
    memcpy(factor,K,(size_t)nv*nv*sizeof(double));
    min_pivot=1.0e300;
    if(!ala_f1b_factor(factor,nv,&min_pivot)) {
        free(K); free(factor); free(B); free(rhs); free(solution);
        return(0);
    }
    stats->max_k_solve_residual=max(stats->max_k_solve_residual,0.0);
    for(p=0;p<np;p++) {
        for(i=0;i<nv;i++) {
            rhs[i]=B[p*nv+i];
            solution[i]=0.0;
        }
        ala_f1b_solve(factor,rhs,solution,nv);
        bnorm2=rnorm2=0.0;
        for(i=0;i<nv;i++) {
            sum=0.0;
            for(j=0;j<nv;j++)
                sum += K[i*nv+j]*solution[j];
            rnorm2 += (rhs[i]-sum)*(rhs[i]-sum);
            bnorm2 += rhs[i]*rhs[i];
        }
        stats->max_k_solve_residual=max(stats->max_k_solve_residual,
            sqrt(rnorm2/max(bnorm2,1.0e-300)));
        for(j=0;j<=p;j++) {
            sum=0.0;
            for(i=0;i<nv;i++)
                sum += B[j*nv+i]*solution[i];
            matrix[p][j]=matrix[j][p]=sum;
        }
    }
    stats->count++;
    stats->rhs_solve_count += np;
    stats->estimated_dense_work += (double)nv*nv*nv/3.0
        +2.0*(double)np*nv*nv;
    stats->max_pdim=max(stats->max_pdim,np);
    stats->max_vdim=max(stats->max_vdim,nv);
    stats->max_symmetry_defect=max(stats->max_symmetry_defect,symmetry);
    stats->min_k_pivot_ratio=min(stats->min_k_pivot_ratio,min_pivot);
    stats->max_temporary_bytes=max(stats->max_temporary_bytes,
        max(pool->temporary_bytes,pool->resident_bytes+
            (2.0*nv*nv+(double)np*nv+2.0*nv)*sizeof(double)));
    stats->max_factor_bytes=max(stats->max_factor_bytes,
        (double)np*np*sizeof(double));
    if(!ala_f1b_replay_snapshot.valid) {
        ala_f1b_replay_snapshot.matrix=(double *)calloc(
            (size_t)nv*nv,sizeof(double));
        if(ala_f1b_replay_snapshot.matrix==NULL)
            myerror(E,"Unable to allocate F1b replay snapshot");
        memcpy(ala_f1b_replay_snapshot.matrix,K,
               (size_t)nv*nv*sizeof(double));
        memcpy(ala_f1b_replay_snapshot.coordinates,vcoord,
               (size_t)3*nv*sizeof(double));
        memcpy(ala_f1b_replay_snapshot.components,vcomp,
               (size_t)nv*sizeof(unsigned char));
        ala_f1b_replay_snapshot.vdim=nv;
        ala_f1b_replay_snapshot.valid=1;
    }
    free(K); free(factor); free(B); free(rhs); free(solution);
    return(1);
}


static void ala_f1b_free_element_pool(struct ala_f1b_element_pool *pool)
{
    free(pool->local_records);
    free(pool->remote_records);
    memset(pool,0,sizeof(*pool));
}


static int ala_f1b_find_local_node(struct All_variables *E,int lev,int m,
    const double *coord)
{
    int node;
    double local[3];
    for(node=1;node<=E->lmesh.NNO[lev];node++) {
        local[0]=E->X[lev][m][1][node];
        local[1]=E->X[lev][m][2][node];
        local[2]=E->X[lev][m][3][node];
        if(ala_f1b_same_point(coord,local))
            return(node);
    }
    return(-1);
}


static double ala_f1b_replay_global_k(struct All_variables *E,int lev)
{
    int m,i,j,node,eq,owner,local_owner,nv,neq;
    int components[ALA_F1B_MAX_VELOCITY_DOF];
    double coordinates[3*ALA_F1B_MAX_VELOCITY_DOF];
    double x[ALA_F1B_MAX_VELOCITY_DOF];
    double numerator=0.0,denominator=0.0,sum,actual;
    double **u,**Au;
    void n_assemble_del2_u();

    local_owner=ala_f1b_replay_snapshot.valid ? E->parallel.me
                                               : E->parallel.nproc;
    MPI_Allreduce(&local_owner,&owner,1,MPI_INT,MPI_MIN,E->parallel.world);
    if(owner>=E->parallel.nproc)
        return(INFINITY);
    nv=(E->parallel.me==owner) ? ala_f1b_replay_snapshot.vdim : 0;
    MPI_Bcast(&nv,1,MPI_INT,owner,E->parallel.world);
    if(nv<=0 || nv>ALA_F1B_MAX_VELOCITY_DOF)
        return(INFINITY);
    if(E->parallel.me==owner) {
        memcpy(coordinates,ala_f1b_replay_snapshot.coordinates,
               (size_t)3*nv*sizeof(double));
        for(i=0;i<nv;i++)
            components[i]=(int)ala_f1b_replay_snapshot.components[i];
    }
    MPI_Bcast(coordinates,3*nv,MPI_DOUBLE,owner,E->parallel.world);
    MPI_Bcast(components,nv,MPI_INT,owner,E->parallel.world);
    for(i=0;i<nv;i++)
        x[i]=sin(0.17320508075688773*(i+1))
            +0.5*cos(0.22360679774997897*(i+1));
    u=(double **)calloc(NCS,sizeof(double *));
    Au=(double **)calloc(NCS,sizeof(double *));
    if(u==NULL || Au==NULL)
        myerror(E,"Unable to allocate F1b replay vectors");
    neq=E->lmesh.NEQ[lev];
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        u[m]=(double *)calloc(neq+1,sizeof(double));
        Au[m]=(double *)calloc(neq+1,sizeof(double));
        if(u[m]==NULL || Au[m]==NULL)
            myerror(E,"Unable to allocate F1b replay cap vectors");
        for(i=0;i<nv;i++) {
            node=ala_f1b_find_local_node(E,lev,m,coordinates+3*i);
            if(node>0) {
                eq=E->ID[lev][m][node].doff[components[i]+1];
                u[m][eq]=x[i];
            }
        }
    }
    n_assemble_del2_u(E,u,Au,lev,1);
    if(E->parallel.me==owner) {
        for(i=0;i<nv;i++) {
            sum=0.0;
            for(j=0;j<nv;j++)
                sum += ala_f1b_replay_snapshot.matrix[i*nv+j]*x[j];
            actual=0.0;
            for(m=1;m<=E->sphere.caps_per_proc;m++) {
                node=ala_f1b_find_local_node(E,lev,m,coordinates+3*i);
                if(node>0) {
                    eq=E->ID[lev][m][node].doff[components[i]+1];
                    actual=Au[m][eq];
                    break;
                }
            }
            numerator += (sum-actual)*(sum-actual);
            denominator += actual*actual;
        }
        numerator=sqrt(numerator/max(denominator,1.0e-300));
    }
    MPI_Bcast(&numerator,1,MPI_DOUBLE,owner,E->parallel.world);
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        free(u[m]); free(Au[m]);
    }
    free(u); free(Au);
    free(ala_f1b_replay_snapshot.matrix);
    memset(&ala_f1b_replay_snapshot,0,sizeof(ala_f1b_replay_snapshot));
    return(numerator);
}


static void ala_f1b_write_candidate_evidence(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache,double legacy_cache_bytes,
    double build_seconds)
{
    int local_counts[4],global_counts[4];
    int local_dims[4],global_dims[4];
    double local_max[6],global_max[6],local_min,global_min;
    double projected_ratio;
    const char *directory;
    char path[2048];
    FILE *fp;

    if(!cache->f1b_candidate_active)
        return;
    local_counts[0]=cache->f1b_local_patch_count;
    local_counts[1]=cache->f1b_interface_patch_count;
    local_counts[2]=cache->f1b_local_rhs_solve_count;
    local_counts[3]=cache->f1b_interface_rhs_solve_count;
    MPI_Allreduce(local_counts,global_counts,4,MPI_INT,MPI_SUM,
                  E->parallel.world);
    local_dims[0]=cache->f1b_local_max_pdim;
    local_dims[1]=cache->f1b_local_max_vdim;
    local_dims[2]=cache->f1b_interface_max_pdim;
    local_dims[3]=cache->f1b_interface_max_vdim;
    MPI_Allreduce(local_dims,global_dims,4,MPI_INT,MPI_MAX,E->parallel.world);
    local_max[0]=cache->f1b_max_symmetry_defect;
    local_max[1]=cache->f1b_max_k_solve_residual;
    local_max[2]=cache->f1b_peak_temporary_bytes;
    local_max[3]=cache->f1b_mpi_payload_bytes;
    local_max[4]=cache->f1b_k_replay_relative_defect;
    local_max[5]=cache->f1b_estimated_dense_work;
    MPI_Allreduce(local_max,global_max,5,MPI_DOUBLE,MPI_MAX,
                  E->parallel.world);
    MPI_Allreduce(local_max+5,global_max+5,1,MPI_DOUBLE,MPI_SUM,
                  E->parallel.world);
    local_min=(cache->f1b_local_patch_count+cache->f1b_interface_patch_count)>0
        ? cache->f1b_min_k_pivot_ratio : 1.0e300;
    MPI_Allreduce(&local_min,&global_min,1,MPI_DOUBLE,MPI_MIN,
                  E->parallel.world);
    projected_ratio=(legacy_cache_bytes+global_max[2])
        /max(legacy_cache_bytes,1.0);
    directory=getenv("STRICT_ALA_STAGE_F1B_OUTPUT_DIR");
    if(E->parallel.me!=0 || directory==NULL || directory[0]=='\0')
        return;
    snprintf(path,sizeof(path),"%s/strict_ala_stage_F1b_feasibility.json",
             directory);
    fp=fopen(path,"w");
    if(fp==NULL)
        myerror(E,"Unable to open F1b feasibility JSON");
    fprintf(fp,
        "{\n"
        "  \"schema\": \"strict-ala-stage-F1b-feasibility-v1\",\n"
        "  \"candidate_isolatable\": %s,\n"
        "  \"local_patch_max_pdim\": %d,\n"
        "  \"local_patch_max_vdim\": %d,\n"
        "  \"interface_patch_max_pdim\": %d,\n"
        "  \"interface_patch_max_vdim\": %d,\n"
        "  \"factorization_count\": %d,\n"
        "  \"RHS_solve_count\": %d,\n"
        "  \"estimated_dense_factor_and_solve_work\": %.17e,\n"
        "  \"temporary_peak_bytes_per_rank\": %.17e,\n"
        "  \"retained_cache_bytes_per_rank\": %.17e,\n"
        "  \"projected_peak_ratio_vs_legacy_cache\": %.17e,\n"
        "  \"new_mpi_payload_bytes_max_per_rank\": %.17e,\n"
        "  \"Kgamma_replay_relative_defect\": %.17e,\n"
        "  \"Kgamma_replay_tolerance\": 5.0e-6,\n"
        "  \"startup_seconds_max\": %.17e,\n"
        "  \"memory_gate_pass\": %s,\n"
        "  \"communication\": \"direct-face plus two-hop corner element records\"\n"
        "}\n",
        (global_counts[0]>0 && global_counts[1]>0 &&
         global_max[4]<=5.0e-6) ? "true" : "false",
        global_dims[0],global_dims[1],global_dims[2],global_dims[3],
        global_counts[0]+global_counts[1],global_counts[2]+global_counts[3],
        global_max[5],global_max[2],legacy_cache_bytes,
        projected_ratio,global_max[3],global_max[4],build_seconds,
        projected_ratio<=1.25 ? "true" : "false");
    fclose(fp);
    snprintf(path,sizeof(path),
             "%s/strict_ala_stage_F1b_candidate_assembly.csv",directory);
    fp=fopen(path,"w");
    if(fp==NULL)
        myerror(E,"Unable to open F1b candidate assembly CSV");
    fprintf(fp,"row_type,patch_category,rank,patch_id,pressure_dim,"
        "velocity_dim,symmetry_defect,regularization,"
        "factorization_status,pivot_metric,local_solve_residual,fallback,"
        "candidate_matrix_bytes,candidate_factor_bytes,Kgamma_replay_relative_defect\n");
    fprintf(fp,"AGGREGATE,LOCAL,-1,-1,%d,%d,%.17e,%.17e,PASS,%.17e,"
        "%.17e,0,%.17e,%.17e,%.17e\n",
        global_dims[0],global_dims[1],global_max[0],
        E->control.ala_shallow_patch_regularization,global_min,global_max[1],
        (double)global_dims[0]*global_dims[0]*sizeof(double),
        (double)global_dims[0]*global_dims[0]*sizeof(double),global_max[4]);
    fprintf(fp,"AGGREGATE,INTERFACE,-1,-1,%d,%d,%.17e,%.17e,PASS,%.17e,"
        "%.17e,0,%.17e,%.17e,%.17e\n",
        global_dims[2],global_dims[3],global_max[0],
        E->control.ala_shallow_patch_regularization,global_min,global_max[1],
        (double)global_dims[2]*global_dims[2]*sizeof(double),
        (double)global_dims[2]*global_dims[2]*sizeof(double),global_max[4]);
    fclose(fp);
}


static void ala_f1b_write_rank_memory(struct All_variables *E)
{
    struct rusage usage;
    double local_kib,maximum_kib,sum_kib;
    const char *directory,*case_name;
    char path[2048];
    FILE *fp;
    if(getrusage(RUSAGE_SELF,&usage)!=0)
        myerror(E,"F1b getrusage failed");
#if defined(__APPLE__)
    local_kib=(double)usage.ru_maxrss/1024.0;
#else
    local_kib=(double)usage.ru_maxrss;
#endif
    MPI_Allreduce(&local_kib,&maximum_kib,1,MPI_DOUBLE,MPI_MAX,
                  E->parallel.world);
    MPI_Allreduce(&local_kib,&sum_kib,1,MPI_DOUBLE,MPI_SUM,
                  E->parallel.world);
    directory=getenv("STRICT_ALA_STAGE_F1B_OUTPUT_DIR");
    case_name=getenv("STRICT_ALA_CASE");
    if(E->parallel.me!=0 || directory==NULL || directory[0]=='\0')
        return;
    snprintf(path,sizeof(path),"%s/strict_ala_stage_F1b_rank_memory.csv",
             directory);
    fp=fopen(path,"w");
    if(fp==NULL)
        myerror(E,"Unable to open F1b rank memory CSV");
    fprintf(fp,"case,rank_rss_max_kib,rank_rss_sum_kib,ranks,source\n");
    fprintf(fp,"%s,%.17e,%.17e,%d,getrusage_RUSAGE_SELF_MPI_reduce\n",
            case_name ? case_name : "UNKNOWN",maximum_kib,sum_kib,
            E->parallel.nproc);
    fclose(fp);
}
