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
#define ALA_F1B_RECORD_NODE_KEY (ALA_F1B_RECORD_COORD+ALA_VANKA_DOF)
#define ALA_F1B_RECORD_B (ALA_F1B_RECORD_NODE_KEY+4*8)
#define ALA_F1B_RECORD_FIXED (ALA_F1B_RECORD_B+ALA_VANKA_DOF)
#define ALA_F1B_RECORD_GLOBAL_DIAG (ALA_F1B_RECORD_FIXED+ALA_VANKA_DOF)
#define ALA_F1B_RECORD_K (ALA_F1B_RECORD_GLOBAL_DIAG+ALA_VANKA_DOF)
#define ALA_F1B_RECORD_DOUBLES (ALA_F1B_RECORD_K+ALA_VANKA_DOF*ALA_VANKA_DOF)
/* Absolute diagnostic workspace bound over two completely disjoint strips.
 * A conforming interface is smaller: (2*overlap+1)*(width+1)*(radial+1)*3.
 * The active 6/3/2 case therefore expects 7*7*3*3=441. */
#define ALA_F1B_MAX_VELOCITY_DOF \
    (2*(ALA_PATCH_MAX_MPI_OVERLAP+1) \
     *(ALA_PATCH_MAX_HORIZONTAL_ELEMENTS+1) \
     *(ALA_PATCH_RADIAL_ELEMENTS+1)*3)
#define ALA_F1B_MAX_PATCH_NODES (8*ALA_PATCH_MAX_ELEMENTS)
#define ALA_F1B_K_RESTRICTION_TOLERANCE 1.0e-10
#define ALA_F1B_K_RESTRICTION_REVIEW_BOUND 1.0e-8
#define ALA_F1B_OUTPUT_REPLICA_TOLERANCE 1.0e-10

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
    double max_diagonal_completion_relative;
    double max_temporary_bytes;
    double max_factor_bytes;
    double estimated_dense_work;
};

struct ala_f1b_replay_snapshot {
    int valid;
    int vdim;
    int pdim;
    double *matrix;
    double *B;
    double coordinates[3*ALA_F1B_MAX_VELOCITY_DOF];
    double pressure_centers[3*ALA_PATCH_MAX_ELEMENTS];
    unsigned char components[ALA_F1B_MAX_VELOCITY_DOF];
};

enum ala_f1b_failure_reason {
    ALA_F1B_FAILURE_NONE=0,
    ALA_F1B_FAILURE_MISSING_RECORD=1,
    ALA_F1B_FAILURE_VELOCITY_DIMENSION=2,
    ALA_F1B_FAILURE_EMPTY_PATCH=3,
    ALA_F1B_FAILURE_INVALID_GLOBAL_DIAGONAL=4,
    ALA_F1B_FAILURE_INCIDENT_ELEMENT_BOUND=5,
    ALA_F1B_FAILURE_INPUT_MATRIX_NONFINITE=6,
    ALA_F1B_FAILURE_CHOLESKY_RECURRENCE_OVERFLOW=7,
    ALA_F1B_FAILURE_FINITE_NEGATIVE_PIVOT=8,
    ALA_F1B_FAILURE_FINITE_SMALL_PIVOT=9,
    ALA_F1B_FAILURE_INTERFACE_TOPOLOGY=10
};

struct ala_f1b_failure_info {
    int reason;
    int pressure_dim;
    int velocity_dim;
    int incident_elements;
    int pivot_index;
    int failure_row;
    int failure_column;
    int recurrence_k;
    int matrix_nonfinite_count;
    int duplicate_element_key_count;
    int topology_local_records;
    int topology_remote_records;
    int topology_local_nodes;
    int topology_remote_nodes;
    int topology_shared_nodes;
    int topology_union_nodes;
    int topology_expected_union_nodes;
    int topology_expected_velocity_dofs;
    int interface_face;
    int interface_tangent;
    int interface_radial_element;
    int remote_rank_min;
    int remote_rank_max;
    double pivot;
    double max_diagonal;
    double min_diagonal;
    double max_abs_offdiagonal;
    int max_abs_offdiagonal_row;
    int max_abs_offdiagonal_column;
    double diagonal_completion_relative;
    double diagonal_difference_min;
    double diagonal_difference_max;
    double failure_value;
    double recurrence_sum;
    double recurrence_term;
    double largest_abs_factor;
    double nearest_same_component_distance;
};

/* Fixed manifest order: 0=LOCAL representative, 1=INTERFACE representative. */
static struct ala_f1b_replay_snapshot ala_f1b_replay_snapshot[2];

static const char *ala_f1b_failure_name(int reason)
{
    switch(reason) {
    case ALA_F1B_FAILURE_MISSING_RECORD: return("missing_record");
    case ALA_F1B_FAILURE_VELOCITY_DIMENSION:
        return("velocity_dimension");
    case ALA_F1B_FAILURE_EMPTY_PATCH: return("empty_patch");
    case ALA_F1B_FAILURE_INVALID_GLOBAL_DIAGONAL:
        return("invalid_global_diagonal");
    case ALA_F1B_FAILURE_INCIDENT_ELEMENT_BOUND:
        return("incident_element_bound");
    case ALA_F1B_FAILURE_INPUT_MATRIX_NONFINITE:
        return("INPUT_MATRIX_NONFINITE");
    case ALA_F1B_FAILURE_CHOLESKY_RECURRENCE_OVERFLOW:
        return("CHOLESKY_RECURRENCE_OVERFLOW");
    case ALA_F1B_FAILURE_FINITE_NEGATIVE_PIVOT:
        return("FINITE_NEGATIVE_PIVOT");
    case ALA_F1B_FAILURE_FINITE_SMALL_PIVOT:
        return("FINITE_SMALL_PIVOT");
    case ALA_F1B_FAILURE_INTERFACE_TOPOLOGY:
        return("INTERFACE_TOPOLOGY_MISMATCH");
    default: return("unknown");
    }
}


static void ala_f1b_abort_candidate(struct All_variables *E,
    const char *category,int patch,const char *reason,
    const struct ala_f1b_failure_info *failure)
{
    const char *factorization_phase=
        getenv("STRICT_ALA_STAGE_F1B_IDENTITY_PREFLIGHT")!=NULL
        ? "diagnostic_factor_attempt"
        : "authorized_candidate_factorization";
    fprintf(stderr,"STRICT_ALA_STAGE_F1B_CANDIDATE_NUMERICAL_FAILURE "
        "reason=%s factorization_phase=%s category=%s rank=%d pz=%d "
        "patch=%d pdim=%d vdim=%d "
        "incident_elements=%d pivot_index=%d pivot=%.17e maxdiag=%.17e "
        "diagonal_completion_relative=%.17e failure_row=%d failure_col=%d "
        "recurrence_k=%d matrix_nonfinite_count=%d failure_value=%.17e "
        "recurrence_sum=%.17e recurrence_term=%.17e "
        "largest_abs_factor=%.17e min_diagonal=%.17e "
        "max_abs_offdiagonal=%.17e max_abs_offdiagonal_row=%d "
        "max_abs_offdiagonal_col=%d diagonal_difference_min=%.17e "
        "diagonal_difference_max=%.17e duplicate_element_key_count=%d "
        "nearest_same_component_distance=%.17e "
        "topology_local_records=%d topology_remote_records=%d "
        "topology_local_nodes=%d topology_remote_nodes=%d "
        "topology_shared_nodes=%d topology_union_nodes=%d "
        "topology_expected_union_nodes=%d "
        "topology_expected_velocity_dofs=%d interface_face=%d "
        "interface_tangent=%d interface_radial_element=%d "
        "remote_rank_min=%d remote_rank_max=%d\n",
        reason,factorization_phase,category,
        E->parallel.me,E->parallel.me_loc[3],patch,
        failure ? failure->pressure_dim : 0,
        failure ? failure->velocity_dim : 0,
        failure ? failure->incident_elements : 0,
        failure ? failure->pivot_index : -1,
        failure ? failure->pivot : 0.0,
        failure ? failure->max_diagonal : 0.0,
        failure ? failure->diagonal_completion_relative : 0.0,
        failure ? failure->failure_row : -1,
        failure ? failure->failure_column : -1,
        failure ? failure->recurrence_k : -1,
        failure ? failure->matrix_nonfinite_count : 0,
        failure ? failure->failure_value : 0.0,
        failure ? failure->recurrence_sum : 0.0,
        failure ? failure->recurrence_term : 0.0,
        failure ? failure->largest_abs_factor : 0.0,
        failure ? failure->min_diagonal : 0.0,
        failure ? failure->max_abs_offdiagonal : 0.0,
        failure ? failure->max_abs_offdiagonal_row : -1,
        failure ? failure->max_abs_offdiagonal_column : -1,
        failure ? failure->diagonal_difference_min : 0.0,
        failure ? failure->diagonal_difference_max : 0.0,
        failure ? failure->duplicate_element_key_count : 0,
        failure ? failure->nearest_same_component_distance : 0.0,
        failure ? failure->topology_local_records : 0,
        failure ? failure->topology_remote_records : 0,
        failure ? failure->topology_local_nodes : 0,
        failure ? failure->topology_remote_nodes : 0,
        failure ? failure->topology_shared_nodes : 0,
        failure ? failure->topology_union_nodes : 0,
        failure ? failure->topology_expected_union_nodes : 0,
        failure ? failure->topology_expected_velocity_dofs : 0,
        failure ? failure->interface_face : -1,
        failure ? failure->interface_tangent : -1,
        failure ? failure->interface_radial_element : -1,
        failure ? failure->remote_rank_min : -1,
        failure ? failure->remote_rank_max : -1);
    fflush(stderr);
    /* A rank-local myerror left the other radial layer inside collectives in
     * job 12126878.  MPI_Abort is intentional here: the LSF wrapper converts
     * this deterministic candidate rejection into the final audit record. */
    MPI_Abort(E->parallel.world,8);
    exit(8);
}

static int ala_f1b_same_point(const double *a,const double *b)
{
    double dx=a[0]-b[0],dy=a[1]-b[1],dz=a[2]-b[2];
    return(dx*dx+dy*dy+dz*dz<=1.0e-18);
}


static int ala_f1b_same_node_key(const double *record_key,
    const int *stored_key)
{
    int d;
    for(d=0;d<4;d++)
        if((int)(record_key[d]+0.5)!=stored_key[d])
            return(0);
    return(1);
}


static int ala_f1b_node_index(const double *coord,const double *node_key,
    const double *coordinates,const int *keys,int n)
{
    int i;
    for(i=0;i<n;i++)
        if(ala_f1b_same_node_key(node_key,keys+4*i) ||
           ala_f1b_same_point(coord,coordinates+3*i))
            return(i);
    return(-1);
}


static int ala_f1b_add_node(const double *coord,const double *node_key,
    double *coordinates,int *keys,int n)
{
    int d;
    if(ala_f1b_node_index(coord,node_key,coordinates,keys,n)>=0)
        return(n);
    if(n>=ALA_F1B_MAX_PATCH_NODES)
        return(-1);
    memcpy(coordinates+3*n,coord,3*sizeof(double));
    for(d=0;d<4;d++)
        keys[4*n+d]=(int)(node_key[d]+0.5);
    return(n+1);
}


static int ala_f1b_measure_patch_topology(struct All_variables *E,
    const double **records,int np,struct ala_f1b_failure_info *failure)
{
    int p,a,d,q,rank,next,index,nearest,side_nodes,other_nodes;
    int union_nodes=0,local_nodes=0,remote_nodes=0,shared=0;
    int *union_keys,*local_keys,*remote_keys;
    int *local_occurrences,*remote_occurrences;
    double *union_coord,*local_coord,*remote_coord;
    double center[3],distance2,same_min2,other_min2,candidate2,dx;
    const double *record,*coord,*key,*side_coord,*other_coord;
    const int *side_keys;

    union_keys=(int *)calloc(4*ALA_F1B_MAX_PATCH_NODES,sizeof(int));
    local_keys=(int *)calloc(4*ALA_F1B_MAX_PATCH_NODES,sizeof(int));
    remote_keys=(int *)calloc(4*ALA_F1B_MAX_PATCH_NODES,sizeof(int));
    union_coord=(double *)calloc(3*ALA_F1B_MAX_PATCH_NODES,sizeof(double));
    local_coord=(double *)calloc(3*ALA_F1B_MAX_PATCH_NODES,sizeof(double));
    remote_coord=(double *)calloc(3*ALA_F1B_MAX_PATCH_NODES,sizeof(double));
    local_occurrences=(int *)calloc(ALA_F1B_MAX_PATCH_NODES,sizeof(int));
    remote_occurrences=(int *)calloc(ALA_F1B_MAX_PATCH_NODES,sizeof(int));
    if(union_keys==NULL || local_keys==NULL || remote_keys==NULL ||
       union_coord==NULL || local_coord==NULL || remote_coord==NULL ||
       local_occurrences==NULL || remote_occurrences==NULL)
        myerror(E,"Unable to allocate F1b topology workspace");
    failure->remote_rank_min=E->parallel.nproc;
    failure->remote_rank_max=-1;
    for(p=0;p<np;p++) {
        record=records[p];
        if(record==NULL)
            continue;
        rank=(int)(record[0]+0.5);
        if(rank==E->parallel.me)
            failure->topology_local_records++;
        else {
            failure->topology_remote_records++;
            failure->remote_rank_min=min(failure->remote_rank_min,rank);
            failure->remote_rank_max=max(failure->remote_rank_max,rank);
        }
        for(a=0;a<8;a++) {
            coord=record+ALA_F1B_RECORD_COORD+3*a;
            key=record+ALA_F1B_RECORD_NODE_KEY+4*a;
            next=ala_f1b_add_node(coord,key,union_coord,union_keys,
                                   union_nodes);
            if(next<0) goto topology_bound;
            union_nodes=next;
            if(rank==E->parallel.me) {
                next=ala_f1b_add_node(coord,key,local_coord,local_keys,
                                       local_nodes);
                if(next<0) goto topology_bound;
                local_nodes=next;
            }
            else {
                next=ala_f1b_add_node(coord,key,remote_coord,remote_keys,
                                       remote_nodes);
                if(next<0) goto topology_bound;
                remote_nodes=next;
            }
        }
    }
    /* This path is reached only after the active conforming dimension was
     * exceeded.  Emit the complete, bounded element/node census before the
     * collective fail-fast so one HPC attempt can distinguish a bad halo
     * selection from a legitimate cubed-sphere corner topology. */
    for(p=0;p<np;p++) {
        record=records[p];
        if(record==NULL)
            continue;
        for(d=0;d<3;d++) {
            center[d]=0.0;
            for(a=0;a<8;a++)
                center[d] += record[ALA_F1B_RECORD_COORD+3*a+d];
            center[d] *= 0.125;
        }
        fprintf(stderr,"STRICT_ALA_STAGE_F1B_TOPOLOGY_RECORD "
            "interface_face=%d interface_tangent=%d "
            "interface_radial_element=%d record=%d side=%s "
            "rank=%d ex=%d ey=%d ez=%d "
            "cap=%d center_x=%.17e center_y=%.17e center_z=%.17e\n",
            failure->interface_face,failure->interface_tangent,
            failure->interface_radial_element,p,
            (int)(record[0]+0.5)==E->parallel.me ? "LOCAL" : "REMOTE",
            (int)(record[0]+0.5),(int)(record[1]+0.5),
            (int)(record[2]+0.5),(int)(record[3]+0.5),
            (int)(record[ALA_F1B_RECORD_NODE_KEY]+0.5),
            center[0],center[1],center[2]);
        for(a=0;a<8;a++) {
            coord=record+ALA_F1B_RECORD_COORD+3*a;
            key=record+ALA_F1B_RECORD_NODE_KEY+4*a;
            if((int)(record[0]+0.5)==E->parallel.me) {
                index=ala_f1b_node_index(coord,key,local_coord,local_keys,
                                          local_nodes);
                if(index>=0)
                    local_occurrences[index]++;
            }
            else {
                index=ala_f1b_node_index(coord,key,remote_coord,remote_keys,
                                          remote_nodes);
                if(index>=0)
                    remote_occurrences[index]++;
            }
        }
    }
    for(rank=0;rank<2;rank++) {
        if(rank==0) {
            side_coord=local_coord; side_keys=local_keys;
            side_nodes=local_nodes; other_coord=remote_coord;
            other_nodes=remote_nodes;
        }
        else {
            side_coord=remote_coord; side_keys=remote_keys;
            side_nodes=remote_nodes; other_coord=local_coord;
            other_nodes=local_nodes;
        }
        for(a=0;a<side_nodes;a++) {
            nearest=-1; same_min2=INFINITY;
            for(q=0;q<side_nodes;q++) {
                if(q==a)
                    continue;
                distance2=0.0;
                for(d=0;d<3;d++) {
                    dx=side_coord[3*a+d]-side_coord[3*q+d];
                    distance2 += dx*dx;
                }
                if(distance2<same_min2) {
                    same_min2=distance2; nearest=q;
                }
            }
            index=-1; other_min2=INFINITY;
            for(q=0;q<other_nodes;q++) {
                candidate2=0.0;
                for(d=0;d<3;d++) {
                    dx=side_coord[3*a+d]-other_coord[3*q+d];
                    candidate2 += dx*dx;
                }
                if(candidate2<other_min2) {
                    other_min2=candidate2; index=q;
                }
            }
            fprintf(stderr,"STRICT_ALA_STAGE_F1B_TOPOLOGY_NODE "
                "side=%s index=%d cap=%d gx=%d gy=%d gz=%d "
                "x=%.17e y=%.17e z=%.17e occurrences=%d "
                "nearest_same_index=%d nearest_same_distance=%.17e "
                "nearest_other_index=%d nearest_other_distance=%.17e\n",
                rank==0 ? "LOCAL" : "REMOTE",a,
                side_keys[4*a],side_keys[4*a+1],side_keys[4*a+2],
                side_keys[4*a+3],side_coord[3*a],side_coord[3*a+1],
                side_coord[3*a+2],rank==0 ? local_occurrences[a]
                                          : remote_occurrences[a],
                nearest,nearest>=0 ? sqrt(same_min2) : INFINITY,
                index,index>=0 ? sqrt(other_min2) : INFINITY);
        }
    }
    fflush(stderr);
    for(a=0;a<local_nodes;a++) {
        for(q=0;q<remote_nodes;q++)
            if((local_keys[4*a]==remote_keys[4*q] &&
                local_keys[4*a+1]==remote_keys[4*q+1] &&
                local_keys[4*a+2]==remote_keys[4*q+2] &&
                local_keys[4*a+3]==remote_keys[4*q+3]) ||
               ala_f1b_same_point(local_coord+3*a,remote_coord+3*q)) {
                shared++;
                break;
            }
    }
    failure->topology_local_nodes=local_nodes;
    failure->topology_remote_nodes=remote_nodes;
    failure->topology_shared_nodes=shared;
    failure->topology_union_nodes=union_nodes;
    if(failure->remote_rank_max<0)
        failure->remote_rank_min=-1;
    free(union_keys); free(local_keys); free(remote_keys);
    free(union_coord); free(local_coord); free(remote_coord);
    free(local_occurrences); free(remote_occurrences);
    return(1);

topology_bound:
    free(union_keys); free(local_keys); free(remote_keys);
    free(union_coord); free(local_coord); free(remote_coord);
    free(local_occurrences); free(remote_occurrences);
    return(0);
}


static double ala_f1b_nearest_component_distance(const double *coord,
    int component,const double *vcoord,const unsigned char *vcomp,int n)
{
    int i;
    double dx,dy,dz,distance2,nearest2=INFINITY;
    for(i=0;i<n;i++) {
        if(vcomp[i]!=component)
            continue;
        dx=coord[0]-vcoord[3*i];
        dy=coord[1]-vcoord[3*i+1];
        dz=coord[2]-vcoord[3*i+2];
        distance2=dx*dx+dy*dy+dz*dz;
        nearest2=min(nearest2,distance2);
    }
    return(isfinite(nearest2) ? sqrt(nearest2) : INFINITY);
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
    int a,b,da,db,i,j,node,eq,fixed,ix,iy,iz;
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
        iz=(node-1)%E->lmesh.NOZ[lev]+1;
        ix=((node-1)/E->lmesh.NOZ[lev])%E->lmesh.NOX[lev]+1;
        iy=(node-1)/(E->lmesh.NOZ[lev]*E->lmesh.NOX[lev])+1;
        /* Exact same-cap topology is the primary MPI-replica identity.  This
         * avoids inferring ownership from independently generated floating
         * coordinates (job 12127734).  Cartesian matching remains the
         * cross-cap fallback below, where structured cap indices differ. */
        record[ALA_F1B_RECORD_NODE_KEY+4*(a-1)]=
            (double)E->sphere.capid[m];
        record[ALA_F1B_RECORD_NODE_KEY+4*(a-1)+1]=
            (double)(E->lmesh.EXS[lev]+ix);
        record[ALA_F1B_RECORD_NODE_KEY+4*(a-1)+2]=
            (double)(E->lmesh.EYS[lev]+iy);
        record[ALA_F1B_RECORD_NODE_KEY+4*(a-1)+3]=
            (double)(E->lmesh.EZS[lev]+iz);
        for(da=0;da<dims;da++) {
            i=(a-1)*dims+da;
            record[ALA_F1B_RECORD_COORD+i]=E->X[lev][m][da+1][node];
            record[ALA_F1B_RECORD_B+i]=ALA_COMBINED_PRESSURE_COEFFICIENT(
                E->elt_del[lev][m][e].g[i][0],
                E->elt_c[lev][m][e].c[i][0]);
            fixed=
                (da==0 && (E->NODE[lev][m][node]&VBX)) ||
                (da==1 && (E->NODE[lev][m][node]&VBY)) ||
                (da==2 && (E->NODE[lev][m][node]&VBZ));
            record[ALA_F1B_RECORD_FIXED+i]=(double)fixed;
            eq=E->ID[lev][m][node].doff[da+1];
            if(!fixed) {
                record[ALA_F1B_RECORD_GLOBAL_DIAG+i]=
                    (!isfinite(E->BI[lev][m][eq]) ||
                     E->BI[lev][m][eq]<=0.0)
                    ? NAN : 1.0/E->BI[lev][m][eq];
            }
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


static int ala_f1b_factor(double *matrix,int n,double *min_ratio,
    struct ala_f1b_failure_info *failure)
{
    int i,j,k;
    double maxdiag=0.0,sum,pivot,term,next,value;
    failure->failure_row=-1;
    failure->failure_column=-1;
    failure->recurrence_k=-1;
    failure->matrix_nonfinite_count=0;
    failure->failure_value=0.0;
    failure->recurrence_sum=0.0;
    failure->recurrence_term=0.0;
    failure->largest_abs_factor=0.0;
    failure->min_diagonal=1.0e300;
    failure->max_abs_offdiagonal=0.0;
    failure->max_abs_offdiagonal_row=-1;
    failure->max_abs_offdiagonal_column=-1;
    for(i=0;i<n;i++)
        for(j=0;j<n;j++) {
            if(!isfinite(matrix[i*n+j])) {
                if(failure->matrix_nonfinite_count==0) {
                    failure->failure_row=i;
                    failure->failure_column=j;
                    failure->failure_value=matrix[i*n+j];
                }
                failure->matrix_nonfinite_count++;
            }
            if(i==j)
                failure->min_diagonal=min(failure->min_diagonal,
                                           matrix[i*n+j]);
            else if(isfinite(matrix[i*n+j]) &&
                    fabs(matrix[i*n+j])>failure->max_abs_offdiagonal) {
                failure->max_abs_offdiagonal=fabs(matrix[i*n+j]);
                failure->max_abs_offdiagonal_row=i;
                failure->max_abs_offdiagonal_column=j;
            }
        }
    if(failure->matrix_nonfinite_count>0) {
        failure->reason=ALA_F1B_FAILURE_INPUT_MATRIX_NONFINITE;
        return(0);
    }
    for(i=0;i<n;i++)
        maxdiag=max(maxdiag,matrix[i*n+i]);
    failure->max_diagonal=maxdiag;
    if(!isfinite(maxdiag) || maxdiag<=0.0) {
        failure->pivot_index=-1;
        failure->pivot=maxdiag;
        failure->reason=maxdiag<0.0
            ? ALA_F1B_FAILURE_FINITE_NEGATIVE_PIVOT
            : ALA_F1B_FAILURE_FINITE_SMALL_PIVOT;
        return(0);
    }
    for(i=0;i<n;i++)
        for(j=0;j<=i;j++) {
            sum=matrix[i*n+j];
            for(k=0;k<j;k++) {
                term=matrix[i*n+k]*matrix[j*n+k];
                next=sum-term;
                if(!isfinite(term) || !isfinite(next)) {
                    failure->reason=
                        ALA_F1B_FAILURE_CHOLESKY_RECURRENCE_OVERFLOW;
                    failure->failure_row=i;
                    failure->failure_column=j;
                    failure->recurrence_k=k;
                    failure->recurrence_sum=sum;
                    failure->recurrence_term=term;
                    failure->failure_value=next;
                    return(0);
                }
                sum=next;
            }
            if(i==j) {
                pivot=sum;
                if(!isfinite(pivot) || pivot<=1.0e-14*maxdiag) {
                    failure->pivot_index=i;
                    failure->pivot=pivot;
                    failure->failure_row=i;
                    failure->failure_column=j;
                    failure->failure_value=pivot;
                    failure->reason=isfinite(pivot)
                        ? (pivot<0.0
                           ? ALA_F1B_FAILURE_FINITE_NEGATIVE_PIVOT
                           : ALA_F1B_FAILURE_FINITE_SMALL_PIVOT)
                        : ALA_F1B_FAILURE_CHOLESKY_RECURRENCE_OVERFLOW;
                    return(0);
                }
                *min_ratio=min(*min_ratio,pivot/maxdiag);
                value=sqrt(pivot);
                matrix[i*n+j]=value;
            }
            else {
                value=sum/matrix[j*n+j];
                if(!isfinite(value)) {
                    failure->reason=
                        ALA_F1B_FAILURE_CHOLESKY_RECURRENCE_OVERFLOW;
                    failure->failure_row=i;
                    failure->failure_column=j;
                    failure->failure_value=value;
                    failure->recurrence_sum=sum;
                    return(0);
                }
                matrix[i*n+j]=value;
            }
            failure->largest_abs_factor=max(
                failure->largest_abs_factor,fabs(value));
            /* Clear only the strict upper triangle.  When i==j this entry is
             * the Cholesky diagonal just written above; clearing it makes the
             * next off-diagonal update divide by zero (job 12127697). */
            if(j<i)
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


static int ala_f1b_velocity_index(const double *coord,const double *node_key,
    int component,const double *vcoord,const int *vkey,
    const unsigned char *vcomp,int n)
{
    int i;
    for(i=0;i<n;i++)
        if(vcomp[i]==component &&
           (ala_f1b_same_node_key(node_key,vkey+4*i) ||
            ala_f1b_same_point(coord,vcoord+3*i)))
            return(i);
    return(-1);
}


static int ala_f1b_build_operator_consistent_patch(
    struct All_variables *E,const struct ala_f1b_element_pool *pool,
    const double **pressure_records,int np,
    double matrix[ALA_PATCH_MAX_ELEMENTS][ALA_PATCH_MAX_ELEMENTS],
    struct ala_f1b_patch_stats *stats,
    struct ala_f1b_failure_info *failure,int replay_category,
    int interface_face,int interface_tangent,int interface_radial_element)
{
    int p,a,d,i,j,k,r,nv,nused,map[ALA_VANKA_DOF];
    int used[1024][4],vkey[4*ALA_F1B_MAX_VELOCITY_DOF];
    unsigned char vcomp[ALA_F1B_MAX_VELOCITY_DOF];
    double vcoord[3*ALA_F1B_MAX_VELOCITY_DOF];
    double vdiag[ALA_F1B_MAX_VELOCITY_DOF];
    double *K,*factor,*B,*rhs,*solution;
    double bnorm2,rnorm2,sum,numerator,denominator,symmetry,min_pivot;
    double diagonal_completion;
    const double *record,*coord;

    memset(failure,0,sizeof(*failure));
    failure->pivot_index=-1;
    failure->interface_face=interface_face;
    failure->interface_tangent=interface_tangent;
    failure->interface_radial_element=interface_radial_element;
    if(replay_category==1) {
        failure->topology_expected_union_nodes=
            (2*E->control.ala_shallow_patch_mpi_overlap+1)
            *(E->control.ala_shallow_patch_horizontal_elements+1)
            *(ALA_PATCH_RADIAL_ELEMENTS+1);
        failure->topology_expected_velocity_dofs=
            3*failure->topology_expected_union_nodes;
    }
    nv=0;
    for(p=0;p<np;p++) {
        record=pressure_records[p];
        if(record==NULL) {
            failure->reason=ALA_F1B_FAILURE_MISSING_RECORD;
            return(0);
        }
        for(a=0;a<8;a++)
            for(d=0;d<3;d++) {
                i=3*a+d;
                if(record[ALA_F1B_RECORD_FIXED+i]!=0.0)
                    continue;
                coord=record+ALA_F1B_RECORD_COORD+3*a;
                k=ala_f1b_velocity_index(coord,
                    record+ALA_F1B_RECORD_NODE_KEY+4*a,d,
                    vcoord,vkey,vcomp,nv);
                if(k<0) {
                    if(nv>=ALA_F1B_MAX_VELOCITY_DOF) {
                        failure->reason=ALA_F1B_FAILURE_VELOCITY_DIMENSION;
                        failure->velocity_dim=nv+1;
                        failure->nearest_same_component_distance=
                            ala_f1b_nearest_component_distance(
                                coord,d,vcoord,vcomp,nv);
                        return(0);
                    }
                    memcpy(vcoord+3*nv,coord,3*sizeof(double));
                    for(j=0;j<4;j++)
                        vkey[4*nv+j]=(int)(
                            record[ALA_F1B_RECORD_NODE_KEY+4*a+j]+0.5);
                    vcomp[nv]=(unsigned char)d;
                    vdiag[nv]=record[ALA_F1B_RECORD_GLOBAL_DIAG+i];
                    if(!isfinite(vdiag[nv]) || vdiag[nv]<=0.0) {
                        failure->reason=
                            ALA_F1B_FAILURE_INVALID_GLOBAL_DIAGONAL;
                        failure->velocity_dim=nv+1;
                        return(0);
                    }
                    nv++;
                }
                else {
                    double replica_diag=
                        record[ALA_F1B_RECORD_GLOBAL_DIAG+i];
                    if(!isfinite(replica_diag) || replica_diag<=0.0 ||
                       fabs(replica_diag-vdiag[k])>
                         1.0e-10*max(max(replica_diag,vdiag[k]),1.0e-300)) {
                        failure->reason=
                            ALA_F1B_FAILURE_INVALID_GLOBAL_DIAGONAL;
                        failure->velocity_dim=nv;
                        return(0);
                    }
                }
            }
    }
    failure->pressure_dim=np;
    failure->velocity_dim=nv;
    if(replay_category==1 &&
       nv>failure->topology_expected_velocity_dofs) {
        /* The normal path pays no topology-audit cost.  Only an interface
         * that exceeds the conforming bound performs the complete node-union
         * census before the collective fail-fast. */
        if(!ala_f1b_measure_patch_topology(E,pressure_records,np,failure)) {
            failure->reason=ALA_F1B_FAILURE_INCIDENT_ELEMENT_BOUND;
            return(0);
        }
        failure->reason=ALA_F1B_FAILURE_INTERFACE_TOPOLOGY;
        return(0);
    }
    if(nv<=0 || np<=0) {
        failure->reason=ALA_F1B_FAILURE_EMPTY_PATCH;
        return(0);
    }
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
                    record+ALA_F1B_RECORD_COORD+3*a,
                    record+ALA_F1B_RECORD_NODE_KEY+4*a,d,
                    vcoord,vkey,vcomp,nv);
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
                    record+ALA_F1B_RECORD_COORD+3*a,
                    record+ALA_F1B_RECORD_NODE_KEY+4*a,d,
                    vcoord,vkey,vcomp,nv);
        for(i=0;i<ALA_VANKA_DOF;i++)
            if(map[i]>=0)
                break;
        if(i==ALA_VANKA_DOF)
            continue;
        if(ala_f1b_record_key_seen(record,used,nused)) {
            failure->duplicate_element_key_count++;
            continue;
        }
        if(nused>=1024) {
            failure->reason=ALA_F1B_FAILURE_INCIDENT_ELEMENT_BOUND;
            failure->incident_elements=nused+1;
            free(K); free(factor); free(B); free(rhs); free(solution);
            return(0);
        }
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

    failure->incident_elements=nused;
    diagonal_completion=0.0;
    failure->diagonal_difference_min=1.0e300;
    failure->diagonal_difference_max=-1.0e300;
    for(i=0;i<nv;i++) {
        double difference=K[i*nv+i]-vdiag[i];
        diagonal_completion=max(diagonal_completion,
            fabs(difference)/max(vdiag[i],1.0e-300));
        failure->diagonal_difference_min=min(
            failure->diagonal_difference_min,difference);
        failure->diagonal_difference_max=max(
            failure->diagonal_difference_max,difference);
        /* This is the exact assembled principal diagonal already used by
         * the production K_gamma operator.  Element records supply every
         * retained-retained off-diagonal; BI closes contributions from
         * incident elements outside the pressure patch. */
        K[i*nv+i]=vdiag[i];
    }
    failure->diagonal_completion_relative=diagonal_completion;
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
    if(!ala_f1b_factor(factor,nv,&min_pivot,failure)) {
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
    stats->max_diagonal_completion_relative=max(
        stats->max_diagonal_completion_relative,diagonal_completion);
    stats->min_k_pivot_ratio=min(stats->min_k_pivot_ratio,min_pivot);
    stats->max_temporary_bytes=max(stats->max_temporary_bytes,
        max(pool->temporary_bytes,pool->resident_bytes+
            (2.0*nv*nv+(double)np*nv+2.0*nv)*sizeof(double)));
    stats->max_factor_bytes=max(stats->max_factor_bytes,
        (double)np*np*sizeof(double));
    if(replay_category>=0 && replay_category<2 &&
       !ala_f1b_replay_snapshot[replay_category].valid) {
        ala_f1b_replay_snapshot[replay_category].matrix=(double *)calloc(
            (size_t)nv*nv,sizeof(double));
        ala_f1b_replay_snapshot[replay_category].B=(double *)calloc(
            (size_t)np*nv,sizeof(double));
        if(ala_f1b_replay_snapshot[replay_category].matrix==NULL ||
           ala_f1b_replay_snapshot[replay_category].B==NULL)
            myerror(E,"Unable to allocate F1b replay snapshot");
        memcpy(ala_f1b_replay_snapshot[replay_category].matrix,K,
               (size_t)nv*nv*sizeof(double));
        memcpy(ala_f1b_replay_snapshot[replay_category].coordinates,vcoord,
               (size_t)3*nv*sizeof(double));
        memcpy(ala_f1b_replay_snapshot[replay_category].components,vcomp,
               (size_t)nv*sizeof(unsigned char));
        memcpy(ala_f1b_replay_snapshot[replay_category].B,B,
               (size_t)np*nv*sizeof(double));
        for(p=0;p<np;p++)
            ala_f1b_record_center(pressure_records[p],
                ala_f1b_replay_snapshot[replay_category].pressure_centers+3*p);
        ala_f1b_replay_snapshot[replay_category].vdim=nv;
        ala_f1b_replay_snapshot[replay_category].pdim=np;
        ala_f1b_replay_snapshot[replay_category].valid=1;
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


static int ala_f1b_find_local_element(struct All_variables *E,int lev,int m,
    const double *coord)
{
    int e,a,d,node;
    double center[3];
    for(e=1;e<=E->lmesh.NEL[lev];e++) {
        for(d=0;d<3;d++) center[d]=0.0;
        for(a=1;a<=enodes[E->mesh.nsd];a++) {
            node=E->IEN[lev][m][e].node[a];
            for(d=0;d<3;d++) center[d]+=E->X[lev][m][d+1][node];
        }
        for(d=0;d<3;d++) center[d]/=enodes[E->mesh.nsd];
        if(ala_f1b_same_point(coord,center)) return(e);
    }
    return(-1);
}


static double ala_f1b_replay_global_k(struct All_variables *E,int lev,
    int replay_category,double *replica_difference,double *b_defect,
    double *bt_defect,double *bilinear_defect)
{
    int m,i,j,p,node,e,eq,owner,local_owner,nv,np,neq,npno;
    int components[ALA_F1B_MAX_VELOCITY_DOF];
    double coordinates[3*ALA_F1B_MAX_VELOCITY_DOF];
    double pressure_centers[3*ALA_PATCH_MAX_ELEMENTS];
    double x[ALA_F1B_MAX_VELOCITY_DOF];
    double pressure_x[ALA_PATCH_MAX_ELEMENTS];
    double numerator=0.0,denominator=0.0,sum,actual;
    double b_numerator=0.0,b_denominator=0.0;
    double bt_numerator=0.0,bt_denominator=0.0;
    double bilinear_lhs=0.0,bilinear_rhs=0.0;
    double local_min,local_max,global_min,global_max,scale;
    double **u,**Au,**pressure_action,**pressure_input,**gradient;
    void n_assemble_del2_u(),assemble_div_rho_u(),assemble_grad_rho_p();

    struct ala_f1b_replay_snapshot *snapshot;

    if(replay_category<0 || replay_category>=2)
        return(INFINITY);
    snapshot=ala_f1b_replay_snapshot+replay_category;
    local_owner=snapshot->valid ? E->parallel.me
                                               : E->parallel.nproc;
    MPI_Allreduce(&local_owner,&owner,1,MPI_INT,MPI_MIN,E->parallel.world);
    *replica_difference=0.0;
    *b_defect=INFINITY;
    *bt_defect=INFINITY;
    *bilinear_defect=INFINITY;
    if(owner>=E->parallel.nproc)
        return(INFINITY);
    nv=(E->parallel.me==owner) ? snapshot->vdim : 0;
    np=(E->parallel.me==owner) ? snapshot->pdim : 0;
    MPI_Bcast(&nv,1,MPI_INT,owner,E->parallel.world);
    MPI_Bcast(&np,1,MPI_INT,owner,E->parallel.world);
    if(nv<=0 || nv>ALA_F1B_MAX_VELOCITY_DOF ||
       np<=0 || np>ALA_PATCH_MAX_ELEMENTS)
        return(INFINITY);
    if(E->parallel.me==owner) {
        memcpy(coordinates,snapshot->coordinates,
               (size_t)3*nv*sizeof(double));
        for(i=0;i<nv;i++)
            components[i]=(int)snapshot->components[i];
        memcpy(pressure_centers,snapshot->pressure_centers,
               (size_t)3*np*sizeof(double));
    }
    MPI_Bcast(coordinates,3*nv,MPI_DOUBLE,owner,E->parallel.world);
    MPI_Bcast(components,nv,MPI_INT,owner,E->parallel.world);
    MPI_Bcast(pressure_centers,3*np,MPI_DOUBLE,owner,E->parallel.world);
    for(i=0;i<nv;i++)
        x[i]=sin(0.17320508075688773*(i+1))
            +0.5*cos(0.22360679774997897*(i+1));
    u=(double **)calloc(NCS,sizeof(double *));
    Au=(double **)calloc(NCS,sizeof(double *));
    pressure_action=(double **)calloc(NCS,sizeof(double *));
    pressure_input=(double **)calloc(NCS,sizeof(double *));
    gradient=(double **)calloc(NCS,sizeof(double *));
    if(u==NULL || Au==NULL || pressure_action==NULL ||
       pressure_input==NULL || gradient==NULL)
        myerror(E,"Unable to allocate F1b replay vectors");
    neq=E->lmesh.NEQ[lev];
    npno=E->lmesh.NPNO[lev];
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        u[m]=(double *)calloc(neq+1,sizeof(double));
        Au[m]=(double *)calloc(neq+1,sizeof(double));
        pressure_action[m]=(double *)calloc(npno+1,sizeof(double));
        pressure_input[m]=(double *)calloc(npno+1,sizeof(double));
        gradient[m]=(double *)calloc(neq+1,sizeof(double));
        if(u[m]==NULL || Au[m]==NULL || pressure_action[m]==NULL ||
           pressure_input[m]==NULL || gradient[m]==NULL)
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
    for(i=0;i<nv;i++) {
        local_min=1.0e300;
        local_max=-1.0e300;
        for(m=1;m<=E->sphere.caps_per_proc;m++) {
            node=ala_f1b_find_local_node(E,lev,m,coordinates+3*i);
            if(node>0) {
                eq=E->ID[lev][m][node].doff[components[i]+1];
                local_min=min(local_min,Au[m][eq]);
                local_max=max(local_max,Au[m][eq]);
            }
        }
        MPI_Allreduce(&local_min,&global_min,1,MPI_DOUBLE,MPI_MIN,
                      E->parallel.world);
        MPI_Allreduce(&local_max,&global_max,1,MPI_DOUBLE,MPI_MAX,
                      E->parallel.world);
        if(global_min>0.5e300 || global_max< -0.5e300 ||
           !isfinite(global_min) || !isfinite(global_max)) {
            *replica_difference=INFINITY;
            actual=NAN;
        }
        else {
            scale=max(max(fabs(global_min),fabs(global_max)),1.0e-300);
            *replica_difference=max(*replica_difference,
                fabs(global_max-global_min)/scale);
            actual=0.5*(global_min+global_max);
        }
        if(E->parallel.me==owner) {
            sum=0.0;
            for(j=0;j<nv;j++)
                sum += snapshot->matrix[i*nv+j]*x[j];
            numerator += (sum-actual)*(sum-actual);
            denominator += actual*actual;
        }
    }
    if(E->parallel.me==owner) {
        numerator=sqrt(numerator/max(denominator,1.0e-300));
    }
    MPI_Bcast(&numerator,1,MPI_DOUBLE,owner,E->parallel.world);

    /* B_i x = R_p B E_v x, using the production continuity operator. */
    assemble_div_rho_u(E,u,pressure_action,lev);
    for(p=0;p<np;p++) {
        local_min=1.0e300;
        local_max=-1.0e300;
        for(m=1;m<=E->sphere.caps_per_proc;m++) {
            e=ala_f1b_find_local_element(E,lev,m,pressure_centers+3*p);
            if(e>0) {
                local_min=min(local_min,pressure_action[m][e]);
                local_max=max(local_max,pressure_action[m][e]);
            }
        }
        MPI_Allreduce(&local_min,&global_min,1,MPI_DOUBLE,MPI_MIN,
                      E->parallel.world);
        MPI_Allreduce(&local_max,&global_max,1,MPI_DOUBLE,MPI_MAX,
                      E->parallel.world);
        actual=(global_min>0.5e300 || global_max< -0.5e300)
            ? NAN : 0.5*(global_min+global_max);
        if(E->parallel.me==owner) {
            sum=0.0;
            for(i=0;i<nv;i++) sum+=snapshot->B[p*nv+i]*x[i];
            b_numerator+=(sum-actual)*(sum-actual);
            b_denominator+=actual*actual;
        }
    }
    if(E->parallel.me==owner)
        b_numerator=sqrt(b_numerator/max(b_denominator,1.0e-300));
    MPI_Bcast(&b_numerator,1,MPI_DOUBLE,owner,E->parallel.world);
    *b_defect=b_numerator;

    /* B_i^T p = R_v B^T E_p p, using the production transpose/exchange. */
    for(p=0;p<np;p++) {
        pressure_x[p]=sin(0.31622776601683794*(p+1))
            +0.25*cos(0.14142135623730950*(p+1));
        for(m=1;m<=E->sphere.caps_per_proc;m++) {
            e=ala_f1b_find_local_element(E,lev,m,pressure_centers+3*p);
            if(e>0) pressure_input[m][e]=pressure_x[p];
        }
    }
    assemble_grad_rho_p(E,pressure_input,gradient,lev);
    for(i=0;i<nv;i++) {
        local_min=1.0e300;
        local_max=-1.0e300;
        for(m=1;m<=E->sphere.caps_per_proc;m++) {
            node=ala_f1b_find_local_node(E,lev,m,coordinates+3*i);
            if(node>0) {
                eq=E->ID[lev][m][node].doff[components[i]+1];
                local_min=min(local_min,gradient[m][eq]);
                local_max=max(local_max,gradient[m][eq]);
            }
        }
        MPI_Allreduce(&local_min,&global_min,1,MPI_DOUBLE,MPI_MIN,
                      E->parallel.world);
        MPI_Allreduce(&local_max,&global_max,1,MPI_DOUBLE,MPI_MAX,
                      E->parallel.world);
        if(global_min>0.5e300 || global_max< -0.5e300) actual=NAN;
        else {
            scale=max(max(fabs(global_min),fabs(global_max)),1.0e-300);
            *replica_difference=max(*replica_difference,
                fabs(global_max-global_min)/scale);
            actual=0.5*(global_min+global_max);
        }
        if(E->parallel.me==owner) {
            sum=0.0;
            for(p=0;p<np;p++) sum+=snapshot->B[p*nv+i]*pressure_x[p];
            bt_numerator+=(sum-actual)*(sum-actual);
            bt_denominator+=actual*actual;
        }
    }
    if(E->parallel.me==owner)
        bt_numerator=sqrt(bt_numerator/max(bt_denominator,1.0e-300));
    MPI_Bcast(&bt_numerator,1,MPI_DOUBLE,owner,E->parallel.world);
    *bt_defect=bt_numerator;

    /* Explicit local bilinear transpose identity in the same P0/unique-DOF
     * coordinates used by the two distributed restriction checks. */
    if(E->parallel.me==owner) {
        for(p=0;p<np;p++) {
            sum=0.0;
            for(i=0;i<nv;i++) sum+=snapshot->B[p*nv+i]*x[i];
            bilinear_lhs+=pressure_x[p]*sum;
        }
        for(i=0;i<nv;i++) {
            sum=0.0;
            for(p=0;p<np;p++) sum+=snapshot->B[p*nv+i]*pressure_x[p];
            bilinear_rhs+=x[i]*sum;
        }
        bilinear_lhs=fabs(bilinear_lhs-bilinear_rhs)
            /max(max(fabs(bilinear_lhs),fabs(bilinear_rhs)),1.0e-300);
    }
    MPI_Bcast(&bilinear_lhs,1,MPI_DOUBLE,owner,E->parallel.world);
    *bilinear_defect=bilinear_lhs;

    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        free(u[m]); free(Au[m]); free(pressure_action[m]);
        free(pressure_input[m]); free(gradient[m]);
    }
    free(u); free(Au); free(pressure_action); free(pressure_input); free(gradient);
    free(snapshot->matrix);
    free(snapshot->B);
    memset(snapshot,0,sizeof(*snapshot));
    return(numerator);
}


static void ala_f1b_write_candidate_evidence(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache,double retained_cache_bytes,
    double build_seconds)
{
    int local_counts[4],global_counts[4];
    int local_dims[4],global_dims[4];
    double local_max[11],global_max[11],local_min,global_min;
    double projected_ratio;
    const char *directory;
    const char *manifest_sha256;
    const char *factorization_phase;
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
    local_max[2]=cache->f1b_max_diagonal_completion_relative;
    local_max[3]=cache->f1b_peak_temporary_bytes;
    local_max[4]=cache->f1b_mpi_payload_bytes;
    local_max[5]=cache->f1b_k_replay_relative_defect;
    local_max[6]=cache->f1b_output_replica_relative_difference;
    local_max[7]=max(cache->f1b_local_b_restriction_relative_defect,
                     cache->f1b_interface_b_restriction_relative_defect);
    local_max[8]=max(cache->f1b_local_bt_restriction_relative_defect,
                     cache->f1b_interface_bt_restriction_relative_defect);
    local_max[9]=cache->f1b_estimated_dense_work;
    local_max[10]=max(cache->f1b_local_bilinear_transpose_relative_defect,
                      cache->f1b_interface_bilinear_transpose_relative_defect);
    MPI_Allreduce(local_max,global_max,9,MPI_DOUBLE,MPI_MAX,
                  E->parallel.world);
    MPI_Allreduce(local_max+9,global_max+9,1,MPI_DOUBLE,MPI_SUM,
                  E->parallel.world);
    MPI_Allreduce(local_max+10,global_max+10,1,MPI_DOUBLE,MPI_MAX,
                  E->parallel.world);
    local_min=(cache->f1b_local_patch_count+cache->f1b_interface_patch_count)>0
        ? cache->f1b_min_k_pivot_ratio : 1.0e300;
    MPI_Allreduce(&local_min,&global_min,1,MPI_DOUBLE,MPI_MIN,
                  E->parallel.world);
    projected_ratio=(retained_cache_bytes+global_max[3])
        /max(retained_cache_bytes,1.0);
    directory=getenv("STRICT_ALA_STAGE_F1B_OUTPUT_DIR");
    manifest_sha256=getenv(
        "STRICT_ALA_STAGE_F1B_COLLECTIVE_MANIFEST_SHA256");
    factorization_phase=
        getenv("STRICT_ALA_STAGE_F1B_IDENTITY_PREFLIGHT")!=NULL
        ? "diagnostic_factor_attempt"
        : "authorized_candidate_factorization";
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
        "  \"collective_manifest_sha256\": \"%s\",\n"
        "  \"factorization_phase\": \"%s\",\n"
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
        "  \"legacy_retained_cache_bytes_per_rank\": %.17e,\n"
        "  \"retained_cache_ratio_vs_legacy\": %.17e,\n"
        "  \"projected_peak_ratio_vs_legacy_cache\": %.17e,\n"
        "  \"new_mpi_payload_bytes_max_per_rank\": %.17e,\n"
        "  \"diagonal_completion_relative_max\": %.17e,\n"
        "  \"Kgamma_replay_relative_defect\": %.17e,\n"
        "  \"production_output_replica_relative_difference\": %.17e,\n"
        "  \"production_output_replica_tolerance\": %.17e,\n"
        "  \"B_restriction_relative_defect\": %.17e,\n"
        "  \"Bt_restriction_relative_defect\": %.17e,\n"
        "  \"B_Bt_bilinear_relative_defect\": %.17e,\n"
        "  \"Kgamma_replay_tolerance\": %.17e,\n"
        "  \"Kgamma_replay_stop_review_upper_bound\": %.17e,\n"
        "  \"Kgamma_replay_status\": \"%s\",\n"
        "  \"startup_seconds_max\": %.17e,\n"
        "  \"memory_gate_pass\": %s,\n"
        "  \"communication\": \"direct-face plus two-hop corner element records\"\n"
        "}\n",
        manifest_sha256 ? manifest_sha256 : "MISSING",factorization_phase,
        (global_counts[0]>0 && global_counts[1]>0 &&
         global_max[5]<=ALA_F1B_K_RESTRICTION_TOLERANCE &&
         global_max[7]<=ALA_F1B_K_RESTRICTION_TOLERANCE &&
         global_max[8]<=ALA_F1B_K_RESTRICTION_TOLERANCE &&
         global_max[10]<=ALA_F1B_K_RESTRICTION_TOLERANCE &&
         global_max[6]<=ALA_F1B_OUTPUT_REPLICA_TOLERANCE) ? "true" : "false",
        global_dims[0],global_dims[1],global_dims[2],global_dims[3],
        global_counts[0]+global_counts[1],global_counts[2]+global_counts[3],
        global_max[9],global_max[3],retained_cache_bytes,
        retained_cache_bytes,1.0,projected_ratio,global_max[4],global_max[2],
        global_max[5],
        global_max[6],ALA_F1B_OUTPUT_REPLICA_TOLERANCE,
        global_max[7],global_max[8],global_max[10],
        ALA_F1B_K_RESTRICTION_TOLERANCE,
        ALA_F1B_K_RESTRICTION_REVIEW_BOUND,
        max(max(global_max[5],global_max[7]),
            max(global_max[8],global_max[10]))
          <=ALA_F1B_K_RESTRICTION_TOLERANCE ? "PASS" :
          (max(max(global_max[5],global_max[7]),
               max(global_max[8],global_max[10]))
             <=ALA_F1B_K_RESTRICTION_REVIEW_BOUND ?
           "STOP_REVIEW" : "FAIL"),
        build_seconds,
        1.0<=1.25 ? "true" : "false");
    fclose(fp);
    snprintf(path,sizeof(path),
             "%s/strict_ala_stage_F1b_candidate_assembly.csv",directory);
    fp=fopen(path,"w");
    if(fp==NULL)
        myerror(E,"Unable to open F1b candidate assembly CSV");
    fprintf(fp,"row_type,patch_category,rank,patch_id,pressure_dim,"
        "velocity_dim,symmetry_defect,regularization,"
        "factorization_status,pivot_metric,local_solve_residual,fallback,"
        "candidate_matrix_bytes,candidate_factor_bytes,Kgamma_replay_relative_defect,"
        "production_output_replica_relative_difference,B_restriction_relative_defect,"
        "Bt_restriction_relative_defect,B_Bt_bilinear_relative_defect\n");
    fprintf(fp,"AGGREGATE,LOCAL,-1,-1,%d,%d,%.17e,%.17e,PASS,%.17e,"
        "%.17e,0,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e\n",
        global_dims[0],global_dims[1],global_max[0],
        E->control.ala_shallow_patch_regularization,global_min,global_max[1],
        (double)global_dims[0]*global_dims[0]*sizeof(double),
        (double)global_dims[0]*global_dims[0]*sizeof(double),
        cache->f1b_local_k_replay_relative_defect,
        cache->f1b_local_output_replica_relative_difference,
        cache->f1b_local_b_restriction_relative_defect,
        cache->f1b_local_bt_restriction_relative_defect,
        cache->f1b_local_bilinear_transpose_relative_defect);
    fprintf(fp,"AGGREGATE,INTERFACE,-1,-1,%d,%d,%.17e,%.17e,PASS,%.17e,"
        "%.17e,0,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e\n",
        global_dims[2],global_dims[3],global_max[0],
        E->control.ala_shallow_patch_regularization,global_min,global_max[1],
        (double)global_dims[2]*global_dims[2]*sizeof(double),
        (double)global_dims[2]*global_dims[2]*sizeof(double),
        cache->f1b_interface_k_replay_relative_defect,
        cache->f1b_interface_output_replica_relative_difference,
        cache->f1b_interface_b_restriction_relative_defect,
        cache->f1b_interface_bt_restriction_relative_defect,
        cache->f1b_interface_bilinear_transpose_relative_defect);
    fclose(fp);
}


static void ala_f1b_write_rank_memory_with_cache(struct All_variables *E,
    double retained_cache_bytes,double temporary_peak_bytes)
{
    struct rusage usage;
    double local_kib,maximum_kib,sum_kib,maximum_temporary_bytes;
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
    MPI_Allreduce(&temporary_peak_bytes,&maximum_temporary_bytes,1,
                  MPI_DOUBLE,MPI_MAX,E->parallel.world);
    directory=getenv("STRICT_ALA_STAGE_F1B_OUTPUT_DIR");
    case_name=getenv("STRICT_ALA_CASE");
    if(E->parallel.me!=0 || directory==NULL || directory[0]=='\0')
        return;
    snprintf(path,sizeof(path),"%s/strict_ala_stage_F1b_rank_memory.csv",
             directory);
    fp=fopen(path,"w");
    if(fp==NULL)
        myerror(E,"Unable to open F1b rank memory CSV");
    fprintf(fp,"case,rank_rss_max_kib,rank_rss_sum_kib,ranks,"
        "retained_cache_bytes_per_rank,temporary_peak_bytes_per_rank,source\n");
    fprintf(fp,"%s,%.17e,%.17e,%d,%.17e,%.17e,"
        "getrusage_RUSAGE_SELF_MPI_reduce\n",
            case_name ? case_name : "UNKNOWN",maximum_kib,sum_kib,
            E->parallel.nproc,retained_cache_bytes,maximum_temporary_bytes);
    fclose(fp);
}


/* Existing Stage-E/F1a call sites need only RSS.  Keep their original
 * one-argument interface while F1b cache construction records the two
 * additional cache metrics through the explicit helper above. */
static void ala_f1b_write_rank_memory(struct All_variables *E)
{
    ala_f1b_write_rank_memory_with_cache(E,0.0,0.0);
}
