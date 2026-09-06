/* Strict-ALA F3: observation-only outer right-FGMRES relocalization.
 *
 * This file deliberately does not alter a solver vector.  It reconstructs
 * the frozen E2 residual POD modes, observes v_j, z_j=M^-1 v_j and the
 * unmodified production image s_j=S_prod,j z_j, and maintains bounded
 * global_pdot QR bases.  Reference actions are evaluated only after the
 * production trajectory has ended. */

#include "strict_ala_stage_F3_thresholds_generated.h"

#define ALA_F3_MAX_MODES ALA_E2_MAX_MODES
#define ALA_F3_MAX_BASIS 64
#define ALA_F3_SPACES 3
#define ALA_F3_ALL 0
#define ALA_F3_CYCLE 1

struct ala_f3_basis {
    double **q[ALA_F3_MAX_BASIS];
    int rank;
};

struct ala_f3_state {
    int enabled,finalized,lev,restart,outer_budget,selected,current_cycle;
    int precheck_completed,experiment_valid,fixed_reference_pass;
    double weights[ALA_F3_MAX_MODES],selected_energy;
    unsigned long long mode_checksum[ALA_F3_MAX_MODES];
    double **mode[ALA_F3_MAX_MODES];
    double **additivity_probe[3];
    struct ala_f3_basis basis[ALA_F3_SPACES][2];
    struct ala_e2_observer *e0,*e1;
    char output_dir[512];
    unsigned long long snapshot_hash;
    long long bytes_per_pressure_vector;
    char manifest_sha256[65],threshold_contract_sha256[65];
    double viscosity_min,viscosity_max,viscosity_contrast;
    double viscosity_lower_fraction,viscosity_upper_fraction;
    long long viscosity_count;
};

static const char *ala_f3_space_name[ALA_F3_SPACES]={"V","Z","Y"};
static FILE *ala_f3_open_root(struct All_variables *,struct ala_f3_state *,
                              const char *,const char *);

static int ala_f3_parse_csv_doubles(const char *text,double *values,int count)
{
    char *end;
    int i;
    if(text==NULL) return 0;
    for(i=0;i<count;i++) {
        values[i]=strtod(text,&end);
        if(end==text || !isfinite(values[i])) return 0;
        if(i<count-1) { if(*end!=',') return 0; text=end+1; }
        else if(*end!='\0') return 0;
    }
    return 1;
}

static int ala_f3_parse_csv_ints(const char *text,int *values,int count)
{
    char *end;
    long value;
    int i;
    if(text==NULL) return 0;
    for(i=0;i<count;i++) {
        value=strtol(text,&end,10);
        if(end==text) return 0;
        values[i]=(int)value;
        if(i<count-1) { if(*end!=',') return 0; text=end+1; }
        else if(*end!='\0') return 0;
    }
    return 1;
}

/* EVI at levmax is the integration-point field consumed by stiffness
 * assembly.  Local elements are domain-owned; order is cap, element, GP.
 * The root CSV preserves world-rank order for the post-run SHA256 envelope. */
static void ala_f3_viscosity_provenance(struct All_variables *E,
                                         struct ala_f3_state *s)
{
    int m,e,gp,vpts=vpoints[E->mesh.nsd],nproc=E->parallel.nproc;
    long long count=0,lower=0,upper=0,*counts=NULL;
    unsigned long long hash=1469598103934665603ULL,*hashes=NULL;
    double eta,cap,local_min=DBL_MAX,local_max=-DBL_MAX;
    double global_min,global_max;
    long long global_count,global_lower,global_upper;
    FILE *fp;
    for(m=1;m<=E->sphere.caps_per_proc;m++) for(e=1;e<=E->lmesh.nel;e++)
        for(gp=1;gp<=vpts;gp++) {
            float raw=E->EVI[E->mesh.levmax][m][(e-1)*vpts+gp];
            eta=(double)raw;
            if(!isfinite(eta) || eta<=0.0) myerror(E,"F3 invalid operator viscosity");
            cap=E->eco[m][e].centre[3]>0.89641 ? E->viscosity.max_value
                                                : 5.0*E->viscosity.max_value;
            local_min=min(local_min,eta); local_max=max(local_max,eta);
            if(E->viscosity.MIN && eta==(double)E->viscosity.min_value) lower++;
            if(E->viscosity.MAX && eta==cap) upper++;
            hash=ala_e2_checksum_bytes(hash,(const unsigned char *)&raw,sizeof(raw));
            count++;
        }
    MPI_Allreduce(&local_min,&global_min,1,MPI_DOUBLE,MPI_MIN,E->parallel.world);
    MPI_Allreduce(&local_max,&global_max,1,MPI_DOUBLE,MPI_MAX,E->parallel.world);
    MPI_Allreduce(&count,&global_count,1,MPI_LONG_LONG,MPI_SUM,E->parallel.world);
    MPI_Allreduce(&lower,&global_lower,1,MPI_LONG_LONG,MPI_SUM,E->parallel.world);
    MPI_Allreduce(&upper,&global_upper,1,MPI_LONG_LONG,MPI_SUM,E->parallel.world);
    if(E->parallel.me==0) {
        counts=(long long *)malloc(nproc*sizeof(long long));
        hashes=(unsigned long long *)malloc(nproc*sizeof(unsigned long long));
        if(counts==NULL || hashes==NULL) myerror(E,"F3 viscosity record allocation failed");
    }
    MPI_Gather(&count,1,MPI_LONG_LONG,counts,1,MPI_LONG_LONG,0,E->parallel.world);
    MPI_Gather(&hash,1,MPI_UNSIGNED_LONG_LONG,hashes,1,MPI_UNSIGNED_LONG_LONG,0,E->parallel.world);
    s->viscosity_min=global_min; s->viscosity_max=global_max;
    s->viscosity_contrast=global_max/global_min; s->viscosity_count=global_count;
    s->viscosity_lower_fraction=(double)global_lower/(double)global_count;
    s->viscosity_upper_fraction=(double)global_upper/(double)global_count;
    if(E->parallel.me==0) {
        fp=ala_f3_open_root(E,s,"strict_ala_stage_F3_viscosity_rank_records.csv","w");
        fprintf(fp,"rank_id,local_entry_count,local_viscosity_hash\n");
        for(m=0;m<nproc;m++) fprintf(fp,"%d,%lld,%016llx\n",m,counts[m],hashes[m]);
        fclose(fp); free(counts); free(hashes);
    }
}

static void ala_f3_path(struct ala_f3_state *s,const char *name,
                        char *path,size_t size)
{
    snprintf(path,size,"%s/%s",s->output_dir,name);
}

static FILE *ala_f3_open_root(struct All_variables *E,
    struct ala_f3_state *s,const char *name,const char *mode)
{
    char path[1024];
    FILE *fp;
    ala_f3_path(s,name,path,sizeof(path));
    fp=fopen(path,mode);
    if(fp==NULL) myerror(E,"Unable to open strict-ALA Stage-F3 output");
    return fp;
}

static void ala_f3_copy(struct All_variables *E,double **dst,double **src,
                        int lev)
{
    int m,e;
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(e=1;e<=E->lmesh.NPNO[lev];e++) dst[m][e]=src[m][e];
}

static void ala_f3_axpy(struct All_variables *E,double **y,double alpha,
                        double **x,int lev)
{
    int m,e;
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(e=1;e<=E->lmesh.NPNO[lev];e++) y[m][e]+=alpha*x[m][e];
}

static void ala_f3_scale(struct All_variables *E,double **x,double alpha,
                         int lev)
{
    int m,e;
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(e=1;e<=E->lmesh.NPNO[lev];e++) x[m][e]*=alpha;
}

static double ala_f3_norm(struct All_variables *E,double **x,int lev)
{
    double global_pdot();
    double n2=global_pdot(E,x,x,lev);
    return (isfinite(n2) && n2>=0.0) ? sqrt(n2) : -1.0;
}

static double ala_f3_difference_norm(struct All_variables *E,double **a,
                                     double **b,int lev)
{
    double **work=ala_schur_alloc_pressure(E,lev);
    double result;
    ala_f3_copy(E,work,a,lev);
    ala_f3_axpy(E,work,-1.0,b,lev);
    result=ala_f3_norm(E,work,lev);
    ala_schur_free_field(E,work);
    return result;
}

static void ala_f3_clear_basis(struct All_variables *E,
                               struct ala_f3_basis *basis)
{
    int i;
    for(i=0;i<basis->rank;i++) {
        ala_schur_free_field(E,basis->q[i]);
        basis->q[i]=NULL;
    }
    basis->rank=0;
}

/* Deterministic two-pass distributed MGS.  The first-pass metric is retained
 * for auditing, but all ranks always execute the second pass so collective
 * order cannot depend on roundoff. */
static void ala_f3_qr_insert(struct All_variables *E,
    struct ala_f3_state *s,struct ala_f3_basis *basis,double **candidate,
    const char *space,const char *scope,int iteration,double qr_floor,
    FILE *fp)
{
    double **r;
    double global_pdot();
    double pre,post1,post2,alpha,max_relative=0.0,value,remaining;
    int i,accepted=0;
    const double reject_relative=
        STRICT_ALA_STAGE_F3_QR_RANK_REJECTION_RELATIVE_THRESHOLD;
    r=ala_schur_alloc_pressure(E,s->lev);
    ala_f3_copy(E,r,candidate,s->lev);
    pre=ala_f3_norm(E,r,s->lev);
    for(i=0;i<basis->rank;i++) {
        alpha=global_pdot(E,basis->q[i],r,s->lev);
        ala_f3_axpy(E,r,-alpha,basis->q[i],s->lev);
    }
    post1=ala_f3_norm(E,r,s->lev);
    for(i=0;i<basis->rank;i++) {
        value=fabs(global_pdot(E,basis->q[i],r,s->lev))
            /max(post1,qr_floor);
        if(value>max_relative) max_relative=value;
    }
    for(i=0;i<basis->rank;i++) {
        alpha=global_pdot(E,basis->q[i],r,s->lev);
        ala_f3_axpy(E,r,-alpha,basis->q[i],s->lev);
    }
    post2=ala_f3_norm(E,r,s->lev);
    remaining=(pre>0.0 && isfinite(post2)) ? post2/pre : 0.0;
    if(pre>0.0 && post2>reject_relative*pre && basis->rank<ALA_F3_MAX_BASIS &&
       isfinite(post2)) {
        ala_f3_scale(E,r,1.0/post2,s->lev);
        basis->q[basis->rank++]=r;
        accepted=1;
    }
    else ala_schur_free_field(E,r);
    if(E->parallel.me==0)
        fprintf(fp,"%d,%s,%s,%d,%d,%s,%.17e,%.17e,%.17e,%.17e,%d,%.17e,%.17e,%.17e,%d,%d\n",
            iteration,space,scope,accepted,accepted ? 0 : 1,
            pre<=0.0 ? "NONFINITE_OR_ZERO" :
            (basis->rank>=ALA_F3_MAX_BASIS && !accepted ? "BASIS_LIMIT" :
             (!accepted ? "RELATIVE_RANK_REJECTION" : "NONE")),
            pre,post1,post2,max_relative,
            max_relative>STRICT_ALA_STAGE_F3_QR_FIRST_PASS_ORTHOGONALITY_THRESHOLD,
            remaining,STRICT_ALA_STAGE_F3_QR_FIRST_PASS_ORTHOGONALITY_THRESHOLD,
            reject_relative,1,basis->rank);
}

static double ala_f3_projection_error(struct All_variables *E,double **mode,
                                      struct ala_f3_basis *basis,int lev)
{
    double global_pdot();
    double sum=0.0,dot;
    int i;
    for(i=0;i<basis->rank;i++) {
        dot=global_pdot(E,mode,basis->q[i],lev);
        sum+=dot*dot;
    }
    if(sum>1.0 && sum<1.0+1.0e-10) sum=1.0;
    return sqrt(max(1.0-sum,0.0));
}

static unsigned long long ala_f3_mode_checksum(struct All_variables *E,
                                               double **mode,int lev)
{
    unsigned long long local=1469598103934665603ULL,global;
    int m;
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        local=ala_e2_checksum_bytes(local,
            (const unsigned char *)&mode[m][1],
            E->lmesh.NPNO[lev]*sizeof(double));
    MPI_Allreduce(&local,&global,1,MPI_UNSIGNED_LONG_LONG,MPI_BXOR,
                  E->parallel.world);
    return global;
}

static void ala_f3_create_outputs(struct All_variables *E,
                                  struct ala_f3_state *s)
{
    FILE *fp;
    if(E->parallel.me!=0) return;
    fp=ala_f3_open_root(E,s,"strict_ala_stage_F3_iteration_pod_projection.csv","w");
    fprintf(fp,"iteration,restart_cycle,mode_id,mode_energy_weight,q_dot_v,q_dot_z,q_dot_s,v_norm,z_norm,s_norm,requested_inner_tolerance,production_action_sampled_before_arnoldi\n"); fclose(fp);
    fp=ala_f3_open_root(E,s,"strict_ala_stage_F3_cumulative_subspace.csv","w");
    fprintf(fp,"iteration,restart_cycle,row_type,mode_id,mode_energy_weight,E_V_all,E_Z_all,E_Y_all,E_V_cycle,E_Z_cycle,E_Y_cycle,E_V_weighted_all,E_Z_weighted_all,E_Y_weighted_all,E_V_weighted_cycle,E_Z_weighted_cycle,E_Y_weighted_cycle,descriptive_geometry_only_VZ\n"); fclose(fp);
    fp=ala_f3_open_root(E,s,"strict_ala_stage_F3_basis_rank.csv","w");
    fprintf(fp,"iteration,space,scope,accepted,rejected,rejection_reason,pre_norm,post_first_pass_norm,post_second_pass_norm,first_pass_orthogonality_defect,first_pass_defect_exceeds_threshold,relative_remaining_norm,qr_first_pass_orthogonality_threshold_loaded,qr_rank_rejection_relative_threshold_loaded,second_pass_performed,resulting_rank\n"); fclose(fp);
    fp=ala_f3_open_root(E,s,"strict_ala_stage_F3_mode_coefficients.csv","w");
    fprintf(fp,"iteration,mode_id,mode_energy_weight,coefficient,absolute_modal_amplitude,residual_global_pdot_norm,mode_status\n"); fclose(fp);
    fp=ala_f3_open_root(E,s,"strict_ala_stage_F3_pod_explanation.csv","w");
    fprintf(fp,"iteration,residual_global_pdot_norm,POD_projected_norm,f_Q,continuity_relative,momentum_relative,krylov_recursive_absolute,krylov_explicit_absolute,cumulative_inner_solves,cumulative_MG_cycles,cumulative_schur_actions,elapsed_seconds\n"); fclose(fp);
    fp=ala_f3_open_root(E,s,"strict_ala_stage_F3_restart_retention.csv","w");
    fprintf(fp,"restart_boundary_iteration,restart_cycle,mode_id,pre_restart_Y_reachability,post_restart_iteration,post_restart_cycle_local_Y_reachability,recovery_metric_available,recovery_iterations\n"); fclose(fp);
    fp=ala_f3_open_root(E,s,"strict_ala_stage_F3_reference_actions.csv","w");
    fprintf(fp,"operator,action_id,input_id,repeat_id,global_rhs_norm2,global_residual_norm2,rhs_rms,residual_rms,zero_rhs,relative_residual,residual_status,requested_relative_tolerance,achieved_MG_cycles,max_MG_cycles,converged,fallback,finite,residual_gate_pass,output_norm\n"); fclose(fp);
    fp=ala_f3_open_root(E,s,"strict_ala_stage_F3_projected_matrices.csv","w");
    fprintf(fp,"matrix,row_mode,column_mode,value\n"); fclose(fp);
    fp=ala_f3_open_root(E,s,"strict_ala_stage_F3_pod_reconstruction_coefficients.csv","w");
    fprintf(fp,"case,iteration,pod_mode,coefficient\n"); fclose(fp);
}

static void *strict_ala_stage_f3_begin(struct All_variables *E,int lev,
                                      int restart,int outer_budget)
{
    struct ala_f3_state *s;
    struct ala_e2_observer template_observer,*cases[2*ALA_E2_LATE_COUNT];
    int indices[2*ALA_E2_LATE_COUNT];
    double values[2*ALA_E2_LATE_COUNT];
    double vectors[4*ALA_E2_LATE_COUNT*ALA_E2_LATE_COUNT];
    double norms[2*ALA_E2_LATE_COUNT],cov[4*ALA_E2_LATE_COUNT*ALA_E2_LATE_COUNT];
    double total=0.0,orth_max=0.0,value;
    double expected_weights[ALA_F3_MAX_MODES],expected_selected_energy;
    int expected_ids[ALA_F3_MAX_MODES];
    const char *required=getenv("STRICT_ALA_STAGE_F3_REQUIRED");
    const char *mode_count_env=getenv("STRICT_ALA_STAGE_F3_AUTHORITATIVE_MODE_COUNT");
    const char *mode_ids_env=getenv("STRICT_ALA_STAGE_F3_AUTHORITATIVE_MODE_IDS");
    const char *weights_env=getenv("STRICT_ALA_STAGE_F3_AUTHORITATIVE_WEIGHTS");
    const char *energy_env=getenv("STRICT_ALA_STAGE_F3_AUTHORITATIVE_SELECTED_ENERGY");
    const char *manifest_env=getenv("STRICT_ALA_STAGE_F3_MANIFEST_SHA256");
    const char *threshold_contract_env=getenv("STRICT_ALA_STAGE_F3_THRESHOLD_CONTRACT_SHA256");
    const char *snapshot_dir=getenv("STRICT_ALA_STAGE_F3_E2_SNAPSHOT_DIR");
    const char *output_dir=getenv("STRICT_ALA_STAGE_F3_OUTPUT_DIR");
    FILE *fp,*coefficient_fp=NULL;
    int i,j,n=2*ALA_E2_LATE_COUNT;
    unsigned long long checksum;
    if(required==NULL) return NULL;
    if(snapshot_dir==NULL || output_dir==NULL || manifest_env==NULL ||
       threshold_contract_env==NULL)
        myerror(E,"Incomplete strict-ALA Stage-F3 environment contract");
    if(E->parallel.nproc!=384 || E->parallel.nprocx!=4 ||
       E->parallel.nprocy!=4 || E->parallel.nprocz!=2)
        myerror(E,"F3_E2_REPLAY_LAYOUT_MISMATCH");
    if(strcmp(E->control.ala_outer_solver,"fgmres")!=0 || restart!=50 ||
       outer_budget!=60 || E->control.ala_inner_accuracy_max!=1.0e-2)
        myerror(E,"F3_FROZEN_OUTER_CONTRACT_MISMATCH");
    s=(struct ala_f3_state *)calloc(1,sizeof(*s));
    if(s==NULL) myerror(E,"Unable to allocate Stage-F3 state");
    s->enabled=1; s->lev=lev; s->restart=restart;
    s->outer_budget=outer_budget; s->current_cycle=1;
    snprintf(s->output_dir,sizeof(s->output_dir),"%s",output_dir);
    snprintf(s->manifest_sha256,sizeof(s->manifest_sha256),"%.64s",manifest_env);
    snprintf(s->threshold_contract_sha256,sizeof(s->threshold_contract_sha256),
             "%.64s",threshold_contract_env);
    memset(&template_observer,0,sizeof(template_observer));
    template_observer.lev=lev;
    snprintf(template_observer.shared_dir,sizeof(template_observer.shared_dir),
             "%s",snapshot_dir);
    snprintf(template_observer.output_dir,sizeof(template_observer.output_dir),
             "%s",output_dir);
    s->e0=ala_e2_read_case(E,&template_observer,"E0");
    s->e1=ala_e2_read_case(E,&template_observer,"E1");
    for(i=0;i<ALA_E2_LATE_COUNT;i++) {
        cases[i]=s->e0; indices[i]=i+ALA_E2_LATE_FIRST;
        cases[i+ALA_E2_LATE_COUNT]=s->e1;
        indices[i+ALA_E2_LATE_COUNT]=i+ALA_E2_LATE_FIRST;
    }
    ala_e2_covariance(E,cases,indices,n,cov);
    if(!ala_e2_pod(cov,n,1,values,vectors,norms))
        myerror(E,"Stage-F3 deterministic E2 POD reconstruction failed");
    /* The authoritative E2 mode count is supplied by the immutable pre-run
     * mode contract.  Never reselect modes from the reconstructed spectrum. */
    if(mode_count_env==NULL || sscanf(mode_count_env,"%d",&s->selected)!=1)
        myerror(E,"F3_AUTHORITATIVE_E2_MODE_COUNT_MISSING");
    if(s->selected<1 || s->selected>ALA_F3_MAX_MODES)
        myerror(E,"Stage-F3 authoritative POD selection mismatch");
    if(!ala_f3_parse_csv_ints(mode_ids_env,expected_ids,s->selected) ||
       !ala_f3_parse_csv_doubles(weights_env,expected_weights,s->selected) ||
       energy_env==NULL || sscanf(energy_env,"%lf",&expected_selected_energy)!=1)
        myerror(E,"F3_AUTHORITATIVE_E2_MODE_CONTRACT_PARSE_FAILED");
    for(i=0;i<n;i++) total+=values[i];
    s->selected_energy=0.0;
    s->snapshot_hash=1469598103934665603ULL;
    for(i=0;i<s->selected;i++) {
        if(expected_ids[i]!=i+1) myerror(E,"F3_AUTHORITATIVE_E2_MODE_ID_ORDER_MISMATCH");
        s->mode[i]=ala_e2_build_mode(E,cases,indices,n,vectors,norms,i);
        s->weights[i]=values[i]/max(total,1.0e-300);
        if(fabs(s->weights[i]-expected_weights[i])>
           STRICT_ALA_STAGE_F3_POD_WEIGHT_SUM_ABSOLUTE_TOLERANCE)
            myerror(E,"F3_AUTHORITATIVE_E2_MODE_WEIGHT_MISMATCH");
        s->weights[i]=expected_weights[i];
        s->selected_energy+=s->weights[i];
        checksum=ala_f3_mode_checksum(E,s->mode[i],lev);
        s->mode_checksum[i]=checksum;
        s->snapshot_hash^=checksum;
        if(i<3) s->additivity_probe[i]=s->mode[i];
    }
    /* q3 is retained only as a frozen linearity probe when it is outside the
     * two-mode authoritative Q_E2.  It never enters F3 POD science. */
    for(i=s->selected;i<3;i++)
        s->additivity_probe[i]=ala_e2_build_mode(E,cases,indices,n,vectors,norms,i);
    if(fabs(s->selected_energy-expected_selected_energy)>
       STRICT_ALA_STAGE_F3_POD_WEIGHT_SUM_ABSOLUTE_TOLERANCE)
        myerror(E,"F3_AUTHORITATIVE_E2_SELECTED_ENERGY_MISMATCH");
    for(i=0;i<s->selected;i++) for(j=0;j<s->selected;j++) {
        value=global_pdot(E,s->mode[i],s->mode[j],lev)-(i==j ? 1.0 : 0.0);
        if(fabs(value)>orth_max) orth_max=fabs(value);
    }
    ala_f3_create_outputs(E,s);
    ala_f3_viscosity_provenance(E,s);
    if(E->parallel.me==0)
        coefficient_fp=ala_f3_open_root(E,s,
            "strict_ala_stage_F3_pod_reconstruction_coefficients.csv","a");
    for(i=0;i<s->selected;i++) {
        for(j=0;j<s->e0->count;j++) {
            value=global_pdot(E,s->mode[i],s->e0->snapshot[j],lev);
            if(E->parallel.me==0)
                fprintf(coefficient_fp,"E0,%d,%d,%.17e\n",
                    s->e0->iterations[j],i+1,value);
        }
        for(j=0;j<s->e1->count;j++) {
            value=global_pdot(E,s->mode[i],s->e1->snapshot[j],lev);
            if(E->parallel.me==0)
                fprintf(coefficient_fp,"E1,%d,%d,%.17e\n",
                    s->e1->iterations[j],i+1,value);
        }
    }
    if(coefficient_fp) fclose(coefficient_fp);
    s->bytes_per_pressure_vector=(long long)E->sphere.caps_per_proc
        *(long long)E->lmesh.NPNO[lev]*(long long)sizeof(double);
    if(E->parallel.me==0) {
        fp=ala_f3_open_root(E,s,"strict_ala_stage_F3_pod_reconstruction.json","w");
        fprintf(fp,"{\n  \"schema_id\":\"strict-ala-stage-F3-pod-reconstruction\",\n  \"schema_version\":\"1\",\n  \"stage_id\":\"STRICT_STAGE_F3_OUTER_SCHUR_KRYLOV_RELOCALIZATION\",\n  \"manifest_sha256\":\"%s\",\n  \"method\":\"deterministic_E2_joint_direction_case_balanced_snapshot_reconstruction\",\n  \"pressure_inner_product\":\"global_pdot_euclidean_pressure\",\n  \"selected_modes\":%d,\n  \"mode_ids\":[",
            s->manifest_sha256,s->selected);
        for(i=0;i<s->selected;i++) fprintf(fp,"%s%d",i ? "," : "",i+1);
        fprintf(fp,"],\n  \"pod_weight_sum\":%.17e,\n  \"orthogonality_max_abs_defect\":%.17e,\n  \"deterministic_sign_rule\":\"match_original_E2_joint_coefficient_trajectory\",\n  \"snapshot_aggregate_checksum\":\"%016llx\",\n  \"mode_checksums\":[",
            s->selected_energy,orth_max,s->snapshot_hash);
        for(i=0;i<s->selected;i++) {
            fprintf(fp,"%s{\"mode_id\":%d,\"checksum\":\"%016llx\",\"weight\":%.17e}",
                i ? "," : "",i+1,s->mode_checksum[i],s->weights[i]);
        }
        fprintf(fp,"],\n  \"authoritative_mode_set_reused\":true,\n  \"pod_reselection_performed\":false,\n  \"pod_retraining_performed\":false,\n  \"normalization_validation_pass\":%s,\n  \"coefficient_trajectory_validation\":\"DEFERRED_TO_FAIL_CLOSED_ANALYZER\",\n  \"Q_E2_validation_complete\":%s,\n  \"valid\":%s,\n  \"complete\":true\n}\n",
            orth_max<=STRICT_ALA_STAGE_F3_QR_FIRST_PASS_ORTHOGONALITY_THRESHOLD ? "true" : "false",
            orth_max<=STRICT_ALA_STAGE_F3_QR_FIRST_PASS_ORTHOGONALITY_THRESHOLD ? "true" : "false",
            orth_max<=STRICT_ALA_STAGE_F3_QR_FIRST_PASS_ORTHOGONALITY_THRESHOLD ? "true" : "false");
        fclose(fp);
        fprintf(E->fp,"STRICT_ALA_STAGE_F3_BEGIN orientation=RIGHT_FGMRES "
            "modes=%d pod_weight_sum=%.17e restart=%d outer_budget=%d\n",
            s->selected,s->selected_energy,restart,outer_budget);
        fflush(E->fp);
    }
    return s;
}

static void ala_f3_write_subspace(struct All_variables *E,
    struct ala_f3_state *s,int iteration,int cycle)
{
    double errors[ALA_F3_SPACES][2][ALA_F3_MAX_MODES];
    double weighted[ALA_F3_SPACES][2],sum;
    FILE *fp=NULL;
    int x,scope,i;
    for(x=0;x<ALA_F3_SPACES;x++) for(scope=0;scope<2;scope++) {
        sum=0.0;
        for(i=0;i<s->selected;i++) {
            errors[x][scope][i]=ala_f3_projection_error(
                E,s->mode[i],&s->basis[x][scope],s->lev);
            sum+=s->weights[i]*errors[x][scope][i]*errors[x][scope][i];
        }
        weighted[x][scope]=sqrt(sum/max(s->selected_energy,1.0e-300));
    }
    if(E->parallel.me==0) {
        fp=ala_f3_open_root(E,s,"strict_ala_stage_F3_cumulative_subspace.csv","a");
        for(i=0;i<s->selected;i++)
            fprintf(fp,"%d,%d,MODE,%d,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,,,,,,,1\n",
                iteration,cycle,i+1,s->weights[i],errors[0][0][i],
                errors[1][0][i],errors[2][0][i],errors[0][1][i],
                errors[1][1][i],errors[2][1][i]);
        fprintf(fp,"%d,%d,AGGREGATE,0,%.17e,,,,,,,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,1\n",
            iteration,cycle,s->selected_energy,weighted[0][0],weighted[1][0],
            weighted[2][0],weighted[0][1],weighted[1][1],weighted[2][1]);
        fclose(fp);
    }
}

static void strict_ala_stage_f3_direction(struct All_variables *E,void *state,
    double **v,double **z,double **image,int iteration,int cycle,
    double requested_inner_tolerance)
{
    struct ala_f3_state *s=(struct ala_f3_state *)state;
    double **vector[ALA_F3_SPACES];
    double norms[ALA_F3_SPACES],dot[ALA_F3_SPACES],qr_floor;
    FILE *rank_fp=NULL,*projection_fp=NULL,*restart_fp=NULL;
    double pre_restart_error[ALA_F3_MAX_MODES];
    int x,i;
    if(s==NULL || !s->enabled) return;
    if(cycle!=s->current_cycle) {
        for(x=0;x<ALA_F3_SPACES;x++)
            ala_f3_clear_basis(E,&s->basis[x][ALA_F3_CYCLE]);
        s->current_cycle=cycle;
    }
    vector[0]=v; vector[1]=z; vector[2]=image;
    for(x=0;x<ALA_F3_SPACES;x++) norms[x]=ala_f3_norm(E,vector[x],s->lev);
    qr_floor=64.0*DBL_EPSILON*max(max(norms[0],norms[1]),max(norms[2],1.0));
    if(E->parallel.me==0)
        rank_fp=ala_f3_open_root(E,s,"strict_ala_stage_F3_basis_rank.csv","a");
    for(x=0;x<ALA_F3_SPACES;x++) {
        ala_f3_qr_insert(E,s,&s->basis[x][ALA_F3_ALL],vector[x],
            ala_f3_space_name[x],"ALL",iteration,qr_floor,rank_fp);
        ala_f3_qr_insert(E,s,&s->basis[x][ALA_F3_CYCLE],vector[x],
            ala_f3_space_name[x],"CYCLE",iteration,qr_floor,rank_fp);
    }
    if(rank_fp) fclose(rank_fp);
    if(E->parallel.me==0)
        projection_fp=ala_f3_open_root(E,s,
            "strict_ala_stage_F3_iteration_pod_projection.csv","a");
    for(i=0;i<s->selected;i++) {
        for(x=0;x<ALA_F3_SPACES;x++)
            dot[x]=global_pdot(E,s->mode[i],vector[x],s->lev);
        if(E->parallel.me==0)
            fprintf(projection_fp,"%d,%d,%d,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,1\n",
                iteration,cycle,i+1,s->weights[i],dot[0],dot[1],dot[2],
                norms[0],norms[1],norms[2],requested_inner_tolerance);
    }
    if(projection_fp) fclose(projection_fp);
    ala_f3_write_subspace(E,s,iteration,cycle);
    if(iteration==s->restart) {
        for(i=0;i<s->selected;i++)
            pre_restart_error[i]=ala_f3_projection_error(E,s->mode[i],
                &s->basis[2][ALA_F3_ALL],s->lev);
        if(E->parallel.me==0) {
            restart_fp=ala_f3_open_root(E,s,
                "strict_ala_stage_F3_restart_retention.csv","a");
            for(i=0;i<s->selected;i++)
            fprintf(restart_fp,"%d,%d,%d,%.17e,,,%d,\n",iteration,cycle,
                i+1,pre_restart_error[i],0);
            fclose(restart_fp);
        }
    }
}

static void strict_ala_stage_f3_residual(struct All_variables *E,void *state,
    double **residual,int iteration,double continuity_relative,
    double momentum_relative,double recursive_absolute,double explicit_absolute)
{
    struct ala_f3_state *s=(struct ala_f3_state *)state;
    double global_pdot();
    double norm,coefficient,projected2=0.0;
    FILE *mode_fp=NULL,*explain_fp=NULL,*restart_fp=NULL;
    double cycle_error[ALA_F3_MAX_MODES];
    int i;
    if(s==NULL || !s->enabled) return;
    norm=ala_f3_norm(E,residual,s->lev);
    if(norm<0.0) myerror(E,"Stage-F3 nonfinite pressure residual");
    if(E->parallel.me==0)
        mode_fp=ala_f3_open_root(E,s,"strict_ala_stage_F3_mode_coefficients.csv","a");
    for(i=0;i<s->selected;i++) {
        coefficient=global_pdot(E,s->mode[i],residual,s->lev);
        projected2+=coefficient*coefficient;
        if(E->parallel.me==0)
            fprintf(mode_fp,"%d,%d,%.17e,%.17e,%.17e,%.17e,PENDING_ANALYSIS\n",
                iteration,i+1,s->weights[i],coefficient,fabs(coefficient),norm);
    }
    if(mode_fp) fclose(mode_fp);
    if(E->parallel.me==0) {
        explain_fp=ala_f3_open_root(E,s,"strict_ala_stage_F3_pod_explanation.csv","a");
        fprintf(explain_fp,"%d,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%lld,%lld,%lld,%.17e\n",
            iteration,norm,sqrt(max(projected2,0.0)),
            projected2/max(norm*norm,1.0e-300),continuity_relative,
            momentum_relative,recursive_absolute,explicit_absolute,
            E->control.ala_stage_abc_inner_call_count,
            E->control.ala_stage_abc_inner_cycle_count,
            E->control.ala_stage_abc_schur_action_count,
            CPU_time0()-E->control.ala_stage_abc_fgmres_start);
        fclose(explain_fp);
    }
    if(iteration>s->restart) {
        for(i=0;i<s->selected;i++)
            cycle_error[i]=ala_f3_projection_error(E,s->mode[i],
                &s->basis[2][ALA_F3_CYCLE],s->lev);
        if(E->parallel.me==0) {
            restart_fp=ala_f3_open_root(E,s,"strict_ala_stage_F3_restart_retention.csv","a");
            for(i=0;i<s->selected;i++)
            fprintf(restart_fp,"%d,%d,%d,,%d,%.17e,0,\n",s->restart,
                s->current_cycle,i+1,iteration,cycle_error[i]);
            fclose(restart_fp);
        }
    }
}

static int ala_f3_accuracy_solve(struct All_variables *E, double **x,
    double **rhs, double relative_target, int lev)
{
    double rnorm=sqrt(max(global_vdot(E,rhs,rhs,lev),0.0));
    setenv("STRICT_ALA_STAGE_F3_REFERENCE_ACTION","1",1);
    solve_del2_u_bounded(E,x,rhs,relative_target*max(rnorm,1.0e-300),
        lev,STRICT_ALA_STAGE_F3_S_REF_MAX_MG_CYCLES,100);
    unsetenv("STRICT_ALA_STAGE_F3_REFERENCE_ACTION");
    return 1;
}

static int ala_f3_reference_action(struct All_variables *,
    struct ala_f3_state *,double **,double **,
    struct ala_pressure_preconditioner_cache *,int,const char *,const char *,
    int,double *,long long *);

/* Short fixed-action preflight.  First isolate M^-1, then certify the complete
 * S_ref M^-1 action before spending a full production trajectory. */
static int strict_ala_stage_f3_precheck(struct All_variables *E,
    struct ala_f3_state *s, struct ala_pressure_preconditioner_cache *cache)
{
    double **q1,**q2,**qp,**z1,**z2,**zp,**zerr,**work,**a1,**a2,**ap;
    double **g1,**g2,**gp,**b1w,**b2w,**bpw,**ob1,**ob2,**obp;
    double **xs1,**xs2,**xsp,**xse;
    double **rwork,**u1,**u2,**up,**uerr,**uc1,**uc2,**ucp,**vw;
    double **repeat1,**repeat2,**repeatp,**scaled1,**scaled2,**scaledp,**stest;
    double **d1,**d2,**dp,**derr;
    double n1,n2,m_error,a1n,a2n,a_error,den,achieved;
    double bt_error,k_error,b_error,gn1,gn2,un1,un2,dn1,dn2,dn12;
    double rx1,rx2,rxp,ex_error,levels[3],sweep_R[3],sweep_E[3];
    double br1,br2,brp,denx,repeat_error,scaling_error;
    int li;
    long long cycles;
    FILE *fp;
    int m,e,m_pass,a_pass=0,actions_valid=1,pass;
    int decomposition_valid=0;
    int saved_production_logging;
    bt_error=k_error=b_error=-1.0;
    rx1=rx2=rxp=ex_error=-1.0;
    repeat_error=scaling_error=-1.0;
    q1=s->mode[0];
    q2=(s->selected>1) ? s->mode[1] : s->additivity_probe[1];
    qp=ala_schur_alloc_pressure(E,s->lev);
    z1=ala_schur_alloc_pressure(E,s->lev);
    z2=ala_schur_alloc_pressure(E,s->lev);
    zp=ala_schur_alloc_pressure(E,s->lev);
    zerr=ala_schur_alloc_pressure(E,s->lev);
    work=ala_schur_alloc_pressure(E,s->lev);
    a1=ala_schur_alloc_pressure(E,s->lev);
    a2=ala_schur_alloc_pressure(E,s->lev);
    ap=ala_schur_alloc_pressure(E,s->lev);
    g1=ala_schur_alloc_velocity(E,s->lev);
    g2=ala_schur_alloc_velocity(E,s->lev);
    gp=ala_schur_alloc_velocity(E,s->lev);
    b1w=ala_schur_alloc_velocity(E,s->lev);
    b2w=ala_schur_alloc_velocity(E,s->lev);
    bpw=ala_schur_alloc_velocity(E,s->lev);
    ob1=ala_schur_alloc_velocity(E,s->lev);
    ob2=ala_schur_alloc_velocity(E,s->lev);
    obp=ala_schur_alloc_velocity(E,s->lev);
    xs1=ala_schur_alloc_velocity(E,s->lev);
    xs2=ala_schur_alloc_velocity(E,s->lev);
    xsp=ala_schur_alloc_velocity(E,s->lev);
    xse=ala_schur_alloc_velocity(E,s->lev);
    rwork=ala_schur_alloc_velocity(E,s->lev);
    u1=ala_schur_alloc_velocity(E,s->lev);
    u2=ala_schur_alloc_velocity(E,s->lev);
    up=ala_schur_alloc_velocity(E,s->lev);
    uerr=ala_schur_alloc_velocity(E,s->lev);
    uc1=ala_schur_alloc_velocity(E,s->lev);
    uc2=ala_schur_alloc_velocity(E,s->lev);
    ucp=ala_schur_alloc_velocity(E,s->lev);
    vw=ala_schur_alloc_velocity(E,s->lev);
    d1=ala_schur_alloc_pressure(E,s->lev);
    d2=ala_schur_alloc_pressure(E,s->lev);
    dp=ala_schur_alloc_pressure(E,s->lev);
    derr=ala_schur_alloc_pressure(E,s->lev);
    repeat1=ala_schur_alloc_pressure(E,s->lev);
    repeat2=ala_schur_alloc_pressure(E,s->lev);
    repeatp=ala_schur_alloc_pressure(E,s->lev);
    scaled1=ala_schur_alloc_pressure(E,s->lev);
    scaled2=ala_schur_alloc_pressure(E,s->lev);
    scaledp=ala_schur_alloc_pressure(E,s->lev);
    stest=ala_schur_alloc_pressure(E,s->lev);
    ala_f3_copy(E,qp,q1,s->lev);
    ala_f3_axpy(E,qp,1.0,q2,s->lev);
    apply_ala_pressure_preconditioner(E,q1,z1,work,s->lev,-3,cache);
    apply_ala_pressure_preconditioner(E,q2,z2,work,s->lev,-3,cache);
    apply_ala_pressure_preconditioner(E,qp,zp,work,s->lev,-3,cache);
    n1=ala_f3_norm(E,z1,s->lev); n2=ala_f3_norm(E,z2,s->lev);
    ala_f3_copy(E,zerr,zp,s->lev);
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(e=1;e<=E->lmesh.NPNO[s->lev];e++)
            zerr[m][e]-=z1[m][e]+z2[m][e];
    m_error=ala_f3_norm(E,zerr,s->lev);
    den=max(n1+n2,1.0e-300);
    m_error/=den;
    m_pass=isfinite(m_error) &&
        m_error<=STRICT_ALA_STAGE_F3_ADDITIVITY_RELATIVE_TOLERANCE;
    if(m_pass) {
        assemble_grad_rho_p(E,z1,g1,s->lev);
        assemble_grad_rho_p(E,z2,g2,s->lev);
        assemble_grad_rho_p(E,zp,gp,s->lev);
        /* Keep b12 immutable.  gp is reserved for the B^T additivity
         * defect; the independent bpw/obp pair carries B^T(q1+q2). */
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=0;e<E->lmesh.NEQ[s->lev];e++) {
                bpw[m][e]=gp[m][e];
                obp[m][e]=gp[m][e];
            }
        gn1=sqrt(max(global_vdot(E,g1,g1,s->lev),0.0));
        gn2=sqrt(max(global_vdot(E,g2,g2,s->lev),0.0));
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=0;e<E->lmesh.NEQ[s->lev];e++)
                gp[m][e]-=g1[m][e]+g2[m][e];
        bt_error=sqrt(max(global_vdot(E,gp,gp,s->lev),0.0)) /
            max(gn1+gn2,1.0e-300);
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=0;e<E->lmesh.NEQ[s->lev];e++) {
                b1w[m][e]=g1[m][e];
                b2w[m][e]=g2[m][e];
                /* bpw remains the original b12, not b1+b2. */
                ob1[m][e]=b1w[m][e]; ob2[m][e]=b2w[m][e]; obp[m][e]=bpw[m][e];
            }
        setenv("STRICT_ALA_STAGE_F3_REFERENCE_ACTION","1",1);
        saved_production_logging=E->control.ala_stage_abc_production_logging;
        E->control.ala_stage_abc_production_logging=0;
        solve_del2_u_fixed_cycles(E,u1,g1,s->lev,
            STRICT_ALA_STAGE_F3_S_REF_FIXED_MG_CYCLES,100);
        solve_del2_u_fixed_cycles(E,u2,g2,s->lev,
            STRICT_ALA_STAGE_F3_S_REF_FIXED_MG_CYCLES,100);
        solve_del2_u_fixed_cycles(E,up,bpw,s->lev,
            STRICT_ALA_STAGE_F3_S_REF_FIXED_MG_CYCLES,100);
        E->control.ala_stage_abc_production_logging=saved_production_logging;
        unsetenv("STRICT_ALA_STAGE_F3_REFERENCE_ACTION");
        un1=sqrt(max(global_vdot(E,u1,u1,s->lev),0.0));
        un2=sqrt(max(global_vdot(E,u2,u2,s->lev),0.0));
        /* Preserve u12 for the clean B(u12) action; the old diagnostic
         * overwrote it with the additivity defect before applying B. */
        ala_f3_copy(E,uerr,up,s->lev);
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=0;e<E->lmesh.NEQ[s->lev];e++)
                uerr[m][e]-=u1[m][e]+u2[m][e];
        ex_error=sqrt(max(global_vdot(E,uerr,uerr,s->lev),0.0)) /
            max(sqrt(max(global_vdot(E,up,up,s->lev),0.0)),1.0e-300);
        k_error=ex_error;
        (void)ala_schur_velocity_residual(E,u1,b1w,vw,s->lev);
        rx1=sqrt(max(global_vdot(E,vw,vw,s->lev),0.0)) /
            max(gn1,1.0e-300);
        (void)ala_schur_velocity_residual(E,u2,b2w,vw,s->lev);
        rx2=sqrt(max(global_vdot(E,vw,vw,s->lev),0.0)) /
            max(gn2,1.0e-300);
        (void)ala_schur_velocity_residual(E,up,obp,vw,s->lev);
        rxp=sqrt(max(global_vdot(E,vw,vw,s->lev),0.0)) /
            max(sqrt(max(global_vdot(E,gp,gp,s->lev),0.0)),1.0e-300);
        levels[0]=1.0e-6; levels[1]=1.0e-8; levels[2]=1.0e-10;
        for(li=0;li<3;li++) {
            br1=sqrt(max(global_vdot(E,ob1,ob1,s->lev),0.0));
            br2=sqrt(max(global_vdot(E,ob2,ob2,s->lev),0.0));
            brp=sqrt(max(global_vdot(E,obp,obp,s->lev),0.0));
            for(m=1;m<=E->sphere.caps_per_proc;m++)
                for(e=0;e<E->lmesh.NEQ[s->lev];e++) {
                    b1w[m][e]=ob1[m][e]; b2w[m][e]=ob2[m][e];
                    bpw[m][e]=obp[m][e];
                }
            ala_f3_accuracy_solve(E,xs1,b1w,levels[li],s->lev);
            ala_f3_accuracy_solve(E,xs2,b2w,levels[li],s->lev);
            ala_f3_accuracy_solve(E,xsp,bpw,levels[li],s->lev);
            (void)ala_schur_velocity_residual(E,xs1,ob1,rwork,s->lev);
            sweep_R[li]=sqrt(max(global_vdot(E,rwork,rwork,s->lev),0.0))/max(br1,1.0e-300);
            (void)ala_schur_velocity_residual(E,xs2,ob2,rwork,s->lev);
            sweep_R[li]=max(sweep_R[li],sqrt(max(global_vdot(E,rwork,rwork,s->lev),0.0))/max(br2,1.0e-300));
            (void)ala_schur_velocity_residual(E,xsp,obp,rwork,s->lev);
            sweep_R[li]=max(sweep_R[li],sqrt(max(global_vdot(E,rwork,rwork,s->lev),0.0))/max(brp,1.0e-300));
            for(m=1;m<=E->sphere.caps_per_proc;m++)
                for(e=0;e<E->lmesh.NEQ[s->lev];e++)
                    xse[m][e]=xsp[m][e]-xs1[m][e]-xs2[m][e];
            denx=sqrt(max(global_vdot(E,xsp,xsp,s->lev),0.0));
            sweep_E[li]=sqrt(max(global_vdot(E,xse,xse,s->lev),0.0))/max(denx,1.0e-300);
        }
        ala_f3_copy(E,uc1,u1,s->lev);
        ala_f3_copy(E,uc2,u2,s->lev);
        ala_f3_copy(E,ucp,up,s->lev);
        strip_bcs_from_residual(E,uc1,s->lev);
        strip_bcs_from_residual(E,uc2,s->lev);
        strip_bcs_from_residual(E,ucp,s->lev);
        assemble_div_rho_u(E,uc1,d1,s->lev);
        assemble_div_rho_u(E,uc2,d2,s->lev);
        assemble_div_rho_u(E,ucp,dp,s->lev);
        dn1=ala_f3_norm(E,d1,s->lev); dn2=ala_f3_norm(E,d2,s->lev);
        dn12=ala_f3_norm(E,dp,s->lev);
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=1;e<=E->lmesh.NPNO[s->lev];e++)
                derr[m][e]=dp[m][e]-d1[m][e]-d2[m][e];
        b_error=ala_f3_norm(E,derr,s->lev)/max(dn12,1.0e-300);
        decomposition_valid=isfinite(bt_error) && isfinite(k_error) &&
            isfinite(b_error);
        if(!ala_f3_reference_action(E,s,q1,a1,cache,1,
            "PRECHECK_H_RIGHT","q1",1,&achieved,&cycles)) actions_valid=0;
        if(!ala_f3_reference_action(E,s,q2,a2,cache,1,
            "PRECHECK_H_RIGHT","q2",1,&achieved,&cycles)) actions_valid=0;
        if(!ala_f3_reference_action(E,s,qp,ap,cache,1,
            "PRECHECK_H_RIGHT","q1_plus_q2",1,&achieved,&cycles)) actions_valid=0;
        a1n=ala_f3_norm(E,a1,s->lev); a2n=ala_f3_norm(E,a2,s->lev);
        if(!ala_f3_reference_action(E,s,q1,repeat1,cache,1,
            "PRECHECK_REPEAT","q1",2,&achieved,&cycles)) actions_valid=0;
        if(!ala_f3_reference_action(E,s,q2,repeat2,cache,1,
            "PRECHECK_REPEAT","q2",2,&achieved,&cycles)) actions_valid=0;
        if(!ala_f3_reference_action(E,s,qp,repeatp,cache,1,
            "PRECHECK_REPEAT","q1_plus_q2",2,&achieved,&cycles)) actions_valid=0;
        ala_f3_copy(E,stest,q1,s->lev); ala_f3_scale(E,stest,0.5,s->lev);
        if(!ala_f3_reference_action(E,s,stest,scaled1,cache,1,
            "PRECHECK_SCALING","0.5q1",1,&achieved,&cycles)) actions_valid=0;
        ala_f3_copy(E,stest,q2,s->lev); ala_f3_scale(E,stest,0.5,s->lev);
        if(!ala_f3_reference_action(E,s,stest,scaled2,cache,1,
            "PRECHECK_SCALING","0.5q2",1,&achieved,&cycles)) actions_valid=0;
        ala_f3_copy(E,stest,qp,s->lev); ala_f3_scale(E,stest,0.5,s->lev);
        if(!ala_f3_reference_action(E,s,stest,scaledp,cache,1,
            "PRECHECK_SCALING","0.5q1_plus_q2",1,&achieved,&cycles)) actions_valid=0;
        repeat_error=max(ala_f3_difference_norm(E,a1,repeat1,s->lev)/max(a1n,1.0e-300),
                         ala_f3_difference_norm(E,a2,repeat2,s->lev)/max(a2n,1.0e-300));
        repeat_error=max(repeat_error,ala_f3_difference_norm(E,ap,repeatp,s->lev)/max(ala_f3_norm(E,ap,s->lev),1.0e-300));
        ala_f3_scale(E,scaled1,2.0,s->lev); ala_f3_scale(E,scaled2,2.0,s->lev);
        ala_f3_scale(E,scaledp,2.0,s->lev);
        scaling_error=max(ala_f3_difference_norm(E,a1,scaled1,s->lev)/max(a1n,1.0e-300),
                          ala_f3_difference_norm(E,a2,scaled2,s->lev)/max(a2n,1.0e-300));
        scaling_error=max(scaling_error,ala_f3_difference_norm(E,ap,scaledp,s->lev)/max(ala_f3_norm(E,ap,s->lev),1.0e-300));
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=1;e<=E->lmesh.NPNO[s->lev];e++)
                ap[m][e]-=a1[m][e]+a2[m][e];
        a_error=ala_f3_norm(E,ap,s->lev)/max(a1n+a2n,1.0e-300);
        a_pass=actions_valid && isfinite(a_error) &&
            a_error<=STRICT_ALA_STAGE_F3_ADDITIVITY_RELATIVE_TOLERANCE &&
            repeat_error<=STRICT_ALA_STAGE_F3_REPEAT_ACTION_RELATIVE_TOLERANCE &&
            scaling_error<=STRICT_ALA_STAGE_F3_SCALING_RELATIVE_TOLERANCE;
    }
    else {
        a1n=a2n=0.0; a_error=-1.0;
    }
    pass=m_pass && a_pass;
    s->fixed_reference_pass=pass;
    s->experiment_valid=actions_valid && isfinite(m_error) && isfinite(a_error);
    s->precheck_completed=1;
    if(E->parallel.me==0) {
        fp=ala_f3_open_root(E,s,"strict_ala_stage_F3_precheck.json","w");
        fprintf(fp,"{\n  \"schema_id\":\"strict-ala-stage-F3-precheck\",\n"
                "  \"M_inverse_definition\":\"||M^-1(q1+q2)-M^-1(q1)-M^-1(q2)||/(||M^-1(q1)||+||M^-1(q2)||)\",\n"
                "  \"M_inverse_q1_norm\":%.17e,\n  \"M_inverse_q2_norm\":%.17e,\n"
                "  \"M_inverse_relative_additivity_error\":%.17e,\n"
                "  \"M_inverse_pass\":%s,\n"
                "  \"B_transpose_relative_additivity_error\":%.17e,\n"
                "  \"K_gamma_path_dependence_diagnostic\":%.17e,\n"
                "  \"B_decomposed_path_dependence_diagnostic\":%.17e,\n"
                "  \"decomposed_diagnostics_authoritative\":false,\n"
                "  \"reference_solution_additivity_error\":%.17e,\n"
                "  \"reference_residual_b1\":%.17e,\n"
                "  \"reference_residual_b2\":%.17e,\n"
                "  \"reference_residual_b12\":%.17e,\n"
                "  \"reference_accuracy_levels\":[1e-6,1e-8,1e-10],\n"
                "  \"reference_accuracy_sweep_R\":[%.17e,%.17e,%.17e],\n"
                "  \"reference_accuracy_sweep_E_x\":[%.17e,%.17e,%.17e],\n"
                "  \"H_right_definition\":\"S_ref M^-1\",\n"
                "  \"H_right_q1_norm\":%.17e,\n  \"H_right_q2_norm\":%.17e,\n"
                "  \"H_right_relative_additivity_error\":%.17e,\n"
                "  \"repeat_action_relative_error\":%.17e,\n"
                "  \"scaling_relative_error\":%.17e,\n"
                "  \"H_right_actions_valid\":%s,\n"
                "  \"H_right_pass\":%s,\n"
                "  \"fixed_reference_pass\":%s,\n"
                "  \"experiment_valid\":%s,\n"
                "  \"threshold\":%.17e,\n  \"pass\":%s\n}\n",
                n1,n2,m_error,m_pass ? "true" : "false",
                bt_error,k_error,b_error,ex_error,rx1,rx2,rxp,
                sweep_R[0],sweep_R[1],sweep_R[2],sweep_E[0],sweep_E[1],sweep_E[2],
                a1n,a2n,a_error,repeat_error,scaling_error,
                actions_valid ? "true" : "false",
                a_pass ? "true" : "false", pass ? "true" : "false",
                s->experiment_valid ? "true" : "false",
                (double)STRICT_ALA_STAGE_F3_ADDITIVITY_RELATIVE_TOLERANCE,
                pass ? "true" : "false");
        fclose(fp);
    }
    ala_schur_free_field(E,qp); ala_schur_free_field(E,z1);
    ala_schur_free_field(E,z2); ala_schur_free_field(E,zp);
    ala_schur_free_field(E,zerr);
    ala_schur_free_field(E,work);
    ala_schur_free_field(E,a1); ala_schur_free_field(E,a2);
    ala_schur_free_field(E,ap);
    ala_schur_free_field(E,g1); ala_schur_free_field(E,g2);
    ala_schur_free_field(E,gp); ala_schur_free_field(E,u1);
    ala_schur_free_field(E,b1w); ala_schur_free_field(E,b2w);
    ala_schur_free_field(E,bpw); ala_schur_free_field(E,xs1);
    ala_schur_free_field(E,ob1); ala_schur_free_field(E,ob2);
    ala_schur_free_field(E,obp);
    ala_schur_free_field(E,xs2); ala_schur_free_field(E,xsp);
    ala_schur_free_field(E,xse); ala_schur_free_field(E,rwork);
    ala_schur_free_field(E,u2); ala_schur_free_field(E,up);
    ala_schur_free_field(E,uerr);
    ala_schur_free_field(E,uc1); ala_schur_free_field(E,uc2);
    ala_schur_free_field(E,ucp);
    ala_schur_free_field(E,vw);
    ala_schur_free_field(E,d1);
    ala_schur_free_field(E,d2); ala_schur_free_field(E,dp);
    ala_schur_free_field(E,derr);
    ala_schur_free_field(E,repeat1); ala_schur_free_field(E,repeat2);
    ala_schur_free_field(E,repeatp); ala_schur_free_field(E,scaled1);
    ala_schur_free_field(E,scaled2); ala_schur_free_field(E,scaledp);
    ala_schur_free_field(E,stest);
    MPI_Barrier(E->parallel.world);
    return s->experiment_valid;
}

static int ala_f3_reference_action(struct All_variables *E,
    struct ala_f3_state *s,double **input,double **output,
    struct ala_pressure_preconditioner_cache *cache,int preconditioned,
    const char *action_id,const char *input_id,int repeat_id,
    double *achieved,long long *cycles)
{
    double **pwork=0,**rhs=0,**saved=0,**u=0,**vwork=0;
    double rhs2,residual2,rhs_rms,residual_rms,relative,output_norm;
    long long total_cycles;
    int m,e,valid,zero_rhs,finite_values;
    int saved_production_logging;
    FILE *fp=NULL;
    (void)cache;
    pwork=ala_schur_alloc_pressure(E,s->lev);
    rhs=ala_schur_alloc_velocity(E,s->lev);
    saved=ala_schur_alloc_velocity(E,s->lev);
    u=ala_schur_alloc_velocity(E,s->lev);
    vwork=ala_schur_alloc_velocity(E,s->lev);
    if(preconditioned)
        apply_ala_pressure_preconditioner(E,input,pwork,output,s->lev,-3,cache);
    else ala_f3_copy(E,pwork,input,s->lev);
    assemble_grad_rho_p(E,pwork,rhs,s->lev);
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(e=0;e<E->lmesh.NEQ[s->lev];e++) saved[m][e]=rhs[m][e];
    rhs2=global_vdot(E,rhs,rhs,s->lev);
    zero_rhs=(rhs2==0.0);
    total_cycles=0;
    if(zero_rhs) {
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=0;e<E->lmesh.NEQ[s->lev];e++) u[m][e]=0.0;
        valid=1;
    }
    else {
        setenv("STRICT_ALA_STAGE_F3_REFERENCE_ACTION","1",1);
        /* A fixed operator requires the same arithmetic path for every RHS.
         * Residual-based stopping made the former solve/refinement split
         * input-dependent and produced O(1e-7) additivity defects even when
         * every final action residual was tiny.  Run one fixed-length solve;
         * the explicit full-action residual below remains authoritative. */
        saved_production_logging=E->control.ala_stage_abc_production_logging;
        E->control.ala_stage_abc_production_logging=0;
        valid=solve_del2_u_fixed_cycles(E,u,rhs,s->lev,
            STRICT_ALA_STAGE_F3_S_REF_FIXED_MG_CYCLES,100);
        E->control.ala_stage_abc_production_logging=saved_production_logging;
        total_cycles=E->control.ala_stage_f3_last_inner_cycles;
        unsetenv("STRICT_ALA_STAGE_F3_REFERENCE_ACTION");
    }
    *cycles=total_cycles;
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(e=0;e<E->lmesh.NEQ[s->lev];e++) rhs[m][e]=saved[m][e];
    (void)ala_schur_velocity_residual(E,u,rhs,vwork,s->lev);
    residual2=global_vdot(E,vwork,vwork,s->lev);
    rhs_rms=sqrt(rhs2/(double)E->mesh.neq);
    residual_rms=sqrt(residual2/(double)E->mesh.neq);
    relative=zero_rhs ? (residual2==0.0 ? 0.0 : -1.0)
                      : residual_rms/rhs_rms;
    *achieved=relative;
    strip_bcs_from_residual(E,u,s->lev);
    assemble_div_rho_u(E,u,output,s->lev);
    output_norm=ala_f3_norm(E,output,s->lev);
    finite_values=isfinite(rhs2) && isfinite(residual2) &&
        isfinite(rhs_rms) && isfinite(residual_rms) && output_norm>=0.0;
    valid=valid && finite_values && relative>=0.0 &&
        relative<=STRICT_ALA_STAGE_F3_S_REF_INNER_RELATIVE_TOLERANCE &&
        *cycles<=STRICT_ALA_STAGE_F3_S_REF_MAX_MG_CYCLES;
    if(E->parallel.me==0) {
        fp=ala_f3_open_root(E,s,"strict_ala_stage_F3_reference_actions.csv","a");
        fprintf(fp,"%s,%s,%s,%d,%.17e,%.17e,%.17e,%.17e,%d,",
            preconditioned ? "S_ref_M_inverse" : "S_ref",action_id,input_id,
            repeat_id,rhs2,residual2,rhs_rms,residual_rms,zero_rhs);
        if(relative>=0.0) fprintf(fp,"%.17e",relative);
        fprintf(fp,",%s,%.17e,%lld,%d,%d,0,%d,%d,%.17e\n",
            zero_rhs ? (residual2==0.0 ? "ZERO_RHS_EXACT" :
                         "ZERO_RHS_NONZERO_RESIDUAL") : "NONZERO_RHS",
            STRICT_ALA_STAGE_F3_S_REF_INNER_RELATIVE_TOLERANCE,*cycles,
            STRICT_ALA_STAGE_F3_S_REF_MAX_MG_CYCLES,valid,finite_values,valid,
            output_norm);
        fclose(fp);
    }
    ala_schur_free_field(E,pwork); ala_schur_free_field(E,rhs);
    ala_schur_free_field(E,saved); ala_schur_free_field(E,u);
    ala_schur_free_field(E,vwork);
    return valid;
}

static void strict_ala_stage_f3_finalize(struct All_variables *E,void *state,
    struct ala_pressure_preconditioner_cache *cache,int iterations,
    int joint_target_reached,double best_continuity,double final_continuity,
    double continuity_target)
{
    struct ala_f3_state *s=(struct ala_f3_state *)state;
    double **Sref[ALA_F3_MAX_MODES],**Aq[ALA_F3_MAX_MODES];
    double **Aprobe[3];
    double **test=0,**out=0;
    double Htrue[ALA_F3_MAX_MODES*ALA_F3_MAX_MODES];
    double Hright[ALA_F3_MAX_MODES*ALA_F3_MAX_MODES];
    double repeat_error[ALA_F3_MAX_MODES],scaling_error[ALA_F3_MAX_MODES];
    double add_error[3],leakage[ALA_F3_MAX_MODES];
    double achieved,norm,action_scale=0.0,numerical_floor,den,sum,orth;
    double symmetry_num=0.0,symmetry_den=0.0,weighted_leakage=0.0;
    double eig[ALA_F3_MAX_MODES],eigvec[ALA_F3_MAX_MODES*ALA_F3_MAX_MODES];
    double sym[ALA_F3_MAX_MODES*ALA_F3_MAX_MODES],maxeig=0.0,mineig;
    long long cycles;
    int i,j,pair_a[3]={0,0,1},pair_b[3]={1,2,2};
    int sref_valid=1,hright_valid=1,htrue_valid=1;
    FILE *fp=NULL,*matrix_fp=NULL;
    if(s==NULL || !s->enabled || s->finalized) return;
    s->finalized=1;
    for(i=0;i<s->selected;i++) { Sref[i]=NULL; Aq[i]=NULL; }
    for(i=0;i<3;i++) Aprobe[i]=NULL;
    for(i=0;i<s->selected;i++) {
        Sref[i]=ala_schur_alloc_pressure(E,s->lev);
        if(!ala_f3_reference_action(E,s,s->mode[i],Sref[i],cache,0,
            "H_TRUE","q",i+1,&achieved,&cycles)) sref_valid=0;
        Aq[i]=ala_schur_alloc_pressure(E,s->lev);
        if(!ala_f3_reference_action(E,s,s->mode[i],Aq[i],cache,1,
            "H_RIGHT","q",i+1,&achieved,&cycles)) hright_valid=0;
        norm=ala_f3_norm(E,Aq[i],s->lev);
        action_scale=max(action_scale,norm);
        if(i<3) Aprobe[i]=Aq[i];
    }
    for(i=s->selected;i<3;i++) {
        Aprobe[i]=ala_schur_alloc_pressure(E,s->lev);
        if(!ala_f3_reference_action(E,s,s->additivity_probe[i],Aprobe[i],cache,1,
            "ADDITIVITY_PROBE","q",i+1,&achieved,&cycles)) hright_valid=0;
    }
    numerical_floor=64.0*DBL_EPSILON*max(action_scale,1.0);
    for(j=0;j<s->selected;j++) for(i=0;i<s->selected;i++) {
        Htrue[i*s->selected+j]=global_pdot(E,s->mode[i],Sref[j],s->lev);
        Hright[i*s->selected+j]=global_pdot(E,s->mode[i],Aq[j],s->lev);
    }
    for(i=0;i<s->selected;i++) {
        out=ala_schur_alloc_pressure(E,s->lev);
        if(!ala_f3_reference_action(E,s,s->mode[i],out,cache,1,
            "REPEAT","q",i+1,&achieved,&cycles)) hright_valid=0;
        repeat_error[i]=ala_f3_difference_norm(E,Aq[i],out,s->lev)
            /max(ala_f3_norm(E,Aq[i],s->lev),numerical_floor);
        ala_schur_free_field(E,out); out=NULL;
        test=ala_schur_alloc_pressure(E,s->lev);
        ala_f3_copy(E,test,s->mode[i],s->lev);
        ala_f3_scale(E,test,STRICT_ALA_STAGE_F3_SCALING_SCALAR_C,s->lev);
        out=ala_schur_alloc_pressure(E,s->lev);
        if(!ala_f3_reference_action(E,s,test,out,cache,1,
            "SCALING","0.5q",i+1,&achieved,&cycles)) hright_valid=0;
        ala_f3_scale(E,out,1.0/STRICT_ALA_STAGE_F3_SCALING_SCALAR_C,s->lev);
        scaling_error[i]=ala_f3_difference_norm(E,Aq[i],out,s->lev)
            /max(ala_f3_norm(E,Aq[i],s->lev),numerical_floor);
        ala_schur_free_field(E,test); ala_schur_free_field(E,out);
        test=out=NULL;
        if(repeat_error[i]>STRICT_ALA_STAGE_F3_REPEAT_ACTION_RELATIVE_TOLERANCE ||
           scaling_error[i]>STRICT_ALA_STAGE_F3_SCALING_RELATIVE_TOLERANCE)
            hright_valid=0;
    }
    for(i=0;i<3;i++) {
        test=ala_schur_alloc_pressure(E,s->lev);
        ala_f3_copy(E,test,s->additivity_probe[pair_a[i]],s->lev);
        ala_f3_axpy(E,test,1.0,s->additivity_probe[pair_b[i]],s->lev);
        out=ala_schur_alloc_pressure(E,s->lev);
        if(!ala_f3_reference_action(E,s,test,out,cache,1,
            "ADDITIVITY","q_pair",i+1,&achieved,&cycles)) hright_valid=0;
        ala_f3_axpy(E,out,-1.0,Aprobe[pair_a[i]],s->lev);
        ala_f3_axpy(E,out,-1.0,Aprobe[pair_b[i]],s->lev);
        den=ala_f3_norm(E,Aprobe[pair_a[i]],s->lev)
            +ala_f3_norm(E,Aprobe[pair_b[i]],s->lev);
        add_error[i]=ala_f3_norm(E,out,s->lev)/max(den,numerical_floor);
        if(add_error[i]>STRICT_ALA_STAGE_F3_ADDITIVITY_RELATIVE_TOLERANCE)
            hright_valid=0;
        ala_schur_free_field(E,test); ala_schur_free_field(E,out);
        test=out=NULL;
    }
    for(j=0;j<s->selected;j++) {
        sum=0.0;
        for(i=0;i<s->selected;i++) {
            orth=global_pdot(E,s->mode[i],Aq[j],s->lev);
            sum+=orth*orth;
        }
        norm=ala_f3_norm(E,Aq[j],s->lev);
        leakage[j]=sqrt(max(norm*norm-sum,0.0))/max(norm,numerical_floor);
        weighted_leakage+=s->weights[j]*leakage[j]*leakage[j];
    }
    weighted_leakage=sqrt(weighted_leakage/max(s->selected_energy,1.0e-300));
    for(i=0;i<s->selected;i++) for(j=0;j<s->selected;j++) {
        double d=Htrue[i*s->selected+j]-Htrue[j*s->selected+i];
        symmetry_num+=d*d;
        symmetry_den+=Htrue[i*s->selected+j]*Htrue[i*s->selected+j];
        sym[i*s->selected+j]=0.5*(Htrue[i*s->selected+j]
                                  +Htrue[j*s->selected+i]);
    }
    htrue_valid=sref_valid && isfinite(symmetry_num) && isfinite(symmetry_den)
        && sqrt(symmetry_num)/max(sqrt(symmetry_den),numerical_floor)<=
           STRICT_ALA_STAGE_F3_H_TRUE_SYMMETRY_RELATIVE_TOLERANCE;
    ala_geneo_jacobi_eigensolve(sym,eigvec,eig,s->selected);
    mineig=eig[0];
    for(i=0;i<s->selected;i++) {
        mineig=min(mineig,eig[i]); maxeig=max(maxeig,fabs(eig[i]));
    }
    if(mineig < -STRICT_ALA_STAGE_F3_H_TRUE_NEGATIVE_CURVATURE_RELATIVE_TOLERANCE*
       max(maxeig,numerical_floor)) htrue_valid=0;
    /* A trajectory-only run may still report descriptive H_right numbers,
     * but they are never authoritative unless the short fixed-reference
     * gate passed. */
    hright_valid=hright_valid && sref_valid && s->fixed_reference_pass;
    if(E->parallel.me==0) {
        matrix_fp=ala_f3_open_root(E,s,"strict_ala_stage_F3_projected_matrices.csv","a");
        for(i=0;i<s->selected;i++) for(j=0;j<s->selected;j++) {
            fprintf(matrix_fp,"H_true,%d,%d,%.17e\n",i+1,j+1,
                Htrue[i*s->selected+j]);
            fprintf(matrix_fp,"H_right,%d,%d,%.17e\n",i+1,j+1,
                Hright[i*s->selected+j]);
        }
        fclose(matrix_fp);
        fp=ala_f3_open_root(E,s,
            "strict_ala_stage_F3_reference_operator_validation_raw.json","w");
        fprintf(fp,"{\n  \"schema\":\"strict-ala-stage-F3-reference-validation-raw-v1\",\n  \"S_ref_validation_pass\":%s,\n  \"H_true_validation_pass\":%s,\n  \"H_right_fixed_operator_valid\":%s,\n  \"H_true_symmetry_relative_defect\":%.17e,\n  \"H_true_minimum_symmetric_eigenvalue\":%.17e,\n  \"reference_action_norm_scale\":%.17e,\n  \"numerical_floor\":%.17e,\n  \"numerical_floor_definition\":\"64*DBL_EPSILON*max(reference_action_norm_scale,1)\",\n  \"repeat_errors\":[",
            sref_valid ? "true" : "false",htrue_valid ? "true" : "false",
            hright_valid ? "true" : "false",
            sqrt(symmetry_num)/max(sqrt(symmetry_den),numerical_floor),mineig,
            action_scale,numerical_floor);
        for(i=0;i<s->selected;i++) fprintf(fp,"%s%.17e",i ? "," : "",repeat_error[i]);
        fprintf(fp,"],\n  \"scaling_errors\":[");
        for(i=0;i<s->selected;i++) fprintf(fp,"%s%.17e",i ? "," : "",scaling_error[i]);
        fprintf(fp,"],\n  \"additivity_errors\":[%.17e,%.17e,%.17e],\n  \"additivity_mode_pairs\":[[1,2],[1,3],[2,3]],\n  \"leakage_per_mode\":[",add_error[0],add_error[1],add_error[2]);
        for(i=0;i<s->selected;i++) fprintf(fp,"%s%.17e",i ? "," : "",leakage[i]);
        fprintf(fp,"],\n  \"leakage_weighted\":%.17e,\n  \"normalized_weighted_rms\":true,\n  \"raw_vector_action_nonfinite\":false\n}\n",weighted_leakage);
        fclose(fp);
        fp=ala_f3_open_root(E,s,"strict_ala_stage_F3_runtime.json","w");
        fprintf(fp,"{\n  \"schema_id\":\"strict-ala-stage-F3-runtime-raw\",\n  \"schema_version\":\"1\",\n  \"stage_id\":\"STRICT_STAGE_F3_OUTER_SCHUR_KRYLOV_RELOCALIZATION\",\n  \"manifest_sha256\":\"%s\",\n  \"canonical_threshold_sha256\":\"%s\",\n  \"generated_threshold_contract_sha256\":\"%s\",\n  \"actual_mpi_world_size\":%d,\n  \"actual_restart\":%d,\n  \"actual_outer_budget\":%d,\n  \"runtime_threshold_values\":{\"S_ref_inner_relative_tolerance\":%.17e,\"S_ref_fixed_MG_cycles\":%d,\"S_ref_max_MG_cycles\":%d,\"repeat_action_relative_tolerance\":%.17e,\"scaling_scalar_c\":%.17e,\"scaling_relative_tolerance\":%.17e,\"additivity_relative_tolerance\":%.17e,\"qr_first_pass_orthogonality_threshold\":%.17e,\"qr_rank_rejection_relative_threshold\":%.17e,\"H_true_symmetry_relative_tolerance\":%.17e,\"H_true_negative_curvature_relative_tolerance\":%.17e},\n  \"outer_preconditioning_orientation\":\"RIGHT_FGMRES\",\n  \"v_j\":\"pressure_space_Arnoldi_basis\",\n  \"z_j\":\"configured_M_inverse_v_j_pressure_correction\",\n  \"s_j\":\"S_prod_j_z_j_sampled_before_Arnoldi_in_place_orthogonalization\",\n  \"pressure_inner_product\":\"global_pdot_euclidean_pressure\",\n  \"iterations\":%d,\n  \"joint_target_reached\":%s,\n  \"best_R_cont\":%.17e,\n  \"final_R_cont\":%.17e,\n  \"R_cont_target\":%.17e,\n  \"base_plateau_confirmed\":%s,\n  \"restart_retention_status\":\"INCONCLUSIVE_SINGLE_BOUNDARY\",\n  \"restart_boundary_observed\":%s,\n  \"viscosity_provenance\":{\"storage\":\"EVI[levmax][cap][(element-1)*vpts+gp]\",\"entity_semantics\":\"owned_element_integration_points\",\"local_order\":\"cap_element_gauss_point\",\"configured_visc_max\":%.17e,\"eta_min\":%.17e,\"eta_max\":%.17e,\"eta_contrast\":%.17e,\"lower_clamp_fraction\":%.17e,\"upper_clamp_fraction\":%.17e,\"global_owned_entry_count\":%lld},\n  \"bytes_per_pressure_vector_per_rank\":%lld,\n  \"maximum_resident_diagnostic_pressure_vector_equivalents_per_rank\":%d,\n  \"estimated_peak_diagnostic_bytes_per_rank\":%lld,\n  \"production_default_change_authorized\":false,\n  \"complete\":true\n}\n",
            s->manifest_sha256,STRICT_ALA_STAGE_F3_THRESHOLD_SOURCE_SHA256,
            s->threshold_contract_sha256,E->parallel.nproc,s->restart,s->outer_budget,
            STRICT_ALA_STAGE_F3_S_REF_INNER_RELATIVE_TOLERANCE,
            STRICT_ALA_STAGE_F3_S_REF_FIXED_MG_CYCLES,
            STRICT_ALA_STAGE_F3_S_REF_MAX_MG_CYCLES,
            STRICT_ALA_STAGE_F3_REPEAT_ACTION_RELATIVE_TOLERANCE,
            STRICT_ALA_STAGE_F3_SCALING_SCALAR_C,
            STRICT_ALA_STAGE_F3_SCALING_RELATIVE_TOLERANCE,
            STRICT_ALA_STAGE_F3_ADDITIVITY_RELATIVE_TOLERANCE,
            STRICT_ALA_STAGE_F3_QR_FIRST_PASS_ORTHOGONALITY_THRESHOLD,
            STRICT_ALA_STAGE_F3_QR_RANK_REJECTION_RELATIVE_THRESHOLD,
            STRICT_ALA_STAGE_F3_H_TRUE_SYMMETRY_RELATIVE_TOLERANCE,
            STRICT_ALA_STAGE_F3_H_TRUE_NEGATIVE_CURVATURE_RELATIVE_TOLERANCE,
            iterations,joint_target_reached ? "true" : "false",best_continuity,
            final_continuity,continuity_target,
            (!joint_target_reached && iterations==60 &&
             best_continuity>continuity_target) ? "true" : "false",
            iterations>=s->restart ? "true" : "false",
            E->viscosity.max_value,s->viscosity_min,s->viscosity_max,
            s->viscosity_contrast,s->viscosity_lower_fraction,
            s->viscosity_upper_fraction,s->viscosity_count,
            s->bytes_per_pressure_vector,
            6*ALA_F3_MAX_BASIS+2*ALA_E2_LATE_COUNT+2*s->selected+20,
            s->bytes_per_pressure_vector*(long long)(6*ALA_F3_MAX_BASIS+
                2*ALA_E2_LATE_COUNT+2*s->selected+20));
        fclose(fp);
        fprintf(E->fp,"STRICT_ALA_STAGE_F3_NUMERICAL_COMPLETE iterations=%d "
            "S_ref_valid=%d H_true_valid=%d H_right_valid=%d\n",iterations,
            sref_valid,htrue_valid,hright_valid);
        fprintf(stderr,"STRICT_ALA_STAGE_F3_NUMERICAL_COMPLETE iterations=%d "
            "S_ref_valid=%d H_true_valid=%d H_right_valid=%d\n",iterations,
            sref_valid,htrue_valid,hright_valid);
        fflush(E->fp);
        fflush(stderr);
    }
    for(i=0;i<s->selected;i++) {
        ala_schur_free_field(E,Sref[i]); ala_schur_free_field(E,Aq[i]);
        ala_schur_free_field(E,s->mode[i]);
    }
    for(i=s->selected;i<3;i++) {
        ala_schur_free_field(E,Aprobe[i]);
        ala_schur_free_field(E,s->additivity_probe[i]);
    }
    for(i=0;i<ALA_F3_SPACES;i++) {
        ala_f3_clear_basis(E,&s->basis[i][ALA_F3_ALL]);
        ala_f3_clear_basis(E,&s->basis[i][ALA_F3_CYCLE]);
    }
    ala_e2_free_observer(E,s->e0); ala_e2_free_observer(E,s->e1);
    free(s);
}
