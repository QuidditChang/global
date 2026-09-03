/* Strict-ALA F3: observation-only outer right-FGMRES relocalization.
 *
 * This file deliberately does not alter a solver vector.  It reconstructs
 * the frozen E2 residual POD modes, observes v_j, z_j=M^-1 v_j and the
 * unmodified production image s_j=S_prod,j z_j, and maintains bounded
 * global_pdot QR bases.  Reference actions are evaluated only after the
 * production trajectory has ended. */

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
    double weights[ALA_F3_MAX_MODES],selected_energy;
    unsigned long long mode_checksum[ALA_F3_MAX_MODES];
    double **mode[ALA_F3_MAX_MODES];
    struct ala_f3_basis basis[ALA_F3_SPACES][2];
    struct ala_e2_observer *e0,*e1;
    char output_dir[512];
    unsigned long long snapshot_hash;
    long long bytes_per_pressure_vector;
};

static const char *ala_f3_space_name[ALA_F3_SPACES]={"V","Z","Y"};

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
    double pre,post1,post2,alpha,max_relative=0.0,value;
    int i,accepted=0;
    const double reject_relative=1.0e-12;
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
    if(pre>0.0 && post2>reject_relative*pre && basis->rank<ALA_F3_MAX_BASIS &&
       isfinite(post2)) {
        ala_f3_scale(E,r,1.0/post2,s->lev);
        basis->q[basis->rank++]=r;
        accepted=1;
    }
    else ala_schur_free_field(E,r);
    if(E->parallel.me==0)
        fprintf(fp,"%d,%s,%s,%d,%d,%s,%.17e,%.17e,%.17e,%.17e,%d,%d\n",
            iteration,space,scope,accepted,accepted ? 0 : 1,
            pre<=0.0 ? "NONFINITE_OR_ZERO" :
            (basis->rank>=ALA_F3_MAX_BASIS && !accepted ? "BASIS_LIMIT" :
             (!accepted ? "RELATIVE_RANK_REJECTION" : "NONE")),
            pre,post1,post2,max_relative,1,basis->rank);
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
    fprintf(fp,"iteration,space,scope,accepted,rejected,rejection_reason,pre_norm,post_first_pass_norm,post_second_pass_norm,first_pass_relative_orthogonality,second_pass_performed,resulting_rank\n"); fclose(fp);
    fp=ala_f3_open_root(E,s,"strict_ala_stage_F3_mode_coefficients.csv","w");
    fprintf(fp,"iteration,mode_id,mode_energy_weight,coefficient,absolute_modal_amplitude,residual_global_pdot_norm,mode_status\n"); fclose(fp);
    fp=ala_f3_open_root(E,s,"strict_ala_stage_F3_pod_explanation.csv","w");
    fprintf(fp,"iteration,residual_global_pdot_norm,POD_projected_norm,f_Q,continuity_relative,momentum_relative,krylov_recursive_absolute,krylov_explicit_absolute,cumulative_inner_solves,cumulative_MG_cycles,cumulative_schur_actions,elapsed_seconds\n"); fclose(fp);
    fp=ala_f3_open_root(E,s,"strict_ala_stage_F3_restart_retention.csv","w");
    fprintf(fp,"restart_boundary_iteration,restart_cycle,mode_id,pre_restart_Y_reachability,post_restart_iteration,post_restart_cycle_local_Y_reachability,recovery_metric_available,recovery_iterations\n"); fclose(fp);
    fp=ala_f3_open_root(E,s,"strict_ala_stage_F3_reference_actions.csv","w");
    fprintf(fp,"operator,action_id,input_id,repeat_id,requested_relative_tolerance,max_MG_cycles,achieved_relative,MG_cycles,status,output_norm\n"); fclose(fp);
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
    const char *required=getenv("STRICT_ALA_STAGE_F3_REQUIRED");
    const char *snapshot_dir=getenv("STRICT_ALA_STAGE_F3_E2_SNAPSHOT_DIR");
    const char *output_dir=getenv("STRICT_ALA_STAGE_F3_OUTPUT_DIR");
    FILE *fp,*coefficient_fp=NULL;
    int i,j,n=2*ALA_E2_LATE_COUNT;
    unsigned long long checksum;
    if(required==NULL) return NULL;
    if(snapshot_dir==NULL || output_dir==NULL)
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
    s->selected=ala_f1a_selected_count(values,n,&s->selected_energy);
    if(s->selected<3 || s->selected>ALA_F3_MAX_MODES ||
       s->selected_energy<0.95-1.0e-12)
        myerror(E,"Stage-F3 authoritative POD selection mismatch");
    for(i=0;i<n;i++) total+=values[i];
    s->snapshot_hash=1469598103934665603ULL;
    for(i=0;i<s->selected;i++) {
        s->mode[i]=ala_e2_build_mode(E,cases,indices,n,vectors,norms,i);
        s->weights[i]=values[i]/max(total,1.0e-300);
        checksum=ala_f3_mode_checksum(E,s->mode[i],lev);
        s->mode_checksum[i]=checksum;
        s->snapshot_hash^=checksum;
    }
    for(i=0;i<s->selected;i++) for(j=0;j<s->selected;j++) {
        value=global_pdot(E,s->mode[i],s->mode[j],lev)-(i==j ? 1.0 : 0.0);
        if(fabs(value)>orth_max) orth_max=fabs(value);
    }
    ala_f3_create_outputs(E,s);
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
        fprintf(fp,"{\n  \"schema\":\"strict-ala-stage-F3-pod-reconstruction-v1\",\n  \"method\":\"deterministic_E2_joint_direction_case_balanced_snapshot_reconstruction\",\n  \"pressure_inner_product\":\"global_pdot_euclidean_pressure\",\n  \"selected_modes\":%d,\n  \"pod_weight_sum\":%.17e,\n  \"orthogonality_max_abs_defect\":%.17e,\n  \"deterministic_sign_rule\":\"match_original_E2_eigensolver_order_then_validate_joint_coefficients\",\n  \"snapshot_aggregate_checksum\":\"%016llx\",\n  \"mode_checksums\":[",
            s->selected,s->selected_energy,orth_max,s->snapshot_hash);
        for(i=0;i<s->selected;i++) {
            fprintf(fp,"%s{\"mode_id\":%d,\"checksum\":\"%016llx\",\"weight\":%.17e}",
                i ? "," : "",i+1,s->mode_checksum[i],s->weights[i]);
        }
        fprintf(fp,"],\n  \"coefficient_trajectory_validation\":\"DEFERRED_TO_FAIL_CLOSED_ANALYZER\",\n  \"valid\":%s\n}\n",
            orth_max<=1.0e-10 ? "true" : "false");
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

static int ala_f3_reference_action(struct All_variables *E,
    struct ala_f3_state *s,double **input,double **output,
    struct ala_pressure_preconditioner_cache *cache,int preconditioned,
    const char *action_id,const char *input_id,int repeat_id,
    double *achieved,long long *cycles)
{
    double **pwork=0,**rhs=0,**saved=0,**u=0,**vwork=0;
    double rhs2,target,output_norm;
    int m,e,valid;
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
    target=1.0e-10*sqrt(rhs2/max((double)E->mesh.neq,1.0));
    setenv("STRICT_ALA_STAGE_F3_REFERENCE_ACTION","1",1);
    valid=solve_del2_u_bounded(E,u,rhs,target,s->lev,1024,100);
    unsetenv("STRICT_ALA_STAGE_F3_REFERENCE_ACTION");
    *cycles=E->control.ala_stage_f3_last_inner_cycles;
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(e=0;e<E->lmesh.NEQ[s->lev];e++) rhs[m][e]=saved[m][e];
    *achieved=ala_schur_velocity_residual(E,u,rhs,vwork,s->lev);
    strip_bcs_from_residual(E,u,s->lev);
    assemble_div_rho_u(E,u,output,s->lev);
    output_norm=ala_f3_norm(E,output,s->lev);
    valid=valid && isfinite(*achieved) && *achieved<=1.0e-10 &&
        *cycles<=1024 && output_norm>=0.0;
    if(E->parallel.me==0) {
        fp=ala_f3_open_root(E,s,"strict_ala_stage_F3_reference_actions.csv","a");
        fprintf(fp,"%s,%s,%s,%d,1.00000000000000004e-10,1024,%.17e,%lld,%s,%.17e\n",
            preconditioned ? "S_ref_M_inverse" : "S_ref",action_id,input_id,
            repeat_id,*achieved,*cycles,valid ? "CONVERGED" : "FAILED",
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
    for(i=0;i<s->selected;i++) {
        Sref[i]=ala_schur_alloc_pressure(E,s->lev);
        if(!ala_f3_reference_action(E,s,s->mode[i],Sref[i],cache,0,
            "H_TRUE","q",i+1,&achieved,&cycles)) sref_valid=0;
        Aq[i]=ala_schur_alloc_pressure(E,s->lev);
        if(!ala_f3_reference_action(E,s,s->mode[i],Aq[i],cache,1,
            "H_RIGHT","q",i+1,&achieved,&cycles)) hright_valid=0;
        norm=ala_f3_norm(E,Aq[i],s->lev);
        action_scale=max(action_scale,norm);
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
        ala_f3_copy(E,test,s->mode[i],s->lev); ala_f3_scale(E,test,0.5,s->lev);
        out=ala_schur_alloc_pressure(E,s->lev);
        if(!ala_f3_reference_action(E,s,test,out,cache,1,
            "SCALING","0.5q",i+1,&achieved,&cycles)) hright_valid=0;
        ala_f3_scale(E,out,2.0,s->lev);
        scaling_error[i]=ala_f3_difference_norm(E,Aq[i],out,s->lev)
            /max(ala_f3_norm(E,Aq[i],s->lev),numerical_floor);
        ala_schur_free_field(E,test); ala_schur_free_field(E,out);
        test=out=NULL;
        if(repeat_error[i]>1.0e-10 || scaling_error[i]>1.0e-8)
            hright_valid=0;
    }
    for(i=0;i<3;i++) {
        test=ala_schur_alloc_pressure(E,s->lev);
        ala_f3_copy(E,test,s->mode[pair_a[i]],s->lev);
        ala_f3_axpy(E,test,1.0,s->mode[pair_b[i]],s->lev);
        out=ala_schur_alloc_pressure(E,s->lev);
        if(!ala_f3_reference_action(E,s,test,out,cache,1,
            "ADDITIVITY","q_pair",i+1,&achieved,&cycles)) hright_valid=0;
        ala_f3_axpy(E,out,-1.0,Aq[pair_a[i]],s->lev);
        ala_f3_axpy(E,out,-1.0,Aq[pair_b[i]],s->lev);
        den=ala_f3_norm(E,Aq[pair_a[i]],s->lev)
            +ala_f3_norm(E,Aq[pair_b[i]],s->lev);
        add_error[i]=ala_f3_norm(E,out,s->lev)/max(den,numerical_floor);
        if(add_error[i]>1.0e-8) hright_valid=0;
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
        && sqrt(symmetry_num)/max(sqrt(symmetry_den),numerical_floor)<=1.0e-8;
    ala_geneo_jacobi_eigensolve(sym,eigvec,eig,s->selected);
    mineig=eig[0];
    for(i=0;i<s->selected;i++) {
        mineig=min(mineig,eig[i]); maxeig=max(maxeig,fabs(eig[i]));
    }
    if(mineig < -1.0e-12*max(maxeig,numerical_floor)) htrue_valid=0;
    hright_valid=hright_valid && sref_valid;
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
        fprintf(fp,"{\n  \"schema\":\"strict-ala-stage-F3-runtime-v1\",\n  \"outer_preconditioning_orientation\":\"RIGHT_FGMRES\",\n  \"v_j\":\"pressure_space_Arnoldi_basis\",\n  \"z_j\":\"configured_M_inverse_v_j_pressure_correction\",\n  \"s_j\":\"S_prod_j_z_j_sampled_before_Arnoldi_in_place_orthogonalization\",\n  \"pressure_inner_product\":\"global_pdot_euclidean_pressure\",\n  \"iterations\":%d,\n  \"joint_target_reached\":%s,\n  \"best_R_cont\":%.17e,\n  \"final_R_cont\":%.17e,\n  \"R_cont_target\":%.17e,\n  \"base_plateau_confirmed\":%s,\n  \"restart_retention_status\":\"INCONCLUSIVE_SINGLE_BOUNDARY\",\n  \"restart_boundary_observed\":%s,\n  \"bytes_per_pressure_vector_per_rank\":%lld,\n  \"maximum_resident_diagnostic_vectors_per_rank\":%d,\n  \"estimated_peak_diagnostic_bytes_per_rank\":%lld,\n  \"production_default_change_authorized\":false\n}\n",
            iterations,joint_target_reached ? "true" : "false",best_continuity,
            final_continuity,continuity_target,
            (!joint_target_reached && iterations==60 &&
             best_continuity>continuity_target) ? "true" : "false",
            iterations>=s->restart ? "true" : "false",
            s->bytes_per_pressure_vector,
            6*ALA_F3_MAX_BASIS+2*s->selected+8,
            s->bytes_per_pressure_vector*(long long)(6*ALA_F3_MAX_BASIS+
                                                     2*s->selected+8));
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
    for(i=0;i<ALA_F3_SPACES;i++) {
        ala_f3_clear_basis(E,&s->basis[i][ALA_F3_ALL]);
        ala_f3_clear_basis(E,&s->basis[i][ALA_F3_CYCLE]);
    }
    ala_e2_free_observer(E,s->e0); ala_e2_free_observer(E,s->e1);
    free(s);
}
