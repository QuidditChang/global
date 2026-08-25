/* Strict-ALA Stage F1a: observation-only, actual-mode Schwarz forensics.
 *
 * This file is included by Stokes_flow_Incomp.c so it can call the exact
 * pressure-map and tight Schur helpers without widening the production API.
 * The opt-in environment gate is mandatory; the default solver path remains
 * bit-for-bit on strict_ala_f1a_patch_stage == -1. */

#define ALA_F1A_MAX_MODES 6

static FILE *ala_f1a_open_root(struct All_variables *E,const char *directory,
                               const char *name,const char *mode)
{
    char path[1024];
    FILE *fp;
    snprintf(path,sizeof(path),"%s/%s",directory,name);
    fp=fopen(path,mode);
    if(fp==NULL) myerror(E,"Unable to open strict-ALA Stage-F1a output");
    return fp;
}

static void ala_f1a_zero_pressure(struct All_variables *E,double **q,int lev)
{
    int m,e;
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(e=1;e<=E->lmesh.NPNO[lev];e++) q[m][e]=0.0;
}

static void ala_f1a_copy_pressure(struct All_variables *E,double **dst,
                                  double **src,int lev)
{
    int m,e;
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(e=1;e<=E->lmesh.NPNO[lev];e++) dst[m][e]=src[m][e];
}

static void ala_f1a_patch_action(struct All_variables *E,double **q,double **z,
    double **work,int lev,int stage,struct ala_pressure_preconditioner_cache *cache)
{
    if(stage<0 || stage>3 || strict_ala_f1a_patch_stage!=-1)
        myerror(E,"Invalid nested strict-ALA Stage-F1a patch action");
    strict_ala_f1a_patch_stage=stage;
    apply_ala_pressure_preconditioner(E,q,z,work,lev,-3,cache);
    strict_ala_f1a_patch_stage=-1;
}

static double ala_f1a_true_action(struct All_variables *E,double **q,
    double **action,double **rhs,double **u,double **velocity_work,int lev,
    FILE *tight,int *call_id,int mode,const char *role)
{
    double rhs2,achieved,target;
    long long before=E->control.ala_stage_abc_inner_cycle_count;
    double global_vdot();
    void assemble_div_rho_u();
    assemble_grad_rho_p(E,q,rhs,lev);
    rhs2=global_vdot(E,rhs,rhs,lev);
    target=E->control.ala_schur_diagnostic_tight_tolerance
        *sqrt(max(rhs2,0.0)/max((double)E->mesh.neq,1.0));
    achieved=ala_schur_tight_solve(E,u,rhs,velocity_work,lev,
        E->control.ala_schur_diagnostic_tight_tolerance);
    assemble_div_rho_u(E,u,action,lev);
    (*call_id)++;
    if(E->parallel.me==0) {
        fprintf(tight,"%d,%s,%d,%.17e,%.17e,%.17e,%.17e,%lld,%d,1\n",
            *call_id,role,mode,sqrt(max(rhs2,0.0)),
            E->control.ala_schur_diagnostic_tight_tolerance,target,achieved,
            E->control.ala_stage_abc_inner_cycle_count-before,
            E->control.ala_schur_diagnostic_max_cycles);
        fflush(tight);
    }
    return achieved;
}

static int ala_f1a_selected_count(double *values,int n,double *energy)
{
    int i,selected=0;
    double total=0.0,cumulative=0.0,tolerance;
    for(i=0;i<n;i++) total+=values[i];
    tolerance=n*DBL_EPSILON*max(values[0],1.0);
    for(i=0;i<n && selected<ALA_F1A_MAX_MODES;i++) {
        if(values[i]<=tolerance) break;
        cumulative+=values[i]/max(total,1.0e-300);
        selected++;
        if(cumulative>=0.95) break;
    }
    *energy=cumulative;
    return selected;
}

static void ala_f1a_write_row(FILE *fp,int mode,const char *stage,double weight,
    double xnorm,double ynorm,double residual,double standalone,double ep,
    double cosine,double qx,double qy,double eby,double yby,double self,
    double cross,double delta,double closure,double achieved)
{
    fprintf(fp,"%d,%s,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,"
        "%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,1\n",
        mode,stage,weight,xnorm,ynorm,residual,standalone,ep,cosine,qx,qy,
        eby,yby,self,cross,delta,closure,achieved);
    fflush(fp);
}

static void ala_f1a_patch_audit(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache,double ***modes,int selected,
    double **source,double **action,double **rhs,double **u,double **vwork,
    int lev,FILE *tight,int *call_id,FILE *patch,FILE *local,FILE *nonlocal)
{
    int mode,m,b,i,j,e,n,owner_m=0,owner_b=0,owner_n=0;
    int owner_elements[ALA_PATCH_MAX_ELEMENTS];
    double best_depth=0.0;
    for(mode=0;mode<selected;mode++) {
        struct { double value; int rank; } local_best,global_best;
        double metrics[9]={0,0,0,0,0,0,0,0,0};
        double distance_local[5]={0,0,0,0,0},distance_global[5];
        memset(owner_elements,0,sizeof(owner_elements));
        ala_f1a_zero_pressure(E,source,lev);
        ala_f1a_zero_pressure(E,action,lev);
        local_best.value=-1.0; local_best.rank=E->parallel.me;
        owner_m=owner_b=owner_n=0; best_depth=0.0;
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(b=0;b<cache->blocks[m];b++) {
                double score=0.0;
                n=cache->size[m][b];
                for(i=0;i<n;i++) {
                    e=cache->elements[m][b*cache->patch_capacity+i];
                    score+=modes[mode][m][e]*modes[mode][m][e];
                }
                if(score>local_best.value) {
                    local_best.value=score; owner_m=m; owner_b=b; owner_n=n;
                }
            }
        MPI_Allreduce(&local_best,&global_best,1,MPI_DOUBLE_INT,MPI_MAXLOC,
                      E->parallel.world);
        if(global_best.value<=1.0e-300)
            myerror(E,"Stage-F1a representative patch has zero contribution");
        if(E->parallel.me==global_best.rank) {
            double *L=cache->chol[owner_m]
                +(size_t)owner_b*cache->patch_capacity*cache->patch_capacity;
            double bvec[ALA_PATCH_MAX_ELEMENTS],sol[ALA_PATCH_MAX_ELEMENTS];
            double temp[ALA_PATCH_MAX_ELEMENTS],avec[ALA_PATCH_MAX_ELEMENTS];
            double residual2=0.0,b2=0.0,sol2=0.0,ai2=0.0,depth=0.0;
            for(i=0;i<owner_n;i++) {
                owner_elements[i]=cache->elements[owner_m]
                    [owner_b*cache->patch_capacity+i];
                source[owner_m][owner_elements[i]]=modes[mode][owner_m]
                    [owner_elements[i]]/sqrt(global_best.value);
                bvec[i]=source[owner_m][owner_elements[i]];
                depth+=(1.0-E->ECO[lev][owner_m][owner_elements[i]].centre[3])
                    *E->data.radius_km;
            }
            best_depth=depth/max((double)owner_n,1.0);
            for(i=0;i<owner_n;i++) {
                double sum=bvec[i];
                for(j=0;j<i;j++) sum-=L[i*cache->patch_capacity+j]*sol[j];
                sol[i]=sum/L[i*cache->patch_capacity+i];
            }
            for(i=owner_n-1;i>=0;i--) {
                double sum=sol[i];
                for(j=i+1;j<owner_n;j++) sum-=L[j*cache->patch_capacity+i]*sol[j];
                sol[i]=sum/L[i*cache->patch_capacity+i];
            }
            for(i=0;i<owner_n;i++) {
                temp[i]=0.0;
                for(j=i;j<owner_n;j++) temp[i]+=L[j*cache->patch_capacity+i]*sol[j];
            }
            for(i=0;i<owner_n;i++) {
                double az=0.0;
                for(j=0;j<=i;j++) az+=L[i*cache->patch_capacity+j]*temp[j];
                residual2+=(az-bvec[i])*(az-bvec[i]);
                b2+=bvec[i]*bvec[i]; sol2+=sol[i]*sol[i];
            }
            for(i=0;i<owner_n;i++) {
                temp[i]=0.0;
                for(j=i;j<owner_n;j++) temp[i]+=L[j*cache->patch_capacity+i]*bvec[j];
            }
            for(i=0;i<owner_n;i++) {
                avec[i]=0.0;
                for(j=0;j<=i;j++) avec[i]+=L[i*cache->patch_capacity+j]*temp[j];
                ai2+=avec[i]*avec[i];
            }
            metrics[0]=sqrt(b2); metrics[1]=sqrt(sol2);
            metrics[2]=sqrt(residual2/max(b2,1.0e-300));
            metrics[3]=sqrt(ai2); metrics[8]=best_depth;
        }
        MPI_Bcast(&owner_m,1,MPI_INT,global_best.rank,E->parallel.world);
        MPI_Bcast(&owner_b,1,MPI_INT,global_best.rank,E->parallel.world);
        MPI_Bcast(&owner_n,1,MPI_INT,global_best.rank,E->parallel.world);
        MPI_Bcast(owner_elements,ALA_PATCH_MAX_ELEMENTS,MPI_INT,
                  global_best.rank,E->parallel.world);
        MPI_Bcast(&best_depth,1,MPI_DOUBLE,global_best.rank,E->parallel.world);
        ala_f1a_true_action(E,source,action,rhs,u,vwork,lev,tight,call_id,
                           mode+1,"REPRESENTATIVE_PATCH_TRUE_ACTION");
        if(E->parallel.me==global_best.rank) {
            double ti2=0.0,diff2=0.0,dot=0.0;
            double *L=cache->chol[owner_m]
                +(size_t)owner_b*cache->patch_capacity*cache->patch_capacity;
            double bvec[ALA_PATCH_MAX_ELEMENTS],temp[ALA_PATCH_MAX_ELEMENTS];
            double avec[ALA_PATCH_MAX_ELEMENTS];
            for(i=0;i<owner_n;i++) bvec[i]=source[owner_m][owner_elements[i]];
            for(i=0;i<owner_n;i++) {
                temp[i]=0.0;
                for(j=i;j<owner_n;j++) temp[i]+=L[j*cache->patch_capacity+i]*bvec[j];
            }
            for(i=0;i<owner_n;i++) {
                avec[i]=0.0;
                for(j=0;j<=i;j++) avec[i]+=L[i*cache->patch_capacity+j]*temp[j];
                ti2+=action[owner_m][owner_elements[i]]*action[owner_m][owner_elements[i]];
                diff2+=(avec[i]-action[owner_m][owner_elements[i]])
                    *(avec[i]-action[owner_m][owner_elements[i]]);
                dot+=avec[i]*action[owner_m][owner_elements[i]];
            }
            metrics[4]=sqrt(ti2);
            metrics[5]=sqrt(diff2/max(ti2,1.0e-300));
            metrics[6]=dot/max(metrics[3]*metrics[4],1.0e-300);
            metrics[7]=global_best.value;
        }
        MPI_Bcast(metrics,9,MPI_DOUBLE,global_best.rank,E->parallel.world);
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=1;e<=E->lmesh.NPNO[lev];e++) {
                int distance=4;
                if(E->parallel.me==global_best.rank && m==owner_m) {
                    int ex=(e-1)/E->lmesh.ELZ[lev]%E->lmesh.ELX[lev];
                    int ey=(e-1)/(E->lmesh.ELZ[lev]*E->lmesh.ELX[lev]);
                    int ez=(e-1)%E->lmesh.ELZ[lev];
                    for(i=0;i<owner_n;i++) {
                        int pe=owner_elements[i];
                        int px=(pe-1)/E->lmesh.ELZ[lev]%E->lmesh.ELX[lev];
                        int py=(pe-1)/(E->lmesh.ELZ[lev]*E->lmesh.ELX[lev]);
                        int pz=(pe-1)%E->lmesh.ELZ[lev];
                        distance=min(distance,abs(ex-px)+abs(ey-py)+abs(ez-pz));
                    }
                }
                distance=min(distance,4);
                distance_local[distance]+=action[m][e]*action[m][e];
            }
        MPI_Allreduce(distance_local,distance_global,5,MPI_DOUBLE,MPI_SUM,
                      E->parallel.world);
        if(E->parallel.me==0) {
            double sum=distance_global[0]+distance_global[1]+distance_global[2]
                +distance_global[3]+distance_global[4];
            long long id=(long long)global_best.rank*1000000LL+owner_b;
            const char *depth=best_depth<200.0 ? "0_200" :
                (best_depth<410.0 ? "200_410" : "410_660");
            fprintf(patch,"%lld,TOP_CONTRIBUTION_PER_MODE,%d,%s,top,%d,0\n",
                    id,mode+1,depth,mode+1);
            fprintf(local,"%lld,%d,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%s,HIGH_FREQUENCY,1\n",
                    id,mode+1,metrics[0],metrics[1],metrics[2],metrics[3],metrics[4],
                    metrics[5],metrics[6],metrics[7],depth);
            fprintf(nonlocal,"%lld_mode_%d,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,1\n",
                    id,mode+1,distance_global[0],distance_global[1],
                    distance_global[2],distance_global[3],distance_global[4],
                    (distance_global[2]+distance_global[3]+distance_global[4])
                        /max(sum,1.0e-300));
            fflush(patch); fflush(local); fflush(nonlocal);
        }
        ala_f1a_zero_pressure(E,source,lev);
        ala_f1a_zero_pressure(E,action,lev);
    }
}

static int strict_ala_stage_f1a_run(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache,int lev)
{
    static int executed=0;
    const char *required=getenv("STRICT_ALA_STAGE_F1A_REQUIRED");
    const char *snapshot_dir=getenv("STRICT_ALA_STAGE_F1A_E2_SNAPSHOT_DIR");
    const char *output_dir=getenv("STRICT_ALA_STAGE_F1A_OUTPUT_DIR");
    struct ala_e2_observer template_observer,*e0,*e1;
    struct ala_e2_observer *cases[2*ALA_E2_LATE_COUNT];
    int indices[2*ALA_E2_LATE_COUNT];
    double values[2*ALA_E2_LATE_COUNT];
    double vectors[4*ALA_E2_LATE_COUNT*ALA_E2_LATE_COUNT];
    double norms[2*ALA_E2_LATE_COUNT],total=0.0,selected_energy=0.0;
    double **q[ALA_F1A_MAX_MODES],**xB[ALA_F1A_MAX_MODES];
    double **yB[ALA_F1A_MAX_MODES],**xS[ALA_F1A_MAX_MODES];
    double **yS[ALA_F1A_MAX_MODES];
    double **x=0,**y=0,**work=0,**rhs=0,**u=0,**vwork=0,**combined=0;
    FILE *tel=NULL,*tight=NULL,*matrix=NULL,*patch=NULL,*local=NULL,*nonlocal=NULL;
    FILE *structure=NULL,*support=NULL,*linearity=NULL;
    int i,j,m,e,stage,n=2*ALA_E2_LATE_COUNT,selected,call_id=0;
    double weight,achieved=0.0,achievedB[ALA_F1A_MAX_MODES];
    double achievedS[ALA_F1A_MAX_MODES];
    static const char *stage_names[4]={"S0_RAW_LOCAL","S1_RIGHT_SCALING",
        "S2_PARTITION_WEIGHTING","S3_LEFT_SCALING"};
    if(required==NULL) return 0;
    if(executed) return 1;
    executed=1;
    if(snapshot_dir==NULL || output_dir==NULL)
        myerror(E,"Incomplete strict-ALA Stage-F1a environment contract");
    if(E->parallel.nproc!=384 || E->parallel.nprocx!=4 ||
       E->parallel.nprocy!=4 || E->parallel.nprocz!=2)
        myerror(E,"E2_REPLAY_LAYOUT_MISMATCH");
    if(!E->control.ala_shallow_patch_preconditioner ||
       strcmp(E->control.ala_shallow_patch_velocity_solver,"element_vanka")!=0)
        myerror(E,"Stage-F1a requires the configured additive patch map");
    memset(&template_observer,0,sizeof(template_observer));
    template_observer.lev=lev;
    snprintf(template_observer.shared_dir,sizeof(template_observer.shared_dir),
             "%s",snapshot_dir);
    snprintf(template_observer.output_dir,sizeof(template_observer.output_dir),
             "%s",output_dir);
    e0=ala_e2_read_case(E,&template_observer,"E0");
    e1=ala_e2_read_case(E,&template_observer,"E1");
    for(i=0;i<ALA_E2_LATE_COUNT;i++) {
        cases[i]=e0; indices[i]=i+ALA_E2_LATE_FIRST;
        cases[i+ALA_E2_LATE_COUNT]=e1;
        indices[i+ALA_E2_LATE_COUNT]=i+ALA_E2_LATE_FIRST;
    }
    {
        double covariance[4*ALA_E2_LATE_COUNT*ALA_E2_LATE_COUNT];
        ala_e2_covariance(E,cases,indices,n,covariance);
        if(!ala_e2_pod(covariance,n,1,values,vectors,norms))
            myerror(E,"Stage-F1a POD reconstruction failed");
    }
    selected=ala_f1a_selected_count(values,n,&selected_energy);
    if(selected<=0 || selected_energy<0.95-1.0e-12)
        myerror(E,"Stage-F1a 95%-or-6 mode selection failed");
    for(i=0;i<n;i++) total+=values[i];
    x=ala_schur_alloc_pressure(E,lev); y=ala_schur_alloc_pressure(E,lev);
    work=ala_schur_alloc_pressure(E,lev); combined=ala_schur_alloc_pressure(E,lev);
    rhs=ala_schur_alloc_velocity(E,lev); u=ala_schur_alloc_velocity(E,lev);
    vwork=ala_schur_alloc_velocity(E,lev);
    if(E->parallel.me==0) {
        tel=ala_f1a_open_root(E,output_dir,"strict_ala_stage_F1a_telescoping.csv","w");
        fprintf(tel,"POD_mode,stage,POD_energy_weight,x_norm,Sx_norm,residual_energy,standalone_residual_energy,E_P,cosine,qTPq,qTSPq,eB_yS,yB_yS,self_term,useful_cross_term,delta_energy,delta_identity_error,tight_solve_achieved,valid\n");
        tight=ala_f1a_open_root(E,output_dir,"strict_ala_stage_F1a_tight_solves.csv","w");
        fprintf(tight,"call_id,role,POD_mode,rhs_norm,requested_relative_tolerance,target_residual,achieved_relative_residual,cycles,max_cycles,converged\n");
        matrix=ala_f1a_open_root(E,output_dir,"strict_ala_stage_F1a_projected_matrices.csv","w");
        fprintf(matrix,"matrix,row_mode,column_mode,value\n");
        patch=ala_f1a_open_root(E,output_dir,"strict_ala_stage_F1a_patch_selection.csv","w");
        fprintf(patch,"patch_ID,selection_rule,POD_mode,depth_stratum,contribution_quantile,deterministic_order,outlier_guard\n");
        local=ala_f1a_open_root(E,output_dir,"strict_ala_stage_F1a_local_operator.csv","w");
        fprintf(local,"patch_or_bin_ID,POD_mode,local_RHS_norm,local_solution_norm,local_solve_relative_residual,Ai_v_norm,Ti_v_norm,local_action_error,local_action_cosine,contribution_weight,depth_category,scale_category,valid\n");
        nonlocal=ala_f1a_open_root(E,output_dir,"strict_ala_stage_F1a_nonlocality.csv","w");
        fprintf(nonlocal,"source_vector_or_bin,inside_patch_energy,distance_1_energy,distance_2_energy,distance_3_energy,farther_energy,outside_patch_plus_one_fraction,valid\n");
        structure=ala_f1a_open_root(E,output_dir,"strict_ala_stage_F1a_patch_structure.csv","w");
        fprintf(structure,"supported_DOFs,multiplicity_min,multiplicity_max,invalid_multiplicity,top_W,mid_W,transition_W,partition_formula,valid\n");
        support=ala_f1a_open_root(E,output_dir,"strict_ala_stage_F1a_support.csv","w");
        fprintf(support,"POD_mode,q_energy_in_support_fraction,eB_energy_in_support_fraction,valid\n");
        linearity=ala_f1a_open_root(E,output_dir,"strict_ala_stage_F1a_linearity.csv","w");
        fprintf(linearity,"role,difference_norm,combined_norm,sum_norm,relative_defect,solve_floor,warn_limit,hard_limit,status,valid\n");
        fflush(tel); fflush(tight); fflush(matrix); fflush(patch); fflush(local);
        fflush(nonlocal); fflush(structure); fflush(support); fflush(linearity);
    }
    {
        double local_stats[3]={0.0,1.0e300,0.0},global_min,global_max;
        long long local_supported=0,global_supported,local_invalid=0,global_invalid;
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=1;e<=E->lmesh.NPNO[lev];e++)
                if(cache->multiplicity[m][e]>0) {
                    double mult=(double)cache->multiplicity[m][e];
                    local_supported++; local_stats[1]=min(local_stats[1],mult);
                    local_stats[2]=max(local_stats[2],mult);
                    if(!isfinite(mult) || mult<=0.0) local_invalid++;
                }
        MPI_Allreduce(&local_supported,&global_supported,1,MPI_LONG_LONG_INT,MPI_SUM,E->parallel.world);
        MPI_Allreduce(&local_invalid,&global_invalid,1,MPI_LONG_LONG_INT,MPI_SUM,E->parallel.world);
        MPI_Allreduce(&local_stats[1],&global_min,1,MPI_DOUBLE,MPI_MIN,E->parallel.world);
        MPI_Allreduce(&local_stats[2],&global_max,1,MPI_DOUBLE,MPI_MAX,E->parallel.world);
        if(E->parallel.me==0) {
            fprintf(structure,"%lld,%.17e,%.17e,%lld,1.0,%.17e,%.17e,"
                "symmetric_inverse_sqrt_multiplicity_with_weighted_local_gauge,1\n",
                global_supported,global_min,global_max,global_invalid,
                sqrt(E->control.ala_shallow_patch_mid_action_scale),
                sqrt(E->control.ala_shallow_patch_transition_action_scale));
            fflush(structure);
        }
    }
    for(i=0;i<selected;i++) {
        double qq,qx,qy,xx,yy,residual,standalone,cosine;
        q[i]=ala_e2_build_mode(E,cases,indices,n,vectors,norms,i);
        xB[i]=ala_schur_alloc_pressure(E,lev); yB[i]=ala_schur_alloc_pressure(E,lev);
        xS[i]=ala_schur_alloc_pressure(E,lev); yS[i]=ala_schur_alloc_pressure(E,lev);
        weight=values[i]/max(total,1.0e-300);
        ala_schur_preconditioner_action(E,q[i],xB[i],work,combined,lev,0,cache);
        achieved=ala_f1a_true_action(E,xB[i],yB[i],rhs,u,vwork,lev,tight,
                                     &call_id,i+1,"BPI");
        achievedB[i]=achieved;
        qq=global_pdot(E,q[i],q[i],lev); qx=global_pdot(E,q[i],xB[i],lev);
        qy=global_pdot(E,q[i],yB[i],lev); xx=global_pdot(E,xB[i],xB[i],lev);
        yy=global_pdot(E,yB[i],yB[i],lev); residual=qq+yy-2.0*qy;
        {
            double local_support[4]={0,0,0,0},global_support[4];
            for(m=1;m<=E->sphere.caps_per_proc;m++)
                for(e=1;e<=E->lmesh.NPNO[lev];e++) {
                    double eb=q[i][m][e]-yB[i][m][e];
                    double depth=(1.0-E->ECO[lev][m][e].centre[3])*E->data.radius_km;
                    local_support[0]+=q[i][m][e]*q[i][m][e];
                    local_support[1]+=eb*eb;
                    if(depth<=E->control.ala_shallow_patch_depth_km) {
                        local_support[2]+=q[i][m][e]*q[i][m][e];
                        local_support[3]+=eb*eb;
                    }
                }
            MPI_Allreduce(local_support,global_support,4,MPI_DOUBLE,MPI_SUM,E->parallel.world);
            if(E->parallel.me==0) {
                fprintf(support,"%d,%.17e,%.17e,1\n",i+1,
                    global_support[2]/max(global_support[0],1.0e-300),
                    global_support[3]/max(global_support[1],1.0e-300));
                fflush(support);
            }
        }
        if(E->parallel.me==0)
            ala_f1a_write_row(tel,i+1,"BPI",weight,sqrt(max(xx,0.0)),
                sqrt(max(yy,0.0)),residual,residual,sqrt(max(residual,0.0)),
                qy/max(sqrt(max(qq*yy,0.0)),1.0e-300),qx,qy,0,0,0,0,0,0,achieved);
        for(stage=0;stage<4;stage++) {
            double yby,eby,self,cross,delta,closure;
            ala_f1a_zero_pressure(E,x,lev); ala_f1a_zero_pressure(E,y,lev);
            ala_f1a_patch_action(E,q[i],x,work,lev,stage,cache);
            achieved=ala_f1a_true_action(E,x,y,rhs,u,vwork,lev,tight,&call_id,
                                         i+1,stage_names[stage]);
            qx=global_pdot(E,q[i],x,lev); qy=global_pdot(E,q[i],y,lev);
            xx=global_pdot(E,x,x,lev); yy=global_pdot(E,y,y,lev);
            yby=global_pdot(E,yB[i],y,lev); eby=qy-yby;
            standalone=qq+yy-2.0*qy;
            residual=qq+global_pdot(E,yB[i],yB[i],lev)+yy
                -2.0*global_pdot(E,q[i],yB[i],lev)-2.0*qy+2.0*yby;
            self=yy; cross=2.0*eby; delta=self-cross;
            closure=fabs((residual-(qq+global_pdot(E,yB[i],yB[i],lev)
                -2.0*global_pdot(E,q[i],yB[i],lev)))-delta);
            cosine=qy/max(sqrt(max(qq*yy,0.0)),1.0e-300);
            if(E->parallel.me==0)
                ala_f1a_write_row(tel,i+1,stage_names[stage],weight,
                    sqrt(max(xx,0.0)),sqrt(max(yy,0.0)),residual,standalone,
                    sqrt(max(standalone,0.0)),cosine,qx,qy,eby,yby,self,cross,
                    delta,closure,achieved);
            if(stage==3) {
                ala_f1a_copy_pressure(E,xS[i],x,lev);
                ala_f1a_copy_pressure(E,yS[i],y,lev);
                achievedS[i]=achieved;
            }
        }
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=1;e<=E->lmesh.NPNO[lev];e++) {
                x[m][e]=xB[i][m][e]+xS[i][m][e];
                y[m][e]=yB[i][m][e]+yS[i][m][e];
            }
        qx=global_pdot(E,q[i],x,lev); qy=global_pdot(E,q[i],y,lev);
        xx=global_pdot(E,x,x,lev); yy=global_pdot(E,y,y,lev);
        residual=qq+yy-2.0*qy;
        {
        double configured_yby=global_pdot(E,yB[i],yS[i],lev);
        if(E->parallel.me==0)
            ala_f1a_write_row(tel,i+1,"CONFIGURED",weight,sqrt(max(xx,0.0)),
                sqrt(max(yy,0.0)),residual,residual,sqrt(max(residual,0.0)),
                qy/max(sqrt(max(qq*yy,0.0)),1.0e-300),qx,qy,0,
                configured_yby,0,0,0,0,
                max(achieved,E->control.ala_schur_diagnostic_tight_tolerance));
        }
    }
    ala_f1a_patch_audit(E,cache,q,selected,x,y,rhs,u,vwork,lev,tight,
                        &call_id,patch,local,nonlocal);
    /* One independently solved dominant-mode linearity guard. */
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(e=1;e<=E->lmesh.NPNO[lev];e++)
            combined[m][e]=xB[0][m][e]+xS[0][m][e];
    achieved=ala_f1a_true_action(E,combined,y,rhs,u,vwork,lev,tight,&call_id,1,
                                 "CONFIGURED_LINEARITY_GUARD");
    {
        double combined_norm=sqrt(max(global_pdot(E,y,y,lev),0.0));
        double sum_norm,difference_norm,relative_defect,solve_floor;
        double warn_limit,hard_limit,scale;
        int linearity_valid;
        const char *status;
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=1;e<=E->lmesh.NPNO[lev];e++)
                combined[m][e]=yB[0][m][e]+yS[0][m][e];
        sum_norm=sqrt(max(global_pdot(E,combined,combined,lev),0.0));
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=1;e<=E->lmesh.NPNO[lev];e++) y[m][e]-=combined[m][e];
        difference_norm=sqrt(max(global_pdot(E,y,y,lev),0.0));
        scale=max(max(combined_norm,sum_norm),1.0e-300);
        relative_defect=difference_norm/scale;
        solve_floor=max(achieved,max(achievedB[0],achievedS[0]));
        warn_limit=max(100.0*solve_floor,1.0e-8);
        hard_limit=max(100000.0*solve_floor,1.0e-5);
        linearity_valid=isfinite(difference_norm) && isfinite(relative_defect)
            && relative_defect<=hard_limit;
        status=!linearity_valid ? "FAIL" :
            (relative_defect>warn_limit ? "WARN" : "PASS");
        if(E->parallel.me==0) {
            fprintf(linearity,"CONFIGURED,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%s,%d\n",
                difference_norm,combined_norm,sum_norm,relative_defect,
                solve_floor,warn_limit,hard_limit,status,linearity_valid);
            fflush(linearity); fflush(tel); fflush(tight); fflush(patch);
            fflush(local); fflush(nonlocal); fflush(structure); fflush(support);
            fprintf(E->fp,"STRICT_ALA_STAGE_F1A_LINEARITY status=%s "
                "difference_norm=%.17e combined_norm=%.17e sum_norm=%.17e "
                "relative_defect=%.17e solve_floor=%.17e warn_limit=%.17e "
                "hard_limit=%.17e\n",status,difference_norm,combined_norm,
                sum_norm,relative_defect,solve_floor,warn_limit,hard_limit);
            fflush(E->fp);
        }
        if(!linearity_valid)
            myerror(E,"Stage-F1a normalized S_gamma linearity guard failed");
    }
    if(call_id>40) myerror(E,"Stage-F1a tight solve budget exceeded");
    for(i=0;i<selected;i++) for(j=0;j<selected;j++) {
        double mb=global_pdot(E,q[i],xB[j],lev);
        double ms=global_pdot(E,q[i],xS[j],lev);
        double hb=global_pdot(E,q[i],yB[j],lev);
        double hs=global_pdot(E,q[i],yS[j],lev);
        if(E->parallel.me==0) {
            fprintf(matrix,"M_B,%d,%d,%.17e\n",i+1,j+1,mb);
            fprintf(matrix,"M_S,%d,%d,%.17e\n",i+1,j+1,ms);
            fprintf(matrix,"M_cfg,%d,%d,%.17e\n",i+1,j+1,mb+ms);
            fprintf(matrix,"H_B,%d,%d,%.17e\n",i+1,j+1,hb);
            fprintf(matrix,"H_S,%d,%d,%.17e\n",i+1,j+1,hs);
            fprintf(matrix,"H_cfg,%d,%d,%.17e\n",i+1,j+1,hb+hs);
        }
    }
    if(E->parallel.me==0) fflush(matrix);
    if(E->parallel.me==0) {
        fprintf(E->fp,"STRICT_ALA_STAGE_F1A_COMPLETE selected_modes=%d "
            "selected_energy=%.17e tight_solves=%d K_gamma_build_count=1 "
            "production_solver_modified=0\n",selected,selected_energy,call_id);
        fflush(E->fp);
        fclose(tel); fclose(tight); fclose(matrix); fclose(patch); fclose(local); fclose(nonlocal);
        fclose(structure); fclose(support); fclose(linearity);
    }
    for(i=0;i<selected;i++) {
        ala_schur_free_field(E,q[i]); ala_schur_free_field(E,xB[i]);
        ala_schur_free_field(E,yB[i]); ala_schur_free_field(E,xS[i]);
        ala_schur_free_field(E,yS[i]);
    }
    ala_schur_free_field(E,x); ala_schur_free_field(E,y);
    ala_schur_free_field(E,work); ala_schur_free_field(E,combined);
    ala_schur_free_field(E,rhs); ala_schur_free_field(E,u);
    ala_schur_free_field(E,vwork);
    ala_e2_free_observer(E,e0); ala_e2_free_observer(E,e1);
    return 1;
}
