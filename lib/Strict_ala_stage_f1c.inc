/* Strict-ALA Stage F1c: actual-mode balanced coarse qualification.
 *
 * This opt-in diagnostic reconstructs the frozen E2 Joint-POD modes and
 * changes one pressure-map object only: it wraps the configured F1b map P in
 * the symmetric balanced coarse form
 *
 *   C + (I-C S) P (I-S C),   C = Q (Q^T S Q)^-1 Q^T.
 *
 * Q and S Q are retained only for this diagnostic run.  The production
 * default and the configured B, K_gamma, patch, W, BPI, and velocity solve
 * paths are unchanged.
 */

#define ALA_F1C_SYMMETRY_TOLERANCE 1.0e-8
#define ALA_F1C_PIVOT_TOLERANCE 1.0e-12

static int ala_f1c_factor(double *matrix,int n,double *minimum_ratio)
{
    int i,j,k;
    double diagonal_max=0.0,sum,pivot;
    for(i=0;i<n;i++) diagonal_max=max(diagonal_max,matrix[i*n+i]);
    if(!isfinite(diagonal_max) || diagonal_max<=0.0) return 0;
    *minimum_ratio=1.0e300;
    for(i=0;i<n;i++)
        for(j=0;j<=i;j++) {
            sum=matrix[i*n+j];
            for(k=0;k<j;k++) sum-=matrix[i*n+k]*matrix[j*n+k];
            if(i==j) {
                pivot=sum;
                if(!isfinite(pivot) ||
                   pivot<=ALA_F1C_PIVOT_TOLERANCE*diagonal_max) return 0;
                *minimum_ratio=min(*minimum_ratio,pivot/diagonal_max);
                matrix[i*n+j]=sqrt(pivot);
            }
            else matrix[i*n+j]=sum/matrix[j*n+j];
        }
    for(i=0;i<n;i++) for(j=i+1;j<n;j++) matrix[i*n+j]=0.0;
    return 1;
}

static void ala_f1c_solve(const double *factor,const double *rhs,double *x,int n)
{
    int i,j;
    double y[ALA_F1A_MAX_MODES],sum;
    for(i=0;i<n;i++) {
        sum=rhs[i];
        for(j=0;j<i;j++) sum-=factor[i*n+j]*y[j];
        y[i]=sum/factor[i*n+i];
    }
    for(i=n-1;i>=0;i--) {
        sum=y[i];
        for(j=i+1;j<n;j++) sum-=factor[j*n+i]*x[j];
        x[i]=sum/factor[i*n+i];
    }
}

static int strict_ala_stage_f1c_run(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache,int lev)
{
    static int executed=0;
    const char *snapshot_dir=getenv("STRICT_ALA_STAGE_F1C_E2_SNAPSHOT_DIR");
    const char *output_dir=getenv("STRICT_ALA_STAGE_F1C_OUTPUT_DIR");
    struct ala_e2_observer template_observer,*e0,*e1;
    struct ala_e2_observer *cases[2*ALA_E2_LATE_COUNT];
    int indices[2*ALA_E2_LATE_COUNT];
    double values[2*ALA_E2_LATE_COUNT];
    double vectors[4*ALA_E2_LATE_COUNT*ALA_E2_LATE_COUNT];
    double norms[2*ALA_E2_LATE_COUNT],total=0.0,selected_energy=0.0;
    double **q[ALA_F1A_MAX_MODES],**Sq[ALA_F1A_MAX_MODES];
    double **z[ALA_F1A_MAX_MODES],**Sz[ALA_F1A_MAX_MODES];
    double **fine_rhs=0,**fine_z=0,**work=0,**aux=0;
    double **velocity_rhs=0,**velocity=0,**velocity_work=0;
    double T[ALA_F1A_MAX_MODES*ALA_F1A_MAX_MODES];
    double factor[ALA_F1A_MAX_MODES*ALA_F1A_MAX_MODES];
    double M[ALA_F1A_MAX_MODES*ALA_F1A_MAX_MODES];
    double H[ALA_F1A_MAX_MODES*ALA_F1A_MAX_MODES];
    FILE *metrics=NULL,*projected=NULL,*tight=NULL,*decomposition=NULL;
    FILE *alpha_scan=NULL,*summary=NULL;
    int i,j,m,e,n=2*ALA_E2_LATE_COUNT,selected,call_id=0;
    double symmetry=0.0,scale=0.0,min_pivot=0.0,weight;
    double alpha_a=0.0,alpha_b=0.0,alpha_c=0.0,weight_sum=0.0;

    if(getenv("STRICT_ALA_STAGE_F1C_REQUIRED")==NULL) return 0;
    if(executed) return 1;
    executed=1;
    if(snapshot_dir==NULL || output_dir==NULL)
        myerror(E,"Incomplete strict-ALA Stage-F1c environment contract");
    if(E->parallel.nproc!=384 || E->parallel.nprocx!=4 ||
       E->parallel.nprocy!=4 || E->parallel.nprocz!=2)
        myerror(E,"F1C_E2_REPLAY_LAYOUT_MISMATCH");
    if(!E->control.ala_shallow_patch_preconditioner ||
       strcmp(E->control.ala_shallow_patch_velocity_solver,"element_vanka")!=0 ||
       strcmp(E->control.ala_shallow_patch_local_operator,
              "operator_consistent")!=0)
        myerror(E,"Stage-F1c requires the validated F1b pressure map");

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
            myerror(E,"Stage-F1c POD reconstruction failed");
    }
    selected=ala_f1a_selected_count(values,n,&selected_energy);
    if(selected<=0 || selected>ALA_F1A_MAX_MODES ||
       selected_energy<0.95-1.0e-12)
        myerror(E,"Stage-F1c frozen POD selection failed");
    for(i=0;i<n;i++) total+=values[i];

    fine_rhs=ala_schur_alloc_pressure(E,lev);
    fine_z=ala_schur_alloc_pressure(E,lev);
    work=ala_schur_alloc_pressure(E,lev);
    aux=ala_schur_alloc_pressure(E,lev);
    velocity_rhs=ala_schur_alloc_velocity(E,lev);
    velocity=ala_schur_alloc_velocity(E,lev);
    velocity_work=ala_schur_alloc_velocity(E,lev);
    if(E->parallel.me==0) {
        metrics=ala_f1a_open_root(E,output_dir,
            "strict_ala_stage_F1c_true_mode.csv","w");
        fprintf(metrics,"POD_mode,POD_energy_weight,E_P,cosine,qTPq,qTSPq,"
                "coarse_residual_projection_max,tight_solve_achieved,valid\n");
        projected=ala_f1a_open_root(E,output_dir,
            "strict_ala_stage_F1c_projected_matrices.csv","w");
        fprintf(projected,"matrix,row_mode,column_mode,value\n");
        tight=ala_f1a_open_root(E,output_dir,
            "strict_ala_stage_F1c_tight_solves.csv","w");
        fprintf(tight,"call_id,role,POD_mode,rhs_norm,requested_relative_tolerance,"
                "target_residual,achieved_relative_residual,cycles,max_cycles,converged\n");
        decomposition=ala_f1a_open_root(E,output_dir,
            "strict_ala_stage_F1c_decomposition.csv","w");
        fprintf(decomposition,"POD_mode,POD_energy_weight,E_coarse_only,"
                "E_conditioned_fine,E_right_additive,E_balanced,"
                "conditioned_rhs_relative,conditioned_fine_action_relative,"
                "left_balance_action_relative,coarse_projection_max,"
                "right_projection_max,balanced_projection_max,valid\n");
        alpha_scan=ala_f1a_open_root(E,output_dir,
            "strict_ala_stage_F1c_alpha_scan.csv","w");
        fprintf(alpha_scan,"alpha,E_P_RMS,is_clipped_optimum,valid\n");
    }
    memset(T,0,sizeof(T)); memset(M,0,sizeof(M)); memset(H,0,sizeof(H));
    for(i=0;i<selected;i++) {
        q[i]=ala_e2_build_mode(E,cases,indices,n,vectors,norms,i);
        Sq[i]=ala_schur_alloc_pressure(E,lev);
        z[i]=ala_schur_alloc_pressure(E,lev);
        Sz[i]=ala_schur_alloc_pressure(E,lev);
        ala_f1a_true_action(E,q[i],Sq[i],velocity_rhs,velocity,velocity_work,
                            lev,tight,&call_id,i+1,"COARSE_BASIS_SQ");
    }
    for(i=0;i<selected;i++)
        for(j=0;j<selected;j++) {
            T[i*selected+j]=global_pdot(E,q[i],Sq[j],lev);
            scale=max(scale,fabs(T[i*selected+j]));
        }
    for(i=0;i<selected;i++)
        for(j=0;j<i;j++) {
            symmetry=max(symmetry,fabs(T[i*selected+j]-T[j*selected+i])
                         /max(scale,1.0e-300));
            T[i*selected+j]=T[j*selected+i]=
                0.5*(T[i*selected+j]+T[j*selected+i]);
        }
    memcpy(factor,T,(size_t)selected*selected*sizeof(double));
    if(symmetry>ALA_F1C_SYMMETRY_TOLERANCE ||
       !ala_f1c_factor(factor,selected,&min_pivot)) {
        if(E->parallel.me==0) {
            fprintf(stderr,"STRICT_ALA_STAGE_F1C_COARSE_REJECTION "
                    "symmetry=%.17e pivot_ratio=%.17e\n",symmetry,min_pivot);
            fflush(stderr);
        }
        myerror(E,"Stage-F1c coarse matrix failed SPD contract");
    }

    for(j=0;j<selected;j++) {
        double rhs_c[ALA_F1A_MAX_MODES]={0},c[ALA_F1A_MAX_MODES]={0};
        double rhs_d[ALA_F1A_MAX_MODES]={0},d[ALA_F1A_MAX_MODES]={0};
        double qq,yy,qy,qz,ep,cosine,projection=0.0,achieved;
        double yc2,qyc,rf2,yf2,rfyf,right2,qright,left2;
        double coarse_ep,fine_ep,right_ep,conditioned_rhs_relative;
        double conditioned_action_relative,left_action_relative;
        double coarse_projection=0.0,right_projection=0.0;
        double delta2,rc_delta,rc2;
        for(i=0;i<selected;i++) rhs_c[i]=global_pdot(E,q[i],q[j],lev);
        ala_f1c_solve(factor,rhs_c,c,selected);
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=1;e<=E->lmesh.NPNO[lev];e++) {
                double value=q[j][m][e];
                for(i=0;i<selected;i++) value-=Sq[i][m][e]*c[i];
                fine_rhs[m][e]=value;
            }
        ala_schur_preconditioner_action(E,fine_rhs,fine_z,work,aux,lev,3,cache);
        for(i=0;i<selected;i++) rhs_d[i]=global_pdot(E,Sq[i],fine_z,lev);
        ala_f1c_solve(factor,rhs_d,d,selected);
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=1;e<=E->lmesh.NPNO[lev];e++) {
                double value=fine_z[m][e];
                for(i=0;i<selected;i++) value+=q[i][m][e]*(c[i]-d[i]);
                z[j][m][e]=value;
            }
        achieved=ala_f1a_true_action(E,z[j],Sz[j],velocity_rhs,velocity,
            velocity_work,lev,tight,&call_id,j+1,"BALANCED_CANDIDATE");
        qq=global_pdot(E,q[j],q[j],lev);
        yy=global_pdot(E,Sz[j],Sz[j],lev);
        qy=global_pdot(E,q[j],Sz[j],lev);
        qz=global_pdot(E,q[j],z[j],lev);
        ep=sqrt(max(qq+yy-2.0*qy,0.0)/max(qq,1.0e-300));
        cosine=qy/sqrt(max(qq*yy,1.0e-300));
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=1;e<=E->lmesh.NPNO[lev];e++) {
                double yc=0.0,yd=0.0;
                for(i=0;i<selected;i++) {
                    yc+=Sq[i][m][e]*c[i];
                    yd+=Sq[i][m][e]*d[i];
                }
                work[m][e]=yc;
                aux[m][e]=Sz[j][m][e]-yc+yd;
            }
        yc2=global_pdot(E,work,work,lev);
        qyc=global_pdot(E,q[j],work,lev);
        rf2=global_pdot(E,fine_rhs,fine_rhs,lev);
        yf2=global_pdot(E,aux,aux,lev);
        rfyf=global_pdot(E,fine_rhs,aux,lev);
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=1;e<=E->lmesh.NPNO[lev];e++)
                fine_z[m][e]=work[m][e]+aux[m][e];
        right2=global_pdot(E,fine_z,fine_z,lev);
        qright=global_pdot(E,q[j],fine_z,lev);
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=1;e<=E->lmesh.NPNO[lev];e++) {
                fine_rhs[m][e]=Sz[j][m][e]-work[m][e];
                aux[m][e]=Sz[j][m][e]-fine_z[m][e];
            }
        delta2=global_pdot(E,fine_rhs,fine_rhs,lev);
        rc_delta=global_pdot(E,q[j],fine_rhs,lev)
            -global_pdot(E,work,fine_rhs,lev);
        left2=global_pdot(E,aux,aux,lev);
        rc2=max(qq+yc2-2.0*qyc,0.0);
        coarse_ep=sqrt(rc2/max(qq,1.0e-300));
        fine_ep=sqrt(max(rf2+yf2-2.0*rfyf,0.0)/max(rf2,1.0e-300));
        right_ep=sqrt(max(qq+right2-2.0*qright,0.0)/max(qq,1.0e-300));
        conditioned_rhs_relative=sqrt(max(rf2,0.0)/max(qq,1.0e-300));
        conditioned_action_relative=sqrt(max(yf2,0.0)/max(qq,1.0e-300));
        left_action_relative=sqrt(max(left2,0.0)/max(qq,1.0e-300));
        for(i=0;i<selected;i++) {
            double residual=global_pdot(E,q[i],q[j],lev)
                -global_pdot(E,q[i],Sz[j],lev);
            double coarse_residual=global_pdot(E,q[i],q[j],lev)
                -global_pdot(E,q[i],work,lev);
            double right_residual=global_pdot(E,q[i],q[j],lev)
                -global_pdot(E,q[i],fine_z,lev);
            projection=max(projection,fabs(residual));
            coarse_projection=max(coarse_projection,fabs(coarse_residual));
            right_projection=max(right_projection,fabs(right_residual));
        }
        weight=values[j]/max(total,1.0e-300);
        alpha_a+=weight*rc2/max(qq,1.0e-300);
        alpha_b+=weight*rc_delta/max(qq,1.0e-300);
        alpha_c+=weight*delta2/max(qq,1.0e-300);
        weight_sum+=weight;
        if(E->parallel.me==0) {
            fprintf(metrics,"%d,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,1\n",
                j+1,weight,ep,cosine,qz,qy,
                projection,achieved);
            fflush(metrics);
            fprintf(decomposition,"%d,%.17e,%.17e,%.17e,%.17e,%.17e,"
                "%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,1\n",
                j+1,weight,coarse_ep,fine_ep,right_ep,ep,
                conditioned_rhs_relative,conditioned_action_relative,
                left_action_relative,coarse_projection,right_projection,
                projection);
            fflush(decomposition);
        }
    }
    for(i=0;i<selected;i++)
        for(j=0;j<selected;j++) {
            M[i*selected+j]=global_pdot(E,q[i],z[j],lev);
            H[i*selected+j]=global_pdot(E,q[i],Sz[j],lev);
            if(E->parallel.me==0) {
                fprintf(projected,"T,%d,%d,%.17e\n",i+1,j+1,
                        T[i*selected+j]);
                fprintf(projected,"M_balanced,%d,%d,%.17e\n",i+1,j+1,
                        M[i*selected+j]);
                fprintf(projected,"H_balanced,%d,%d,%.17e\n",i+1,j+1,
                        H[i*selected+j]);
            }
        }
    if(E->parallel.me==0) {
        static const double fixed_alpha[]={0.0,0.125,0.25,0.5,0.75,1.0};
        double optimum=alpha_c>0.0?alpha_b/alpha_c:0.0;
        int a;
        optimum=max(0.0,min(1.0,optimum));
        for(a=0;a<6;a++) {
            double alpha=fixed_alpha[a];
            double rms=sqrt(max(alpha_a-2.0*alpha*alpha_b+
                                alpha*alpha*alpha_c,0.0)
                            /max(weight_sum,1.0e-300));
            fprintf(alpha_scan,"%.17e,%.17e,0,1\n",alpha,rms);
        }
        fprintf(alpha_scan,"%.17e,%.17e,1,1\n",optimum,
            sqrt(max(alpha_a-2.0*optimum*alpha_b+
                     optimum*optimum*alpha_c,0.0)
                 /max(weight_sum,1.0e-300)));
        fflush(alpha_scan);
        fflush(projected);
        summary=ala_f1a_open_root(E,output_dir,
            "strict_ala_stage_F1c_runtime.json","w");
        fprintf(summary,"{\n  \"schema\": \"strict-ala-stage-F1c-runtime-v1\",\n"
            "  \"selected_modes\": %d,\n  \"selected_energy\": %.17e,\n"
            "  \"coarse_matrix_symmetry_defect\": %.17e,\n"
            "  \"coarse_matrix_minimum_pivot_ratio\": %.17e,\n"
            "  \"tight_solve_count\": %d,\n"
            "  \"candidate\": \"balanced_actual_mode_coarse\",\n"
            "  \"causal_review\": \"coarse_only_conditioned_fine_right_additive_left_balance_alpha_family\",\n"
            "  \"fixed_alpha_count\": 6,\n"
            "  \"production_default_change_authorized\": false\n}\n",
            selected,selected_energy,symmetry,min_pivot,call_id);
        fclose(summary);
        fprintf(E->fp,"STRICT_ALA_STAGE_F1C_CAUSAL_REVIEW_V2 "
                "modes=%d tight_solves=%d\n",selected,call_id);
        fprintf(E->fp,"STRICT_ALA_STAGE_F1C_COMPLETE modes=%d tight_solves=%d "
                "candidate=balanced_actual_mode_coarse\n",selected,call_id);
        fflush(E->fp);
    }
    if(metrics) fclose(metrics); if(projected) fclose(projected);
    if(tight) fclose(tight); if(decomposition) fclose(decomposition);
    if(alpha_scan) fclose(alpha_scan);
    for(i=0;i<selected;i++) {
        ala_schur_free_field(E,q[i]); ala_schur_free_field(E,Sq[i]);
        ala_schur_free_field(E,z[i]); ala_schur_free_field(E,Sz[i]);
    }
    ala_schur_free_field(E,fine_rhs); ala_schur_free_field(E,fine_z);
    ala_schur_free_field(E,work); ala_schur_free_field(E,aux);
    ala_schur_free_field(E,velocity_rhs); ala_schur_free_field(E,velocity);
    ala_schur_free_field(E,velocity_work);
    ala_e2_free_observer(E,e0); ala_e2_free_observer(E,e1);
    return 1;
}
