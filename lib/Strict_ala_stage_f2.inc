/* Strict-ALA Stage F2: frozen actual-mode K_gamma multigrid review.
 *
 * F1d exhausted pressure-side coarse enrichment.  This observation-only
 * experiment keeps the assembled K_gamma operator and frozen E2 pressure
 * modes fixed, and changes only the velocity multigrid smoother:
 *
 *   CONFIGURED_VANKA   current strict-ALA element-Vanka smoother
 *   DIAGONAL_ROLLBACK  historical nodal/diagonal smoother
 *
 * Each trajectory is compared with a tight K_gamma inverse on the identical
 * G^T q right-hand side.  Production defaults are restored before return.
 */

#define ALA_F2_CHECKPOINTS 7

static void ala_f2_copy_velocity(struct All_variables *E,double **dst,
                                 double **src,int lev)
{
    int m,i;
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(i=0;i<=E->lmesh.NEQ[lev];i++) dst[m][i]=src[m][i];
}

static int strict_ala_stage_f2_run(struct All_variables *E,int lev)
{
    static int executed=0;
    const char *snapshot_dir=getenv("STRICT_ALA_STAGE_F2_E2_SNAPSHOT_DIR");
    const char *output_dir=getenv("STRICT_ALA_STAGE_F2_OUTPUT_DIR");
    const int checkpoints[ALA_F2_CHECKPOINTS]={1,2,4,8,16,32,64};
    const char *variant_name[2]={"CONFIGURED_VANKA","DIAGONAL_ROLLBACK"};
    struct ala_e2_observer template_observer,*e0,*e1;
    struct ala_e2_observer *cases[2*ALA_E2_LATE_COUNT];
    int indices[2*ALA_E2_LATE_COUNT];
    double values[2*ALA_E2_LATE_COUNT];
    double vectors[4*ALA_E2_LATE_COUNT*ALA_E2_LATE_COUNT];
    double norms[2*ALA_E2_LATE_COUNT],total=0.0,selected_energy=0.0;
    double **q=0,**rhs=0,**rhs_saved=0,**u_exact=0,**u=0;
    double **velocity_work=0,**Ku=0,**error=0,**Kerror=0;
    double **S_exact=0,**S_approx=0,**S_error=0;
    FILE *trajectory=NULL,*tight=NULL,*summary=NULL;
    int i,j,m,e,n=2*ALA_E2_LATE_COUNT,selected,variant,cycle,checkpoint;
    int original_vanka,call_id=0;
    double rhs2,exact2,exact_energy,target,achieved,weight;
    double global_vdot(),global_pdot(),multi_grid();
    void assemble_grad_rho_p(),assemble_div_rho_u(),assemble_del2_u();

    if(getenv("STRICT_ALA_STAGE_F2_REQUIRED")==NULL) return 0;
    if(executed) return 1;
    executed=1;
    if(snapshot_dir==NULL || output_dir==NULL)
        myerror(E,"Incomplete strict-ALA Stage-F2 environment contract");
    if(E->parallel.nproc!=384 || E->parallel.nprocx!=4 ||
       E->parallel.nprocy!=4 || E->parallel.nprocz!=2)
        myerror(E,"F2_E2_REPLAY_LAYOUT_MISMATCH");
    if(!E->control.NMULTIGRID && !E->control.EMULTIGRID)
        myerror(E,"Stage-F2 requires velocity multigrid");
    if(!E->control.ala_element_vanka_smoother ||
       E->control.ala_augmented_lagrangian_gamma<=0.0)
        myerror(E,"Stage-F2 requires configured K_gamma element-Vanka");

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
            myerror(E,"Stage-F2 POD reconstruction failed");
    }
    selected=ala_f1a_selected_count(values,n,&selected_energy);
    if(selected<=0 || selected>ALA_F1A_MAX_MODES ||
       selected_energy<0.95-1.0e-12)
        myerror(E,"Stage-F2 frozen POD selection failed");
    for(i=0;i<n;i++) total+=values[i];

    rhs=ala_schur_alloc_velocity(E,lev);
    rhs_saved=ala_schur_alloc_velocity(E,lev);
    u_exact=ala_schur_alloc_velocity(E,lev);
    u=ala_schur_alloc_velocity(E,lev);
    velocity_work=ala_schur_alloc_velocity(E,lev);
    Ku=ala_schur_alloc_velocity(E,lev);
    error=ala_schur_alloc_velocity(E,lev);
    Kerror=ala_schur_alloc_velocity(E,lev);
    S_exact=ala_schur_alloc_pressure(E,lev);
    S_approx=ala_schur_alloc_pressure(E,lev);
    S_error=ala_schur_alloc_pressure(E,lev);
    if(E->parallel.me==0) {
        trajectory=ala_f1a_open_root(E,output_dir,
            "strict_ala_stage_F2_velocity_trajectory.csv","w");
        fprintf(trajectory,"variant,POD_mode,POD_energy_weight,MG_cycles,"
            "momentum_residual_relative,velocity_error_relative,"
            "K_energy_error_relative,Schur_action_error_relative,"
            "velocity_cosine,exact_tight_achieved,valid\n");
        tight=ala_f1a_open_root(E,output_dir,
            "strict_ala_stage_F2_tight_solves.csv","w");
        fprintf(tight,"call_id,POD_mode,rhs_norm,requested_relative_tolerance,"
            "target_residual,achieved_relative_residual,cycles,max_cycles,"
            "converged\n");
    }

    original_vanka=E->control.ala_element_vanka_smoother;
    for(j=0;j<selected;j++) {
        long long before=E->control.ala_stage_abc_inner_cycle_count;
        double qnorm;
        /* The preceding mode ends with the rollback trajectory.  Restore the
         * configured smoother before constructing every tight reference. */
        E->control.ala_element_vanka_smoother=original_vanka;
        q=ala_e2_build_mode(E,cases,indices,n,vectors,norms,j);
        qnorm=sqrt(max(global_pdot(E,q,q,lev),0.0));
        if(!isfinite(qnorm) || qnorm<=0.0)
            myerror(E,"Stage-F2 zero POD mode");
        ala_f1d_scale_pressure(E,q,1.0/qnorm,lev);
        assemble_grad_rho_p(E,q,rhs_saved,lev);
        ala_f2_copy_velocity(E,rhs,rhs_saved,lev);
        rhs2=global_vdot(E,rhs_saved,rhs_saved,lev);
        target=E->control.ala_schur_diagnostic_tight_tolerance*
            sqrt(max(rhs2,0.0)/max((double)E->mesh.neq,1.0));
        achieved=ala_schur_tight_solve(E,u_exact,rhs,velocity_work,lev,
            E->control.ala_schur_diagnostic_tight_tolerance);
        assemble_div_rho_u(E,u_exact,S_exact,lev);
        exact2=global_vdot(E,u_exact,u_exact,lev);
        exact_energy=global_vdot(E,u_exact,rhs_saved,lev);
        if(!isfinite(exact2) || exact2<=0.0 ||
           !isfinite(exact_energy) || exact_energy<=0.0)
            myerror(E,"Stage-F2 invalid exact K_gamma solution energy");
        call_id++;
        if(E->parallel.me==0) {
            fprintf(tight,"%d,%d,%.17e,%.17e,%.17e,%.17e,%lld,%d,1\n",
                call_id,j+1,sqrt(max(rhs2,0.0)),
                E->control.ala_schur_diagnostic_tight_tolerance,target,achieved,
                E->control.ala_stage_abc_inner_cycle_count-before,
                E->control.ala_schur_diagnostic_max_cycles);
            fflush(tight);
        }
        weight=values[j]/max(total,1.0e-300);
        for(variant=0;variant<2;variant++) {
            E->control.ala_element_vanka_smoother=variant==0 ?
                original_vanka : 0;
            ala_schur_zero_velocity(E,u,lev);
            ala_f2_copy_velocity(E,rhs,rhs_saved,lev);
            checkpoint=0;
            for(cycle=1;cycle<=checkpoints[ALA_F2_CHECKPOINTS-1];cycle++) {
                double r2,e2,kerr2,s2,serr2,unorm2,dot,cosine;
                multi_grid(E,u,rhs,target,lev);
                if(cycle!=checkpoints[checkpoint]) continue;
                assemble_del2_u(E,u,Ku,lev,1);
                for(m=1;m<=E->sphere.caps_per_proc;m++)
                    for(i=0;i<E->lmesh.NEQ[lev];i++) {
                        rhs[m][i]=rhs_saved[m][i]-Ku[m][i];
                        error[m][i]=u[m][i]-u_exact[m][i];
                    }
                assemble_del2_u(E,error,Kerror,lev,1);
                assemble_div_rho_u(E,u,S_approx,lev);
                for(m=1;m<=E->sphere.caps_per_proc;m++)
                    for(e=1;e<=E->lmesh.NPNO[lev];e++)
                        S_error[m][e]=S_approx[m][e]-S_exact[m][e];
                r2=global_vdot(E,rhs,rhs,lev);
                e2=global_vdot(E,error,error,lev);
                kerr2=global_vdot(E,error,Kerror,lev);
                s2=global_pdot(E,S_exact,S_exact,lev);
                serr2=global_pdot(E,S_error,S_error,lev);
                unorm2=global_vdot(E,u,u,lev);
                dot=global_vdot(E,u,u_exact,lev);
                cosine=dot/sqrt(max(unorm2*exact2,1.0e-300));
                if(!isfinite(r2) || !isfinite(e2) || !isfinite(kerr2) ||
                   !isfinite(serr2) || kerr2 < -1.0e-12*exact_energy)
                    myerror(E,"Stage-F2 non-finite velocity trajectory");
                if(E->parallel.me==0) {
                    fprintf(trajectory,"%s,%d,%.17e,%d,%.17e,%.17e,"
                        "%.17e,%.17e,%.17e,%.17e,1\n",variant_name[variant],
                        j+1,weight,cycle,sqrt(r2/max(rhs2,1.0e-300)),
                        sqrt(e2/max(exact2,1.0e-300)),
                        sqrt(max(kerr2,0.0)/max(exact_energy,1.0e-300)),
                        sqrt(serr2/max(s2,1.0e-300)),cosine,achieved);
                    fflush(trajectory);
                }
                checkpoint++;
            }
        }
        ala_schur_free_field(E,q); q=0;
    }
    E->control.ala_element_vanka_smoother=original_vanka;
    if(E->parallel.me==0) {
        summary=ala_f1a_open_root(E,output_dir,
            "strict_ala_stage_F2_runtime.json","w");
        fprintf(summary,"{\n  \"schema\": \"strict-ala-stage-F2-runtime-v1\",\n"
            "  \"selected_modes\": %d,\n  \"selected_energy\": %.17e,\n"
            "  \"variants\": [\"CONFIGURED_VANKA\", \"DIAGONAL_ROLLBACK\"],\n"
            "  \"maximum_MG_cycles\": 64,\n  \"tight_solve_count\": %d,\n"
            "  \"operator\": \"frozen_current_K_gamma\",\n"
            "  \"changed_object\": \"velocity_multigrid_smoother_only\",\n"
            "  \"production_default_change_authorized\": false\n}\n",
            selected,selected_energy,call_id);
        fclose(summary);
        fprintf(E->fp,"STRICT_ALA_STAGE_F2_AVV_REVIEW_V1 modes=%d "
            "tight_solves=%d maximum_MG_cycles=64\n",selected,call_id);
        fprintf(E->fp,"STRICT_ALA_STAGE_F2_COMPLETE candidate="
            "configured_vanka_vs_diagonal_rollback\n");
        fflush(E->fp);
    }
    if(trajectory) fclose(trajectory);
    if(tight) fclose(tight);
    ala_schur_free_field(E,rhs); ala_schur_free_field(E,rhs_saved);
    ala_schur_free_field(E,u_exact); ala_schur_free_field(E,u);
    ala_schur_free_field(E,velocity_work); ala_schur_free_field(E,Ku);
    ala_schur_free_field(E,error); ala_schur_free_field(E,Kerror);
    ala_schur_free_field(E,S_exact); ala_schur_free_field(E,S_approx);
    ala_schur_free_field(E,S_error);
    ala_e2_free_observer(E,e0); ala_e2_free_observer(E,e1);
    return 1;
}
