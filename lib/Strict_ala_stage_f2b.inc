/* Strict-ALA Stage F2b: actual-mode velocity-MG level localization.
 *
 * F2 proved that element Vanka is stable and essential, while its late tight
 * convergence remains slow.  F2b keeps the configured Vanka K_gamma solve
 * fixed and observes selected cycles inside multi_grid(): down smoothing,
 * restriction, bottom solve, and upward correction. */

#define ALA_F2B_CHECKPOINTS 7

static int ala_f2b_is_checkpoint(int cycle)
{
    const int checkpoints[ALA_F2B_CHECKPOINTS]={1,2,4,8,16,32,64};
    int i;
    for(i=0;i<ALA_F2B_CHECKPOINTS;i++)
        if(cycle==checkpoints[i]) return 1;
    return 0;
}

static int strict_ala_stage_f2b_run(struct All_variables *E,int lev)
{
    static int executed=0;
    const char *snapshot_dir=getenv("STRICT_ALA_STAGE_F2B_E2_SNAPSHOT_DIR");
    const char *output_dir=getenv("STRICT_ALA_STAGE_F2B_OUTPUT_DIR");
    struct ala_e2_observer template_observer,*e0,*e1;
    struct ala_e2_observer *cases[2*ALA_E2_LATE_COUNT];
    int indices[2*ALA_E2_LATE_COUNT];
    double values[2*ALA_E2_LATE_COUNT];
    double vectors[4*ALA_E2_LATE_COUNT*ALA_E2_LATE_COUNT];
    double norms[2*ALA_E2_LATE_COUNT],total=0.0,selected_energy=0.0;
    double **q=0,**rhs=0,**rhs_saved=0,**u=0,**Ku=0;
    FILE *level_log=NULL,*trajectory=NULL,*mode_log=NULL,*summary=NULL;
    int i,j,m,n=2*ALA_E2_LATE_COUNT,selected,cycle;
    double rhs2,target,weight,residual2,qnorm;
    double global_vdot(),global_pdot(),multi_grid();
    void assemble_grad_rho_p(),assemble_del2_u();
    char mode_text[32],weight_text[64],cycle_text[32];

    if(getenv("STRICT_ALA_STAGE_F2B_REQUIRED")==NULL) return 0;
    if(executed) return 1;
    executed=1;
    if(snapshot_dir==NULL || output_dir==NULL)
        myerror(E,"Incomplete strict-ALA Stage-F2b environment contract");
    if(E->parallel.nproc!=384 || E->parallel.nprocx!=4 ||
       E->parallel.nprocy!=4 || E->parallel.nprocz!=2)
        myerror(E,"F2B_E2_REPLAY_LAYOUT_MISMATCH");
    if((!E->control.NMULTIGRID && !E->control.EMULTIGRID) ||
       !E->control.ala_element_vanka_smoother ||
       E->control.ala_augmented_lagrangian_gamma<=0.0)
        myerror(E,"Stage-F2b requires configured K_gamma element-Vanka MG");

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
            myerror(E,"Stage-F2b POD reconstruction failed");
    }
    selected=ala_f1a_selected_count(values,n,&selected_energy);
    if(selected<=0 || selected>ALA_F1A_MAX_MODES ||
       selected_energy<0.95-1.0e-12)
        myerror(E,"Stage-F2b frozen POD selection failed");
    for(i=0;i<n;i++) total+=values[i];

    rhs=ala_schur_alloc_velocity(E,lev);
    rhs_saved=ala_schur_alloc_velocity(E,lev);
    u=ala_schur_alloc_velocity(E,lev);
    Ku=ala_schur_alloc_velocity(E,lev);
    if(E->parallel.me==0) {
        level_log=ala_f1a_open_root(E,output_dir,
            "strict_ala_stage_F2b_level_transfer.csv","w");
        fprintf(level_log,"POD_mode,POD_energy_weight,MG_cycle,top_level,"
            "level,phase,input_rms,output_rms,reduction_ratio,alpha,valid\n");
        fclose(level_log); level_log=NULL;
        trajectory=ala_f1a_open_root(E,output_dir,
            "strict_ala_stage_F2b_velocity_trajectory.csv","w");
        fprintf(trajectory,"POD_mode,POD_energy_weight,MG_cycles,"
            "momentum_residual_relative,valid\n");
        mode_log=ala_f1a_open_root(E,output_dir,
            "strict_ala_stage_F2b_mode_rhs.csv","w");
        fprintf(mode_log,"POD_mode,POD_energy_weight,rhs_norm,valid\n");
    }

    setenv("STRICT_ALA_STAGE_F2B_OUTPUT_DIR",output_dir,1);
    for(j=0;j<selected;j++) {
        q=ala_e2_build_mode(E,cases,indices,n,vectors,norms,j);
        qnorm=sqrt(max(global_pdot(E,q,q,lev),0.0));
        if(!isfinite(qnorm) || qnorm<=0.0)
            myerror(E,"Stage-F2b zero POD mode");
        ala_f1d_scale_pressure(E,q,1.0/qnorm,lev);
        assemble_grad_rho_p(E,q,rhs_saved,lev);
        ala_f2_copy_velocity(E,rhs,rhs_saved,lev);
        rhs2=global_vdot(E,rhs_saved,rhs_saved,lev);
        if(!isfinite(rhs2) || rhs2<=0.0)
            myerror(E,"Stage-F2b invalid mode RHS");
        target=E->control.ala_schur_diagnostic_tight_tolerance*
            sqrt(rhs2/max((double)E->mesh.neq,1.0));
        weight=values[j]/max(total,1.0e-300);
        if(E->parallel.me==0) {
            fprintf(mode_log,"%d,%.17e,%.17e,1\n",j+1,weight,sqrt(rhs2));
            fflush(mode_log);
        }
        snprintf(mode_text,sizeof(mode_text),"%d",j+1);
        snprintf(weight_text,sizeof(weight_text),"%.17e",weight);
        setenv("STRICT_ALA_STAGE_F2B_MODE",mode_text,1);
        setenv("STRICT_ALA_STAGE_F2B_WEIGHT",weight_text,1);
        ala_schur_zero_velocity(E,u,lev);
        for(cycle=1;cycle<=64;cycle++) {
            if(ala_f2b_is_checkpoint(cycle)) {
                snprintf(cycle_text,sizeof(cycle_text),"%d",cycle);
                setenv("STRICT_ALA_STAGE_F2B_CYCLE",cycle_text,1);
                setenv("STRICT_ALA_STAGE_F2B_OBSERVE","1",1);
            }
            else unsetenv("STRICT_ALA_STAGE_F2B_OBSERVE");
            multi_grid(E,u,rhs,target,lev);
            unsetenv("STRICT_ALA_STAGE_F2B_OBSERVE");
            if(!ala_f2b_is_checkpoint(cycle)) continue;
            assemble_del2_u(E,u,Ku,lev,1);
            for(m=1;m<=E->sphere.caps_per_proc;m++)
                for(i=0;i<E->lmesh.NEQ[lev];i++)
                    rhs[m][i]=rhs_saved[m][i]-Ku[m][i];
            residual2=global_vdot(E,rhs,rhs,lev);
            if(!isfinite(residual2) || residual2<0.0)
                myerror(E,"Stage-F2b invalid checkpoint residual");
            if(E->parallel.me==0) {
                fprintf(trajectory,"%d,%.17e,%d,%.17e,1\n",j+1,weight,
                    cycle,sqrt(residual2/max(rhs2,1.0e-300)));
                fflush(trajectory);
            }
        }
        ala_schur_free_field(E,q); q=0;
    }
    unsetenv("STRICT_ALA_STAGE_F2B_OBSERVE");
    unsetenv("STRICT_ALA_STAGE_F2B_MODE");
    unsetenv("STRICT_ALA_STAGE_F2B_WEIGHT");
    unsetenv("STRICT_ALA_STAGE_F2B_CYCLE");
    unsetenv("STRICT_ALA_STAGE_F2B_OUTPUT_DIR");
    if(E->parallel.me==0) {
        summary=ala_f1a_open_root(E,output_dir,
            "strict_ala_stage_F2b_runtime.json","w");
        fprintf(summary,"{\n  \"schema\": \"strict-ala-stage-F2b-runtime-v1\",\n"
            "  \"selected_modes\": %d,\n  \"selected_energy\": %.17e,\n"
            "  \"maximum_MG_cycles\": 64,\n"
            "  \"checkpoints\": [1, 2, 4, 8, 16, 32, 64],\n"
            "  \"operator\": \"frozen_current_K_gamma\",\n"
            "  \"velocity_smoother\": \"configured_element_vanka\",\n"
            "  \"production_default_change_authorized\": false\n}\n",
            selected,selected_energy);
        fclose(summary);
        fprintf(E->fp,"STRICT_ALA_STAGE_F2B_LEVEL_LOCALIZATION_V1 modes=%d "
            "maximum_MG_cycles=64\n",selected);
        fprintf(E->fp,"STRICT_ALA_STAGE_F2B_COMPLETE candidate="
            "configured_vanka_level_localization\n");
        fflush(E->fp);
    }
    if(trajectory) fclose(trajectory);
    if(mode_log) fclose(mode_log);
    ala_schur_free_field(E,rhs); ala_schur_free_field(E,rhs_saved);
    ala_schur_free_field(E,u); ala_schur_free_field(E,Ku);
    ala_e2_free_observer(E,e0); ala_e2_free_observer(E,e1);
    return 1;
}
