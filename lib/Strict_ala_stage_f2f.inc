/* Strict-ALA Stage F2f: frozen-mode overlapping L4 face-patch audit.
 *
 * F2e showed that disjoint two-element face blocks are substantially weaker
 * than configured element Vanka.  F2f changes exactly one structural
 * variable: every rank-local internal face in one orientation contributes a
 * two-element 36x36 K_gamma block, so adjacent blocks overlap by one element.
 * Directions are tested separately and only one cache exists at a time.
 * Lower levels, transfers, K_gamma, D+C, viscosity, sweeps, damping, and
 * production defaults remain frozen. */

#define ALA_F2F_VARIANTS 4

static int strict_ala_stage_f2f_run(struct All_variables *E,int lev)
{
    static int executed=0;
    const char *snapshot_dir=getenv("STRICT_ALA_STAGE_F2F_E2_SNAPSHOT_DIR");
    const char *output_dir=getenv("STRICT_ALA_STAGE_F2F_OUTPUT_DIR");
    const char *names[ALA_F2F_VARIANTS]={
        "CONFIGURED","L4_OVERLAP_FACE_Z","L4_OVERLAP_FACE_X",
        "L4_OVERLAP_FACE_Y"};
    const char *modes[ALA_F2F_VARIANTS]={
        "configured","overlap_face_z","overlap_face_x","overlap_face_y"};
    struct ala_e2_observer template_observer,*e0,*e1;
    struct ala_e2_observer *cases[2*ALA_E2_LATE_COUNT];
    int indices[2*ALA_E2_LATE_COUNT];
    double values[2*ALA_E2_LATE_COUNT];
    double vectors[4*ALA_E2_LATE_COUNT*ALA_E2_LATE_COUNT];
    double norms[2*ALA_E2_LATE_COUNT],total=0.0,selected_energy=0.0;
    double **q=0,**rhs=0,**rhs_saved=0,**u=0,**Ku=0;
    FILE *level_log=NULL,*trajectory=NULL,*mode_log=NULL,*performance=NULL;
    FILE *cache_log=NULL,*runtime=NULL;
    int i,j,v,m,n=2*ALA_E2_LATE_COUNT,selected,cycle;
    int configured_sweeps,sweeps[ALA_F2F_VARIANTS];
    double configured_damping,damping[ALA_F2F_VARIANTS];
    double rhs2,target,weight,residual2,qnorm;
    double variant_start,variant_seconds,variant_seconds_max;
    double global_vdot(),global_pdot(),multi_grid(),CPU_time0();
    void assemble_grad_rho_p(),assemble_del2_u();
    char mode_text[32],weight_text[64],cycle_text[32];

    if(getenv("STRICT_ALA_STAGE_F2F_REQUIRED")==NULL) return 0;
    if(executed) return 1;
    executed=1;
    if(snapshot_dir==NULL || output_dir==NULL)
        myerror(E,"Incomplete strict-ALA Stage-F2f environment contract");
    if(E->parallel.nproc!=384 || E->parallel.nprocx!=4 ||
       E->parallel.nprocy!=4 || E->parallel.nprocz!=2)
        myerror(E,"F2F_E2_REPLAY_LAYOUT_MISMATCH");
    if((!E->control.NMULTIGRID && !E->control.EMULTIGRID) ||
       !E->control.ala_element_vanka_smoother ||
       E->control.ala_augmented_lagrangian_gamma<=0.0)
        myerror(E,"Stage-F2f requires configured K_gamma element-Vanka MG");

    configured_sweeps=E->control.v_steps_high;
    configured_damping=E->control.ala_element_vanka_damping;
    if(configured_sweeps<1 || configured_damping<=0.0 ||
       configured_damping>1.0)
        myerror(E,"Invalid Stage-F2f configured smoother baseline");
    for(v=0;v<ALA_F2F_VARIANTS;v++) {
        sweeps[v]=configured_sweeps;
        damping[v]=configured_damping;
    }

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
            myerror(E,"Stage-F2f POD reconstruction failed");
    }
    selected=ala_f1a_selected_count(values,n,&selected_energy);
    if(selected<=0 || selected>ALA_F1A_MAX_MODES ||
       selected_energy<0.95-1.0e-12)
        myerror(E,"Stage-F2f frozen POD selection failed");
    for(i=0;i<n;i++) total+=values[i];

    rhs=ala_schur_alloc_velocity(E,lev);
    rhs_saved=ala_schur_alloc_velocity(E,lev);
    u=ala_schur_alloc_velocity(E,lev);
    Ku=ala_schur_alloc_velocity(E,lev);
    if(E->parallel.me==0) {
        level_log=ala_f1a_open_root(E,output_dir,
            "strict_ala_stage_F2f_level_transfer.csv","w");
        fprintf(level_log,"candidate,POD_mode,POD_energy_weight,MG_cycle,"
            "top_level,level,phase,input_rms,output_rms,reduction_ratio,"
            "alpha,valid\n");
        fclose(level_log); level_log=NULL;
        trajectory=ala_f1a_open_root(E,output_dir,
            "strict_ala_stage_F2f_velocity_trajectory.csv","w");
        fprintf(trajectory,"candidate,smoother_composition,finest_sweeps,"
            "finest_damping,POD_mode,POD_energy_weight,MG_cycles,"
            "momentum_residual_relative,valid\n");
        mode_log=ala_f1a_open_root(E,output_dir,
            "strict_ala_stage_F2f_mode_rhs.csv","w");
        fprintf(mode_log,"POD_mode,POD_energy_weight,rhs_norm,valid\n");
        performance=ala_f1a_open_root(E,output_dir,
            "strict_ala_stage_F2f_performance.csv","w");
        fprintf(performance,"candidate,smoother_block,finest_sweeps,"
            "finest_damping,block_dof,seconds_max,valid\n");
        cache_log=ala_f1a_open_root(E,output_dir,
            "strict_ala_stage_F2f_cache.csv","w");
        fprintf(cache_log,"orientation,patches_per_cap,block_dof,"
            "cache_mb_max,build_seconds_max,min_pivot_ratio,valid\n");
        fclose(cache_log); cache_log=NULL;
    }

    setenv("STRICT_ALA_STAGE_F2B_OUTPUT_DIR",output_dir,1);
    for(v=0;v<ALA_F2F_VARIANTS;v++) {
        setenv("STRICT_ALA_STAGE_F2F_ACTIVE","1",1);
        setenv("STRICT_ALA_STAGE_F2F_CANDIDATE",names[v],1);
        setenv("STRICT_ALA_STAGE_F2F_L4_MODE",modes[v],1);
        MPI_Barrier(E->parallel.world);
        variant_start=CPU_time0();
        if(E->parallel.me==0) {
            fprintf(E->fp,"STRICT_ALA_STAGE_F2F_CANDIDATE_BEGIN name=%s "
                "block=%s sweeps=%d damping=%.17e block_dof=%d\n",
                names[v],modes[v],sweeps[v],damping[v],v==0 ? 24 : 36);
            fflush(E->fp);
        }
        for(j=0;j<selected;j++) {
            q=ala_e2_build_mode(E,cases,indices,n,vectors,norms,j);
            qnorm=sqrt(max(global_pdot(E,q,q,lev),0.0));
            if(!isfinite(qnorm) || qnorm<=0.0)
                myerror(E,"Stage-F2f zero POD mode");
            ala_f1d_scale_pressure(E,q,1.0/qnorm,lev);
            assemble_grad_rho_p(E,q,rhs_saved,lev);
            ala_f2_copy_velocity(E,rhs,rhs_saved,lev);
            rhs2=global_vdot(E,rhs_saved,rhs_saved,lev);
            if(!isfinite(rhs2) || rhs2<=0.0)
                myerror(E,"Stage-F2f invalid mode RHS");
            target=E->control.ala_schur_diagnostic_tight_tolerance*
                sqrt(rhs2/max((double)E->mesh.neq,1.0));
            weight=values[j]/max(total,1.0e-300);
            if(E->parallel.me==0 && v==0) {
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
                    myerror(E,"Stage-F2f invalid checkpoint residual");
                if(E->parallel.me==0) {
                    fprintf(trajectory,
                        "%s,%s,%d,%.17e,%d,%.17e,%d,%.17e,1\n",
                        names[v],modes[v],sweeps[v],damping[v],j+1,weight,
                        cycle,sqrt(residual2/max(rhs2,1.0e-300)));
                    fflush(trajectory);
                }
            }
            ala_schur_free_field(E,q); q=0;
        }
        MPI_Barrier(E->parallel.world);
        variant_seconds=CPU_time0()-variant_start;
        MPI_Allreduce(&variant_seconds,&variant_seconds_max,1,MPI_DOUBLE,
                      MPI_MAX,E->parallel.world);
        if(E->parallel.me==0) {
            fprintf(performance,"%s,%s,%d,%.17e,%d,%.17e,1\n",
                names[v],modes[v],sweeps[v],damping[v],v==0 ? 24 : 36,
                variant_seconds_max);
            fflush(performance);
            fprintf(E->fp,"STRICT_ALA_STAGE_F2F_CANDIDATE_END name=%s "
                "seconds_max=%.17e\n",names[v],variant_seconds_max);
            fflush(E->fp);
        }
    }
    unsetenv("STRICT_ALA_STAGE_F2B_OBSERVE");
    unsetenv("STRICT_ALA_STAGE_F2B_MODE");
    unsetenv("STRICT_ALA_STAGE_F2B_WEIGHT");
    unsetenv("STRICT_ALA_STAGE_F2B_CYCLE");
    unsetenv("STRICT_ALA_STAGE_F2B_OUTPUT_DIR");
    unsetenv("STRICT_ALA_STAGE_F2F_ACTIVE");
    unsetenv("STRICT_ALA_STAGE_F2F_CANDIDATE");
    unsetenv("STRICT_ALA_STAGE_F2F_L4_MODE");
    if(E->parallel.me==0) {
        runtime=ala_f1a_open_root(E,output_dir,
            "strict_ala_stage_F2f_runtime.json","w");
        fprintf(runtime,
            "{\n  \"schema\": \"strict-ala-stage-F2f-runtime-v1\",\n"
            "  \"selected_modes\": %d,\n  \"selected_energy\": %.17e,\n"
            "  \"maximum_MG_cycles\": 64,\n"
            "  \"configured_finest_sweeps\": %d,\n"
            "  \"configured_finest_damping\": %.17e,\n"
            "  \"candidate_count\": %d,\n"
            "  \"scope\": \"finest_level_only\",\n"
            "  \"candidate_structure\": \"overlapping_two_element_face_pairs_all_local_internal_faces\",\n"
            "  \"cached_blocks\": \"one_orientation_at_a_time_36x36_K_gamma\",\n"
            "  \"face_orientations\": [\"z\", \"x\", \"y\"],\n"
            "  \"face_patch_dof\": 36,\n"
            "  \"cross_rank_faces\": false,\n"
            "  \"viscosity\": \"frozen_current\",\n"
            "  \"operator\": \"frozen_current_K_gamma_D_plus_C\",\n"
            "  \"production_default_change_authorized\": false\n}\n",
            selected,selected_energy,configured_sweeps,configured_damping,
            ALA_F2F_VARIANTS);
        fclose(runtime);
        fprintf(E->fp,
            "STRICT_ALA_STAGE_F2F_L4_OVERLAPPING_FACE_PATCH_V1 modes=%d "
            "candidates=%d maximum_MG_cycles=64\n",
            selected,ALA_F2F_VARIANTS);
        fprintf(E->fp,"STRICT_ALA_STAGE_F2F_COMPLETE\n");
        fflush(E->fp);
    }
    if(trajectory) fclose(trajectory);
    if(mode_log) fclose(mode_log);
    if(performance) fclose(performance);
    ala_schur_free_field(E,rhs); ala_schur_free_field(E,rhs_saved);
    ala_schur_free_field(E,u); ala_schur_free_field(E,Ku);
    ala_e2_free_observer(E,e0); ala_e2_free_observer(E,e1);
    return 1;
}
