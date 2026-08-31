/* Strict-ALA Stage F1d: residual-closed actual-mode coarse qualification.
 *
 * F1c showed that the three-mode Galerkin correction leaves a large residual
 * outside span(Q).  This opt-in diagnostic changes exactly one object: the
 * coarse basis is enlarged from Q to [Q,R], where each R column is the
 * orthonormalized coarse-only residual
 *
 *   r_j = q_j - S Q (Q^T S Q)^-1 Q^T q_j.
 *
 * The configured F1b fine map, complete B=D+C, K_gamma, frozen E2 modes, and
 * production defaults remain unchanged.  All new S actions use the same tight
 * frozen-state velocity solve as F1a--F1c.
 */

#define ALA_F1D_MAX_BASIS (2*ALA_F1A_MAX_MODES)
#define ALA_F1D_SYMMETRY_TOLERANCE 1.0e-8
#define ALA_F1D_PIVOT_TOLERANCE 1.0e-12
#define ALA_F1D_BASIS_TOLERANCE 1.0e-8

static int ala_f1d_factor(double *matrix,int n,double *minimum_ratio,
                          double *pivots,int *failure_step)
{
    int i,j,k;
    double diagonal_max=0.0,sum,pivot;
    *failure_step=-1;
    for(i=0;i<n;i++) pivots[i]=0.0;
    for(i=0;i<n;i++) diagonal_max=max(diagonal_max,matrix[i*n+i]);
    if(!isfinite(diagonal_max) || diagonal_max<=0.0) {
        *failure_step=0;
        return 0;
    }
    *minimum_ratio=1.0e300;
    for(i=0;i<n;i++)
        for(j=0;j<=i;j++) {
            sum=matrix[i*n+j];
            for(k=0;k<j;k++) sum-=matrix[i*n+k]*matrix[j*n+k];
            if(i==j) {
                pivot=sum;
                pivots[i]=pivot;
                if(!isfinite(pivot) ||
                   pivot<=ALA_F1D_PIVOT_TOLERANCE*diagonal_max) {
                    *failure_step=i;
                    return 0;
                }
                *minimum_ratio=min(*minimum_ratio,pivot/diagonal_max);
                matrix[i*n+j]=sqrt(pivot);
            }
            else matrix[i*n+j]=sum/matrix[j*n+j];
        }
    for(i=0;i<n;i++) for(j=i+1;j<n;j++) matrix[i*n+j]=0.0;
    return 1;
}

static int ala_f1d_popcount(int value)
{
    int count=0;
    while(value) { count+=value&1; value>>=1; }
    return count;
}

static double ala_f1d_submatrix_symmetry(const double *raw,
    const int *active,int n,int leading)
{
    int i,j;
    double defect=0.0,scale=0.0;
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            scale=max(scale,fabs(raw[active[i]*leading+active[j]]));
    for(i=0;i<n;i++)
        for(j=0;j<i;j++)
            defect=max(defect,fabs(raw[active[i]*leading+active[j]]-
                                  raw[active[j]*leading+active[i]])/
                               max(scale,1.0e-300));
    return defect;
}

static void ala_f1d_extract_symmetric(const double *symmetric,
    const int *active,int n,int leading,double *result)
{
    int i,j;
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            result[i*n+j]=symmetric[active[i]*leading+active[j]];
}

static void ala_f1d_solve(const double *factor,const double *rhs,double *x,int n)
{
    int i,j;
    double y[ALA_F1D_MAX_BASIS],sum;
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

static void ala_f1d_scale_pressure(struct All_variables *E,double **x,
                                   double scale,int lev)
{
    int m,e;
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(e=1;e<=E->lmesh.NPNO[lev];e++) x[m][e]*=scale;
}

static int strict_ala_stage_f1d_run(struct All_variables *E,
    struct ala_pressure_preconditioner_cache *cache,int lev)
{
    static int executed=0;
    const char *snapshot_dir=getenv("STRICT_ALA_STAGE_F1D_E2_SNAPSHOT_DIR");
    const char *output_dir=getenv("STRICT_ALA_STAGE_F1D_OUTPUT_DIR");
    struct ala_e2_observer template_observer,*e0,*e1;
    struct ala_e2_observer *cases[2*ALA_E2_LATE_COUNT];
    int indices[2*ALA_E2_LATE_COUNT];
    double values[2*ALA_E2_LATE_COUNT];
    double vectors[4*ALA_E2_LATE_COUNT*ALA_E2_LATE_COUNT];
    double norms[2*ALA_E2_LATE_COUNT],total=0.0,selected_energy=0.0;
    double **q[ALA_F1A_MAX_MODES],**Sq[ALA_F1A_MAX_MODES];
    double **basis[ALA_F1D_MAX_BASIS],**Sbasis[ALA_F1D_MAX_BASIS];
    double **z[ALA_F1A_MAX_MODES],**Sz[ALA_F1A_MAX_MODES];
    double **fine_rhs=0,**fine_z=0,**work=0,**aux=0;
    double **velocity_rhs=0,**velocity=0,**velocity_work=0;
    double original[ALA_F1A_MAX_MODES*ALA_F1A_MAX_MODES];
    double original_factor[ALA_F1A_MAX_MODES*ALA_F1A_MAX_MODES];
    double original_pivots[ALA_F1A_MAX_MODES];
    int original_failure_step=-1;
    double coarse_raw[ALA_F1D_MAX_BASIS*ALA_F1D_MAX_BASIS];
    double coarse[ALA_F1D_MAX_BASIS*ALA_F1D_MAX_BASIS];
    double factor[ALA_F1D_MAX_BASIS*ALA_F1D_MAX_BASIS];
    double full_factor[ALA_F1D_MAX_BASIS*ALA_F1D_MAX_BASIS];
    double full_pivots[ALA_F1D_MAX_BASIS],selected_pivots[ALA_F1D_MAX_BASIS];
    double eigen_work[ALA_F1D_MAX_BASIS*ALA_F1D_MAX_BASIS];
    double eigenvectors[ALA_F1D_MAX_BASIS*ALA_F1D_MAX_BASIS];
    double eigenvalues[ALA_F1D_MAX_BASIS];
    int active_index[ALA_F1D_MAX_BASIS];
    int basis_source[ALA_F1D_MAX_BASIS];
    FILE *metrics=NULL,*basis_log=NULL,*projected=NULL,*tight=NULL,*summary=NULL;
    FILE *matrix_log=NULL,*pivot_log=NULL,*spectrum_log=NULL,*subset_log=NULL;
    int i,j,m,e,pass,n=2*ALA_E2_LATE_COUNT,selected,coarse_dim,accepted=0;
    int call_id=0,full_factor_pass=0,full_failure_step=-1;
    int selected_factor_pass=0,selected_failure_step=-1;
    int selected_mask=0,retained=0,active_dim=0,mask,target_count;
    double original_symmetry=0.0,coarse_symmetry=0.0,scale=0.0;
    double original_pivot=0.0,full_pivot=0.0,selected_pivot=0.0,weight;
    double selected_symmetry=0.0,eigen_min=0.0,eigen_max=0.0,condition=0.0;
    double full_diagonal_scale=0.0,selected_diagonal_scale=0.0;

    if(getenv("STRICT_ALA_STAGE_F1D_REQUIRED")==NULL) return 0;
    if(executed) return 1;
    executed=1;
    if(snapshot_dir==NULL || output_dir==NULL)
        myerror(E,"Incomplete strict-ALA Stage-F1d environment contract");
    if(E->parallel.nproc!=384 || E->parallel.nprocx!=4 ||
       E->parallel.nprocy!=4 || E->parallel.nprocz!=2)
        myerror(E,"F1D_E2_REPLAY_LAYOUT_MISMATCH");
    if(!E->control.ala_shallow_patch_preconditioner ||
       strcmp(E->control.ala_shallow_patch_velocity_solver,"element_vanka")!=0 ||
       strcmp(E->control.ala_shallow_patch_local_operator,
              "operator_consistent")!=0)
        myerror(E,"Stage-F1d requires the validated F1b pressure map");

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
            myerror(E,"Stage-F1d POD reconstruction failed");
    }
    selected=ala_f1a_selected_count(values,n,&selected_energy);
    if(selected<=0 || selected>ALA_F1A_MAX_MODES || selected_energy<0.95-1.0e-12)
        myerror(E,"Stage-F1d frozen POD selection failed");
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
            "strict_ala_stage_F1d_true_mode.csv","w");
        fprintf(metrics,"POD_mode,POD_energy_weight,E_P,cosine,qTPq,qTSPq,"
                "original_projection_max,enriched_projection_max,"
                "tight_solve_achieved,valid\n");
        basis_log=ala_f1a_open_root(E,output_dir,
            "strict_ala_stage_F1d_basis.csv","w");
        fprintf(basis_log,"basis_index,basis_type,source_mode,"
                "preorthogonal_norm,postorthogonal_norm,accepted,"
                "maximum_orthogonality,S_energy,valid\n");
        projected=ala_f1a_open_root(E,output_dir,
            "strict_ala_stage_F1d_projected_matrices.csv","w");
        fprintf(projected,"matrix,row_mode,column_mode,value\n");
        tight=ala_f1a_open_root(E,output_dir,
            "strict_ala_stage_F1d_tight_solves.csv","w");
        fprintf(tight,"call_id,role,POD_mode,rhs_norm,requested_relative_tolerance,"
                "target_residual,achieved_relative_residual,cycles,max_cycles,converged\n");
        matrix_log=ala_f1a_open_root(E,output_dir,
            "strict_ala_stage_F1d_coarse_matrix.csv","w");
        fprintf(matrix_log,"scope,row_basis,column_basis,value\n");
        pivot_log=ala_f1a_open_root(E,output_dir,
            "strict_ala_stage_F1d_cholesky_pivots.csv","w");
        fprintf(pivot_log,"scope,step,pivot,pivot_ratio,evaluated,passed\n");
        spectrum_log=ala_f1a_open_root(E,output_dir,
            "strict_ala_stage_F1d_coarse_spectrum.csv","w");
        fprintf(spectrum_log,"eigen_index,eigenvalue,relative_to_max,positive\n");
        subset_log=ala_f1a_open_root(E,output_dir,
            "strict_ala_stage_F1d_selected_basis.csv","w");
        fprintf(subset_log,"selected_index,full_basis_index,basis_type,"
                "source_mode,retained\n");
    }

    memset(original,0,sizeof(original));
    for(i=0;i<selected;i++) {
        double qnorm;
        q[i]=ala_e2_build_mode(E,cases,indices,n,vectors,norms,i);
        Sq[i]=ala_schur_alloc_pressure(E,lev);
        z[i]=ala_schur_alloc_pressure(E,lev);
        Sz[i]=ala_schur_alloc_pressure(E,lev);
        ala_f1a_true_action(E,q[i],Sq[i],velocity_rhs,velocity,velocity_work,
                           lev,tight,&call_id,i+1,"ORIGINAL_BASIS_SQ");
        qnorm=sqrt(max(global_pdot(E,q[i],q[i],lev),0.0));
        if(!isfinite(qnorm) || qnorm<=0.0)
            myerror(E,"Stage-F1d zero original basis vector");
        ala_f1d_scale_pressure(E,q[i],1.0/qnorm,lev);
        ala_f1d_scale_pressure(E,Sq[i],1.0/qnorm,lev);
        basis[i]=q[i]; Sbasis[i]=Sq[i]; basis_source[i]=i+1;
    }
    for(i=0;i<selected;i++)
        for(j=0;j<selected;j++) {
            original[i*selected+j]=global_pdot(E,q[i],Sq[j],lev);
            scale=max(scale,fabs(original[i*selected+j]));
        }
    for(i=0;i<selected;i++)
        for(j=0;j<i;j++) {
            original_symmetry=max(original_symmetry,
                fabs(original[i*selected+j]-original[j*selected+i])/
                max(scale,1.0e-300));
            original[i*selected+j]=original[j*selected+i]=
                0.5*(original[i*selected+j]+original[j*selected+i]);
        }
    memcpy(original_factor,original,(size_t)selected*selected*sizeof(double));
    if(original_symmetry>ALA_F1D_SYMMETRY_TOLERANCE ||
       !ala_f1d_factor(original_factor,selected,&original_pivot,
                       original_pivots,&original_failure_step))
        myerror(E,"Stage-F1d original coarse matrix failed SPD contract");

    coarse_dim=selected;
    for(j=0;j<selected;j++) {
        double rhs[ALA_F1A_MAX_MODES]={0},c[ALA_F1A_MAX_MODES]={0};
        double pre2,post2,postnorm,max_orth=0.0,senergy=0.0;
        int accepted_this=0,basis_index=0;
        for(i=0;i<selected;i++) rhs[i]=global_pdot(E,q[i],q[j],lev);
        ala_f1d_solve(original_factor,rhs,c,selected);
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=1;e<=E->lmesh.NPNO[lev];e++) {
                double value=q[j][m][e];
                for(i=0;i<selected;i++) value-=Sq[i][m][e]*c[i];
                fine_rhs[m][e]=value;
            }
        pre2=global_pdot(E,fine_rhs,fine_rhs,lev);
        for(pass=0;pass<2;pass++)
            for(i=0;i<coarse_dim;i++) {
                double projection=global_pdot(E,basis[i],fine_rhs,lev);
                for(m=1;m<=E->sphere.caps_per_proc;m++)
                    for(e=1;e<=E->lmesh.NPNO[lev];e++)
                        fine_rhs[m][e]-=projection*basis[i][m][e];
            }
        post2=global_pdot(E,fine_rhs,fine_rhs,lev);
        postnorm=sqrt(max(post2,0.0));
        if(isfinite(postnorm) && postnorm>
           ALA_F1D_BASIS_TOLERANCE*sqrt(max(pre2,1.0e-300)) &&
           coarse_dim<ALA_F1D_MAX_BASIS) {
            basis_index=coarse_dim+1;
            basis[coarse_dim]=ala_schur_alloc_pressure(E,lev);
            Sbasis[coarse_dim]=ala_schur_alloc_pressure(E,lev);
            basis_source[coarse_dim]=j+1;
            ala_f1a_copy_pressure(E,basis[coarse_dim],fine_rhs,lev);
            ala_f1d_scale_pressure(E,basis[coarse_dim],1.0/postnorm,lev);
            for(i=0;i<coarse_dim;i++)
                max_orth=max(max_orth,fabs(global_pdot(
                    E,basis[i],basis[coarse_dim],lev)));
            ala_f1a_true_action(E,basis[coarse_dim],Sbasis[coarse_dim],
                velocity_rhs,velocity,velocity_work,lev,tight,&call_id,j+1,
                "RESIDUAL_CLOSURE_SR");
            senergy=global_pdot(E,basis[coarse_dim],Sbasis[coarse_dim],lev);
            if(!isfinite(senergy) || senergy<=0.0 ||
               max_orth>ALA_F1D_BASIS_TOLERANCE)
                myerror(E,"Stage-F1d residual basis validation failed");
            coarse_dim++; accepted++; accepted_this=1;
        }
        if(E->parallel.me==0) {
            fprintf(basis_log,"%d,coarse_residual,%d,%.17e,%.17e,%d,"
                    "%.17e,%.17e,1\n",basis_index,j+1,
                    sqrt(max(pre2,0.0)),postnorm,accepted_this,max_orth,senergy);
            fflush(basis_log);
        }
    }
    memset(coarse_raw,0,sizeof(coarse_raw));
    memset(coarse,0,sizeof(coarse)); scale=0.0;
    for(i=0;i<coarse_dim;i++)
        for(j=0;j<coarse_dim;j++) {
            coarse_raw[i*coarse_dim+j]=
                global_pdot(E,basis[i],Sbasis[j],lev);
            coarse[i*coarse_dim+j]=coarse_raw[i*coarse_dim+j];
            scale=max(scale,fabs(coarse_raw[i*coarse_dim+j]));
        }
    for(i=0;i<coarse_dim;i++)
        for(j=0;j<i;j++) {
            coarse_symmetry=max(coarse_symmetry,
                fabs(coarse[i*coarse_dim+j]-coarse[j*coarse_dim+i])/
                max(scale,1.0e-300));
            coarse[i*coarse_dim+j]=coarse[j*coarse_dim+i]=
                0.5*(coarse[i*coarse_dim+j]+coarse[j*coarse_dim+i]);
        }
    if(E->parallel.me==0) {
        for(i=0;i<coarse_dim;i++)
            for(j=0;j<coarse_dim;j++) {
                fprintf(matrix_log,"FULL_RAW,%d,%d,%.17e\n",i+1,j+1,
                        coarse_raw[i*coarse_dim+j]);
                fprintf(matrix_log,"FULL_SYMMETRIZED,%d,%d,%.17e\n",i+1,j+1,
                        coarse[i*coarse_dim+j]);
            }
        fflush(matrix_log);
    }

    memcpy(full_factor,coarse,(size_t)coarse_dim*coarse_dim*sizeof(double));
    for(i=0;i<coarse_dim;i++)
        full_diagonal_scale=max(full_diagonal_scale,
                                coarse[i*coarse_dim+i]);
    full_factor_pass=ala_f1d_factor(full_factor,coarse_dim,&full_pivot,
                                    full_pivots,&full_failure_step);
    if(scale>0.0 && isfinite(scale)) {
        for(i=0;i<coarse_dim*coarse_dim;i++) eigen_work[i]=coarse[i]/scale;
        ala_geneo_jacobi_eigensolve(eigen_work,eigenvectors,eigenvalues,
                                    coarse_dim);
        for(i=0;i<coarse_dim;i++) eigenvalues[i]*=scale;
        for(i=0;i<coarse_dim;i++)
            for(j=i+1;j<coarse_dim;j++)
                if(eigenvalues[j]<eigenvalues[i]) {
                    double swap=eigenvalues[i];
                    eigenvalues[i]=eigenvalues[j]; eigenvalues[j]=swap;
                }
        eigen_min=eigenvalues[0]; eigen_max=eigenvalues[coarse_dim-1];
        condition=eigen_min>0.0 ? eigen_max/eigen_min : 0.0;
    }
    if(E->parallel.me==0) {
        for(i=0;i<coarse_dim;i++) {
            int evaluated=full_factor_pass || i<=full_failure_step;
            int pivot_pass=evaluated &&
                (full_factor_pass || i<full_failure_step);
            fprintf(pivot_log,"FULL,%d,%.17e,%.17e,%d,%d\n",i+1,
                    full_pivots[i],
                    full_pivots[i]/max(full_diagonal_scale,1.0e-300),
                    evaluated,pivot_pass);
            fprintf(spectrum_log,"%d,%.17e,%.17e,%d\n",i+1,
                    eigenvalues[i],eigenvalues[i]/max(fabs(eigen_max),1.0e-300),
                    eigenvalues[i]>0.0);
        }
        fflush(pivot_log); fflush(spectrum_log);
    }

    /* Retain all original POD modes and choose the largest deterministic
     * residual subset whose raw Galerkin defect and symmetrized Cholesky
     * factor both satisfy the frozen contract.  A full-space failure is a
     * scientific result, not an MPI-fatal instrumentation failure. */
    selected_mask=-1;
    for(target_count=accepted;target_count>=0 && selected_mask<0;target_count--)
        for(mask=0;mask<(1<<accepted);mask++) {
            if(ala_f1d_popcount(mask)!=target_count) continue;
            active_dim=0;
            for(i=0;i<selected;i++) active_index[active_dim++]=i;
            for(i=0;i<accepted;i++)
                if(mask&(1<<i)) active_index[active_dim++]=selected+i;
            selected_symmetry=ala_f1d_submatrix_symmetry(coarse_raw,
                active_index,active_dim,coarse_dim);
            ala_f1d_extract_symmetric(coarse,active_index,active_dim,
                                     coarse_dim,factor);
            selected_factor_pass=ala_f1d_factor(factor,active_dim,
                &selected_pivot,selected_pivots,&selected_failure_step);
            if(selected_symmetry<=ALA_F1D_SYMMETRY_TOLERANCE &&
               selected_factor_pass) { selected_mask=mask; break; }
        }
    if(selected_mask<0)
        myerror(E,"Stage-F1d original POD subspace lost SPD contract");
    retained=ala_f1d_popcount(selected_mask);
    for(i=0;i<active_dim;i++)
        selected_diagonal_scale=max(selected_diagonal_scale,
            coarse[active_index[i]*coarse_dim+active_index[i]]);

    if(E->parallel.me==0) {
        for(i=0;i<active_dim;i++) {
            int full_i=active_index[i];
            fprintf(subset_log,"%d,%d,%s,%d,1\n",i+1,full_i+1,
                    full_i<selected ? "POD" : "coarse_residual",
                    basis_source[full_i]);
            fprintf(pivot_log,"SELECTED,%d,%.17e,%.17e,1,1\n",i+1,
                    selected_pivots[i],
                    selected_pivots[i]/max(selected_diagonal_scale,1.0e-300));
            for(j=0;j<active_dim;j++) {
                int full_j=active_index[j];
                fprintf(matrix_log,"SELECTED_RAW,%d,%d,%.17e\n",i+1,j+1,
                        coarse_raw[full_i*coarse_dim+full_j]);
                fprintf(matrix_log,"SELECTED_SYMMETRIZED,%d,%d,%.17e\n",
                        i+1,j+1,coarse[full_i*coarse_dim+full_j]);
            }
        }
        fflush(subset_log); fflush(pivot_log); fflush(matrix_log);
    }

    for(j=0;j<selected;j++) {
        double rhs_c[ALA_F1D_MAX_BASIS]={0},c[ALA_F1D_MAX_BASIS]={0};
        double rhs_d[ALA_F1D_MAX_BASIS]={0},d[ALA_F1D_MAX_BASIS]={0};
        double qq,yy,qy,qz,ep,cosine,achieved;
        double original_projection=0.0,enriched_projection=0.0;
        for(i=0;i<active_dim;i++)
            rhs_c[i]=global_pdot(E,basis[active_index[i]],q[j],lev);
        ala_f1d_solve(factor,rhs_c,c,active_dim);
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=1;e<=E->lmesh.NPNO[lev];e++) {
                double value=q[j][m][e];
                for(i=0;i<active_dim;i++)
                    value-=Sbasis[active_index[i]][m][e]*c[i];
                fine_rhs[m][e]=value;
            }
        ala_schur_preconditioner_action(E,fine_rhs,fine_z,work,aux,lev,3,cache);
        for(i=0;i<active_dim;i++)
            rhs_d[i]=global_pdot(E,Sbasis[active_index[i]],fine_z,lev);
        ala_f1d_solve(factor,rhs_d,d,active_dim);
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=1;e<=E->lmesh.NPNO[lev];e++) {
                double value=fine_z[m][e];
                for(i=0;i<active_dim;i++)
                    value+=basis[active_index[i]][m][e]*(c[i]-d[i]);
                z[j][m][e]=value;
            }
        achieved=ala_f1a_true_action(E,z[j],Sz[j],velocity_rhs,velocity,
            velocity_work,lev,tight,&call_id,j+1,"RESIDUAL_CLOSED_CANDIDATE");
        qq=global_pdot(E,q[j],q[j],lev);
        yy=global_pdot(E,Sz[j],Sz[j],lev);
        qy=global_pdot(E,q[j],Sz[j],lev);
        qz=global_pdot(E,q[j],z[j],lev);
        ep=sqrt(max(qq+yy-2.0*qy,0.0)/max(qq,1.0e-300));
        cosine=qy/sqrt(max(qq*yy,1.0e-300));
        for(i=0;i<selected;i++) {
            double residual=global_pdot(E,q[i],q[j],lev)-
                            global_pdot(E,q[i],Sz[j],lev);
            original_projection=max(original_projection,fabs(residual));
        }
        for(i=0;i<active_dim;i++) {
            double residual=global_pdot(E,basis[active_index[i]],q[j],lev)-
                global_pdot(E,basis[active_index[i]],Sz[j],lev);
            enriched_projection=max(enriched_projection,fabs(residual));
        }
        weight=values[j]/max(total,1.0e-300);
        if(E->parallel.me==0) {
            fprintf(metrics,"%d,%.17e,%.17e,%.17e,%.17e,%.17e,"
                    "%.17e,%.17e,%.17e,1\n",j+1,weight,ep,cosine,qz,qy,
                    original_projection,enriched_projection,achieved);
            fflush(metrics);
        }
    }
    for(i=0;i<active_dim;i++)
        for(j=0;j<active_dim;j++)
            if(E->parallel.me==0)
                fprintf(projected,"T_enriched,%d,%d,%.17e\n",i+1,j+1,
                        coarse[active_index[i]*coarse_dim+active_index[j]]);
    for(i=0;i<selected;i++)
        for(j=0;j<selected;j++) {
            double map=global_pdot(E,q[i],z[j],lev);
            double action=global_pdot(E,q[i],Sz[j],lev);
            if(E->parallel.me==0) {
                fprintf(projected,"M_balanced,%d,%d,%.17e\n",i+1,j+1,map);
                fprintf(projected,"H_balanced,%d,%d,%.17e\n",i+1,j+1,action);
            }
        }
    if(E->parallel.me==0) {
        fflush(projected);
        summary=ala_f1a_open_root(E,output_dir,
            "strict_ala_stage_F1d_runtime.json","w");
        fprintf(summary,"{\n  \"schema\": \"strict-ala-stage-F1d-runtime-v2\",\n"
            "  \"selected_modes\": %d,\n  \"selected_energy\": %.17e,\n"
            "  \"original_basis_dimension\": %d,\n"
            "  \"accepted_residual_modes\": %d,\n"
            "  \"retained_residual_modes\": %d,\n"
            "  \"full_enriched_basis_dimension\": %d,\n"
            "  \"enriched_basis_dimension\": %d,\n"
            "  \"original_matrix_symmetry_defect\": %.17e,\n"
            "  \"original_matrix_minimum_pivot_ratio\": %.17e,\n"
            "  \"full_matrix_symmetry_defect\": %.17e,\n"
            "  \"full_matrix_cholesky_pass\": %s,\n"
            "  \"full_matrix_failure_step\": %d,\n"
            "  \"full_matrix_minimum_pivot_ratio\": %.17e,\n"
            "  \"full_matrix_eigenvalue_minimum\": %.17e,\n"
            "  \"full_matrix_eigenvalue_maximum\": %.17e,\n"
            "  \"full_matrix_spectral_condition\": %.17e,\n"
            "  \"selected_residual_mask\": %d,\n"
            "  \"enriched_matrix_symmetry_defect\": %.17e,\n"
            "  \"enriched_matrix_minimum_pivot_ratio\": %.17e,\n"
            "  \"structural_fallback\": %s,\n"
            "  \"tight_solve_count\": %d,\n"
            "  \"candidate\": \"maximal_spd_residual_closed_balanced_actual_mode_coarse\",\n"
            "  \"fine_map\": \"F1b_operator_consistent_frozen\",\n"
            "  \"production_default_change_authorized\": false\n}\n",
            selected,selected_energy,selected,accepted,retained,coarse_dim,
            active_dim,original_symmetry,original_pivot,coarse_symmetry,
            full_factor_pass ? "true" : "false",full_failure_step,
            full_pivot,eigen_min,eigen_max,condition,selected_mask,
            selected_symmetry,selected_pivot,
            retained<accepted ? "true" : "false",call_id);
        fclose(summary);
        fprintf(E->fp,"STRICT_ALA_STAGE_F1D_MAXIMAL_SPD_SUBSPACE_V2 "
                "modes=%d accepted=%d retained=%d full_dim=%d selected_dim=%d "
                "tight_solves=%d\n",selected,accepted,retained,coarse_dim,
                active_dim,call_id);
        fprintf(E->fp,"STRICT_ALA_STAGE_F1D_COMPLETE candidate="
                "maximal_spd_residual_closed_balanced_actual_mode_coarse\n");
        fflush(E->fp);
    }
    if(metrics) fclose(metrics); if(basis_log) fclose(basis_log);
    if(projected) fclose(projected); if(tight) fclose(tight);
    if(matrix_log) fclose(matrix_log); if(pivot_log) fclose(pivot_log);
    if(spectrum_log) fclose(spectrum_log); if(subset_log) fclose(subset_log);
    for(i=selected;i<coarse_dim;i++) {
        ala_schur_free_field(E,basis[i]); ala_schur_free_field(E,Sbasis[i]);
    }
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
