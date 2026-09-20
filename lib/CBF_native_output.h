/* Native Q1 CBF data: no interpolation or NetCDF dependency. */
#ifndef CITCOMS_CBF_OUTPUT_H
#define CITCOMS_CBF_OUTPUT_H
#include <errno.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>

static void CBF_io_check(struct All_variables *E, int bad)
{
    int global_bad;
    void parallel_process_termination();
    MPI_Allreduce(&bad,&global_bad,1,MPI_INT,MPI_MAX,E->parallel.world);
    if(global_bad) {
        if(E->parallel.me==0) fprintf(stderr,"CBF native output failed; step is incomplete\n");
        parallel_process_termination();
    }
}

static void CBF_native_boundary(struct All_variables *E,int top,
        double *rhs[NCS],double *mass[NCS],double *q[NCS],double totals[2])
{
    int m,i,node,e,a,d,bad=0,step=E->monitor.solution_cycles;
    int active=E->parallel.me_loc[3]==(top ? E->parallel.nprocz-1:0);
    int lev=E->mesh.levmax,elz=E->lmesh.elz,side=top ? SIDE_TOP:SIDE_BOTTOM;
    double length=E->data.radius_km*1000.,x[4][3],dm[4];
    char path[512],tmp[520];
    FILE *fp=NULL;
    /* data_dir already contains this rank's expanded cfg datadir. */
    if(active) {
        snprintf(path,sizeof(path),"%s/q.%s.%d.%d",E->control.data_dir,
                 top ? "surf":"botm",E->parallel.me,step);
        snprintf(tmp,sizeof(tmp),"%s.tmp",path);
        fp=fopen(tmp,"w");
        if(!fp) bad=1;
        else {
            fprintf(fp,"# CBF_NATIVE_Q1_V1 boundary=%s step=%d rank=%d\n",top ? "top":"bottom",step,E->parallel.me);
            fprintf(fp,"# time_nd=%.17g length_scale_m=%.17g k0_W_m_K=%.17g deltaT_K=%.17g\n",
                (double)E->monitor.elapsed_time,length,E->data.k0,E->data.ref_temperature);
            fprintf(fp,"# positive=%s global_heat_W=%.17g global_area_m2=%.17g\n",
                top ? "mantle_to_surface":"core_to_mantle",totals[0]*length*length,totals[1]*length*length);
            fprintf(fp,"# state=output_T_and_solver_Tdot filter=%d lith_age=%d initial_or_restart_state=%d\n",
                E->advection.filter_temperature,E->control.lith_age,step==E->monitor.solution_cycles_init);
            fprintf(fp,"# RHS is nondimensional outward Galerkin energy residual; mass is assembled GLL area_nd.\n");
            fprintf(fp,"# N: cap_global node_local x_nd y_nd z_nd theta_rad phi_rad r_nd q_W_m2 rhs_nd mass_nd\n");
            for(m=1;m<=E->sphere.caps_per_proc;++m)
                for(i=1;i<=E->lmesh.nsf;++i) {
                    node=top ? E->surf_node[m][i]:E->surf_node[m][i]-E->lmesh.noz+1;
                    fprintf(fp,"N %d %d %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
                        E->sphere.capid[m],node,(double)E->X[lev][m][1][node],(double)E->X[lev][m][2][node],(double)E->X[lev][m][3][node],
                        (double)E->SX[lev][m][1][node],(double)E->SX[lev][m][2][node],(double)E->SX[lev][m][3][node],
                        q[m][node],rhs[m][node],mass[m][node]);
                }
            fprintf(fp,"# F: cap_global element_local node1 node2 node3 node4 weight1_nd weight2_nd weight3_nd weight4_nd\n");
            for(m=1;m<=E->sphere.caps_per_proc;++m)
                for(e=top ? elz:1;e<=E->lmesh.nel;e+=elz) {
                    for(a=0;a<4;++a) {
                        node=E->ien[m][e].node[sidenodes[side][a+1]];
                        for(d=0;d<3;++d)x[a][d]=E->X[lev][m][d+1][node];
                    }
                    CBF_face_gll_mass(x,dm);
                    fprintf(fp,"F %d %d",E->sphere.capid[m],e);
                    for(a=0;a<4;++a)fprintf(fp," %d",E->ien[m][e].node[sidenodes[side][a+1]]);
                    for(a=0;a<4;++a)fprintf(fp," %.17g",dm[a]);
                    fprintf(fp,"\n");
                }
            if(ferror(fp))bad=1;
            if(fclose(fp)!=0)bad=1;
        }
    }
    CBF_io_check(E,bad);
    if(active && rename(tmp,path)!=0)bad=1;
    CBF_io_check(E,bad);
}
#endif
