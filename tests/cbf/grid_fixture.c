/* Synthetic closed six-face mesh for the production mapper/NetCDF writer.
 * Not a physical CitcomS solution or a spherical-shell accuracy benchmark. */
#include <string.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "cbf_output.h"
void parallel_process_termination(void) {MPI_Abort(MPI_COMM_WORLD,2);}
int main(int argc,char **argv)
{
    static const double face[6][4][3]={
      {{1,-1,-1},{1,1,-1},{1,1,1},{1,-1,1}},
      {{-1,1,-1},{-1,-1,-1},{-1,-1,1},{-1,1,1}},
      {{1,1,-1},{-1,1,-1},{-1,1,1},{1,1,1}},
      {{-1,-1,-1},{1,-1,-1},{1,-1,1},{-1,-1,1}},
      {{-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1}},
      {{-1,1,-1},{1,1,-1},{1,-1,-1},{-1,-1,-1}}};
    struct All_variables *E=calloc(1,sizeof(*E));
    int rank,size,nface,e,a,d,node,f;double sum;
    MPI_Init(&argc,&argv);MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
    if(6%size)MPI_Abort(MPI_COMM_WORLD,3);
    E->parallel.world=MPI_COMM_WORLD;E->parallel.me=rank;E->parallel.nprocz=1;
    E->parallel.me_loc[3]=0;E->sphere.caps_per_proc=1;E->sphere.ro=1;E->sphere.ri=0.5;
    E->data.kappa0=1.e-6;E->data.k0=4;E->data.ref_temperature=3400;
    E->data.radius_km=1;E->monitor.elapsed_time=0.125;
    E->output.cbf_output_shflux=E->output.cbf_output_bhflux=1;
    nface=6/size;E->lmesh.nel=nface;E->lmesh.elz=1;E->lmesh.noz=2;
    E->lmesh.nno=nface*8;E->lmesh.nsf=nface*4;E->mesh.levmax=1;
    E->ien[1]=calloc(nface+1,sizeof(struct IEN));
    E->surf_node[1]=calloc(nface*4+1,sizeof(int));
    E->slice.shflux_CBF[1]=calloc(nface*4+1,sizeof(double));
    E->slice.bhflux_CBF[1]=calloc(nface*4+1,sizeof(double));
    for(d=1;d<=3;++d)E->X[1][1][d]=calloc(nface*8+1,sizeof(double));
    for(e=1;e<=nface;++e) {
        f=rank*nface+e-1;
        for(a=0;a<4;++a) {
            node=(e-1)*8+2*a+1;
            E->ien[1][e].node[a+1]=node;E->ien[1][e].node[a+5]=node+1;
            E->surf_node[1][(e-1)*4+a+1]=node+1;sum=0;
            for(d=0;d<3;++d) {
                E->X[1][1][d+1][node]=0.5*face[f][a][d];
                E->X[1][1][d+1][node+1]=face[f][a][d];sum+=face[f][a][d];
            }
            E->slice.shflux_CBF[1][(e-1)*4+a+1]=3+0.1*sum;
            E->slice.bhflux_CBF[1][(e-1)*4+a+1]=7-0.2*sum;
        }
    }
    cbf_output_grids(E,50);
    MPI_Finalize();return 0;
}
