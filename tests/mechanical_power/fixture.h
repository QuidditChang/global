#include <math.h>
#include <float.h>
#include <string.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "qvis_limiter.h"
#include "drive_solvers.h"
#include "npz_writer.h"

static int fixture_size=1;
static void exchange(struct All_variables *E,double **v,int level) {
    MPI_Allreduce(MPI_IN_PLACE,v[1],24,MPI_DOUBLE,MPI_SUM,E->parallel.world);
}
/* Deliberately unusable: production diagnostics must bypass the BC-masked
 * nodal action. Coupled get_elt_k below is the full matrix authority. */
void assemble_del2_u(struct All_variables *E,double **u,double **v,int lev,int strip) {
    MPI_Abort(E->parallel.world,99);
}
void get_global_shape_fn(struct All_variables *E,int e,struct Shape_function *n,
    struct Shape_function_dx *dx,struct Shape_function_dA *da,int p,int sphere,
    double rtf[4][9],int lev,int m) {
    int i; for(i=1;i<=8;i++) da->vpt[i]=0.125/fixture_size;
}
void construct_c3x3matrix_el(struct All_variables *E,int e,struct CC *c,
    struct CCX *cx,int lev,int m,int p) {
    int a,i; memset(c,0,sizeof(*c));
    for(a=1;a<=8;a++) for(i=1;i<=8;i++) {
        c->vpt[BVINDEX(3,1,a,i)]=0.6;
        c->vpt[BVINDEX(3,3,a,i)]=0.8;
    }
}
void get_elt_k(struct All_variables *E,int e,double *k,int lev,int m,int flag) {
    int i; memset(k,0,24*24*sizeof(double));
    for(i=0;i<24;i++) k[i*24+i]=(2+i*0.1)/fixture_size;
    k[1]=k[24]=-0.5/fixture_size; /* prescribed/free coupling */
}
void get_aug_k(struct All_variables *E,int e,double *k,int lev,int m) {
    int i; for(i=0;i<24;i++) k[i*24+i]+=0.2/fixture_size;
}
void get_elt_tr(struct All_variables *E,int i,int side,double *f,int m) {
    if(side==SIDE_BEGIN) f[2]+=0.7/fixture_size;
}
void remove_horiz_ave2(struct All_variables *E,double **f) {
    double avg=0; int i; for(i=1;i<=8;i++) avg+=f[1][i]/8;
    for(i=1;i<=8;i++) f[1][i]-=avg;
}
float phase_change_reference_fraction(struct All_variables *E,int p,int m,int i) {
    return 0.25f;
}
void strain_rate_2_inv(struct All_variables *E,int m,float *s,int root) { s[1]=E->U[1][0]==0 ? 0 : 3; }
