#include <math.h>
#include <float.h>
#include <string.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "drive_solvers.h"
#include "npz_writer.h"

static void exchange(struct All_variables *E,double **v,int level) {}
double global_vdot(struct All_variables *E,double **x,double **y,int lev) {
    double v=0; int i; for(i=0;i<24;i++) v+=x[1][i]*y[1][i]; return v;
}
void assemble_del2_u(struct All_variables *E,double **u,double **v,int lev,int strip) {
    int i; for(i=0;i<24;i++) v[1][i]=(2+i*0.1)*u[1][i];
}
void get_global_shape_fn(struct All_variables *E,int e,struct Shape_function *n,
    struct Shape_function_dx *dx,struct Shape_function_dA *da,int p,int sphere,
    double rtf[4][9],int lev,int m) {
    int i; for(i=1;i<=8;i++) da->vpt[i]=0.125;
}
void construct_c3x3matrix_el(struct All_variables *E,int e,struct CC *c,
    struct CCX *cx,int lev,int m,int p) {
    int a,i; memset(c,0,sizeof(*c));
    for(a=1;a<=8;a++) for(i=1;i<=8;i++) {
        c->vpt[BVINDEX(3,1,a,i)]=0.6;
        c->vpt[BVINDEX(3,3,a,i)]=0.8;
    }
}
void get_elt_k(struct All_variables *E,int e,double *k,int lev,int m,int flag) {}
void get_elt_tr(struct All_variables *E,int i,int side,double *f,int m) {
    if(side==SIDE_BEGIN) f[2]+=0.7;
}
void remove_horiz_ave2(struct All_variables *E,double **f) {
    double avg=0; int i; for(i=1;i<=8;i++) avg+=f[1][i]/8;
    for(i=1;i<=8;i++) f[1][i]-=avg;
}
float phase_change_reference_fraction(struct All_variables *E,int p,int m,int i) {
    return 0.25f;
}
void strain_rate_2_inv(struct All_variables *E,int m,float *s,int root) { s[1]=E->U[1][0]==0 ? 0 : 3; }
