#ifndef CITCOMS_CBF_GEOMETRY_H
#define CITCOMS_CBF_GEOMETRY_H
#include <math.h>

/* Cyclic Q1 face nodes: (-1,-1),(1,-1),(1,1),(-1,1).
 * The area is that of the actual trilinear FE face, not its planar projection.
 * GLL weights at the four Q1 support points are all one. */
static double cbf_face_jacobian(const double x[4][3], double u, double v)
{
    static const double su[4]={-1,1,1,-1}, sv[4]={-1,-1,1,1};
    double du[3]={0,0,0}, dv[3]={0,0,0}, cross[3];
    int a,d;
    for(a=0;a<4;++a)
        for(d=0;d<3;++d) {
            du[d] += 0.25*su[a]*(1+sv[a]*v)*x[a][d];
            dv[d] += 0.25*sv[a]*(1+su[a]*u)*x[a][d];
        }
    cross[0]=du[1]*dv[2]-du[2]*dv[1];
    cross[1]=du[2]*dv[0]-du[0]*dv[2];
    cross[2]=du[0]*dv[1]-du[1]*dv[0];
    return sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
}
static void cbf_face_gll_mass(const double x[4][3], double mass[4])
{
    mass[0]=cbf_face_jacobian(x,-1,-1);
    mass[1]=cbf_face_jacobian(x, 1,-1);
    mass[2]=cbf_face_jacobian(x, 1, 1);
    mass[3]=cbf_face_jacobian(x,-1, 1);
}
/* Intersect a ray with the Q1 face. Work in a tangent basis to the ray;
 * solve its two projected coordinates for (u,v), then interpolate nodal q.
 * Cartesian geometry makes the operation periodic and pole independent. */
static int cbf_face_ray_value(const double x[4][3], const double q[4],
                              const double ray[3], double *value)
{
    static const double su[4]={-1,1,1,-1},sv[4]={-1,-1,1,1};
    double t[3],b[3],p[4][2],u=0,v=0,fu,fv,uu,uv,vu,vv;
    double det,du,dv,norm,weight,radial;
    int a,d,it;
    if(fabs(ray[2])<0.9) {t[0]=-ray[1];t[1]=ray[0];t[2]=0;}
    else {t[0]=ray[2];t[1]=0;t[2]=-ray[0];}
    norm=sqrt(t[0]*t[0]+t[1]*t[1]+t[2]*t[2]);
    for(d=0;d<3;++d)t[d]/=norm;
    b[0]=ray[1]*t[2]-ray[2]*t[1];
    b[1]=ray[2]*t[0]-ray[0]*t[2];
    b[2]=ray[0]*t[1]-ray[1]*t[0];
    for(a=0;a<4;++a) {
        p[a][0]=p[a][1]=0;
        for(d=0;d<3;++d) {p[a][0]+=t[d]*x[a][d];p[a][1]+=b[d]*x[a][d];}
    }
    for(it=0;it<20;++it) {
        fu=fv=uu=uv=vu=vv=0;
        for(a=0;a<4;++a) {
            weight=0.25*(1+su[a]*u)*(1+sv[a]*v);
            fu+=weight*p[a][0];fv+=weight*p[a][1];
            uu+=0.25*su[a]*(1+sv[a]*v)*p[a][0];
            vu+=0.25*su[a]*(1+sv[a]*v)*p[a][1];
            uv+=0.25*sv[a]*(1+su[a]*u)*p[a][0];
            vv+=0.25*sv[a]*(1+su[a]*u)*p[a][1];
        }
        det=uu*vv-uv*vu;
        if(!isfinite(det) || fabs(det)<1.e-30)return 0;
        du=(fu*vv-fv*uv)/det;dv=(fv*uu-fu*vu)/det;
        u-=du;v-=dv;
        if(!isfinite(u)||!isfinite(v)||fabs(u)>4||fabs(v)>4)return 0;
        if(fabs(du)+fabs(dv)<1.e-12)break;
    }
    if(it==20 || fabs(u)>1+1.e-9 || fabs(v)>1+1.e-9)return 0;
    *value=0;radial=0;
    for(a=0;a<4;++a) {
        weight=0.25*(1+su[a]*u)*(1+sv[a]*v);
        *value+=weight*q[a];
        for(d=0;d<3;++d)radial+=weight*x[a][d]*ray[d];
    }
    return radial>0 && isfinite(*value);
}
#endif
