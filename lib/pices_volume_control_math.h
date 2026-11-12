#ifndef CITCOMS_PICES_VOLUME_CONTROL_MATH_H
#define CITCOMS_PICES_VOLUME_CONTROL_MATH_H
#include <math.h>
/* A steady compact vortex. q.u=0 keeps |q| and angular speed constant.
 * div(omega(|q|) ez cross q)=0, so the exact map preserves volume, even
 * though its centre differs from the centre of the spherical shell. */
struct Pices_control_flow { double center[3], radius, omega; };
static double pices_control_omega(const struct Pices_control_flow *f,const double x[3])
{
    double s=0,q;int d;
    for(d=0;d<3;d++){q=x[d]-f->center[d];s+=q*q;}
    s/=f->radius*f->radius;
    return s<1 ? f->omega*(1-s)*(1-s) : 0;
}
static void pices_control_velocity(const struct Pices_control_flow *f,const double x[3],double u[3])
{
    double w=pices_control_omega(f,x);
    u[0]=-w*(x[1]-f->center[1]);u[1]=w*(x[0]-f->center[0]);u[2]=0;
}
static void pices_control_exact(const struct Pices_control_flow *f,const double x[3],double t,double y[3])
{
    double w=pices_control_omega(f,x),a=w*t,c,s,qx,qy;
    if(w==0 || t==0){y[0]=x[0];y[1]=x[1];y[2]=x[2];return;}
    c=cos(a);s=sin(a);qx=x[0]-f->center[0];qy=x[1]-f->center[1];
    y[0]=f->center[0]+c*qx-s*qy;
    y[1]=f->center[1]+s*qx+c*qy;y[2]=x[2];
}
#endif
