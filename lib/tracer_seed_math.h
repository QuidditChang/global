#ifndef CITCOMS_TRACER_SEED_MATH_H
#define CITCOMS_TRACER_SEED_MATH_H
#include <math.h>
#include <stdint.h>
/* A private stream: enabling a new seeder must not consume legacy drand48. */
static uint64_t tracer_seed_bits(uint64_t *state) {
    uint64_t z=(*state+=UINT64_C(0x9e3779b97f4a7c15));
    z=(z^(z>>30))*UINT64_C(0xbf58476d1ce4e5b9);
    z=(z^(z>>27))*UINT64_C(0x94d049bb133111eb);
    return z^(z>>31);
}
static double tracer_seed_unit(uint64_t *state) {
    /* Exactly representable midpoints, strictly inside (0,1). */
    return ((double)(tracer_seed_bits(state)>>12)+0.5)*0x1p-52;
}
static double tracer_seed_area(const double a[2],const double b[2],const double c[2]) {
    return 0.5*fabs((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]));
}
static double tracer_seed_solid_angle(const double a[3],const double b[3],const double c[3]) {
    double det=a[0]*(b[1]*c[2]-b[2]*c[1])+a[1]*(b[2]*c[0]-b[0]*c[2])+a[2]*(b[0]*c[1]-b[1]*c[0]);
    double den=1+a[0]*b[0]+a[1]*b[1]+a[2]*b[2]+a[0]*c[0]+a[1]*c[1]+a[2]*c[2]+b[0]*c[0]+b[1]*c[1]+b[2]*c[2];
    return 2*atan2(fabs(det),den);
}
/* Proposal is uniform in the two planar triangles; acceptance converts its
 * area measure to the unit sphere. The universal upper bound on J is one. */
static int tracer_seed_direction(uint64_t *state,const double uv[4][2],
                                 const double area[2],double *u,double *v) {
    static const int corner[2][3]={{0,1,2},{0,2,3}};
    int trial;
    for(trial=0;trial<100000;trial++) {
        int t=tracer_seed_unit(state)*(area[0]+area[1])<area[0]?0:1;
        double s=sqrt(tracer_seed_unit(state)),q=tracer_seed_unit(state);
        double w[3]={1-s,s*(1-q),s*q},r2;
        *u=w[0]*uv[corner[t][0]][0]+w[1]*uv[corner[t][1]][0]+w[2]*uv[corner[t][2]][0];
        *v=w[0]*uv[corner[t][0]][1]+w[1]*uv[corner[t][1]][1]+w[2]*uv[corner[t][2]][1];
        r2=1+(*u)*(*u)+(*v)*(*v);
        if(tracer_seed_unit(state)<=1/(r2*sqrt(r2)))return trial+1;
    }
    return 0;
}
static double tracer_seed_radius(uint64_t *state,double lo,double hi) {
    double r;int trial;
    /* Reject the measure-zero floating-point endpoints, never clamp them. */
    for(trial=0;trial<128;trial++) {
        r=cbrt(lo*lo*lo+tracer_seed_unit(state)*(hi*hi*hi-lo*lo*lo));
        if(r>lo&&r<hi)return r;
    }
    return NAN; /* The caller's unchanged host predicate fails closed. */
}
#endif
