#ifndef CITCOMS_QVIS_LIMITER_H
#define CITCOMS_QVIS_LIMITER_H
#include <math.h>
/* S = 2 eps:eps = 4 eps_II^2, tau_II = eta sqrt(S).
 * Tensor radial scaling, never componentwise clipping. Pressure is SI Pa.
 * sigma_y uses the sqrt(J2) convention of the documented DP formula.
 * This changes thermal conversion only; EVi and the Stokes operator stay intact. */
static double qvis_cap_factor(struct All_variables *E, int element,
                              double viscosity, double strain_sqr)
{
    int z;
    double pressure, phi, yield_pa, stress_scale, stress_pa;
    if(E->control.qvis_mode == 0 || strain_sqr <= 0.0 || viscosity <= 0.0)
        return 1.0;
    z = (element-1) % E->lmesh.elz + 1;
    pressure = 0.5*(E->refstate.lithostatic_pressure_pa[z]
                   + E->refstate.lithostatic_pressure_pa[z+1]);
    phi = E->control.qvis_friction_angle_rad;
    yield_pa = 6.0*(E->control.qvis_cohesion_pa*cos(phi)+pressure*sin(phi))
                  /(sqrt(3.0)*(3.0+sin(phi)));
    stress_scale = (double)E->data.ref_viscosity * E->data.kappa0
                   /pow(E->data.radius_km*1000.0, 2.0);
    stress_pa = viscosity*sqrt(strain_sqr)*stress_scale;
    return stress_pa > yield_pa ? yield_pa/stress_pa : 1.0;
}
#endif
