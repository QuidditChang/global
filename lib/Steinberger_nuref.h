/* Frozen Figure 1 -> adiabatic Figure 4 -> quadratic M2-A reference.
 * Provenance and conventions: ../doc/rheol7_nuref.md.
 * Dimensional inputs: depth in km, reference temperature in K, eta0 in Pa s.
 * Returns eta_ref/eta0 (dynamic viscosity ratio, NOT eta/rho).
 * Temperature correction uses Az, with a configurable cold-side multiplier.
 */
#ifndef CITCOMS_STEINBERGER_NUREF_H
#define CITCOMS_STEINBERGER_NUREF_H
#include <math.h>
/* Raw activation-temperature profile (K); normalization offsets excluded. */
static double steinberger_Az(double depth)
{
    static const double h[6] = {525.0, 267.7486860312165, -151.850572928359, 557.680505165775, -757.6458479851975, 318.24791988551436};
    static const double tm[3] = {3386.046041109414, 3012.4459371739035, -1201.4200029742149};
    double x, value;
    int j;
    if(depth < 660.) {
        x = depth/660.;
        value = h[5];
        for(j=4;j>=0;j--) value = value*x+h[j];
        return 1000.*value/(3.5*8.3144);
    }
    x = (depth-660.)/2231.;
    return 12.*(tm[0]+x*(tm[1]+x*tm[2]));
}

static double steinberger_nuref(double depth, double tref_K, double eta0)
{
    static const double factor[4][3] = {
        {5.067917152408204, -8.15474461356726, 4.799095524666932},
        {0.22571234789343927, 0.0, 0.0},
        {0.3127706714585134, 0.0, 0.0},
        {2.073829265582211, -0.7521860487837174, 0.09721318357053432}
    };
    static const double boundary[5] = {0.,410.,520.,660.,2891.};
    double value, exponent, u, f;
    int segment;
    if (!isfinite(depth) || !isfinite(tref_K) || !isfinite(eta0) ||
        depth < -1.e-5 || depth > 2891.+1.e-5 || tref_K <= 0. || eta0 <= 0.)
        return -1.;
    /* Roundoff tolerance at physical endpoints only. */
    if(depth < 0.) depth = 0.;
    if(depth > 2891.) depth = 2891.;
    exponent = steinberger_Az(depth)/tref_K;
    if(depth >= 660.) exponent -= 7.466790490109812;
    exponent -= 13.163586986394256;
    segment = depth < 410. ? 0 : (depth < 520. ? 1 : (depth < 660. ? 2 : 3));
    u = (depth-boundary[segment])/(boundary[segment+1]-boundary[segment]);
    f = factor[segment][0]+u*(factor[segment][1]+u*factor[segment][2]);
    value = (1.e21/eta0)*f*exp(exponent);
    return isfinite(value) && value > 0. ? value : -1.;
}
/* Exact Kelvin form of rheol=3 reference-relative Arrhenius structure.
 * Cold scaling affects only the anomaly exponent, never nuref itself. */
static double steinberger_viscosity(double depth, double tref_K,
                                   double temperature_K, double eta0, double cold_scale)
{
    double reference, side_scale, value;
    if(!isfinite(cold_scale) || cold_scale < 0.) return -1.;
    if(!isfinite(temperature_K) || temperature_K <= 0.) return -1.;
    reference = steinberger_nuref(depth,tref_K,eta0);
    if(reference <= 0.) return -1.;
    if(depth < 0.) depth = 0.;
    if(depth > 2891.) depth = 2891.;
    side_scale = temperature_K < tref_K ? cold_scale : 1.;
    value = reference * exp(side_scale * steinberger_Az(depth)
                           * (1./temperature_K - 1./tref_K));
    return isfinite(value) && value > 0. ? value : -1.;
}
#endif
