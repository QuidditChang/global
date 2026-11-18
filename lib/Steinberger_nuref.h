/* Figure 1 + EBA Tref reference family: 7 baseline, 71 layer factors,
 * 72 layer polynomials (the historical M2-A coefficients by default).
 * Provenance and conventions: ../doc/rheol7_nuref.md.
 * Dimensional inputs: depth in km, reference temperature in K, eta0 in Pa s.
 * Returns eta_ref/eta0 (dynamic viscosity ratio, NOT eta/rho).
 * Temperature correction uses Az, with a configurable cold-side multiplier.
 */
#ifndef CITCOMS_STEINBERGER_NUREF_H
#define CITCOMS_STEINBERGER_NUREF_H
#include <math.h>
#include <stddef.h>

#define STEINBERGER_LAYERS 4

static inline int steinberger_is_family(int rheol)
{
    return rheol == 7 || rheol == 71 || rheol == 72;
}

/* Frozen M2-A coefficients are defaults, not part of the physical baseline. */
static inline void steinberger_defaults(double *q, double *a, double *b, double *c)
{
    static const double coefficients[4][3] = {
        {5.067917152408204, -8.15474461356726, 4.799095524666932},
        {0.22571234789343927, 0.0, 0.0},
        {0.3127706714585134, 0.0, 0.0},
        {2.073829265582211, -0.7521860487837174, 0.09721318357053432}
    };
    int i;
    for(i=0;i<STEINBERGER_LAYERS;i++) {
        q[i] = 1.0;
        a[i] = coefficients[i][0];
        b[i] = coefficients[i][1];
        c[i] = coefficients[i][2];
    }
}

/* Validate the whole interval u in [0,1], including an interior minimum. */
static inline int steinberger_validate_correction(int rheol, const double *q,
                                                 const double *a, const double *b,
                                                 const double *c)
{
    int i;
    if(rheol == 7) return 1;
    if(rheol == 71) {
        if(!q) return 0;
        for(i=0;i<STEINBERGER_LAYERS;i++)
            if(!isfinite(q[i]) || q[i] <= 0.0) return 0;
        return 1;
    }
    if(rheol != 72 || !a || !b || !c) return 0;
    for(i=0;i<STEINBERGER_LAYERS;i++) {
        double end, u, minimum;
        if(!isfinite(a[i]) || !isfinite(b[i]) || !isfinite(c[i])) return 0;
        end = a[i]+b[i]+c[i];
        if(a[i] <= 0.0 || !isfinite(end) || end <= 0.0) return 0;
        if(c[i] > 0.0 && b[i] < 0.0) {
            u = -0.5*(b[i]/c[i]);
            if(u > 0.0 && u < 1.0) {
                minimum = a[i]+u*(b[i]+u*c[i]);
                if(!isfinite(minimum) || minimum <= 0.0) return 0;
            }
        }
    }
    return 1;
}

static inline double steinberger_factor(int rheol, double depth, const double *q,
                                        const double *a, const double *b,
                                        const double *c)
{
    static const double boundary[5] = {0.,410.,520.,660.,2891.};
    int segment;
    double u, f;
    if(!isfinite(depth) || depth < -1.e-5 || depth > 2891.+1.e-5) return -1.;
    if(depth < 0.) depth = 0.;
    if(depth > 2891.) depth = 2891.;
    if(rheol == 7) return 1.;
    segment = depth < 410. ? 0 : (depth < 520. ? 1 : (depth < 660. ? 2 : 3));
    if(rheol == 71 && q) f = q[segment];
    else if(rheol == 72 && a && b && c) {
        u = (depth-boundary[segment])/(boundary[segment+1]-boundary[segment]);
        f = a[segment]+u*(b[segment]+u*c[segment]);
    } else return -1.;
    return isfinite(f) && f > 0. ? f : -1.;
}

/* Raw activation-temperature profile (K); normalization offsets excluded. */
static inline double steinberger_Az(double depth)
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

static inline double steinberger_nuref_mode(int rheol, double depth, double tref_K,
                                           double eta0, const double *q,
                                           const double *a, const double *b,
                                           const double *c)
{
    double value, exponent, f;
    if (!isfinite(depth) || !isfinite(tref_K) || !isfinite(eta0) ||
        depth < -1.e-5 || depth > 2891.+1.e-5 || tref_K <= 0. || eta0 <= 0.)
        return -1.;
    /* Roundoff tolerance at physical endpoints only. */
    if(depth < 0.) depth = 0.;
    if(depth > 2891.) depth = 2891.;
    exponent = steinberger_Az(depth)/tref_K;
    if(depth >= 660.) exponent -= 7.466790490109812;
    exponent -= 13.163586986394256;
    f = steinberger_factor(rheol,depth,q,a,b,c);
    if(f <= 0.) return -1.;
    value = (1.e21/eta0)*f*exp(exponent);
    return isfinite(value) && value > 0. ? value : -1.;
}
/* Exact Kelvin form of rheol=3 reference-relative Arrhenius structure.
 * Cold scaling affects only the anomaly exponent, never nuref itself. */
static inline double steinberger_viscosity_mode(int rheol, double depth, double tref_K,
                                               double temperature_K, double eta0,
                                               double cold_scale, const double *q,
                                               const double *a, const double *b,
                                               const double *c)
{
    double reference, side_scale, value;
    if(!isfinite(cold_scale) || cold_scale < 0.) return -1.;
    if(!isfinite(temperature_K) || temperature_K <= 0.) return -1.;
    reference = steinberger_nuref_mode(rheol,depth,tref_K,eta0,q,a,b,c);
    if(reference <= 0.) return -1.;
    if(depth < 0.) depth = 0.;
    if(depth > 2891.) depth = 2891.;
    side_scale = temperature_K < tref_K ? cold_scale : 1.;
    value = reference * exp(side_scale * steinberger_Az(depth)
                           * (1./temperature_K - 1./tref_K));
    return isfinite(value) && value > 0. ? value : -1.;
}

/* Unqualified helpers now explicitly mean the rheol=7 physical baseline. */
static inline double steinberger_nuref(double depth, double tref_K, double eta0)
{
    return steinberger_nuref_mode(7,depth,tref_K,eta0,NULL,NULL,NULL,NULL);
}
static inline double steinberger_viscosity(double depth, double tref_K,
                                          double temperature_K, double eta0,
                                          double cold_scale)
{
    return steinberger_viscosity_mode(7,depth,tref_K,temperature_K,eta0,cold_scale,
                                     NULL,NULL,NULL,NULL);
}
#endif
