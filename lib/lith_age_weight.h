#ifndef CITCOMS_LITH_AGE_WEIGHT_H
#define CITCOMS_LITH_AGE_WEIGHT_H

#include <math.h>

static double lith_age_softplus(double x)
{
    return x > 0.0 ? x + log1p(exp(-x)) : log1p(exp(x));
}

/* Integral of sigmoid(b+t)-sigmoid(b), t=0..h.  The expansion avoids
 * subtracting two nearly equal O(h) terms at the zero-weight bottom. */
static double lith_age_logistic_integral(double b, double h)
{
    double eb=exp(b), p=eb/(1.0+eb);
    if(h <= 0.0) return 0.0;
    if(h < 1.e-4) {
        double q=1.0-p;
        return p*q*h*h*(0.5 + h*((1.0-2.0*p)/6.0
                  + h*(1.0-6.0*p+6.0*p*p)/24.0));
    }
    if(h < 0.5)
        return log1p(p*expm1(h))-p*h;
    return lith_age_softplus(b+h)-lith_age_softplus(b)-p*h;
}

/* c and s are fractions of the local trench-dependent H.  With
 * L(x)=sigmoid((c-x)/s), w(x)=integral_x^1[L(t)-L(1)]dt / integral_0^1.
 * Thus w(0)=1, w(1)=w'(1)=0, w'<0 and w''>0 in the interior. */
static double lith_age_logistic_rate_weight(double fraction,
                                           double center, double width)
{
    double b, numerator, denominator, value;
    if(fraction <= 0.0) return 1.0;
    if(fraction >= 1.0) return 0.0;
    b=(center-1.0)/width;
    numerator=lith_age_logistic_integral(b,(1.0-fraction)/width);
    denominator=lith_age_logistic_integral(b,1.0/width);
    value=numerator/denominator;
    /* Remove floating-point endpoint roundoff, not a temperature clip. */
    if(value < 0.0) value=0.0;
    if(value > 1.0) value=1.0;
    return value;
}

#endif
