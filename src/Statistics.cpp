#include "oemmpa/Statistics.h"

#include <cmath>
#include <limits>

namespace OEMMPA {
namespace {

// Continued-fraction evaluation for the regularized incomplete beta function
// (Numerical Recipes betacf, modified Lentz). Used to compute the Student's t
// two-sided p-value without a third-party dependency, matching scipy's
// stats.t.sf(|t|, df) * 2 to machine precision.
double incomplete_beta_cf(double a, double b, double x) {
    const int MAX_ITERATIONS = 200;
    const double EPSILON = 3.0e-16;
    const double FLOOR = 1.0e-300;
    const double qab = a + b;
    const double qap = a + 1.0;
    const double qam = a - 1.0;
    double c = 1.0;
    double d = 1.0 - qab * x / qap;
    if (std::fabs(d) < FLOOR) {
        d = FLOOR;
    }
    d = 1.0 / d;
    double h = d;
    for (int m = 1; m <= MAX_ITERATIONS; ++m) {
        const double m2 = 2.0 * m;
        double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (std::fabs(d) < FLOOR) {
            d = FLOOR;
        }
        c = 1.0 + aa / c;
        if (std::fabs(c) < FLOOR) {
            c = FLOOR;
        }
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (std::fabs(d) < FLOOR) {
            d = FLOOR;
        }
        c = 1.0 + aa / c;
        if (std::fabs(c) < FLOOR) {
            c = FLOOR;
        }
        d = 1.0 / d;
        const double delta = d * c;
        h *= delta;
        if (std::fabs(delta - 1.0) < EPSILON) {
            break;
        }
    }
    return h;
}

}  // namespace

// Regularized incomplete beta function I_x(a, b).
double regularized_incomplete_beta(double a, double b, double x) {
    if (x <= 0.0) {
        return 0.0;
    }
    if (x >= 1.0) {
        return 1.0;
    }
    const double front = std::exp(
        std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b) +
        a * std::log(x) + b * std::log1p(-x)
    );
    if (x < (a + 1.0) / (a + b + 2.0)) {
        return front * incomplete_beta_cf(a, b, x) / a;
    }
    return 1.0 - front * incomplete_beta_cf(b, a, 1.0 - x) / b;
}

// Two-sided Student's t p-value for a t statistic with the given degrees of
// freedom. Equals scipy.stats.t.sf(|t|, df) * 2 via the identity
// 2 * sf(|t|, df) == I_{df/(df+t^2)}(df/2, 1/2).
double two_sided_t_p_value(double t, double df) {
    if (df <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double x = df / (df + t * t);
    return regularized_incomplete_beta(df / 2.0, 0.5, x);
}

}  // namespace OEMMPA
