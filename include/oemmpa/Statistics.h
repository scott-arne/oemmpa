#pragma once

namespace OEMMPA {

/// \brief Regularized incomplete beta function I_x(a, b).
///
/// Continued-fraction evaluation. Used to obtain Student's t tail probability
/// without a SciPy dependency.
///
/// \param a First shape parameter, positive.
/// \param b Second shape parameter, positive.
/// \param x Upper limit, in [0, 1].
/// \returns I_x(a, b).
double regularized_incomplete_beta(double a, double b, double x);

/// \brief Two-sided p-value for a Student's t statistic.
///
/// Equivalent to ``scipy.stats.t.sf(abs(t), df) * 2`` via the identity
/// ``2 * sf(|t|, df) == I_{df/(df+t^2)}(df/2, 1/2)``.
///
/// \param t_statistic The t statistic; the sign is immaterial.
/// \param degrees_of_freedom Degrees of freedom.
/// \returns The two-sided p-value in [0, 1], or NaN when
///     ``degrees_of_freedom`` is not positive.
double two_sided_t_p_value(double t_statistic, double degrees_of_freedom);

}  // namespace OEMMPA
