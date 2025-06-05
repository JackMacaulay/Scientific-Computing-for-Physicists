#ifndef FRICTIONH
#define FRICTIONH

#include <rarray>

/**
 * @brief Graduate Assignment PHY1610 - Friction Rate Calculator
 *
 * Computes an estimate of the friction rate from a series of velocity samples.
 * The method uses the average ratio of successive accelerations (obtained via numerical differentiation)
 * to determine the exponential decay constant.
 *
 * @param dt Time interval between velocity measurements.
 * @param v Vector of velocity samples.
 * @return double Estimated friction rate.
 */
double frictionrate(double dt, const rvector<double>& v);

/**
 * @brief Graduate Assignment PHY1610 - Numerical Differentiation Utility
 *
 * Estimates velocities by computing finite differences on a series of position samples.
 *
 * @param dt Time interval between position samples.
 * @param z Vector of position samples.
 * @return rvector<double> Vector of estimated velocities.
 */
rvector<double> numdiff(double dt, const rvector<double>& z);

#endif
