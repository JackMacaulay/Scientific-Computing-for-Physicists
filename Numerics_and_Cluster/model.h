#ifndef MODELH
#define MODELH

#include <rarray>

/**
 * @brief Graduate Assignment PHY1610 - Physical Model Parameters
 *
 * This structure encapsulates the parameters for modeling the vertical motion 
 * of a marble in a viscous fluid.
 */
struct ModelParameters
{
    double alpha;  /**< Friction constant specific to the experimental setup. */
    double g;      /**< Gravitational acceleration (in m/s^2) used in the model. */
    double v0;     /**< Initial vertical velocity (in m/s) of the marble. */
    double z0;     /**< Initial vertical position (in m) of the marble. */
};

/**
 * @brief Graduate Assignment PHY1610 - Compute Vertical Position
 *
 * Computes the vertical position of the marble at time t using the provided model parameters.
 *
 * @param t Time (in seconds) at which to compute the position.
 * @param p Model parameters.
 * @return double Vertical position at time t.
 */
double z(double t, const ModelParameters& p);

/**
 * @brief Graduate Assignment PHY1610 - Compute Vertical Velocity
 *
 * Computes the vertical velocity of the marble at time t based on the model parameters.
 *
 * @param t Time (in seconds) at which to compute the velocity.
 * @param p Model parameters.
 * @return double Vertical velocity at time t.
 */
double v(double t, const ModelParameters& p);

/**
 * @brief Graduate Assignment PHY1610 - Generate Model Velocity Data
 *
 * Generates a vector of vertical velocity values for the time interval [t1, t2] using a fixed time step dt.
 *
 * @param t1 Starting time (in seconds).
 * @param t2 Ending time (in seconds).
 * @param dt Time step between samples.
 * @param p Model parameters.
 * @return rvector<double> Vector of computed vertical velocities.
 */
rvector<double> compute_model_v(double t1, double t2, double dt,
                                const ModelParameters& p);

/**
 * @brief Graduate Assignment PHY1610 - Generate Model Position Data
 *
 * Generates a vector of vertical position values for the time interval [t1, t2] using a fixed time step dt.
 *
 * @param t1 Starting time (in seconds).
 * @param t2 Ending time (in seconds).
 * @param dt Time step between samples.
 * @param p Model parameters.
 * @return rvector<double> Vector of computed vertical positions.
 */
rvector<double> compute_model_z(double t1, double t2, double dt,
                                const ModelParameters& p);

#endif
