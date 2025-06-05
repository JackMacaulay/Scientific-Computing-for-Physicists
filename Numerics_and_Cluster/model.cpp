#include "model.h"
#include <cmath>

double v(double t, const ModelParameters& p)
{
    const double v0 = p.v0;
    const double g = p.g;
    const double alpha = p.alpha;
    const double ia = 1.0f/alpha;
    return (v0 + g*ia)*exp(-alpha*t) - g*ia;
}

double z(double t, const ModelParameters& p)
{
    const double z0 = p.z0;
    const double v0 = p.v0;
    const double g = p.g;
    const double alpha = p.alpha;
    const double ia = 1.0f/alpha;
    return z0 + ia*(v0 + g*ia)*(1 - exp(-alpha*t)) - g*ia*t;
}

rvector<double> compute_model_z(double t1, double t2, double dt,
                               const ModelParameters& p)
{
    const int n = int((t2-t1)/dt)+1;
    const rvector<double> t = linspace(t1,t2,n);
    rvector<double> zdata(n);
    for (int i = 0; i < n; i++)
        zdata[i] = z(t[i], p);
    return zdata;
}

rvector<double> compute_model_v(double t1, double t2, double dt,
                               const ModelParameters& p)
{
    const int n = int((t2-t1)/dt)+1;
    const rvector<double> t = linspace(t1,t2,n);
    rvector<double> vdata(n);
    for (int i = 0; i < n; i++)
        vdata[i] = v(t[i], p);
    return vdata;
}
