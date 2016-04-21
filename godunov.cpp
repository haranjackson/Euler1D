#include <vector>

#include "common.h"
#include "hllc.h"
#include "variables.h"


void godunov_stepper(std::vector<state> & X, const std::vector<state> & X0,
                     double dx, double dt)
{
    for (unsigned int i=0; i<X.size(); i++)
    {
        std::vector<double> U = X[i].U();
        std::vector<double> FL = hllc_flux(X0[i],X0[i+1]);
        std::vector<double> FR = hllc_flux(X0[i+1],X0[i+2]);
        X[i].setU(U + (dt/dx) * (FL-FR));
    }
}


std::vector<state> godunov(const std::vector<state> & X, int BCL, int BCR,
                           double dx, double C, double t)
{
    double t0 = 0;
    int iter = 0;
    std::vector<state> ret = X;

    while (t0 < t)
    {
        std::vector<state> X0 = add_boundary_conditions(ret, BCL, BCR, 1);
        double dt = timestep(X0, dx, C, t0, t, iter);
        godunov_stepper(ret, X0, dx, dt);
        t0 += dt;
        iter++;
    }
    return ret;
}
