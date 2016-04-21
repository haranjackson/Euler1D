#include <vector>

#include "common.h"
#include "force.h"
#include "variables.h"


std::vector<double> slic_flux(state CLL, state CL, state CR, state CRR,
                              double dx, double dt, double w, int l)
{
    state CRi = evolve_state(CLL, CL, CR, dx, dt, w, l, 1);
    state CLi_1 = evolve_state(CL, CR, CRR, dx, dt, w, l, 0);
    return force_flux(CRi, CLi_1, dx, dt);
}


void slic_stepper(std::vector<state> & X, const std::vector<state> & X0,
                  double dx, double dt, double w, int l)
{
    for (unsigned int i=0; i<X.size(); i++)
    {
        std::vector<double> U = X[i].U();
        std::vector<double> FL = slic_flux(X0[i], X0[i+1], X0[i+2], X0[i+3],
                                           dx, dt, w, l);
        std::vector<double> FR = slic_flux(X0[i+1], X0[i+2], X0[i+3], X0[i+4],
                                           dx, dt, w, l);
        X[i].setU(U + (dt/dx) * (FL-FR));
    }
}


std::vector<state> slic(const std::vector<state> & X, int BCL, int BCR,
                        double dx, double C, double t, double w, int l)
{
    double t0 = 0;
    int iter = 0;
    std::vector<state> ret = X;

    while (t0 < t)
    {
        std::vector<state> X0 = add_boundary_conditions(ret, BCL, BCR, 2);
        double dt = timestep(X0, dx, C, t0, t, iter);
        slic_stepper(ret, X0, dx, dt, w, l);
        t0 += dt;
        iter++;
    }
    return ret;
}
