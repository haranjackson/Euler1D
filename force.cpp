#include <QDebug>

#include <vector>

#include "common.h"
#include "variables.h"


std::vector<double> lax_friedrichs_flux(state CL, state CR, double dx,
                                        double dt)
{
    // Returns the Lax-Friedrichs flux at the interface of neighbouring cells
    // with states CL and CR

    std::vector<double> FL = CL.F();
    std::vector<double> FR = CR.F();
    std::vector<double> UL = CL.U();
    std::vector<double> UR = CR.U();
    return 0.5 * (FL + FR + (dx/dt) * (UL - UR));
}


std::vector<double> richtmyer_flux(state CL, state CR, double dx, double dt)
{
    // Returns the Richtmyer flux at the interface of neighbouring cells with
    // states CL and CR

    std::vector<double> FL = CL.F();
    std::vector<double> FR = CR.F();
    std::vector<double> UL = CL.U();
    std::vector<double> UR = CR.U();
    // Richtmyer vector of conserved variables
    std::vector<double> URI = 0.5 * (UL + UR + (dt/dx) * (FL-FR));
    return F(URI);
}


std::vector<double> force_flux(state CL, state CR, double dx, double dt)
{
    // Returns the FORCE flux at the interface of neighbouring cells with states
    // CL and CR

    std::vector<double> FLF = lax_friedrichs_flux(CL, CR, dx, dt);
    std::vector<double> FRI = richtmyer_flux(CL, CR, dx, dt);
    return 0.5 * (FLF + FRI);
}


void force_stepper(std::vector<state> & X, const std::vector<state> & X0,
                   double dx, double dt)
{
    for (unsigned int i=0; i<X.size(); i++)
    {
        std::vector<double> U = X[i].U();
        std::vector<double> FL = force_flux(X0[i], X0[i+1], dx, dt);
        std::vector<double> FR = force_flux(X0[i+1], X0[i+2], dx, dt);
        X[i].setU(U + (dt/dx) * (FL-FR));
    }
}


std::vector<state> force(const std::vector<state> & X, int BCL, int BCR,
                         double dx, double C, double t)
{
    double t0 = 0;
    int iter = 0;
    std::vector<state> ret = X;

    while (t0 < t)
    {
        std::vector<state> X0 = add_boundary_conditions(ret, BCL, BCR, 1);
        double dt = timestep(X0, dx, C, t0, t, iter);
        force_stepper(ret, X0, dx, dt);
        t0 += dt;
        iter++;
    }
    return ret;
}
