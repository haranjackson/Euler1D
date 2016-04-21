#include <QDebug>

#include <algorithm>
#include <vector>

#include "common.h"
#include "flux.h"
#include "force.h"
#include "variables.h"


double phiG(double r, double c)
{
    return (r<1) ? 0 : (1-c)/(1+c);
}


double phi_flic(double r, double c, int l)
{
    assert(0<=l && l<=3);
    switch (l)
    {
        case 0:
            return (r<=1) ? psiSB(r) : std::min(phiG(r,c)+(1-phiG(r,c))*r, 2.);
        case 1:
            return phiG(r,c) + (1 - phiG(r,c)) * psiVL(r);
        case 2:
            return phiG(r,c) + (1 - phiG(r,c)) * psiVA(r);
        default:
            return psiMB(r);
    }
}


std::vector<double> flic_flux(state CLL, state CL, state CR, state CRR,
                              double dx, double dt, double C, int l)
{
    // Returns the FLIC flux at the interface of neighbouring cells with states
    // CL and CR

    double TOL = 1.0e-5;
    double q0 = bound_below(CL.r()  - CLL.r(), TOL);
    double q1 = bound_below(CR.r()  - CL.r(),  TOL);
    double q2 = bound_below(CRR.r() - CR.r(),  TOL);
    double rLeft = q0 / q1;
    double rRight = q2 / q1;
    double phi = std::min(phi_flic(rLeft,C,l), phi_flic(rRight,C,l));
    std::vector<double> FFO = force_flux(CL, CR, dx, dt);
    std::vector<double> FRI = richtmyer_flux(CL, CR, dx, dt);
    return FFO + phi * (FRI - FFO);
}


void flic_stepper(std::vector<state> & X, const std::vector<state> & X0,
                  double dx, double dt, double C, int l)
{
    for (unsigned int i=0; i<X.size(); i++)
    {
        std::vector<double> U = X[i].U();
        std::vector<double> FL = flic_flux(X0[i], X0[i+1], X0[i+2], X0[i+3],
                                           dx, dt, C, l);
        std::vector<double> FR = flic_flux(X0[i+1], X0[i+2], X0[i+3], X0[i+4],
                                           dx, dt, C, l);
        X[i].setU(U + (dt/dx) * (FL-FR));
    }
}


std::vector<state> flic(const std::vector<state> & X, int BCL, int BCR,
                        double dx, double C, double t, int l)
{
    double t0 = 0;
    int iter = 0;
    std::vector<state> ret = X;

    while (t0 < t)
    {
        std::vector<state> X0 = add_boundary_conditions(ret, BCL, BCR, 2);
        double dt = timestep(X0, dx, C, t0, t, iter);
        flic_stepper(ret, X0, dx, dt, C, l);
        t0 += dt;
        iter++;
    }
    return ret;
}
