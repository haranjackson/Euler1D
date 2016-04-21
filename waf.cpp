#include <QDebug>

#include <cmath>
#include <vector>

#include "common.h"
#include "flux.h"
#include "hllc.h"
#include "variables.h"


double phi_waf(double r, double c, int l)
{
    assert(0<=l && l<=3);

    switch (l)
    {
        case 0:
            return 1 - (1-std::abs(c)) * psiSB(r);
        case 1:
            return 1 - (1-std::abs(c)) * psiVL(r);
        case 2:
            return 1 - (1-std::abs(c)) * psiVA(r);
        default:
            return 1 - (1-std::abs(c)) * psiMB(r);
    }
}


std::vector<double> density_changes(state CL, state CR)
{
    // Returns the changes in density across the HLLC wavefronts, given
    // neighbouring cells with states CL and CR

    triple S = wavespeeds(CL, CR);
    std::vector<double> ret(3);
    double rL = CL.r();
    double rL_ = r_star(CL.r(), S.L, CL.u(), S.Star);
    double rR_ = r_star(CR.r(), S.R, CR.u(), S.Star);
    double rR = CR.r();
    ret[0] = rL_ - rL;
    ret[1] = rR_ - rL_;
    ret[2] = rR - rR_;
    return ret;
}


std::vector<std::vector<double> > flux_changes(state CL, state CR, triple S)
{
    // Returns the changes in flux across the HLLC wavefronts, given
    // neighbouring cells with states CL and CR, and wavespeeds S

    std::vector<std::vector<double> > ret(3);
    std::vector<double> FL_ = f_star(CL, S.L, S.Star);
    std::vector<double> FR_ = f_star(CR, S.R, S.Star);
    ret[0] = FL_ - CL.F();
    ret[1] = FR_ - FL_;
    ret[2] = CR.F() - FR_;
    return ret;
}


std::vector<double> ratio(state CLL, state CL, state CR, state CRR,
                          std::vector<double> c)
{
    double TOL = 1.0e-6;
    std::vector<double> q_L = density_changes(CLL, CL);
    std::vector<double> q_M = density_changes(CL, CR);
    std::vector<double> q_R = density_changes(CR, CRR);
    std::vector<double> ret(3);
    for (int i=0; i<3; i++)
    {
        double qL = bound_below(q_L[i], TOL);
        double qM = bound_below(q_M[i], TOL);
        double qR = bound_below(q_R[i], TOL);
        ret[i] = (0<=c[i]) ? qL/qM : qR/qM;
    }
    return ret;
}


std::vector<double> waf_flux(state CLL, state CL, state CR, state CRR,
                             double dx, double dt, int l)
{
    // Returns the WAF flux at the interface of neighbouring cells with states
    // CL and CR

    triple S = wavespeeds(CL, CR);
    std::vector<double> c = (dt/dx) * S.toVector();
    std::vector<double> r = ratio(CLL, CL, CR, CRR, c);
    std::vector<std::vector<double> > dF = flux_changes(CL, CR, S);
    std::vector<double> ret = CL.F() + CR.F();
    for (int i=0; i<3; i++)
        ret = ret - sgn(c[i]) * phi_waf(r[i],c[i],l) * dF[i];
    return 0.5 * ret;
}


void waf_stepper(std::vector<state> & X, const std::vector<state> & X0,
                 double dx, double dt, int l)
{
    for (unsigned int i=0; i<X.size(); i++)
    {
        std::vector<double> FL = waf_flux(X0[i], X0[i+1], X0[i+2], X0[i+3],
                                          dx, dt, l);
        std::vector<double> FR = waf_flux(X0[i+1], X0[i+2], X0[i+3], X0[i+4],
                                          dx, dt, l);
        std::vector<double> U = X[i].U();
        X[i].setU(U + (dt/dx) * (FL-FR));
    }
}


std::vector<state> waf(const std::vector<state> & X, int BCL, int BCR,
                       double dx, double C, double t, int l)
{
    double t0 = 0;
    int iter = 0;
    std::vector<state> ret = X;

    while (t0 < t)
    {
        std::vector<state> X0 = add_boundary_conditions(ret, BCL, BCR, 2);
        double dt = timestep(X0, dx, C, t0, t, iter);;
        waf_stepper(ret, X0, dx, dt, l);
        t0 += dt;
        iter++;
    }
    return ret;
}
