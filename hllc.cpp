#include <QDebug>

#include <algorithm>
#include <cmath>
#include <vector>

#include "variables.h"


double r_star(double r, double S, double u, double S_)
{
    // Returns density in the (left/right) side of the star region, given (left/
    // right) density, wavespeed, velocity, and intermediate wavespeed

    return r * (S - u) / (S - S_);
}


std::vector<double> f_star(state C, double S, double S_)
{
    // Returns HLLC flux in (left/right) side of the star region, given the
    // state of the (left/right) cell, (left/right) wavespeed, and intermediate
    // wavespeed

    double r = C.r();
    double u = C.u();
    double p = C.p();
    double E = C.E();

    std::vector<double> UStar(7);
    double r_ = r_star(r, S, u, S_);
    double temp = S_ + p / (r*(S-u));
    temp *= S_ - u;
    UStar[0] = r_;
    UStar[1] = r_ * S_;
    UStar[2] = r_ * C.v();
    UStar[3] = r_ * C.w();
    UStar[4] = r_ * (E/r + temp);
    UStar[5] = r_ * C.y();
    UStar[6] = r_ * C.p_inf();

    std::vector<double> U = C.U();
    std::vector<double> F = C.F();
    return F + S * (UStar - U);
}


double q(double pStar, state C)
{
    // Returns weight used in wavespeed estimate, given pressure in star region
    // and state of (left/right) cell

    double p = C.p();
    double y = C.y();
    double p_inf = C.p_inf();

    if (pStar <= p)
        return 1;
    else
        return sqrt( 1 + (y+1)/(2*y) * ((pStar+p_inf)/(p+p_inf) - 1) );
}


triple wavespeeds(state CL, state CR)
{
    // Returns the left, middle, and right wavespeeds for the Riemann problem
    // for neighbouring cells with states CL and CR

    triple S;
    double QUSER = 2;

    double rL = CL.r();
    double rR = CR.r();
    double aL = CL.a();
    double aR = CR.a();
    double pL = CL.p();
    double pR = CR.p();
    double uL = CL.u();
    double uR = CR.u();
    double yL = CL.y();
    double yR = CR.y();
    double p_infL = CL.p_inf();
    double p_infR = CR.p_inf();

    // Compute guess pressure from PVRS Riemann solver
    double CUP  = 0.25 * (rL + rR) * (aL + aR);
    double pPV  = 0.5 * (pL + pR) + 0.5 * (uL - uR) * CUP;
    pPV  = std::max(0., pPV);
    double pMin = std::min(pL, pR);
    double pMax = std::max(pL, pR);
    double qMax = pMax / pMin;

    double p_;
    double u_;
    double TOL = 1e-6;
    double yDiff = std::abs(yL-yR);
    bool PVRS = (qMax <= QUSER) && (pMin <= pPV) && (pPV <= pMax);
    if (PVRS  ||  yDiff>TOL || p_infL || p_infR)
    {
        // Select PRVS Riemann solver
        p_ = pPV;
        u_ = 0.5 * (uL + uR) + 0.5 * (pL - pR) / CUP;
    }
    else
    {
        double y = (yL+yR)/2;  // (=yR)
        if (pPV < pMin)
        {
            // Select Two-Rarefaction Riemann solver
            double G1 = (y-1) / (2*y);
            double PQ  = pow(pL/pR, G1);
            double G4 = 2 / (y-1);
            u_  = (PQ*uL/aL + uR/aR + G4*(PQ-1)) / (PQ/aL + 1/aR);
            double G7 = (y-1) / 2;
            double PTL = 1 + G7 * (uL-u_) / aL;
            double PTR = 1 + G7 * (u_-uR) / aR;
            double G3 = 2 * y / (y-1);
            p_  = 0.5 * (pL*pow(PTL,G3) + pR*pow(PTR,G3));
        }
        else
        {
            // Use Two-Shock Riemann solver with PVRS as estimate
            double G5 = 2 / (y+1);
            double G6 = (y-1) / (y+1);
            double GEL = sqrt((G5/rL) / (G6*pL + pPV));
            double GER = sqrt((G5/rR) / (G6*pR + pPV));
            p_  = (GEL*pL + GER*pR - (uR - uL)) / (GEL + GER);
            u_  = 0.5 * (uL + uR) + 0.5 * (GER*(p_ - pR) - GEL*(p_ - pL));
        }
    }

    S.Star = u_;
    S.L = uL - q(p_,CL)*aL;
    S.R = uR + q(p_,CR)*aR;
    return S;
}


std::vector<double> hllc_flux(state CL, state CR)
{
    // Returns HLLC flux at x/t=0 for the Riemann problem for neighbouring cells
    // with states CL and CR

    triple S = wavespeeds(CL, CR);

    if (0 <= S.L)
        return CL.F();

    else if (S.R <= 0)
        return CR.F();

    else if (0 <= S.Star)
        return f_star(CL, S.L, S.Star);

    else // S.Star <= 0
        return f_star(CR, S.R, S.Star);
}
