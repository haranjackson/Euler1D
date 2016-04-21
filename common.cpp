#include <cmath>
#include <vector>

#include "slope.h"
#include "variables.h"


std::vector<state> add_boundary_conditions(const std::vector<state> & X,
                                           int BCL, int BCR, int boundarySize)
{
    // Returns X with boundary cells added either side

    int n = X.size();
    std::vector<state> ret(n + 2 * boundarySize);

    if (boundarySize == 1)
    {
        ret[0] = X[0];
        ret[n+1] = X[n-1];

        if (BCL == 1)
            ret[0].setu(-ret[0].u());
        if (BCR == 1)
            ret[n+1].setu(-ret[n+1].u());
    }
    else if (boundarySize == 2)
    {
        ret[0] = X[1];
        ret[1] = X[0];
        ret[n+2] = X[n-1];
        ret[n+3] = X[n-2];

        if (BCL == 1)
        {
            ret[0].setu(-ret[0].u());
            ret[1].setu(-ret[1].u());
        }
        if (BCR == 1)
        {
            ret[n+2].setu(-ret[n+2].u());
            ret[n+3].setu(-ret[n+3].u());
        }
    }
    for (int i=0; i<n; i++)
        ret[i + boundarySize] = X[i];

    return ret;
}


double s_max(const std::vector<state> & X)
{
    // Returns an estimate of the maximum wavespeed across the domain

    double ret = 0;
    for (unsigned int i=0; i<X.size(); i++)
    {
        double temp = std::abs(X[i].u()) + X[i].a();
        if (temp > ret)
            ret = temp;
    }
    return ret;
}


double timestep(const std::vector<state> & X, double dx, double C,
                double tCurrent, double tFinal, int iter)
{
    double dt = C * dx / s_max(X);
    if (iter <= 5)
        dt *= 0.2;
    if (tCurrent+dt > tFinal)
        dt = tFinal - tCurrent;
    return dt;
}


state evolve_state(state CL, state CM, state CR, double dx, double dt, double w,
                   int l, int side)
{
    // Evolves state CM by half a time step according to the gradient across
    // cells with states CL, CM, CR. Returns the state of the middle cell on the
    // specified side.

    std::vector<double> d = slope(CL, CM, CR, w, l);
    std::vector<double> UL = CM.U() - 0.5*d;
    std::vector<double> UR = CM.U() + 0.5*d;
    std::vector<double> dU = 0.5 * (dt/dx) * (F(UL) - F(UR));
    state ret;
    if (side==0)
        ret.setU(UL + dU);
    else
        ret.setU(UR + dU);
    return ret;
}
