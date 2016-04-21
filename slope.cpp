#include <algorithm>
#include <vector>

#include "variables.h"


double xiR(double r, double w)
{
    return 2 / (1-w + (1+w)*r);
}


double xiSB(double r, double w)
{
    if (r < 0)
        return 0;

    else if (r < 0.5)
        return 2*r;

    else if (r < 1)
        return 1;

    else
        return std::min(std::min(r, xiR(r,w)), 2.);
}


double xiVL(double r, double w)
{
    return (r<0) ? 0 : std::min(2*r/(1+r), xiR(r,w));
}


double xiVA(double r, double w)
{
    return (r<0) ? 0 : std::min(r*(1+r)/(1+r*r), xiR(r,w));
}


double xiMB(double r, double w)
{
    if (r < 0)
        return 0;

    else if (r < 1)
        return r;

    else
        return std::min(1., xiR(r,w));
}


double xi(double r, double w, int l)
{
    assert(0<=l && l<=3);

    switch (l)
    {
        case 0:
            return xiSB(r,w);
        case 1:
            return xiVL(r,w);
        case 2:
            return xiVA(r,w);
        default:
            return xiMB(r,w);
    }
}


std::vector<double> slope(state CL, state C, state CR, double w, int l)
{
    std::vector<double> dUL = C.U() - CL.U();
    std::vector<double> dUR = CR.U() - C.U();
    std::vector<double> d = 0.5*(1+w)*dUL + 0.5*(1-w)*dUR;
    for (unsigned int i=0; i<d.size(); i++)
    {
        double TOL = 1.0e-5;
        double num = bound_below(dUL[i], TOL);
        double den = bound_below(dUR[i], TOL);
        double r = num / den;
        d[i] *= xi(r,w,l);
    }
    return d;
}
