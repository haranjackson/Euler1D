#include <cmath>
#include <vector>

#include "variables.h"


double A(double r, double y)
{
    return 2 / (r * (y+1));
}


double B(double p, double y, double p_inf)
{
    return p * (y-1) / (y+1) + (2*y*p_inf) / (y+1);
}


double f0(double z, state C)
{
    double r = C.r();
    double p = C.p();
    double a = C.a();
    double y = C.y();
    double p_inf = C.p_inf();

    if (z > p)
    {
        double temp = sqrt(A(r,y) / (z+B(p,y,p_inf)));
        return (z-p) * temp;
    }
    else
    {
        double temp = pow((z+p_inf)/(p+p_inf), (y-1)/(2*y));
        return 2 * a / (y-1) * (temp-1);
    }
}


double f1(double z, state C)
{
    double r = C.r();
    double p = C.p();
    double a = C.a();
    double y = C.y();
    double p_inf = C.p_inf();

    if (z > p)
    {
        double temp = sqrt( A(r,y) / (z+B(p,y,p_inf)) );
        return (1 - (z-p) / (2*(z+B(p,y,p_inf))) ) * temp;
    }
    else
    {
        double temp = pow((z+p_inf)/(p+p_inf), -(y+1)/(2*y));
        return temp * a / (y*(p+p_inf));
    }
}


double f(double z, state CL, state CR)
{
    return f0(z,CL) + f0(z,CR) + CR.u() - CL.u();
}


double f_deriv(double z, state CL, state CR)
{
    return f1(z,CL) + f1(z,CR);
}


double CHA(double pk, double pk_1)
{
    return 2 * std::abs(pk-pk_1) / std::abs(pk+pk_1);
}


double p_star(state CL, state CR)
{
    double TOL = 1.0e-6;
    double p0 = (CL.p() + CR.p()) / 2;
    double p1 = 3*p0;                   // CHA(p1,p0)=1 on first iteration

    while (CHA(p1,p0) > TOL)
    {
        p0 = p1;
        p1 = p0 - f(p0, CL, CR) / f_deriv(p0, CL, CR);
        if (p1<0)
            p1 = TOL;
    }

    return p1;
}


double u_star(double p_, state CL, state CR)
{
    return (CL.u() + CR.u() + f0(p_,CR) - f0(p_,CL)) / 2;
}


double Q(double p_, state C)
{
    double r = C.r();
    double p = C.p();
    double y = C.y();
    double p_inf = C.p_inf();
    return sqrt( (p_ + B(p,y,p_inf)) / A(r,y) );
}


std::vector<double> Wfan(double S, state C, double a)
{
    double r = C.r();
    double u = C.u();
    double p = C.p();
    double y = C.y();
    double p_inf = C.p_inf();
    double temp = 2/(y+1) + (y-1)*(u-S)/((y+1)*a);

    std::vector<double> ret(7);
    ret[0] = r * pow(temp, 2/(y-1));
    ret[1] = 2 * (a + (y-1)*u/2 + S) / (y+1);
    ret[2] = 0;
    ret[3] = 0;
    ret[4] = (p+p_inf) * pow(temp, 2*y/(y-1)) - p_inf;
    ret[5] = y;
    ret[6] = C.p_inf();
    return ret;
}


double r_star_shock(double p_, state C)
{
    double r = C.r();
    double p = C.p();
    double y = C.y();
    double p_inf = C.p_inf();
    double temp1 = (p_+p_inf)/(p+p_inf) + (y-1)/(y+1);
    double temp2 = (y-1)/(y+1)*(p_+p_inf)/(p+p_inf) + 1;
    return r * temp1 / temp2;
}


double r_star_fan(double p_, state C)
{
    double r = C.r();
    double p = C.p();
    double y = C.y();
    double p_inf = C.p_inf();
    return r * pow((p_+p_inf)/(p+p_inf), 1/y);
}


double a_star(double p_, state C)
{
    double y = C.y();
    double p = C.p();
    double a = C.a();
    double p_inf = C.p_inf();
    return a * pow((p_+p_inf)/(p+p_inf), (y-1)/(2*y));
}


state exact_euler(double x, double t, double x0, state CL, state CR)
{
    // Returns the exact state solution to the Euler equations at (x,t), given
    // initial states CL for x<x0 and CR for x>x0

    double S = (x-x0)/t;
    double rL = CL.r();
    double rR = CR.r();
    double uL = CL.u();
    double uR = CR.u();
    double pL = CL.p();
    double pR = CR.p();
    double aL = CL.a();
    double aR = CR.a();

    double p_ = p_star(CL, CR);
    double u_ = u_star(p_, CL, CR);

    if (S < u_)
    {
        if (p_ < pL)		// Left fan
        {
            if (S < uL-aL)
                return CL;

            else
            {
                double STL = u_ - a_star(p_,CL);

                if (S < STL)
                    return state(Wfan(S,CL,aL));
                else
                    return state(r_star_fan(p_,CL), u_, 0, 0, p_, CL.y(),
                                 CL.p_inf());
            }
        }
        else				// Left shock
        {
            double SL = uL - Q(p_,CL)/rL;

            if (S < SL)
                return CL;
            else
                return state(r_star_shock(p_,CL), u_, 0, 0, p_, CL.y(),
                             CL.p_inf());
        }
    }
    else
    {
        if (p_ < pR)		// Right fan
        {
            if (uR+aR < S)
                return CR;

            else
            {
                double STR = u_ + a_star(p_,CR);

                if (STR < S)
                    return state(Wfan(S,CR,-aR));
                else
                    return state(r_star_fan(p_,CR), u_, 0, 0, p_, CR.y(),
                                 CR.p_inf());
            }
        }
        else				// Right shock
        {
            double SR = uR + Q(p_,CR)/rR;

            if (SR < S)
                return CR;
            else
                return state(r_star_shock(p_,CR), u_, 0, 0, p_, CR.y(),
                             CR.p_inf());
        }
    }
}
