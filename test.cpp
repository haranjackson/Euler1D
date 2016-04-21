#include <vector>

#include "variables.h"


std::vector<state> initial_vector(state C0, state C1, state C2, state C3,
                                  double x0, double x1, double x2, double M,
                                  int nInt)
{
    std::vector<state> ret(M);
    for (int i=0; i<M; i++)
    {
        double x = (i+0.5)/M;

        if (x <= x0)
            ret[i] = C0;

        else
        {
            if (nInt == 1 || (nInt > 1 && x <= x1))
                ret[i] = C1;

            else
            {
                if (nInt == 2 || (nInt == 3 && x <= x2))
                    ret[i] = C2;

                else
                    ret[i] = C3;
            }
        }
    }
    return ret;
}


std::vector<double> density(const std::vector<state> & X)
{
    // Returns a vector of the densities in the cells of X

    int n = X.size();
    std::vector<double> ret(n);
    for (int i=0; i<n; i++)
        ret[i] = X[i].r();
    return ret;
}


std::vector<double> velocity(const std::vector<state> & X)
{
    // Returns a vector of the velocities in the cells of X

    int n = X.size();
    std::vector<double> ret(n);
    for (int i=0; i<n; i++)
        ret[i] = X[i].u();
    return ret;
}


std::vector<double> pressure(const std::vector<state> & X)
{
    // Returns a vector of the pressures in the cells of X

    int n = X.size();
    std::vector<double> ret(n);
    for (int i=0; i<n; i++)
        ret[i] = X[i].p();
    return ret;
}


std::vector<double> internal_energy(const std::vector<state> & X)
{
    // Returns a vector of the energies in the cells of X

    int n = X.size();
    std::vector<double> ret(n);
    for (int i=0; i<n; i++)
        ret[i] = X[i].e();
    return ret;
}
