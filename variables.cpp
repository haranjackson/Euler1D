#define NDEBUG

#include <QDebug>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

#include "variables.h"


double bound_below(double x, double TOL)
{
    return (std::abs(x)<TOL) ? copysign(TOL,x) : x;
}


double e(double r, double p, double y, double p_inf)
{
    // Returns specific internal energy of a gas obeying the stiffened gas EOS.
    // NB p_inf=0 gives the EOS of an ideal gas.

    return (p + y * p_inf) / ((y - 1)*r);
}


double p(double r, double e, double y, double p_inf)
{
    // Returns pressure of a gas obeying the stiffened gas EOS.
    // NB p_inf=0 gives the EOS of an ideal gas.

    return e * (y - 1) * r - y * p_inf;
}


double a(double r, double p, double y, double p_inf)
{
    // Returns speed of sound of a gas obeying the stiffened gas EOS.
    // NB p_inf=0 gives the EOS of an ideal gas.

    return sqrt((y * (p + p_inf)) / r);
}


double S(double r, double p, double y)
{
    // Returns the entropy of an ideal gas, given density, pressure, and gamma
    // constant

    return p / pow(r,y);
}


std::vector<double> F(std::vector<double> U)
{
    // Returns flux vector in the x direction, given vector of conserved
    // variables

    double r = U[0];
    double u = U[1] / r;
    double v = U[2] / r;
    double w = U[3] / r;
    double y = U[5] / r;
    double p_inf = U[6] / r;
    double V2 = pow(u, 2) + pow(v, 2) + pow(w, 2);
    double pressure = p(r, U[4]/r - V2/2, y, p_inf);

    std::vector<double> ret(7);
    ret[0] = U[1];                      // ru
    ret[1] = U[1] * u + pressure;       // ru^2+p
    ret[2] = U[1] * v;                  // ruv
    ret[3] = U[1] * w;                  // ruw
    ret[4] = u * (U[4] + pressure);     // u(E+p)
    ret[5] = U[1] * y;                  // ruy
    ret[6] = U[1] * p_inf;              // rup_inf
    return ret;
}


std::vector<double> G(std::vector<double> U)
{
    // Returns flux vector in the y direction, given vector of conserved
    // variables

    double r = U[0];
    double u = U[1] / r;
    double v = U[2] / r;
    double w = U[3] / r;
    double y = U[5] / r;
    double p_inf = U[6] / r;
    double V2 = pow(u, 2) + pow(v, 2) + pow(w, 2);
    double pressure = p(r, U[4]/r - V2/2, y, p_inf);

    std::vector<double> ret(7);
    ret[0] = U[2];                      // rv
    ret[1] = U[2] * u;                  // ruv
    ret[2] = U[2] * v + pressure;       // rv^2 + p
    ret[3] = U[2] * w;                  // rvw
    ret[4] = v * (U[4] + pressure);     // v(E+p)
    ret[5] = U[2] * y;                  // rvy
    ret[6] = U[2] * p_inf;              // rvp_inf
    return ret;
}


std::vector<double> H(std::vector<double> U)
{
    // Returns flux vector in the z direction, given vector of conserved
    // variables

    double r = U[0];
    double u = U[1] / r;
    double v = U[2] / r;
    double w = U[3] / r;
    double y = U[5] / r;
    double p_inf = U[6] / r;
    double V2 = pow(u, 2) + pow(v, 2) + pow(w, 2);
    double pressure = p(r, U[4]/r - V2/2, y, p_inf);

    std::vector<double> ret(7);
    ret[0] = U[3];                      // rw
    ret[1] = U[3] * u;                  // ruw
    ret[2] = U[3] * v;                  // rvw
    ret[3] = U[3] * w + pressure;       // rw^2 + p
    ret[4] = w * (U[4] + pressure);     // w(E+p)
    ret[5] = U[3] * y;                  // rwy
    ret[6] = U[3] * p_inf;              // rwp_inf
    return ret;
}


state::state()
{
    // Construct a state object with all variables set to 0

    density = 1;
    pressure = 1;
    uVelocity = 0;
    vVelocity = 0;
    wVelocity = 0;
    gamma = 1;
    pressure_inf = 0;
}


state::state(double r, double u, double v, double w, double p, double y,
             double p_inf)
{
    // Construct a state object with primitive variables

    assert(r >= 0);

    density = r;
    uVelocity = u;
    vVelocity = v;
    wVelocity = w;
    pressure = p;
    gamma = y;
    pressure_inf = p_inf;
}


state::state(std::vector<double> W)
{
    // Construct a state object with vector of primitive variables

    assert(W[0] >= 0);

    density = W[0];
    uVelocity = W[1];
    vVelocity = W[2];
    wVelocity = W[3];
    pressure = W[4];
    gamma = W[5];
    pressure_inf = W[6];
}


double state::r() const
{
    return density;
}


double state::u() const
{
    return uVelocity;
}


double state::v() const
{
    return vVelocity;
}


double state::w() const
{
    return wVelocity;
}


double state::p() const
{
    return pressure;
}


double state::y() const
{
    return gamma;
}


double state::p_inf() const
{
    return pressure_inf;
}


double state::e() const
{
    return ::e(density,pressure,gamma,pressure_inf);
}


double state::E() const
{
    double V2 = pow(uVelocity,2) + pow(vVelocity,2) + pow(wVelocity,2);
    return density * (V2/2 + e());
}


double state::a() const
{
    return ::a(density,pressure,gamma,pressure_inf);
}


double state::S() const
{
    return ::S(density,pressure,gamma);
}


void state::setr(double r)
{
    assert(r >= 0);
    density = r;
}


void state::setu(double u)
{
    uVelocity = u;
}


void state::setv(double v)
{
    vVelocity = v;
}


void state::setw(double w)
{
    wVelocity = w;
}


void state::setp(double p)
{
    assert(p>=0);
    pressure = p;
}


void state::sety(double y)
{
    gamma = y;
}


void state::setp_inf(double p_inf)
{
    pressure_inf = p_inf;
}


std::vector<double> state::W() const
{
    // Returns vector of primitive variables for the state

    std::vector<double> ret(7);
    ret[0] = density;
    ret[1] = uVelocity;
    ret[2] = vVelocity;
    ret[3] = wVelocity;
    ret[4] = pressure;
    ret[5] = gamma;
    ret[6] = pressure_inf;
    return ret;
}


std::vector<double> state::U() const
{
    // Returns vector of conserved variables for the state

    std::vector<double> ret(7);
    ret[0] = density;
    ret[1] = density * uVelocity;
    ret[2] = density * vVelocity;
    ret[3] = density * wVelocity;
    ret[4] = E();
    ret[5] = density * gamma;
    ret[6] = density * pressure_inf;
    return ret;
}


std::vector<double> state::F() const
{
    // Returns flux vector for the state in the x direction, F(U)

    return ::F(U());
}


std::vector<double> state::G() const
{
    // Returns flux vector for the state in the y direction, G(U)

    return ::G(U());
}


std::vector<double> state::H() const
{
    // Returns flux vector for the state in the z direction, H(U)

    return ::H(U());
}


void state::setW(std::vector<double> W)
{
    // Sets state parameters given vector of primitive variables

    assert(W[0] >= 0);
    density = W[0];
    uVelocity = W[1];
    vVelocity = W[2];
    wVelocity = W[3];
    pressure = W[4];
    gamma = W[5];
    pressure_inf = W[6];
}


void state::setU(std::vector<double> U)
{
    // Sets state parameters given vector of conserved variables

    assert(U[0] >= 0);
    density = U[0];
    uVelocity = U[1] / density;
    vVelocity = U[2] / density;
    wVelocity = U[3] / density;
    gamma = U[5] / density;
    pressure_inf = U[6] / density;
    double V2 = pow(uVelocity,2) + pow(vVelocity,2) + pow(wVelocity,2);
    pressure = ::p(density, U[4]/density - V2/2, gamma, pressure_inf);
}


bool state::operator ==(const state & other)
{
    return ((density == other.density) && (pressure == other.pressure) &&
            (uVelocity == other.uVelocity) && (vVelocity == other.vVelocity) &&
            (wVelocity == other.wVelocity) && (gamma == other.gamma) &&
            (pressure_inf == other.pressure_inf));
}


bool state::operator !=(const state & other)
{
    return !(*this == other);
}


triple triple::operator* (double c)
{
    L *= c;
    Star *= c;
    R *= c;
    return *this;
}


std::vector<double> triple::toVector() const
{
    std::vector<double> ret(3);
    ret[0] = L;
    ret[1] = Star;
    ret[2] = R;
    return ret;
}
