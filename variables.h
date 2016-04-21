#ifndef __VARIABLES__
#define __VARIABLES__

#include <algorithm>
#include <cassert>
#include <functional>
#include <vector>


double bound_below(double x, double TOL);

std::vector<double> F(std::vector<double> U);


class state
{
    public:
        state();
        state(double r, double u, double v, double w, double p, double y,
              double p_inf);
        state(std::vector<double> W);
        double r() const;
        double u() const;
        double v() const;
        double w() const;
        double p() const;
        double y() const;
        double e() const;
        double E() const;
        double a() const;
        double S() const;
        double p_inf() const;
        void setr(double r);
        void setu(double u);
        void setv(double v);
        void setw(double w);
        void setp(double p);
        void sety(double y);
        void setp_inf(double p_inf);
        std::vector<double> W() const;
        std::vector<double> U() const;
        std::vector<double> F() const;
        std::vector<double> G() const;
        std::vector<double> H() const;
        void setW(std::vector<double> W);
        void setU(std::vector<double> U);
        bool operator ==(const state & other);
        bool operator !=(const state & other);
    private:
        double density;
        double pressure;
        double uVelocity;
        double vVelocity;
        double wVelocity;
        double gamma;
        double pressure_inf;
};


class triple
{
    public:
        double L;
        double Star;
        double R;
        triple operator* (double c);
        std::vector<double> toVector() const;
};


// Vector addition
template <typename T>
std::vector<T> operator+ (const std::vector<T>& v1, const std::vector<T>& v2)
{
    assert(v1.size() == v2.size());

    std::vector<T> ret;
    ret.reserve(v1.size());

    std::transform(v1.begin(), v1.end(), v2.begin(), std::back_inserter(ret),
                   std::plus<T>());
    return ret;
}


// Vector subtraction
template <typename T>
std::vector<T> operator- (const std::vector<T>& v1, const std::vector<T>& v2)
{
    assert(v1.size() == v2.size());

    std::vector<T> ret;
    ret.reserve(v1.size());

    std::transform(v1.begin(), v1.end(), v2.begin(), std::back_inserter(ret),
                   std::minus<T>());
    return ret;
}


// Left-multiplication of vectors by scalars
template <class T, class Q>
std::vector<T> operator* (const Q c, const std::vector<T>& v)
{
    std::vector<T> ret = v;

    std::transform(v.begin(), v.end(), ret.begin(),
                   std::bind1st(std::multiplies<T>(),c));
    return ret;
}


// Sign function (sign(1)=1, sign(0)=1, sign(-1)=-1)
template <typename T> int sgn(T val) {

    assert(!isnan(val));
    return (T(0) <= val) - (val < T(0));
}


#endif
