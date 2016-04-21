#ifndef __HLLC__
#define __HLLC__

#include <vector>

#include "variables.h"


double r_star(double r, double S, double u, double SStar);

std::vector<double> f_star(state C, double S, double S_);

triple wavespeeds(state CL, state CR);

std::vector<double> hllc_flux(state CL, state CR);


#endif