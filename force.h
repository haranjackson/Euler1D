#ifndef __FORCE__
#define __FORCE__

#include <vector>

#include "variables.h"

std::vector<double> richtmyer_flux(state CL, state CR, double dx, double dt);

std::vector<double> force_flux(state CL, state CR, double dx, double dt);

void force_stepper(std::vector<state> & X, const std::vector<state> & X0,
                   double dx, double dt);

std::vector<state> force(const std::vector<state> & X, int BCL, int BCR,
                         double dx, double C, double t);


#endif
