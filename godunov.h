#ifndef __GODUNOV__
#define __GODUNOV__

#include <vector>

#include "variables.h"


void godunov_stepper(std::vector<state> & X, const std::vector<state> & X0,
                     double dx, double dt);

std::vector<state> godunov(const std::vector<state> & X, int BCL, int BCR,
                           double dx, double C, double t);


#endif
