#ifndef __COMMON__
#define __COMMON__

#include <vector>

#include "variables.h"


std::vector<state> add_boundary_conditions(const std::vector<state> & X,
                                           int BCL, int BCR, int size);

double s_max(const std::vector<state> & X);

state evolve_state(state CL, state CM, state CR, double dx, double dt, double w,
                   int l, int side);

double timestep(const std::vector<state> & X, double dx, double C,
                double tCurrent, double tFinal, int iter);


#endif
