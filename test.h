#ifndef __TEST__
#define __TEST__

#include <vector>

#include "variables.h"


std::vector<state> initial_vector(state C0, state C1, state C2, state C3,
                                  double x0, double x1, double x2, double M,
                                  int nInt);

std::vector<double> density(const std::vector<state> & X);

std::vector<double> velocity(const std::vector<state> & X);

std::vector<double> pressure(const std::vector<state> & X);

std::vector<double> internal_energy(const std::vector<state> & X);


#endif
