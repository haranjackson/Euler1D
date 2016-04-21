#ifndef __SLIC__
#define __SLIC__

#include <vector>

#include "variables.h"


void slic_stepper(std::vector<state> & X, const std::vector<state> & X0,
                  double dx, double dt, double w, int l);

std::vector<state> slic(const std::vector<state> & X, int BCL, int BCR,
                        double dx, double C, double t, double w, int l);


#endif
