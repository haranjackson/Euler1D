#ifndef __FLIC__
#define __FLIC__

#include <vector>

#include "variables.h"


void flic_stepper(std::vector<state> & X, const std::vector<state> & X0,
                  double dx, double dt, double C, int l);

std::vector<state> flic(const std::vector<state> & X, int BCL, int BCR,
                        double dx, double C, double t, int l);


#endif
