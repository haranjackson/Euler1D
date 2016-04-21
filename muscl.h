#ifndef __MUSCL__
#define __MUSCL__

#include <vector>

#include "variables.h"


void muscl_stepper(std::vector<state> & X, const std::vector<state> & X0,
                   double dx, double dt, double w, int l);

std::vector<state> muscl(const std::vector<state> & X, int BCL, int BCR,
                         double dx, double C, double t, double w, int l);


#endif
