#ifndef __WAF__
#define __WAF__

#include <vector>

#include "variables.h"


void waf_stepper(std::vector<state> & X, const std::vector<state> & X0,
                 double dx, double dt, int l);

std::vector<state> waf(const std::vector<state> & X, int BCL, int BCR,
                       double dx, double C, double t, int l);


#endif
