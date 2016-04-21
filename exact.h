#ifndef __EXACT__
#define __EXACT__

#include "variables.h"


double p_star(state CL, state CR);

double u_star(double p_, state CL, state CR);

double r_star_shock(double p_, state C);

double r_star_fan(double p_, state C);

state exact_euler(double x, double t, double x0, state CL, state CR);


#endif
