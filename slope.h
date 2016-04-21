#ifndef __SLOPE__
#define __SLOPE__

#include <vector>

#include "variables.h"


double xiR(double r, double w);

double xiSB(double r, double w);

double xiVL(double r, double w);

double xiVA(double r, double w);

double xiMB(double r, double w);

double xi(double r, double w, int l);

std::vector<double> slope(state CL, state C, state CR, double w, int l);


#endif