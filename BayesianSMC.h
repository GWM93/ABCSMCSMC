#include <iostream>
#include <cmath>

// My files
#include "dist.h"
#include "CRN.h"
#include "Gillespie.h"
#include "SIR.h"
#include "Massart.h"


// Function to track and check whether a trace satisfies a given property
double PosteriorProb(double t_0, double t_1,double alpha, double beta, randomVariable& rv);

std::vector<double> BayesianSMC(double delta, double c, CRN srn, std::vector<double> theta, double T, randomVariable& rv,double alpha=1, double beta=1);
