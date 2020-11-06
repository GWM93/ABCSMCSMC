#include <iostream>
#include <cmath>

// My files
#include "dist.h"
#include "CRN.h"
#include "Gillespie.h"
#include "SIR.h"

void conf_int(int m,int k, double alpha,double (&I) [2]);

bool inRange(double x, double low, double high);

int Okamoto(double epsilon, double delta);


int Massart(double gamma, double epsilon, double delta, double alpha);

double h_a(double gamma, double epsilon);
double h_r(double gamma, double epsilon);

// Function to track and check whether a trace satisfies a given property
bool check_trace(std::map<double,std::vector<int>> trace);
bool check_trace(std::map<double,std::vector<double>> trace);

std::vector<double> AbsoluteErrorMassart(double epsilon,                        // Absolute Error parameter
                            double delta,                                       // Confidence parameter delta
                            double alpha,                                       // Confidence parameter alpha s.t. alpha < delta
                            CRN srn,                                            // Chemical Reaction Network structure
                            std::vector<double> theta,                          // Parameters for propensity function
                            double T,                                           // Running time of algorithm, T
                            randomVariable& rv );                                // Random variable



std::vector<double> AbsoluteErrorMassart(double epsilon,                        // Absolute Error parameter
                            double delta,                                       // Confidence parameter delta
                            double alpha,                                       // Confidence parameter alpha s.t. alpha < delta
                            CRN srn,                                            // Chemical Reaction Network structure
                            std::vector<double> theta,                          // Parameters for propensity function
                            double T,                                           // Running time of algorithm, T
                            randomVariable& rv,                                 // Random variable
                            double epsilon_i,                                   // Current threshold
                            std::map<double,std::vector<double>> X);             // True data X


double RelativeErrorMassart(double epsilon,                                     // Absolute Error parameter
                            double delta,                                       // Confidence parameter delta
                            double alpha,                                       // Confidence parameter alpha s.t. alpha < delta
                            CRN srn,                                            // Chemical Reaction Network structure
                            std::vector<double> theta,                          // Parameters for propensity function
                            double T,                                           // Running time of algorithm, T
                            randomVariable& rv );                                // Random variable

std::vector<double> AbsoluteErrorMassart(double epsilon, double delta, double alpha, CRN srn, std::vector<double> theta, double T, randomVariable& rv, bool (property_checker)(std::map<double,std::vector<int>>) );
