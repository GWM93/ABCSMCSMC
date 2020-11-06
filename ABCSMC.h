#ifndef ABCSMC_H_
#define ABCSMC_H_

#include "dist.h"
#include "CRN.h"
#include "Gillespie.h"
#include "Massart.h"
#include "BayesianSMC.h"
// Function to calculate and save data produced from given parameters
std::map<double,std::vector<double>> trueData(CRN srn,std::vector<double> theta, double T,randomVariable& rv, double timestep,std::string directory, int N);

// Function to calculate weighted mean and covariance used in the perturbation kernel
void calculateCovariance(std::vector<std::vector<double>>& covariance, std::vector<std::vector<double>> theta_prev, std::vector<double> weights);

// Perturbation kernel that returns sampled parameters from given mean, covariance and weights
std::vector<double> K_t(std::vector<double> Mean, std::vector<std::vector<double>> covariance, double (*prior_pdf_pt) (std::vector<double>, randomVariable& ), randomVariable& rv);

// Probability density function of perturbation kernel
double K_t_pdf(std::vector<double> theta_star,std::vector<double> theta_old, std::vector<std::vector<double>> covariance, randomVariable& rv);
// Function to normalise the weights
void normaliseWeights(std::vector<double>& weights);

// Effective sample size
double ESS(std::vector<double> weights);
// Functions to resample
double g(double epsilon_new, double epsilon_old, double alpha, std::vector<double> weights_old,std::vector<double> distances);
double bisection( double (*func) (double, double, double, std::vector<double>,std::vector<double>),double* endpoints, double epsilon_old, double alpha, std::vector<double> weights_old,std::vector<double> distances, double tol, int NMAX);

// Sample random parameter from given weights and parameter
std::vector<double> randomParameter(std::vector<double> weights, std::vector<std::vector<double>> theta, randomVariable& rv);



double quantile(std::vector<double> distances, double p);

// Approximate Bayesian Computation Sequential Monte Carlo algorithm
void ABCSMC(std::vector<std::vector<double>>& theta_out,                         // Outputted sampled parameters
            std::vector<double>& weights_out,                                   // corresponding weights
            std::map<double,std::vector<double>> X,                         // True data
            CRN srn,                                                            // Chemical Reaction Network of choice
            std::vector<double> epsilon,                                        // Discrepancy threshold
            int N,                                                              // Number of particles
            int B_t,
            double T,
            randomVariable& rv);
// Approximate Bayesian Computation Sequential Monte Carlo algorithm
void ABCSMCSMC(std::vector<std::vector<double>>& theta_out,                         // Outputted sampled parameters
            std::vector<double>& weights_out,                                   // corresponding weights
            std::map<double,std::vector<double>> X,                         // True data
            CRN srn,                                                            // Chemical Reaction Network of choice
            std::vector<double> epsilon,                                        // Discrepancy threshold
            int N,                                                              // Number of particles
            double T,
            randomVariable& rv);
#endif
