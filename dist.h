#ifndef DIST_H_
#define DIST_H_

#include <iostream>
#include <vector>

// GSL libraries
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_cdf.h>

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

class randomVariable
{
public:
  // Constructor
  randomVariable();
  // Destructor
  ~randomVariable();

  // Return uniform distribution U ~ (a,b)
  double uniform(double a, double b);
  double uniform_pdf(double x, double a, double b);

  // Gaussian distribution and pdf
  double Normal(double mean, double sigma);
  double Normal_pdf(double x, double mean, double sigma);

  double UnitNormalCdf(double x);
  double InvUnitNormalCdf(double x);

  double BetaCdf(double x, double a, double b);

  // Multivariate Gaussian distribution and pdf
  std::vector<double> mNormal(std::vector<double> mean, std::vector<std::vector<double>> covariance);
  double mNormal_pdf(std::vector<double> x, std::vector<double> mean, std::vector<std::vector<double>> covariance);

private:
  gsl_rng * r;
};

#endif
