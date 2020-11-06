#include "dist.h"



// Constructor
randomVariable::randomVariable()
  {
     r = gsl_rng_alloc (gsl_rng_taus);
  };

// Destructor
  randomVariable::~randomVariable()
  {
    gsl_rng_free(r);
  }

double randomVariable::uniform(double a, double b)
  { // Return uniform distribution U ~ (a,b)
    return gsl_ran_flat(r,a,b);
  }

double randomVariable::uniform_pdf(double x, double a, double b)
{
  return gsl_ran_flat_pdf(x, a, b);
}

// Gaussian distribution and pdf
double randomVariable::Normal(double mean, double sigma){
  return mean + gsl_ran_gaussian(r, sigma);
}
double randomVariable::Normal_pdf(double x, double mean, double sigma){
  return gsl_ran_gaussian_pdf(x - mean, sigma);
}

double randomVariable::UnitNormalCdf(double x)
{
  return gsl_cdf_ugaussian_P(x);
}

double randomVariable::InvUnitNormalCdf(double x)
{
  return gsl_cdf_ugaussian_Pinv(x);
}

double randomVariable::BetaCdf(double x, double a, double b){
  return gsl_cdf_beta_P(x,a,b);
}

// Multivariate Gaussian distribution and pdf
std::vector<double> randomVariable::mNormal(std::vector<double> mean, std::vector<std::vector<double>> covariance){
  int n = (int) mean.size();
  // Defining gsl vector and matrix for the mean vector and covariance matrix respectively
  std::vector<double> sample(n,0);
  gsl_vector * gsl_mean = gsl_vector_alloc(n), *gsl_sample = gsl_vector_alloc(n);
  gsl_matrix * gsl_covariance = gsl_matrix_alloc(n,n);

  // Setting elements of vector v, mean vector
  for(int i = 0; i < n; i++){
    gsl_vector_set(gsl_mean,i,mean[i]);
  }
  // Setting elements for matrix M, covariance matrix
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      gsl_matrix_set(gsl_covariance,i,j,covariance[i][j]);
    }
  }

  //Cholesky decomposition for multivariate gaussian sample
  gsl_linalg_cholesky_decomp1(gsl_covariance);

    for(int i = 0; i < n; i++){
      for(int j = 0; j < n; j++){
        if( i < j){
          gsl_matrix_set(gsl_covariance,i,j,0.0);
        }
      }
  }

  // Sample from multivariate normal
  gsl_ran_multivariate_gaussian(this->r, gsl_mean,gsl_covariance,gsl_sample);

  for(int i = 0; i < n; i++){
    sample[i] = gsl_vector_get(gsl_sample,i);
  }
  gsl_vector_free (gsl_mean);
  gsl_vector_free (gsl_sample);
  gsl_matrix_free (gsl_covariance);

  return sample;
}

double randomVariable::mNormal_pdf(std::vector<double> x, std::vector<double> mean, std::vector<std::vector<double>> covariance){
  int n = (int) mean.size();
  // Defining gsl vector and matrix for the mean vector and covariance matrix respectively
  gsl_vector * gsl_x = gsl_vector_alloc(n), * gsl_mean = gsl_vector_alloc(n), *gsl_work = gsl_vector_alloc(n);
  gsl_matrix * gsl_covariance = gsl_matrix_alloc(n,n);

  double result;

  // Setting elements of vector v, mean vector
  for(int i = 0; i < n; i++){
    gsl_vector_set(gsl_mean,i,mean[i]);
    gsl_vector_set(gsl_x, i, mean[i]);
  }
  // Setting elements for matrix M, covariance matrix
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      gsl_matrix_set(gsl_covariance,i,j,covariance[i][j]);
    }
  }

  // Cholesky decomposition for multivariate gaussian sample
  gsl_linalg_cholesky_decomp1(gsl_covariance);
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      if( i < j){
        gsl_matrix_set(gsl_covariance,i,j,0.0);
      }
    }
  }

  gsl_ran_multivariate_gaussian_pdf(gsl_x, gsl_mean, gsl_covariance, &result, gsl_work);
  gsl_vector_free(gsl_x);
  gsl_vector_free (gsl_mean);
  gsl_vector_free (gsl_work);
  gsl_matrix_free (gsl_covariance);

  return result;
}
