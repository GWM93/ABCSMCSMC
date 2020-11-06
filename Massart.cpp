// My files
#include "dist.h"
#include "CRN.h"
#include "Gillespie.h"
#include "SIR.h"
#include <limits>

int Okamoto(double epsilon, double delta)
{
  return ceil((1/(2*pow(epsilon,2)))*log(2/delta));
  //return ceil((-1/(2*pow(epsilon,2)))*log(delta / 2));
}

int Massart(double gamma, double epsilon, double delta, double alpha)
{
  return  ceil((2/(h_a(gamma, epsilon)*pow(epsilon,2)))*log(2/(delta - alpha)));
}

void conf_int(int m,int k, double alpha, double (&I) [2])
{
  randomVariable rv;

  double z = rv.InvUnitNormalCdf(1.0 - alpha), z_half = rv.InvUnitNormalCdf(1.0 - alpha/2.0);

  double mean = (m + (pow(z,2.0)/2.0))/(k + pow(z,2.0)), sigma = z_half * sqrt((1.0/(k + pow(z,2.0)))*mean*(1.0 - mean));
  I[0] = mean - sigma;
  I[1] = mean + sigma;
}

bool inRange(double x, double low, double high)
{
  return (low <= x && x <= high);
}

double h_a(double gamma, double epsilon)
{
  if(gamma < 0.5)
  {
    //return (9.0/2.0)*pow((3.0*gamma + epsilon)*(3.0*(1.0 - gamma) - epsilon),-1.0);
    return (9.0/2.0)*(1/((3.0*gamma + epsilon)*(3.0*(1-gamma) - epsilon)));
  }
  else if (0.5 <= gamma && gamma < 1)
  {
    //return (9.0/2.0)*pow((3.0*(1.0 - gamma) + epsilon)*(3.0*gamma + epsilon),-1.0);
    return (9.0/2.0)*(1/((3*(1-gamma) + epsilon)*(3*gamma + epsilon)));
  }
  else
  {
    return 0.0;
  }

  //return (0 < gamma && gamma < 0.5 ) ? (9/2)*pow((3 * gamma + epsilon)*(3*(1 - gamma) - epsilon),-1) : (9/2)*pow(((3*(1 - gamma) + epsilon)*(3*gamma + epsilon)),-1);

}
double h_r(double gamma, double epsilon)
{
  return (gamma < 0.5) ? 9*gamma / 2 * pow((3 + epsilon)*(3 - gamma*(3 + epsilon)),-1) : 9*gamma/2*pow((3 - epsilon)*(3 - gamma*(3 - epsilon)),-1) ;
}


bool check_trace(std::map<double,std::vector<int>> trace) {

  std::map<double,std::vector<int>>::iterator itr;
  std::map<double,std::vector<int>>::iterator upperBound1,upperBound2;
  upperBound1 = trace.upper_bound(100);
  upperBound2 = trace.upper_bound(150);
  bool prop1, prop2;

  if(std::prev(upperBound1,1)->second[1] > 0){
    prop1 = true;
  }
  else{
    // std::prev(upperBound1,1)->second[i] or I = 0
    prop1 = false;
  }

  if(std::prev(upperBound2,1)->second[1] == 0){
    prop2 = true;
  }
  else{
    //std::prev(upperBound2,1)->second[i] or I > 0
    prop2 = false;
  }

  if(prop1 && prop2){
    //std::cout << "\u03C6 = (I > 0) U[100,150] (I = 0) is true" << std::endl;
    return true;
  }
  else{
    //std::cout << "\u03C6 = (I > 0) U[100,150] (I = 0) is false"<< std::endl;
    return false;
  }
}

bool check_trace(std::map<double,std::vector<double>> trace) {

  std::map<double,std::vector<double>>::iterator itr;
  std::map<double,std::vector<double>>::iterator upperBound1,upperBound2;
  upperBound1 = trace.upper_bound(100);
  upperBound2 = trace.upper_bound(150);
  bool prop1, prop2;

  if(std::prev(upperBound1,1)->second[1] > 0){
    prop1 = true;
  }
  else{
    // std::prev(upperBound1,1)->second[i] or I = 0
    prop1 = false;
  }

  if(std::prev(upperBound2,1)->second[1] == 0){
    prop2 = true;
  }
  else{
    //std::prev(upperBound2,1)->second[i] or I > 0
    prop2 = false;
  }

  if(prop1 && prop2){
    //std::cout << "\u03C6 = (I > 0) U[100,150] (I = 0) is true" << std::endl;
    return true;
  }
  else{
    //std::cout << "\u03C6 = (I > 0) U[100,150] (I = 0) is false"<< std::endl;
    return false;
  }
}



std::vector<double> AbsoluteErrorMassart(double epsilon,                        // Absolute Error parameter
                            double delta,                                       // Confidence parameter delta
                            double alpha,                                       // Confidence parameter alpha s.t. alpha < delta
                            CRN srn,                                            // Chemical Reaction Network structure
                            std::vector<double> theta,                          // Parameters for propensity function
                            double T,                                           // Running time of algorithm, T
                            randomVariable& rv )                                // Random variable
{
  // Worst case sample size determined by Okamoto Bound
  int M = Okamoto(epsilon,delta), n_k = M;
  int m = 0; // Number of successes
  int k = 0;

  // Initiale Confidence interval
  double a_k = 0, b_k = 1;
  double I_k[2];
  I_k[0] = a_k; I_k[1] = b_k;

  std::map<double,std::vector<int>> omega_k;
  bool z;
  while(k < n_k)
  {
    k++;
    //Generate omega(k)
    omega_k = gillespie(srn,theta,T,rv);
    // Evalute whether z(omega(k)) = 1(omega(k) |= phi)
    z = check_trace(omega_k);
    m += z;

    conf_int(m,k,alpha, I_k);
    a_k = I_k[0];
    b_k = I_k[1];
    if(inRange(0.5,a_k, b_k))
    {
      n_k = M;
    }
    else if(b_k < 0.5)
    {
      n_k = Massart(b_k, epsilon, delta, alpha);
    }
    else
    {
      n_k = Massart(a_k, epsilon, delta, alpha);
    }
    n_k = std::min(n_k, M);
  }
  std::cout << "n = " << n_k << std::endl;
  return {(double) m / k,a_k,b_k};
}

std::vector<double> AbsoluteErrorMassart(double epsilon,                        // Absolute Error parameter
                            double delta,                                       // Confidence parameter delta
                            double alpha,                                       // Confidence parameter alpha s.t. alpha < delta
                            CRN srn,                                            // Chemical Reaction Network structure
                            std::vector<double> theta,                          // Parameters for propensity function
                            double T,                                           // Running time of algorithm, T
                            randomVariable& rv,                                 // Random variable
                            double epsilon_i,                                   // Current threshold
                            std::map<double,std::vector<double>> X)             // True data X

{
  // Worst case sample size determined by Okamoto Bound
  int M = Okamoto(epsilon, delta), n_k = Okamoto(epsilon, delta);
  int m = 0; // Number of successes
  int k = 0; // Number of iterations

  // Initial Confidence interval
  double a_k = 0, b_k = 1;
  double I_k[2];
  I_k[0] = a_k; I_k[1] = b_k;

  double timestep = 0;
  // Calculate timestep needed from data
  timestep = std::next(X.begin(),1)->first - X.begin()->first;

  std::map<double,std::vector<int>> omega_k, Y;
  bool z;
  int b_i = 0;                                                                  // Number of traces that satisfy distance
  double distance_i = 0;
  while(k < n_k)
  {
    k++;
    //Generate omega(k)
    omega_k = gillespie(srn,theta,T,rv);
    // Organise trace
    Y = organiseTrace(omega_k,timestep);

    if(d(X,Y) < epsilon_i)
    {// Calculate K_{h_m} or indicator function as to whether |s - s_{obs}| < epsilon_is is true
      b_i++;
    }
    // Calculate distance between simulated and observed data
    distance_i += d(X,Y);
    // Evalute whether z(omega(k)) = 1(omega(k) |= phi)
    z = check_trace(omega_k);
    m += z;

    conf_int(m,k,alpha, I_k);
    a_k = I_k[0];
    b_k = I_k[1];

    if(inRange(0.5,a_k, b_k))
    {
      n_k = M;
    }
    else if(b_k < 0.5)
    {
      n_k = Massart(b_k, epsilon, delta, alpha);
    }
    else
    {
      n_k = Massart(a_k, epsilon, delta, alpha);
    }

    n_k = std::min(n_k, M);

  }
  std::cout << "n = " << n_k << std::endl;
  return {(double) m / k,(double) b_i / k, (double) distance_i/k, a_k, b_k};
}

double RelativeErrorMassart(double epsilon,                                     // Absolute Error parameter
                            double delta,                                       // Confidence parameter delta
                            double alpha,                                       // Confidence parameter alpha s.t. alpha < delta
                            CRN srn,                                            // Chemical Reaction Network structure
                            std::vector<double> theta,                          // Parameters for propensity function
                            double T,                                           // Running time of algorithm, T
                            randomVariable& rv )                                // Random variable
{
  // Worst case sample size determined by Okamoto Bound
  double gamma_min = 0.00001;                    // Machine epsilon
  long int M = ceil((1/(pow(epsilon,2)*h_r(gamma_min,epsilon)))*log(2/delta)), n_k = M;
  M = ceil(log(2/delta)*(1/(pow(epsilon,2)*h_r(gamma_min,epsilon))));
  int m = 0; // Number of successes
  int k = 0;

  // Initiale Confidence interval
  double a_k = 0, b_k = 1;
  double I_k[2];
  I_k[0] = gamma_min; I_k[1] = b_k;


  std::map<double,std::vector<int>> omega_k;
  bool z;
  while(k < n_k)
  {
    k++;
    //Generate omega(k)
    omega_k = gillespie(srn,theta,T,rv);
    // Evalute whether z(omega(k)) = 1(omega(k) |= phi)
    z = check_trace(omega_k);
    m += z;

    conf_int(m,k,alpha, I_k);
    a_k = I_k[0];
    b_k = I_k[1];

    if(gamma_min >= a_k)
    {
      n_k = M;
    }
    else
    {
      n_k = ceil((1/(pow(epsilon,2)*h_r(a_k,epsilon)))*log(2/(delta - alpha)));
    }

    n_k = std::min(n_k, M);
  }
  return (double) m / k;
}

std::vector<double> AbsoluteErrorMassart(double epsilon, double delta, double alpha, CRN srn, std::vector<double> theta, double T, randomVariable& rv, bool (property_checker)(std::map<double,std::vector<int>>) ){
  // Worst case sample size determined by Okamoto Bound
  int M = Okamoto(epsilon,delta), n_k = M;
  int m = 0; // Number of successes
  int k = 0;

  // Initiale Confidence interval
  double a_k = 0, b_k = 1;
  double I_k[2];
  I_k[0] = a_k; I_k[1] = b_k;

  std::map<double,std::vector<int>> omega_k;
  bool z;
  while(k < n_k)
  {
    k++;
    //Generate omega(k)
    omega_k = gillespie(srn,theta,T,rv);
    // Evalute whether z(omega(k)) = 1(omega(k) |= phi)
    z = property_checker(omega_k);
    m += z;

    conf_int(m,k,alpha, I_k);
    a_k = I_k[0];
    b_k = I_k[1];
    if(inRange(0.5,a_k, b_k))
    {
      n_k = M;
    }
    else if(b_k < 0.5)
    {
      n_k = Massart(b_k, epsilon, delta, alpha);
    }
    else
    {
      n_k = Massart(a_k, epsilon, delta, alpha);
    }
    n_k = std::min(n_k, M);

  }

  return {(double) m / k,a_k,b_k, (double) k};
}
