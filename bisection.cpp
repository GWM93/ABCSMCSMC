#include <iostream>
#include <cmath>
#include <cassert>
#include <limits>
#include <vector>

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

double f(double x)
{
  return (pow(x,3.0) - x - 2.0);
}

double ESS(std::vector<double> weights){
  double sum = 0, sq_sum = 0;
  for(int i = 0; i < (int) weights.size(); i++){
    sum += weights[i];
    sq_sum += pow(weights[i], 2.0);
  }

  return pow(sum,2.0)/sq_sum;
}

double g(double epsilon_new, double epsilon_old, double alpha, std::vector<double> weights_old,std::vector<double> distances){
  double epsilon = 0;
  assert(weights_old.size() == distances.size());
  int N = weights_old.size();
  std::vector<double> weights(N, 0);

  for(int i = 0; i < N; i++){
    bool d(distances[i] < epsilon_new);
    weights[i] = weights_old[i]*d;
  }

  return ESS(weights) - alpha*ESS(weights_old);
}



double bisection( double (*func) (double) ,double* endpoints, double tol, int NMAX){
  // Assigning endpoint
  double a = endpoints[0], b = endpoints[1];
  // Assert functions satisfy bisection algorithm conditions
  assert(a < b && ((func(a) < 0 && func(b) > 0) || (func(a) > 0 && func(b) < 0)));
  int N = 1;
  double c = 0;
  while(N <= NMAX){ // limit iterations to prevent infinite loop
    c = (a + b)/2;  // new midpoint
    if(func(c) == 0 || ((b-a)/2) < tol){
      return c;
    }
    N++;
    if( sgn(f(c)) == sgn(f(a))){ // new interval
      a = c;
    }
    else{
      b = c;
    }
  }
  std::cout << "Method failed." << std::endl;
  return std::numeric_limits<double>::quiet_NaN();
}

double bisection( double (*func) (double, double, double, std::vector<double>,std::vector<double>),double* endpoints,double epsilon_new, double epsilon_old, double alpha, std::vector<double> weights_old,std::vector<double> distances, double tol, int NMAX){
  // Assigning endpoint
  double a = endpoints[0], b = endpoints[1];
  // Assert functions satisfy bisection algorithm conditions
  assert(a < b && ((func(a,epsilon_old,alpha, weights_old, distances) < 0 && func(b,epsilon_old,alpha, weights_old, distances)) > 0) || (func(a,epsilon_old,alpha, weights_old, distances) > 0 && func(b,epsilon_old,alpha, weights_old, distances) < 0));
  int N = 1;
  double c = 0;
  while(N <= NMAX){ // limit iterations to prevent infinite loop
    c = (a + b)/2;  // new midpoint
    if(func(c,epsilon_old,alpha, weights_old, distances) == 0 || ((b-a)/2) < tol){
      return c;
    }
    N++;
    if(sgn(func(c,epsilon_old,alpha, weights_old, distances)) == sgn(func(a,epsilon_old,alpha, weights_old, distances))){ // new interval
      a = c;
    }
    else{
      b = c;
    }
  }
  std::cout << "Method failed." << std::endl;
  return std::numeric_limits<double>::quiet_NaN();
}

int main()
{
  double x = 1.0;
  double endpoint[2]; endpoint[0] = 1.0; endpoint[1] = 2.0;

  double c = bisection(f,endpoint,1e-10,10000);

  std::cout << "f(c) = " << f(c) << std::endl;


}
