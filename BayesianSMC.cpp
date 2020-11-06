// My files
#include "dist.h"
#include "CRN.h"
#include "Gillespie.h"
#include "Massart.h"
#include "SIR.h"
#include <limits>


double PosteriorProb(double t_0, double t_1,double alpha, double beta, randomVariable& rv){
return rv.BetaCdf(t_1, alpha, beta) - rv.BetaCdf(t_0, alpha, beta);
}

std::vector<double> BayesianSMC(double delta, double c, CRN srn, std::vector<double> theta, double T, randomVariable& rv, double alpha, double beta){

  int x = 0; // Number of successes
  int n = 0; // Number of simulations
  double t_0 = 0.0, t_1 = 1.0;

  double p_hat = 0;
  double gamma = 0;
  std::map<double,std::vector<int>> omega_k;
  bool z;

  do{
    n++;
    // Simulate trace
    omega_k = gillespie(srn, theta, T, rv);
    // Check trace
    z = check_trace(omega_k);
    x += z;
    p_hat = (x + alpha) / (n + alpha + beta);

    t_0 = p_hat - delta;
    t_1 = p_hat + delta;

    if(t_1 > 1){
      t_0 = 1-2*delta;
      t_1 = 1;
    }
    else if(t_0 < 0){
      t_0 = 0;
      t_1 = 2*delta;
    }
    gamma = PosteriorProb(t_0, t_1,alpha + x,beta + n - x,rv);
}while (gamma < c);

std::cout << "n = " << n << std::endl;
  return {p_hat,t_0,t_1};
}
