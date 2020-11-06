#include "SIR.h"
#include <random>
#include <chrono>


int main()
{
  using clock = std::chrono::system_clock;
  using sec = std::chrono::duration<double>;

  const auto before = clock::now();

  // Running ABCSMC for SIR case study
  randomVariable rv;
  SIR_ABCSMCSMC(rv);
  //SIR_paramSynth(rv);
  //SIR_distance(rv);
  //SIR_Massart(rv);
  //SIR_BayesSMC(rv);
  //SIR_paramSynth_BayesianSMC(rv);


  return 0;
}
