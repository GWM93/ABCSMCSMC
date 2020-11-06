#include "ABCSMC.h"

std::map<double,std::vector<double>> trueData(CRN srn,std::vector<double> theta, double T,randomVariable& rv, double timestep,  std::string directory, int N){
  /****************************************************************************
    GENERATING TRUE DATA
  *****************************************************************************/
  // String for directory

  std::map<double,std::vector<int>> trace;                                      // Ouput of Gillespie / SSA put into trace
  std::map<double,std::vector<double>> meanTrace;                               // Trace is reordered to be observed every timestep seconds of simulation to make it easier to analyse
  std::set<std::map<double,std::vector<int>>> traceSetX;                        // Store traces in set
  std::vector<int> initialState = srn.retInitial();
  std::string tempDirectory = directory;

  int smc_1 = 0, smc_2 = 0;
  for(int i = 0; i < N; i++){
    trace = gillespie(srn,theta,T,rv);
    smc_1 += check_trace(trace);
    std::cout << "check_trace = " << check_trace(trace) << std::endl;
    // printTrace(trace, srn.retNames());
    trace = organiseTrace(trace,timestep);
    smc_2 += check_trace(trace);
    traceSetX.insert(trace);
    tempDirectory = tempDirectory + "traces/trace" + std::to_string(i) + ".csv";
    outputTrace(tempDirectory,trace,srn.retNames());
    tempDirectory = directory;
  }
  meanTrace = meanTraces(traceSetX,timestep);
  tempDirectory = tempDirectory + "meanTrace.csv";
  outputTrace(tempDirectory,meanTrace,srn.retNames());

  std::cout << "SMC for true data = " << (double) smc_1 / (double) N << std::endl;
  std::cout << "SMC for true data observed at timestep = " << timestep << "; " << (double) smc_2 / N << std::endl;
  // Storing true data in map X
  std::map<double,std::vector<double>> X;                                       // Trace is reordered to be observed every 0.5 seconds of simulation to make it easier to analyse
  X = meanTrace;
  std::cout << "Saved true traces and mean in directory: " << directory << std::endl;
  return X;
}

std::vector<double> randomParameter(std::vector<double> weights, std::vector<std::vector<double>> theta, randomVariable& rv){
  // Function returns a random particle i from 0 <= i < particles.size() with the probabilities given according to weights
  double prob =  rv.uniform(0,1);           // Create uniform distribution from (0, T)
  double weightSum = 0, weightCumSum = 0;
  int i = 0;                                                                    // Random integer to return depenedent on distribution of weights
  for(int ii = 0; ii < (int) weights.size(); ii++)
  {
    weightSum += weights[ii];
  }

  for(int ii = 0; ii < (int) weights.size(); ii++)
  {
    if(weightCumSum/weightSum <= prob && prob < (weightCumSum + weights[ii])/weightSum)
    {
      i = ii;
    }
    weightCumSum += weights[ii];
  }
  return theta[i];
}

void calculateCovariance(std::vector<std::vector<double>>& covariance, std::vector<std::vector<double>> theta_prev, std::vector<double> weights){

  // Number of particles
  int N = weights.size();
  // Number of parameters
  int n_theta = theta_prev[0].size();
  // Mean
  std::vector<double> mean(n_theta,0);

  // Temporary covariance function
  std::vector<std::vector<double>> C(covariance.size(),std::vector<double> (covariance[0].size(),0));

  // Calculate mean
  for(int p = 0; p < n_theta; p++){
    for(int i = 0; i < N; i++){
      mean[p] += weights[i]*theta_prev[i][p];
      }
  }
  double weightsSquared = 0;
  for(int n = 0; n < N; n++){
    weightsSquared += pow(weights[n],2.0);
  }


  for(int i = 0; i < n_theta; i++){
    for(int j = 0; j < n_theta; j++){
      for(int n = 0; n < N; n++){
        C[i][j] += weights[n] * (theta_prev[n][j] - mean[j])*(theta_prev[n][i] - mean[i]);
      }
    }
  }

  for(int i = 0; i < n_theta; i++){
    for(int j = 0; j < n_theta; j++){
      C[i][j] += (C[i][j] / (1 - weightsSquared));
    }
  }

  // Ensuring positive definitiveness, alternatively could use covariance matrix s.t. cov(i,j) = 0 for i != j
  for(int i = 0; i < n_theta; i++){
    for(int j = 0; j < n_theta; j++){
      C[i][j] *= 2;
      if(i == j)
      {
        C[i][j] += 1e-10;
      }
    }
  }
  covariance =  C;
}


std::vector<double> K_t(std::vector<double> weightedMean, std::vector<std::vector<double>> covariance, double (*prior_pdf_pt) (std::vector<double>, randomVariable& ), randomVariable& rv)
{
    // Sampled parameters
  std::vector<double> theta_star(weightedMean.size(),0);

  double pdf = 0;
  // Repeat sampling until sample is within the prior probability distribution / parameter bounds
  while(pdf == 0){
    theta_star = rv.mNormal(weightedMean, covariance);
    pdf = prior_pdf_pt(theta_star,rv);
  }

  return theta_star;
}

double K_t_pdf(std::vector<double> theta_star,std::vector<double> theta_old, std::vector<std::vector<double>> covariance, randomVariable& rv){
  double pdf = 0;
  pdf = rv.mNormal_pdf(theta_star, theta_old, covariance);
  return pdf;
}



void normaliseWeights(std::vector<double>& weights){
  double sum = 0;
  std::vector<double> temp = weights;
  for(int i = 0; i < (int) temp.size(); i++){
    sum += temp[i];
  }

  for(int i = 0; i < (int) weights.size(); i++){
    weights[i] /= sum;
  }
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
  assert(weights_old.size() == distances.size());
  int N = weights_old.size();
  std::vector<double> weights(N, 0);
  double num = 0, denom = (double) N;
  for(int i = 0; i < N; i++){
    bool d(distances[i] < epsilon_new);
    num += d;
  }
  double sum = 0;
  for(int i = 0; i < N; i++){
    weights[i] = weights_old[i] * (num / denom);
    sum += weights[i];
  }

  for(int i = 0; i < N; i++){
    weights[i] /= sum;
  }
  return ESS(weights) - alpha*ESS(weights_old);
}

double bisection( double (*func) (double, double, double, std::vector<double>,std::vector<double>),double* endpoints, double epsilon_old, double alpha, std::vector<double> weights_old,std::vector<double> distances, double tol, int NMAX){
  // Assigning endpoint
  double a = endpoints[0], b = endpoints[1];
  // Assert functions satisfy bisection algorithm conditions
  assert(a < b);
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



double quantile(std::vector<double> distances, double p){

  int N = distances.size();
  std::set<double> organisedDistances;
  for(int i = 0; i < N; i++){
    organisedDistances.insert(distances[i]);
  }

  int n = (N * (p / 100)) + 1;
  std::set<double>::iterator it = organisedDistances.begin();
  std::advance(it, n);


  return *it;

}

void ABCSMC(std::vector<std::vector<double>>& theta_out,                         // Outputted sampled parameters
            std::vector<double>& weights_out,                                   // corresponding weights
            std::map<double,std::vector<double>> X,                         // True data
            CRN srn,                                                            // Chemical Reaction Network of choice
            std::vector<double> epsilon,                                        // Discrepancy threshold
            int N,                                                              // Number of particles
            int B_t,
            double T,
            randomVariable& rv){                                                             // Number of stochastic simulations
            /*
                ABC Algorithm
              */
              std::cout << "\n*************************************************\n";
              std::cout << "************** ABC-SMC Algorithm ****************\n";
              std::cout << "*************************************************\n";
              assert(N == (int) theta_out.size() && N == (int) weights_out.size());

              //int N_T = N/2;
              //double epsilon_0 = 0.01, delta_0 = 0.05, alpha_0 = delta_0/50;    // Parameters for Massart Algorithm

              int numThresholds = epsilon.size();
              std::map<double,std::vector<int>> Y;                                  // Trace is reordered to be observed every 0.5 seconds of simulation to make it easier to analyse
              std::map<double,std::vector<int>> Yfull;                                  // Trace is reordered to be observed every 0.5 seconds of simulation to make it easier to analyse


              double timestep = 0;
              // Calculate timestep needed from data
              timestep = std::next(X.begin(),1)->first - X.begin()->first;

              std::vector<int> initial = srn.retInitial();
              std::vector<int> b_t(N,0);                                  // Number of simulations that d(X,Y) < epsilon

              // Storing the sampled parameters in theta with corresponding weights
              std::vector<std::vector<double>> prev_theta(N,std::vector<double>(srn.prior_pt(rv).size(),0));
              std::vector<std::vector<double>> theta(N,std::vector<double>(srn.prior_pt(rv).size(),0));
              std::vector<std::vector<double>> weights( int (numThresholds),std::vector<double> (N,0));

              std::vector<double> distances(N,0), empty(N,0);

              std::vector<double> sample(srn.prior_pt(rv).size(),0);

              // Storing the sampled parameters that were rejected with corresponding SMC values
              std::map<std::vector<double>,double> paramsSMC;


              // Storing weighted mean and covariances:
              std::vector<std::vector<double>> covariance(sample.size(),std::vector<double>(sample.size(),0));

              for(int t = 0; t < numThresholds; t++)
              { // Setting population indicator

                if( t != 0){
                  epsilon[t] = quantile(distances, 0.5);

                  calculateCovariance(covariance, prev_theta, weights[t-1]);
                  std::cout << "Old \u03B5 = " << epsilon[t-1] << std::endl;
                  std::cout << "New \u03B5 = " << epsilon[t] << std::endl;

                  double endpoints[2]; endpoints[0] = 0; endpoints[1] = epsilon[t]*10;
                  double ey = bisection( g,endpoints, epsilon[t-1], 0.5, weights[t-1],distances, 1e-5, 10000);
                  std::cout << "possible new epsilon  = " << bisection( g,endpoints, epsilon[t-1], 0.5, weights[t-1],distances, 1e-5, 10000) << std::endl;
                  // Calculate weights here and return
                  std::cout << "g = " <<  g(ey, epsilon[t-1], 0.5, weights[t-1],distances) << std::endl;

                  // Reset distances
                  distances = empty;


                }

                for(int i = 0; i < N; i++){// Setting particle indicator
                // Clearing set from last particle
                std::cout << "i = " << i << "\r";
                b_t[i] = 0;

                while(b_t[i] == 0){
                  b_t[i] = 0;  // Reset for new set of simulations

                  if(t == 0){ // Sample parameters
                    // Sampling from prior
                    sample = srn.prior_pt(rv);
                  }
                  else{
                    // Sampling from previously accepted particles
                    sample = randomParameter(weights[t-1],prev_theta,rv);
                    sample = K_t(sample, covariance, srn.prior_pdf_pt,rv);
                  }
                /****************************************************************************
                  GENERATING Synthetic DATA
                *****************************************************************************/
                // int num = 0, denom = B_t;
                //double gamma_hat = 0;
                for(int b = 0; b < B_t; b++){
                  // Storing simulated data in map Y

                  Yfull = gillespie(srn,sample,T,rv);
                  Y = organiseTrace(Yfull,timestep);
                  // num += check_trace(Yfull);
                  if(d(X,Y) < epsilon[t]){
                    b_t[i]++;
                    // Save accepted particle i

                    }
                    distances[i] += (double) d(X,Y) / B_t;

                  }
                  //std::vector<double> prob(3,0);

                  //prob = AbsoluteErrorMassart(epsilon_0,delta_0,alpha_0,srn,sample,T,rv);
                  //gamma_hat = prob[0];
                  //paramsSMC.insert(std::pair<std::vector<double>,double>(sample, gamma_hat));
                }

                // Save accepted particle i
                //distances[i] = d(X,Y);
                theta[i] = sample;
                if(t == 0){
                  // calculating weight for particle i when t == 0
                  for(int i = 0; i < N; i++){
                    weights[t][i] = b_t[i];

                  }
                  // Calculating initial covariance
                  normaliseWeights(weights[t]);
                  calculateCovariance(covariance, theta, weights[t]);
                }
                else{
                  // Calculating weight for particle i when t != 0
                  for(int i = 0; i < N; i++){
                    double sum = 0;
                    // Calculating sum present in denominator
                    for(int j = 0; j < N; j++){
                      sum += weights[t-1][j] * K_t_pdf(theta[i],prev_theta[j],covariance ,rv);
                      }

                    weights[t][i] = srn.prior_pdf_pt(theta[i],rv) * b_t[i] / sum;
                  }
                }
                // Print out all sampled particles and corresponding probabilities
                std::ofstream output;
                output.open("simulations/SIR/smc/threshold" + std::to_string(t) + ".csv");

                output << "epsilon = " << epsilon[t] << std::endl;
                for(int kk = 0; kk < (int) srn.retParameterNames().size();kk++){
                  output << srn.retParameterNames()[kk] << ",";
                }
                output << "prob,\n";

              std::map<std::vector<double>,double>::iterator itr;

                for(itr = paramsSMC.begin(); itr != paramsSMC.end(); itr++)
                {
                  for(int ii = 0; ii < (int) itr->first.size(); ii++)
                  {
                    output << itr->first[ii] << ",";
                  }
                  output << itr->second << ",\n";
                }

                fflush(stdout);
                // end print

              } // end of particle loop
                paramsSMC.clear();
                normaliseWeights(weights[t]);

                prev_theta = theta;

                std::vector<double> tempMean(theta[0].size(),0);
                std::cout << "current mean \u03B8 = (";
                for(int p = 0; p < (int)theta[0].size(); p++){
                  for(int n = 0; n < N; n++){
                    tempMean[p] += prev_theta[n][p]*weights[t][n];
                  }
                   std::cout << tempMean[p] << ", ";
                }
                std::cout << ")\n";
                std::cout << "ESS = " << ESS(weights[t]) << std::endl;
              } // end of threshold loop


              double sumWeights = 0;
              std::vector<double> meanTheta(theta[0].size(),0);
              for(int i = 0; i < N; i++){
                sumWeights += weights[numThresholds-1][i];
              }
              for(int j = 0; j < (int) theta[0].size(); j++){
                for(int i = 0; i < N; i++){
                  meanTheta[j] += weights[numThresholds-1][i]*theta[i][j] / sumWeights;
                }
              }

              for(int i = 0; i < N; i++)
              {
                std::cout << "\u03B8[" << i << "] = " << theta[i][0] << ", " << theta[i][1] << "\n";
              }

              std::cout << "Mean \u03B8:\n";
              for(int j = 0; j < (int) theta[0].size(); j++){
                std::cout << meanTheta[j] << ", ";
              }
              std::cout << std::endl;

              // Print out all sampled particles and corresponding probabilities
              std::ofstream outputFinal;
              outputFinal.open("simulations/SIR/outputABCSMC.csv");

              for(int kk = 0; kk < (int) srn.retParameterNames().size();kk++){
                outputFinal << srn.retParameterNames()[kk] << ",";
              }
              outputFinal << "weights,\n";

              for(int i = 0; i < N; i++){
                for(int ii = 0; ii < (int) theta[0].size(); ii++){
                  outputFinal << theta[i][ii] << ",";
                }
                outputFinal << weights[numThresholds-1][i] << ",\n";
              }

            std::map<std::vector<double>,double>::iterator itr;

              for(itr = paramsSMC.begin(); itr != paramsSMC.end(); itr++)
              {
                for(int ii = 0; ii < (int) itr->first.size(); ii++)
                {
                  outputFinal << itr->first[ii] << ",";
                }
                outputFinal << itr->second << ",\n";
              }

              // Outputting sampled parameters and corresponding weights
              theta_out = theta;
              weights_out = weights[numThresholds-1];
};

void ABCSMCSMC(std::vector<std::vector<double>>& theta_out,                         // Outputted sampled parameters
            std::vector<double>& weights_out,                                   // corresponding weights
            std::map<double,std::vector<double>> X,                         // True data
            CRN srn,                                                            // Chemical Reaction Network of choice
            std::vector<double> epsilon,                                        // Discrepancy threshold
            int N,                                                              // Number of particles
            double T,
            randomVariable& rv){                                                             // Number of stochastic simulations

              std::cout << "\n*************************************************\n";
              std::cout << "************** ABC-SMC-SMC Algorithm ****************\n";
              std::cout << "*************************************************\n";
              assert(N == (int) theta_out.size() && N == (int) weights_out.size());

              double epsilon_0 = 0.01, delta_0 = 0.05, alpha_0 = delta_0/50;    // Parameters for Massart Algorithm

              int numThresholds = epsilon.size();
              std::map<double,std::vector<int>> Y;                                  // Trace is reordered to be observed every 0.5 seconds of simulation to make it easier to analyse
              std::map<double,std::vector<int>> Yfull;                                  // Trace is reordered to be observed every 0.5 seconds of simulation to make it easier to analyse


              // Calculate timestep needed from data
              double timestep = 0;
              timestep = std::next(X.begin(),1)->first - X.begin()->first;

              std::vector<int> initial = srn.retInitial();
              std::vector<double> b_t(N,0);                                  // Number of simulations that d(X,Y) < epsilon

              // Storing the sampled parameters in theta with corresponding weights
              std::vector<std::vector<double>> prev_theta(N,std::vector<double>(srn.prior_pt(rv).size(),0));
              std::vector<std::vector<double>> theta(N,std::vector<double>(srn.prior_pt(rv).size(),0));
              std::vector<std::vector<double>> weights( int (numThresholds),std::vector<double> (N,0));

              std::vector<double> distances(N,0), empty(N,0);

              std::vector<double> sample(srn.prior_pt(rv).size(),0);

              // Storing the sampled parameters that were rejected with corresponding SMC values
              std::map<std::vector<double>,std::vector<double>> paramsSMC;


              // Storing covariances:
              std::vector<std::vector<double>> covariance(sample.size(),std::vector<double>(sample.size(),0));

              for(int t = 0; t < numThresholds; t++){ // Setting population indicator

                if( t != 0)
                {
                  epsilon[t] = quantile(distances, 0.5);
                  calculateCovariance(covariance, prev_theta, weights[t-1]);
                  std::cout << "Old \u03B5 = " << epsilon[t-1] << std::endl;
                  std::cout << "New \u03B5 = " << epsilon[t] << std::endl;

                  double endpoints[2]; endpoints[0] = 0; endpoints[1] = epsilon[t]*10;
                  double neweps = bisection(g, endpoints, epsilon[t-1], 0.5, weights[t-1], distances,1e-5, 10000);
                  std::cout << "possible new epsilon = " << neweps << std::endl;

                  // Reset distances
                  distances=empty;
                }

                // Print out all sampled particles and corresponding probabilities
                std::ofstream output;
                output.open("simulations/SIR/smc/threshold" + std::to_string(t) + ".csv");
                output << "epsilon = " << epsilon[t] << std::endl;
                for(int kk = 0; kk < (int) srn.retParameterNames().size();kk++){
                  output << srn.retParameterNames()[kk] << ",";
                }
                output << "prob,lower bound, upper bound," << std::endl;

                // Setting particle indicator
                for(int i = 0; i < N; i++)
                {
                b_t[i] = 0;
                distances[i] = 0;
                while(b_t[i] == 0)
                {
                  b_t[i] = 0;   // Reset for new set of simulations
                  distances[i] = 0;
                  if(t == 0)    // Sample parameters
                  {
                    // Sampling from prior
                    sample = srn.prior_pt(rv);
                  }
                  else
                  { // Sampling from previously accepted particles
                    sample = randomParameter(weights[t-1],prev_theta,rv);
                    sample = K_t(sample, covariance, srn.prior_pdf_pt,rv);
                  }
                /****************************************************************************
                  GENERATING Synthetic DATA
                *****************************************************************************/
                double gamma_hat = 0,a_k = 0, b_k = 0;
                std::cout << "t = " << t << ", i = " << i << std::endl;
                std::vector<double> Mass(5,0);


                //Something funky is going on with the following functions
                Mass = AbsoluteErrorMassart(epsilon_0, delta_0,alpha_0,srn,sample,T,rv,epsilon[t],X);
                // Storing outputs of Modified Massart Algorithm
                gamma_hat = Mass[0]; b_t[i] = Mass[1]; distances[i] = Mass[2]; a_k = Mass[3]; b_k = Mass[4];

                // This seems to work (inferred mean is near true mean) so the number of simulations for each particle doesn't affect the mean.
                // Instead I've tried writing b_t as a ratio of acceptance rates and will run that to see if that will affect the inference technique
                /*
                int B_t = ceil( rv.uniform(0,1800));
                for(int b = 0; b < B_t; b++){
                  // Storing simulated data in map Y
                  Yfull = gillespie(srn,sample,T,rv);
                  Y = organiseTrace(Yfull,timestep);
                  // num += check_trace(Yfull);
                  if(d(X,Y) < epsilon[t]){
                    b_t[i]++;
                    // Save accepted particle i
                    }
                    distances[i] += (double) d(X,Y) / B_t;
                  }
                  b_t[i] /= B_t;
                  */



                std::cout << "distances[" << i << "] = " << distances[i] << std::endl;
                std::vector<double> MassOutput(3,0); MassOutput[0] = gamma_hat; MassOutput[1] = a_k; MassOutput[2] = b_k;

                paramsSMC.insert(std::pair<std::vector<double>,std::vector<double>>(sample, MassOutput));

                // Save sampled parameters and corresponding estimated SMC and CIs
                for(int m = 0; m < (int) sample.size(); m++)
                {
                  output << sample[m] << ",";
                }
                output << MassOutput[0] << "," << MassOutput[1] << "," << MassOutput[2] << "," << std::endl;
                }

                // Save accepted particle i
                theta[i] = sample;
                if(t == 0)
                {
                  // calculating weight for particle i when t == 0
                  for(int i = 0; i < N; i++)
                  {
                    weights[t][i] = b_t[i];
                  }
                }
                else
                {
                  // Calculating weight for particle i when t != 0
                  for(int i = 0; i < N; i++)
                  {
                    double sum = 0;
                    for(int j = 0; j < N; j++)
                    {
                      sum += weights[t-1][j] * K_t_pdf(theta[i],prev_theta[j],covariance ,rv);
                    }
                    weights[t][i] = (srn.prior_pdf_pt(theta[i],rv) * b_t[i]) / sum;
                  }
                }

                } // end of particle loop
                paramsSMC.clear();
                normaliseWeights(weights[t]);

                prev_theta = theta;

                std::ofstream outputinter;
                outputinter.open("simulations/SIR/weights/threshold" + std::to_string(t) + ".csv");
                for(int kk = 0; kk < (int) srn.retParameterNames().size();kk++){
                  outputinter << srn.retParameterNames()[kk] << ",";
                }
                outputinter << std::endl;

                for(int j = 0; j < N; j++){
                  for(int kk = 0; kk < (int) srn.retParameterNames().size(); kk++){
                    outputinter << theta[j][kk] << ",";
                  }
                  outputinter << weights[t][j] << std::endl;
                }

                std::vector<double> tempMean(theta[0].size(),0);
                std::cout << "current mean \u03B8 = (";
                for(int p = 0; p < (int)theta[0].size(); p++){
                  for(int n = 0; n < N; n++){
                    tempMean[p] += prev_theta[n][p]*weights[t][n];
                  }
                   std::cout << tempMean[p] << ", ";
                }
                std::cout << ")\n";
                std::cout << "ESS = " << ESS(weights[t]) << std::endl;
              } // end of threshold loop


              double sumWeights = 0;
              std::vector<double> meanTheta(theta[0].size(),0);
              for(int i = 0; i < N; i++){
                sumWeights += weights[numThresholds-1][i];
              }
              for(int j = 0; j < (int) theta[0].size(); j++){
                for(int i = 0; i < N; i++){
                  meanTheta[j] += weights[numThresholds-1][i]*theta[i][j] / sumWeights;
                }
              }

              for(int i = 0; i < N; i++)
              {
                std::cout << "\u03B8[" << i << "] = " << theta[i][0] << ", " << theta[i][1] << "\n";
              }

              std::cout << "Mean \u03B8:\n";
              for(int j = 0; j < (int) theta[0].size(); j++){
                std::cout << meanTheta[j] << ", ";
              }
              std::cout << std::endl;

              // Print out all sampled particles and corresponding probabilities
              std::ofstream outputFinal;
              outputFinal.open("simulations/SIR/outputABCSMC.csv");

              for(int kk = 0; kk < (int) srn.retParameterNames().size();kk++){
                outputFinal << srn.retParameterNames()[kk] << ",";
              }
              outputFinal << "weights,\n";

              for(int i = 0; i < N; i++){
                for(int ii = 0; ii < (int) theta[0].size(); ii++){
                  outputFinal << theta[i][ii] << ",";
                }
                outputFinal << weights[numThresholds-1][i] << ",\n";
              }

            std::map<std::vector<double>,std::vector<double>>::iterator itr;

              for(itr = paramsSMC.begin(); itr != paramsSMC.end(); itr++)
              {
                for(int ii = 0; ii < (int) itr->first.size(); ii++)
                {
                  outputFinal << itr->first[ii] << ",";
                }
                outputFinal << itr->second[0] << "," << itr->second[1] << "," << itr->second[2] << ",\n";
              }

              // Outputting sampled parameters and corresponding weights
              theta_out = theta;
              weights_out = weights[numThresholds-1];
};
