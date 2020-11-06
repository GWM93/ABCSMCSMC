#include "SIR.h"

/* SIR Epidemiological model model
  S + I -{k_i}-> 2I
  I -{k_r}-> R
*/

// Prior distribution
std::vector<double> prior_SIR(randomVariable& rv){
  std::vector<double> theta(2,0);
  double k_i, k_r;
  k_i = rv.uniform(5e-5,0.003);
  k_r = rv.uniform(0.005,0.2);

  theta[0] = k_i;
  theta[1] = k_r;

  return theta;

}
// Prior pdf
double prior_SIR_pdf(std::vector<double> theta, randomVariable& rv){

double pdf;
pdf = rv.uniform_pdf(theta[0], 5e-5,0.003) * rv.uniform_pdf(theta[1], 0.005,0.2);
return pdf;
}

std::vector<double> propensity_SIR(std::vector<double> theta, std::vector<int> state){

std::vector<double> r(2,0);
double k_i = theta[0], k_r = theta[1];
int S = state[0], I = state[1];

r[0] = k_i * S * I;
r[1] = k_r * I;

return r;
}

void SIR_ABCSMCSMC(randomVariable& rv){

  std::vector<double> sample = prior_SIR(rv);                                     // Initialise size of parameter sample
  std::vector<int> initial = {95,5,0};                                       // Initial molecule count
  std::vector<std::vector<int>> S;                                                // Stoichiometric Matrix
  double T = 150;                                                                  // Runtime of SSA

  std::vector<std::string> speciesNames;
  speciesNames.push_back("S"); speciesNames.push_back("I"); speciesNames.push_back("R");
  std::vector<std::string> paramNames;
  paramNames.push_back("k_i"); paramNames.push_back("k_r");
  /*  SIR model
      R1: S + I -{k_i}-> 2I
      R2: I -{k_r}-> R
  */
    //  R1 R2
  S = {{-1, 0}, // S
       {1, -1}, // I
       {0, 1}};// R

  // True parameter values
  sample[0] = 0.002; sample[1] = 0.125;
  // Defining Chemical Reaction Network using CRN class
  CRN SIR(S);
  SIR.setInitial(initial);
  SIR.inputtingNames(speciesNames);
  SIR.inputtingParameterNames(paramNames);
  SIR.propensity_pt = propensity_SIR;
  SIR.prior_pt = prior_SIR;
  SIR.prior_pdf_pt = prior_SIR_pdf;

  std::string directory;
  directory = "simulations/SIR/";
  /****************************************************************************
      GENERATING TRUE DATA
  *****************************************************************************/
  // Storing true data in map X
  std::map<double, std::vector<double>> X;
  X = trueData(SIR, sample, T, rv, 15.0, directory, 5);
  std::cout << "Property true on data trace?: " << check_trace(X) << std::endl;
  printTrace(X, speciesNames);

  std::map<double,std::vector<double>> Y;

  /*// Discrepancy threshold
  std::vector<double> epsilon(6,0);
  epsilon = {1500,500,300,200,199,75};
  */
  int thresholds = 20;
  std::vector<double> epsilon(thresholds,0);
  epsilon[0] = std::numeric_limits<double>::infinity(); // Set initial threshold
  // Number of stochastic simulations
  // int B_t = 100;
  // Number of particles
  int N = 500;
  // Store sampled parameter and corresponding weights in following vectors
  std::vector<std::vector<double>> theta_out(N, std::vector<double> (4,0));
  std::vector<double> weights_out(N,0);

  //ABCSMC(theta_out, weights_out, X, SIR, epsilon, N,5, T, rv);
  ABCSMCSMC(theta_out, weights_out, X, SIR, epsilon, N, T, rv);
  // Ensure outputs are of equal length
  assert(theta_out.size() == weights_out.size() && (int) theta_out.size() == N);

  // Writing to file
  std::ofstream output;
  output.open("simulations/SIR/output.txt");
  if(output.is_open())
  {
    for(int j = 0; j < (int) theta_out[0].size(); j++)
    {
      std::cout << "theta_" << j << ", ";
      output << "theta_" << j << ", ";
    }
    std::cout << "weight" << std::endl;
    output << "weight" <<  std::endl;

    for(int i = 0; i < (int) theta_out.size(); i++)
    {
      for(int j = 0; j < (int) theta_out[0].size(); j++)
      {
        std::cout << theta_out[i][j] << ", ";
        output << theta_out[i][j] << ", ";
      }
      output << weights_out[i] << std::endl;
      std::cout << weights_out[i] << std::endl;
    }
    // Closing file
    output.close();
  }
  else
  {
    std::cout <<"FILE NOT OPENED\n";
  }

  std::cout << "Mean = " << std::endl;
  std::vector<double> mean(theta_out[0].size(), 0);
  for(int i = 0; i < (int) theta_out[0].size(); i++)
  {
  for(int j = 0; j < N; j++)
    {
      mean[i] += theta_out[j][i]*weights_out[j];
    }
    std::cout << mean[i] << ", ";
  }
  std::cout << std::endl;
}

void SIR_paramSynth(randomVariable& rv)
{
  std::vector<double> sample = prior_SIR(rv);                                     // Initialise size of parameter sample
  std::vector<int> initial = {95,5,0};                                            // Initial molecule count
  std::vector<std::vector<int>> S;                                                // Stoichiometric Matrix
  double T = 150.0;                                                               // Runtime of SSA

  std::vector<std::string> speciesNames;
  speciesNames.push_back("S"); speciesNames.push_back("I"); speciesNames.push_back("R");
  std::vector<std::string> paramNames;
  paramNames.push_back("k_i"); paramNames.push_back("k_r");
  /*  SIR model
      R1: S + I -{k_i}-> 2I
      R2: I -{k_r}-> R
  */
    //  R1 R2
  S = {{-1, 0}, // S
       {1, -1}, // I
       {0, 1}};// R

  // True parameter values
  sample[0] = 0.0025; sample[1] = 0.11;

  // Defining Chemical Reaction Network using CRN class
  CRN SIR(S);
  SIR.setInitial(initial);
  SIR.inputtingNames(speciesNames);
  SIR.inputtingParameterNames(paramNames);
  SIR.propensity_pt = propensity_SIR;
  SIR.prior_pt = prior_SIR;
  SIR.prior_pdf_pt = prior_SIR_pdf;

  // Massart algorithm parameters
  /* Massart parameters
     double epsilon = 0.01, delta = 0.05, alpha = delta/50.0;
    double gamma = 0.0;
  */
  double epsilon = 0.01, delta = 0.05, alpha = delta/50.0;

  std::vector<double> prob(3,0);
  double prob_hat, prob_lb, prob_ub;

  // Writing to file
  std::ofstream output;
  output.open("simulations/SIR/paramSynth.txt");
  if(output.is_open())
  {
    // Print parameters onto top of text file
    for(int j = 0; j < (int) sample.size(); j++)
    {
      std::cout << "theta_" << j << ", ";
      output << "theta_" << j << ", ";
    }
    // print corresponding estimated probabilities
    output << "prob_hat,prob_lb, prob_ub," << std::endl;

  for(int i = 0; i < 1000; i++)
  {
    std::cout << "i = " << i << std::endl;
    sample = prior_SIR(rv);
    prob = AbsoluteErrorMassart(epsilon,delta,alpha,SIR,sample,T,rv );
    for(int j = 0; j < (int) sample.size(); j++)
    {
      output << sample[j] << ", ";
    }
    prob_hat = prob[0]; prob_lb = prob[1]; prob_ub = prob[2];
    std::cout << "(" << prob_hat << ", " << prob_lb << ", " << prob_ub << ")" << std::endl;
    output << prob_hat << ", " << prob_lb << ", " << prob_ub << std::endl;
  }

  // Closing file
  output.close();
  }
}

void SIR_paramSynth_BayesianSMC(randomVariable& rv)
{
  std::vector<double> sample = prior_SIR(rv);                                     // Initialise size of parameter sample
  std::vector<int> initial = {95,5,0};                                            // Initial molecule count
  std::vector<std::vector<int>> S;                                                // Stoichiometric Matrix
  double T = 150.0;                                                               // Runtime of SSA

  std::vector<std::string> speciesNames;
  speciesNames.push_back("S"); speciesNames.push_back("I"); speciesNames.push_back("R");
  std::vector<std::string> paramNames;
  paramNames.push_back("k_i"); paramNames.push_back("k_r");
  /*  SIR model
      R1: S + I -{k_i}-> 2I
      R2: I -{k_r}-> R
  */
    //  R1 R2
  S = {{-1, 0}, // S
       {1, -1}, // I
       {0, 1}};// R

  // True parameter values
  sample[0] = 0.0025; sample[1] = 0.11;

  // Defining Chemical Reaction Network using CRN class
  CRN SIR(S);
  SIR.setInitial(initial);
  SIR.inputtingNames(speciesNames);
  SIR.inputtingParameterNames(paramNames);
  SIR.propensity_pt = propensity_SIR;
  SIR.prior_pt = prior_SIR;
  SIR.prior_pdf_pt = prior_SIR_pdf;

  // Massart algorithm parameters
  /* Massart parameters
     double epsilon = 0.01, delta = 0.05, alpha = delta/50.0;
    double gamma = 0.0;
  */
    double delta = 0.01, c = 0.99;

  std::vector<double> prob(3,0);
  double prob_hat, prob_lb, prob_ub;

  // Writing to file
  std::ofstream output;
  output.open("simulations/SIR/paramSynth_BayesianSMC.txt");
  if(output.is_open())
  {
    // Print parameters onto top of text file
    for(int j = 0; j < (int) sample.size(); j++)
    {
      std::cout << "theta_" << j << ", ";
      output << "theta_" << j << ", ";
    }
    // print corresponding estimated probabilities
    output << "prob_hat,prob_lb, prob_ub," << std::endl;

  for(int i = 0; i < 1000; i++)
  {
    std::cout << "i = " << i << std::endl;
    sample = prior_SIR(rv);
    prob = BayesianSMC(delta, c, SIR, sample,T,rv);
    for(int j = 0; j < (int) sample.size(); j++)
    {
      output << sample[j] << ", ";
    }
    prob_hat = prob[0]; prob_lb = prob[1]; prob_ub = prob[2];
    std::cout << "(" << prob_hat << ", " << prob_lb << ", " << prob_ub << ")" << std::endl;
    output << prob_hat << ", " << prob_lb << ", " << prob_ub << std::endl;
  }

  // Closing file
  output.close();
  }
}


void SIR_distance(randomVariable& rv)
{

    std::vector<double> sample = prior_SIR(rv);                                     // Initialise size of parameter sample
    std::vector<int> initial = {95,5,0};                                       // Initial molecule count
    std::vector<std::vector<int>> S;                                                // Stoichiometric Matrix
    double T = 150;                                                                  // Runtime of SSA

    std::vector<std::string> speciesNames;
    speciesNames.push_back("S"); speciesNames.push_back("I"); speciesNames.push_back("R");
    std::vector<std::string> paramNames;
    paramNames.push_back("k_i"); paramNames.push_back("k_r");
    /*  SIR model
        R1: S + I -{k_i}-> 2I
        R2: I -{k_r}-> R
    */
      //  R1 R2
    S = {{-1, 0}, // S
         {1, -1}, // I
         {0, 1}};// R

    // True parameter values
    sample[0] = 0.0025; sample[1] = 0.11;
    // Defining Chemical Reaction Network using CRN class
    CRN SIR(S);
    SIR.setInitial(initial);
    SIR.inputtingNames(speciesNames);
    SIR.inputtingParameterNames(paramNames);
    SIR.propensity_pt = propensity_SIR;
    SIR.prior_pt = prior_SIR;
    SIR.prior_pdf_pt = prior_SIR_pdf;

    std::string directory;
    directory = "simulations/SIR/";
    /****************************************************************************
        GENERATING TRUE DATA
    *****************************************************************************/
    // Storing true data in map X
    std::map<double, std::vector<double>> X;
    double timestep = 15.0;
    int numTrueSamples = 5;
    X = trueData(SIR, sample, T, rv, timestep, directory, numTrueSamples);
    printTrace(X, speciesNames);
    std::map<double,std::vector<int>> Y,Yfull;
    int n_samp = 15000, B_t = 100;
    std::vector<double> distances(n_samp,0);

  // Writing to file
  std::ofstream output;
  output.open("simulations/SIR/distances.txt");
  if(output.is_open())
  {
    // Print parameters onto top of text file
    for(int j = 0; j < (int) sample.size(); j++)
    {
      std::cout << "theta_" << j << ",";
      output << "theta_" << j << ",";
    }
    // print corresponding estimated probabilities
    output << "distances," << std::endl;

  for(int i = 0; i < n_samp; i++)
  {
    std::cout << "i = " << i;
    sample = prior_SIR(rv);
    for(int b = 0; b < B_t; b++){
      // Storing simulated data in map Y
      Yfull = gillespie(SIR,sample,T,rv);
      Y = organiseTrace(Yfull,timestep);
      distances[i] += (double) d(X,Y) / B_t;
    }
    std::cout << ", theta[" << i << "] = (";
    for(int j = 0; j < (int) sample.size(); j++)
    {
      std::cout << sample[j] << ", ";
      output << sample[j] << ", ";
    }
      std::cout << "), distances[" << i << "] = " << distances[i] << std::endl;
    output << distances[i] << ", " <<  std::endl;
  }

  // Closing file
  output.close();
  }
}

void SIR_Massart(randomVariable& rv)
{
  std::vector<double> sample = prior_SIR(rv);                                     // Initialise size of parameter sample
  std::vector<int> initial = {95,5,0};                                            // Initial molecule count
  std::vector<std::vector<int>> S;                                                // Stoichiometric Matrix
  double T = 150.0;                                                               // Runtime of SSA

  std::vector<std::string> speciesNames;
  speciesNames.push_back("S"); speciesNames.push_back("I"); speciesNames.push_back("R");
  std::vector<std::string> paramNames;
  paramNames.push_back("k_i"); paramNames.push_back("k_r");
  /*  SIR model
      R1: S + I -{k_i}-> 2I
      R2: I -{k_r}-> R
  */
    //  R1 R2
  S = {{-1, 0}, // S
       {1, -1}, // I
       {0, 1}};// R

  // True parameter values
  sample[0] = 0.0025; sample[1] = 0.11;

  // Defining Chemical Reaction Network using CRN class
  CRN SIR(S);
  SIR.setInitial(initial);
  SIR.inputtingNames(speciesNames);
  SIR.inputtingParameterNames(paramNames);
  SIR.propensity_pt = propensity_SIR;
  SIR.prior_pt = prior_SIR;
  SIR.prior_pdf_pt = prior_SIR_pdf;

  // Massart algorithm parameters
  /* Massart parameters
     double epsilon = 0.01, delta = 0.05, alpha = delta/50.0;
    double gamma = 0.0;
  */
  double epsilon = 0.01, delta = 0.05, alpha = delta/50.0;

  std::vector<double> prob(3,0);

  std::vector<double> theta_phi(2,0),theta_negPhi(2,0),theta_U(2,0);
  theta_phi = {0.002, 0.075};
  theta_negPhi = {0.001,0.15};
  theta_U = {0.002,0.125};
  std::cout << "theta_phi:\n";
  prob = AbsoluteErrorMassart(epsilon,delta,alpha,SIR,theta_phi,T,rv );
  std::cout << "\\hat{lambda}_{\\phi}(\\theta) = " << prob[0] << ", [L,U] = [" << prob[1] << ", " << prob[2] << "]\n";
  std::cout << "delta = " << prob[2] - prob[1] << std::endl;
  std::cout << "\n\n";

  std::cout << "theta_negPhi:\n";
  prob = AbsoluteErrorMassart(epsilon,delta,alpha,SIR,theta_negPhi,T,rv );
  std::cout << "\\hat{lambda}_{\\phi}(\\theta) = " << prob[0] << ", [L,U] = [" << prob[1] << ", " << prob[2] << "]\n";
std::cout << "delta = " << prob[2] - prob[1] << std::endl;
  std::cout << "\n\n";

  std::cout << "theta_U:\n";
  prob = AbsoluteErrorMassart(epsilon,delta,alpha,SIR,theta_U,T,rv );
  std::cout << "\\hat{lambda}_{\\phi}(\\theta) = " << prob[0] << ", [L,U] = [" << prob[1] << ", " << prob[2] << "]\n";
  std::cout << "delta = " << prob[2] - prob[1] << std::endl;
  std::cout << "\n\n";

}

void SIR_BayesSMC(randomVariable& rv)
{
  std::vector<double> sample = prior_SIR(rv);                                     // Initialise size of parameter sample
  std::vector<int> initial = {95,5,0};                                            // Initial molecule count
  std::vector<std::vector<int>> S;                                                // Stoichiometric Matrix
  double T = 150.0;                                                               // Runtime of SSA

  std::vector<std::string> speciesNames;
  speciesNames.push_back("S"); speciesNames.push_back("I"); speciesNames.push_back("R");
  std::vector<std::string> paramNames;
  paramNames.push_back("k_i"); paramNames.push_back("k_r");
  /*  SIR model
      R1: S + I -{k_i}-> 2I
      R2: I -{k_r}-> R
  */
    //  R1 R2
  S = {{-1, 0}, // S
       {1, -1}, // I
       {0, 1}};// R

  // True parameter values
  sample[0] = 0.0025; sample[1] = 0.11;

  // Defining Chemical Reaction Network using CRN class
  CRN SIR(S);
  SIR.setInitial(initial);
  SIR.inputtingNames(speciesNames);
  SIR.inputtingParameterNames(paramNames);
  SIR.propensity_pt = propensity_SIR;
  SIR.prior_pt = prior_SIR;
  SIR.prior_pdf_pt = prior_SIR_pdf;

  // Massart algorithm parameters
  /* Massart parameters
     double epsilon = 0.01, delta = 0.05, alpha = delta/50.0;
    double gamma = 0.0;
  */
  double delta = 0.005, c = 0.95;
  std::vector<double> prob(3,0);

  std::vector<double> theta_phi(2,0),theta_negPhi(2,0),theta_U(2,0);
  theta_phi = {0.002, 0.075};
  theta_negPhi = {0.001,0.15};
  theta_U = {0.002,0.125};
  std::cout << "theta_phi:\n";
  prob = BayesianSMC(delta, c, SIR, theta_phi,T,rv);
  std::cout << "\\hat{lambda}_{\\phi}(\\theta) = " << prob[0] << ", [L,U] = [" << prob[1] << ", " << prob[2] << "]\n";
std::cout << "delta = " << prob[2] - prob[1] << std::endl;
  std::cout << "\n\n";

  std::cout << "theta_negPhi:\n";
  prob = BayesianSMC(delta, c, SIR, theta_negPhi,T,rv);
  std::cout << "\\hat{lambda}_{\\phi}(\\theta) = " << prob[0] << ", [L,U] = [" << prob[1] << ", " << prob[2] << "]\n";
  std::cout << "delta = " << prob[2] - prob[1] << std::endl;
  std::cout << "\n\n";

  std::cout << "theta_U:\n";
  prob = BayesianSMC(delta, c, SIR, theta_U,T,rv);
  std::cout << "\\hat{lambda}_{\\phi}(\\theta) = " << prob[0] << ", [L,U] = [" << prob[1] << ", " << prob[2] << "]\n";
  std::cout << "delta = " << prob[2] - prob[1] << std::endl;
  std::cout << "\n\n";

}
