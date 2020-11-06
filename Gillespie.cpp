#include "Gillespie.h"

// Return every reaction in the SSA
std::map<double,std::vector<int>> gillespie(
  CRN srn,
  std::vector<double> theta,                                                    // Parameters for propensity function
  double T,                                                                     // Running time of algorithm, T
  randomVariable& rv)                                                           // Pass initialization of randomVariable class to allow for random variable distributions
{
    // Initialize map that tracks the times and states of the chemical reaction network during the SSA
  std::map<double,std::vector<int>> trace;

  double r1, r2, tau; // Parameters for Gillespie algorithm
  double time = 0;        // Variable tracking time
  double sumOfPropensities; // variable for the sum of propensities
  double cumulativeSum = 0; // Cumulative sum that's required to calculate which reaction occurs

  // Initial state of CTMC
  std::vector<int> initialState = srn.retInitial(), currentState = initialState;                                                // Initial state of CTMC

  // Defining stoichiometricMatrix
  std::vector<std::vector<int>> stoichiometricMatrix;
  stoichiometricMatrix = srn.retStoichiometric();
  // Defining vector prop to reduce calls to function propensity
  std::vector<double> prop = srn.propensity_pt(theta,initialState);

  int numReactions = prop.size();
  // Initialise Gillespie trace
  trace.insert(std::make_pair(0.0,initialState));

  /* Gillespie Algorithm */
  // ********************/
  while (time < T)
  {

    // Usual SSA where we don't care about the type of data we are working with
    trace.insert(std::make_pair(time, currentState));

    // Generate two random variables ~ U(0,1)
    r1 = rv.uniform(0,1);
    r2 = rv.uniform(0,1);;

    sumOfPropensities = 0.0;
    prop = srn.propensity_pt(theta,currentState);
    for(int i = 0; i < numReactions; i++)
    {
      sumOfPropensities += prop[i];
    }

    // Compute time of next reaction
    tau = (1.0 / sumOfPropensities) * log(1 / r1);

    // Checking cases where tau < 0 or propensity is zero which means simulation is finished
    if(tau < 0)
    {
      trace.insert(std::make_pair(time,currentState));
      trace.insert(std::make_pair(T,currentState));
      break;}
    else if(time + tau > T)
    {
      trace.insert(std::make_pair(time,currentState));
      trace.insert(std::make_pair(T,currentState));
      break;
    }

    // Reset cumulativeSum
    cumulativeSum = 0;
    // Calculate which reaction occured and adding corresponding stoichiometric vector
    for(int j = 0; j < numReactions; j++)
    {
      if( ((cumulativeSum/sumOfPropensities) <= r2 ) && (r2 < ((cumulativeSum + prop[j])/sumOfPropensities)))
      {
        for(int i = 0; i < (int) currentState.size(); i++)
        {
          currentState[i] += stoichiometricMatrix[i][j];
        }
      }
      cumulativeSum += prop[j];
    }
    // Updating time
    time += tau;
  }

  return trace;
}


// Return the SSA for each timestep only
std::map<double,std::vector<int>> gillespie(
  CRN srn,
  std::vector<double> theta,                                                    // Parameters for propensity function
  double T,                                                                     // Running time of algorithm, T
  randomVariable& rv,                                                            // Pass initialization of randomVariable class to allow for random variable distributions
  double timestep)                                                              // Only record value at the timesteps
{
    // Initialize map that tracks the times and states of the chemical reaction network during the SSA
  std::map<double,std::vector<int>> trace;

  double r1, r2, tau; // Parameters for Gillespie algorithm
  double time = 0;        // Variable tracking time
  double sumOfPropensities; // variable for the sum of propensities
  double cumulativeSum = 0; // Cumulative sum that's required to calculate which reaction occurs

  // Initial state of CTMC
  std::vector<int> initialState = srn.retInitial(), currentState = initialState;

  // Defining stoichiometricMatrix
  std::vector<std::vector<int>> stoichiometricMatrix;
  stoichiometricMatrix = srn.retStoichiometric();
  // Defining vector prop to reduce calls to function propensity
  std::vector<double> prop = srn.propensity_pt(theta,initialState);

  int numReactions = prop.size();
  // Initialise Gillespie trace
  trace.insert(std::make_pair(0.0,initialState));

  // Define two counters:
  int timestepCounter = 0; // 1) Keep track of which timestep to input into the trace file
  int timestepTracker = 0;// 2) Constrantly track how many timesteps have occured so far in trace

  // timestepTracker should be calculated for every reaction to ensure we haven't passed any timesteps
  // once we have passed a timestep, we update timestepCounter and input the
  // corresponding times (timestepCounter * timestep) and states
  // This is to ensure that if the reaction time is larger than the timestep, we
  // don't miss any of the timesteps



  /* Gillespie Algorithm NOT SAVING EVERY VALUE */
  // ********************/
  while (time < T){
    // Keep track of how many timesteps have passed thus far
    timestepTracker = time / timestep;
    if(timestepCounter < timestepTracker){
      // Another timestep has passed and timestepCounter is less than timestepTracker
      // Inputting times:
      int temp = timestepTracker - timestepCounter;
      for(int i = 0; i < temp; i++){
        timestepCounter++;
        trace.insert(std::make_pair(timestepCounter*timestep,currentState));
      }
      // timestepCounter and timestepTracker should be same
    }

    // Generate two random variables ~ U(0,1)
    r1 = rv.uniform(0,1);
    r2 = rv.uniform(0,1);
    sumOfPropensities = 0.0;
    prop = srn.propensity_pt(theta,currentState);

    for(int i = 0; i < numReactions; i++){
      sumOfPropensities += prop[i];
    }
    // Compute time of next reaction
    tau = (1.0 / sumOfPropensities) * log(1 / r1);

    // Checking cases where tau < 0 or propensity is zero which means simulation is finished
    if(tau < 0){
      timestepTracker = T / timestep;
      int temp = timestepTracker - timestepCounter;
      for(int i = 0; i < temp; i++){
        if(timestepCounter*timestep <= T)
          trace.insert(std::make_pair(timestepCounter*timestep,currentState));
        timestepCounter++;
      }
      trace.insert(std::make_pair(T,currentState));
      break;}
    else if(time + tau > T){
      timestepTracker = T / timestep;
      int temp = timestepTracker - timestepCounter;
      for(int i = 0; i < temp; i++){
        if(timestepCounter*timestep <= T)
          trace.insert(std::make_pair(timestepCounter*timestep,currentState));
        timestepCounter++;
      }
      trace.insert(std::make_pair(T,currentState));
      break;
    }

    // Reset cumulativeSum
    cumulativeSum = 0;
    // Calculate which reaction occured and adding corresponding stoichiometric vector
    for(int j = 0; j < numReactions; j++){
      if( ((cumulativeSum/sumOfPropensities) <= r2 ) && (r2 < ((cumulativeSum + prop[j])/sumOfPropensities))){
        for(int i = 0; i < (int) currentState.size(); i++){
          currentState[i] += stoichiometricMatrix[i][j];
        }
      }
      cumulativeSum += prop[j];
    }
    // Updating time
    time += tau;
  }

  return trace;
}

void printTrace(std::map<double,std::vector<int>> trace,std::vector<std::string> speciesNames)
{
  assert(speciesNames.size() == trace.begin()->second.size());
  for(std::map<double,std::vector<int>>::iterator it = trace.begin(); it != trace.end(); it++){
    std::cout << "(t";
    for(int i = 0; i < (int) speciesNames.size(); i++){
      std::cout << ", " << speciesNames[i];
    }
    std::cout << ") = (" << it->first;
    for(int i = 0; i < (int) it->second.size(); i++){
      std::cout << ", " << it->second[i];
    }
    std::cout << ")" << std::endl;
  }
}

void printTrace(std::map<double,std::vector<double>> trace,std::vector<std::string> speciesNames)
{
    assert(speciesNames.size() == trace.begin()->second.size());
    for(std::map<double,std::vector<double>>::iterator it = trace.begin(); it != trace.end(); it++){
      std::cout << "(t";
      for(int i = 0; i < (int) speciesNames.size(); i++){
        std::cout << ", " << speciesNames[i];
      }
      std::cout << ") = (" << it->first;
      for(int i = 0; i < (int) it->second.size(); i++){
        std::cout << ", " << it->second[i];
      }
      std::cout << ")" << std::endl;
    }
}

std::map<double,std::vector<int>> organiseTrace(std::map<double,std::vector<int>> inputTrace,double timestep)
{
  // Initialise new map to return after sorting
  std::map<double,std::vector<int>> sortedTrace;

  // Iterator to loop through inputted map
  std::map<double,std::vector<int>>::iterator it;
  std::map<double,std::vector<int>>::iterator LB;
  it = std::prev(inputTrace.end(),1);


  // Calculate total number of timesteps
  int numTimesteps = it->first / timestep;

  // Restart iterator from beginning of inputted map
  it = inputTrace.begin();

  sortedTrace.insert(std::make_pair(0,inputTrace.begin()->second));

  for(int i = 0; i <= numTimesteps; i++)
  {
    // Returns iterator pointing to first element that is not less than i*timestep
    LB = inputTrace.lower_bound(i*timestep);

    if(i != 0)
    {
      // Return iterator pointing to prev(LB,1) which is less than i*timestep
      sortedTrace.insert(std::make_pair(i*timestep, std::prev(LB,1)->second));
    }
  }

  sortedTrace.insert(std::make_pair(numTimesteps*timestep,it->second));
  return sortedTrace;
}


std::map<double,std::vector<double>> meanTraces(std::set<std::map<double,std::vector<int>>> tracesSet,double timestep)
{
  int n = tracesSet.size();
  std::map<double,std::vector<int>> trace, initialTrace;
  std::map<double,std::vector<int>>::iterator itTrace;
  std::map<double,std::vector<double>> mean;
  std::set<std::map<double,std::vector<int>>>::iterator itSet;

  // Average using timestep != 0 if the traces are unsorted
  if(timestep != 0)
  {
    std::set<std::map<double,std::vector<int>>> organisedSet;

    for(itSet = tracesSet.begin(); itSet!= tracesSet.end(); itSet++)
    {
      trace = *itSet;
      trace = organiseTrace(trace,timestep);
      organisedSet.insert(trace);
    }
    tracesSet = organisedSet;
  }

  // Use elements of the first inputted trace
  itSet = tracesSet.begin();
  initialTrace = *itSet;

  // Defining number of species and a temporary vector to calculate means
  int numSpecies = initialTrace[0].size();
  std::vector<double> tempVec(numSpecies,0);

  double time = 0;
  // Looping over the times in the traces followed by the individual traces themselves
  for(itTrace = initialTrace.begin(); itTrace != initialTrace.end(); itTrace++)
  {

    time = itTrace->first;
    for(auto x : tracesSet)
    {
      for(int i = 0; i < numSpecies; i++)
      {
        tempVec[i] += (double) (x[time][i])/n;
      }
    }
    mean.insert(std::make_pair(time,tempVec));

    for(int i = 0; i < numSpecies; i++)
    {
      tempVec[i] = 0;
    }
  }
  return mean;
}


void outputTrace(std::string directory,std::map<double,std::vector<double>> trace,std::vector<std::string> speciesNames)                                                     // print simulated trace to screen
{
  std::ofstream output;
  output.open(directory);

  output << "Time,";
  for(int kk = 0; kk < (int) speciesNames.size();kk++)
  {
    output << speciesNames[kk] << ",";
  }
  output << "\n";

std::map<double,std::vector<double>>::iterator itr;

  for(itr = trace.begin(); itr != trace.end(); itr++)
  {
    output << itr->first;
    for(int ii = 0; ii < (int) itr->second.size(); ii++)
    {

      output << "," << itr->second[ii];
    }
    output << "\n";
  }
};

void outputTrace(std::string directory,std::map<double,std::vector<int>> trace,std::vector<std::string> speciesNames)                                                     // print simulated trace to screen
{
  std::ofstream output;
  output.open(directory);

  output << "Time,";
  for(int kk = 0; kk < (int) speciesNames.size();kk++)
  {
    output << speciesNames[kk] << ",";
  }
  output << "\n";

std::map<double,std::vector<int>>::iterator itr;

  for(itr = trace.begin(); itr != trace.end(); itr++)
  {
    output << itr->first;
    for(int ii = 0; ii < (int) itr->second.size(); ii++)
    {

      output << "," << itr->second[ii];
    }
    output << "\n";
  }
};

double d(std::map<double,std::vector<double>> X, std::map<double,std::vector<double>> Y){
  // Function to calculate the L2 Norm distance between the observed data, X and
  // the simulated data, Y.
  int n = X.size(), m = Y.size();
  assert(n == m);
  int d = X[0].size();            // Data dimension

  std::map<double,std::vector<double>>::iterator itX;
  double time = 0;
  double l2norm = 0;
  for(itX = X.begin(); itX != X.end(); itX++){
    time = itX->first;
    for(int i = 0; i < d; i++){
      l2norm += pow((X[time][i] - Y[time][i]),2.0);
    }
  }

  l2norm = sqrt(l2norm);

  return l2norm;
}

double d(std::map<double,std::vector<double>> X, std::map<double,std::vector<int>> Y){
  // Function to calculate the L2 Norm distance between the observed data, X and
  // the simulated data, Y.
  int n = X.size(), m = Y.size();

  assert(n == m);
  int d = X[0].size();            // Data dimension

  std::map<double,std::vector<double>>::iterator itX;
  double time = 0;
  double l2norm = 0;
  for(itX = X.begin(); itX != X.end(); itX++){
    time = itX->first;
    for(int i = 0; i < d; i++){
      l2norm += pow((X[time][i] - Y[time][i]),2.0);
    }
  }

  l2norm = sqrt(l2norm);

  return l2norm;
}
