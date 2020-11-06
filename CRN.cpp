#include "CRN.h"

CRN::CRN(std::vector<std::vector<int>> S){
    this->stoichiometric = S;
    this->numSpecies = S.size();
}

std::vector<std::vector<int>> CRN::retStoichiometric(){
  return this->stoichiometric;
}

void CRN::inputtingNames(std::vector<std::string> names){
  this->speciesNames = names;
}

std::vector<std::string> CRN::retNames(){
  return this->speciesNames;
}

void CRN::inputtingParameterNames(std::vector<std::string> names){
  this->parameterNames = names;
}

std::vector<std::string> CRN::retParameterNames(){
  return this->parameterNames;
}

void CRN::setInitial(std::vector<int> initial){
  this->initialCount = initial;
}

std::vector<int> CRN::retInitial(){
  return this->initialCount;
}

void print_trace(std::string directory, CRN srn, std::map<double,std::vector<int>> trace){
  // Print out all sampled particles and corresponding probabilities
  std::ofstream output;
  output.open(directory + ".csv");

  output << "t,";
  for(int kk = 0; kk < (int) srn.retNames().size();kk++){
    output << srn.retNames()[kk] << ",";
  }
  output << "\n";

  std::map<double,std::vector<int>>::iterator itr;

  for(itr = trace.begin(); itr != trace.end(); itr++){
    output << itr->first << ",";
    for(int i = 0; i < (int) trace[0].size(); i++){
      output << itr->second[i] << ", ";
    }
    output << "\n";
  }

};

void print_trace(std::string directory, CRN srn, std::map<double,std::vector<double>> trace){
  // Print out all sampled particles and corresponding probabilities
  std::ofstream output;
  output.open(directory + ".csv");

  output << "t,";
  for(int kk = 0; kk < (int) srn.retNames().size();kk++){
    output << srn.retNames()[kk] << ",";
  }
  output << "\n";

  std::map<double,std::vector<double>>::iterator itr;

  for(itr = trace.begin(); itr != trace.end(); itr++){
    output << itr->first << ",";
    for(int i = 0; i < (int) trace[0].size(); i++){
      output << itr->second[i] << ", ";
    }
    output << "\n";
  }

};
