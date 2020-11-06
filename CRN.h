#ifndef CRN_H_
#define CRN_H_

#include "dist.h"
#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <filesystem>

class CRN
{
private:
  // Number of species
  int numSpecies;
  // Stoichiometric matrix for Gillespie algorithm updates
  std::vector<std::vector<int>> stoichiometric;
  // Vector of strings to store names of chemical species
  std::vector<std::string> speciesNames;
  // Vector of strings to store parameter names
  std::vector<std::string> parameterNames;
  // Parameter bounds to ensure the perturbed parameter remain in prior
  std::vector<std::vector<double>> parameterBounds();
  // Initial molecule number
  std::vector<int> initialCount;


public:
  // Class constructor
  CRN(std::vector<std::vector<int>> S);

  // Return the stoichiometric matrix, S
  std::vector<std::vector<int>> retStoichiometric();
  // Assign speciesNames
  void inputtingNames(std::vector<std::string> names);
  // Return the species names
  std::vector<std::string> retNames();

  // Assign speciesNames
  void inputtingParameterNames(std::vector<std::string> paramNames);
  // Return the species names
  std::vector<std::string> retParameterNames();

  // Set initial number of molecules
  void setInitial(std::vector<int> initial);

  std::vector<int> retInitial();

  // Prior probability distribution for CRN parameters
  std::vector<double> (*prior_pt) (randomVariable& );
  // Prior probability distribution function for CRN parameters
  double (*prior_pdf_pt) (std::vector<double>, randomVariable& );

  // Propensity function for CRN
  std::vector<double> (*propensity_pt)(std::vector<double>, std::vector<int>);

};


void print_trace(std::string directory, CRN srn,std::map<double,std::vector<int>> trace);
void print_trace(std::string directory, CRN srn,std::map<double,std::vector<double>> trace);
#endif
