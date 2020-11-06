#ifndef GILLESPIE_H_
#define GILLESPIE_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <limits>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "dist.h"
#include "CRN.h"


std::map<double,std::vector<int>> gillespie(CRN srn, std::vector<double> theta, double T,randomVariable& rv);
std::map<double,std::vector<int>> gillespie(CRN srn, std::vector<double> theta, double T,randomVariable& rv, double timestep);


void printTrace(std::map<double,std::vector<int>> trace, std::vector<std::string> speciesNames);
void printTrace(std::map<double,std::vector<double>> trace, std::vector<std::string> speciesNames);
std::map<double,std::vector<int>> organiseTrace(std::map<double,std::vector<int>> trace,double timestep);

std::map<double,std::vector<double>> meanTraces(std::set<std::map<double,std::vector<int>>> tracesSet, double timestep = 0.0);

void outputTrace(std::string directory,std::map<double,std::vector<int>> trace,std::vector<std::string> speciesNames);                                                     // print simulated trace to screen
void outputTrace(std::string directory,std::map<double,std::vector<double>> trace,std::vector<std::string> speciesNames);                                                     // print simulated trace to screen

double d(std::map<double,std::vector<double>> X, std::map<double,std::vector<int>> Y);
double d(std::map<double,std::vector<double>> X, std::map<double,std::vector<double>> Y);
#endif
