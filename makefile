# *****************************************************
# Variables to control Makefile operation

CXX = g++
CXXFLAGS = -Wall -g -ggdb3
DEPS =  dist.h CRN.h Gillespie.h ABCSMC.h   SIR.h  Massart.h BayesianSMC.h
OBJ = main.o dist.o CRN.o Gillespie.o ABCSMC.o  SIR.o Massart.o BayesianSMC.o

# ****************************************************
# Targets needed to bring the executable up to date

%.o: %.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<


main_ABCSMCSMC: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lgsl -lgslcblas -lm -I/usr/local/include
