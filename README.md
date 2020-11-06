# ABCSMCSMC
Code for ABC(SMC)^2: SIR Case study 

Code for our CMSB2020 paper, ABC(SMC)^2: Simultaneous inference and model checking of chemical reaction networks with the SIR case study considered. GSL (https://www.gnu.org/software/gsl/) is needed to compile and run code. 

The main bulk of the code that generated the data used in the paper is within the SIR.cpp source file. In the case study here we initialise the SIR chemical reaction network, define priors and appropriate propensities associated with the CRN (needed in order to simulate using Gillespie/ SSA Algorithm) and generate true data, then run the ABCSMCSMC algorithm to obtain the outputs of the algorithm, namely the parameter samples, associated weights and the mean, upper and lower confidence bound of the satisfaction probabilities. 

Please feel free to email me (gareth.molyneux@cs.ox.ac.uk or garethwynmolyneux@gmail.com) if you have further comments or questions.
