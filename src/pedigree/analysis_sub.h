/*!
\file 
\brief  IBD-analysis of subpopulations using HMM 
\author Martin Boer, Biometris
\date   2006-2010

The calculation of IBD-probabilities in small populations of inbred lines, using Hidden
Markov Models (HMMs). It is assumed that the number non-founders is less than 15. At the moment
we only calculate the transition probabilities on marker positions.

For the calculation of IBD-probabilities on arbitrary positions on the genome we can define a 
function objects of the form:

\code
	std::vector<double> TransProb(mbl::Locus& QTLpos)
\endcode
*/
#ifndef ANALYSIS_SUBPOPULATIONS_HEADER
#define ANALYSIS_SUBPOPULATIONS_HEADER

#include <vector>
#include "matvec.h"
#include "pedigree.h"

mbl::matrix<double> 
calc_IBD_prob_last_ind(const Pedigree& ped,
					   const mbl::matrix<int>& geno,
					   const std::vector<double>& r,double eps);

#endif



