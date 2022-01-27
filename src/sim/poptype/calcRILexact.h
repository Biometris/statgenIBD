#ifndef CALC_RIL_EXACT_HEADER
#define CALC_RIL_EXACT_HEADER

#include <vector>
#include "util_genetics.h"
#include "Genotype.h"
#include "ObsGeno.h"
#include "Loc.h"
#include "poptype.h"

namespace ibd
{

int countbits(unsigned int x);
bool shift_right(unsigned int& x);

typedef unsigned int InhVec; 
typedef unsigned int Group;

// all Inheritance vectors consistent with marker scores
typedef std::vector<InhVec> ObsInhVec; 

std::vector<ibd::ObsInhVec> obs_inh_vectors(const std::vector<ibd::ObsGeno>& g, int Nind);

ibd::matrix<double> calc_prob_left(const std::vector<ObsInhVec>& obs_inh_vec, 
								   int M, const std::vector<double>& r);

ibd::matrix<double> calc_prob_right(const std::vector<ObsInhVec>& obs_inh_vec, 
									int M, const std::vector<double>& r);

ibd::matrix<double> calc_QTL_prob(const ibd::ObsGeno& U,
								  const LinkageMap& markermap,
								  const ibd::matrix3D<double>& Tr_l,
								  const ibd::matrix3D<double>& Tr_r,
							 int nind, int M,
							 const Locus& QTLpos);

int test_calc_prob_RIL();

}

#endif