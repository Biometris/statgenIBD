#ifndef GENOPROB_HEADER
#define GENOPROB_HEADER

#include <vector>
#include "util_genetics.h"
#include "Genotype.h"
#include "ObsGeno.h"
#include "Loc.h"
#include "poptype.h"
#include "calcRILexact.h"

namespace ibd
{

extern bool EXACT_RIL;

class GenoProb
{
public:
	GenoProb(const ibd::matrix<ibd::ObsGeno>& obs_geno,
			const ibd::PopulationType& popt,
			const LinkageMap& markermap);
	ibd::matrix<double> operator()(const Locus& QTLpos) const;
	int Nind() const { return nind; } 
private:
	LinkageMap	      markermap_;
	ibd::PopulationType	  popt_;
	ibd::matrix3D<double> Tr_l, Tr_r;
	ibd::ObsGeno		  U;
	int			          dimU, nind, nloc;
};

ibd::matrix<double> calc_prob_left(const std::vector<ibd::ObsGeno>& obs_geno,
								   const ibd::PopulationType& popt,
								   const std::vector<double>& r);

ibd::matrix<double> calc_prob_right(const std::vector<ibd::ObsGeno>& obs_geno,
									const ibd::PopulationType& popt,
									const std::vector<double>& r);

int IsRIpop(const ibd::PopulationType& popt);
}

#endif