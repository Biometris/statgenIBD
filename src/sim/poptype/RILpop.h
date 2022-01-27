#ifndef RILPOP_HEADER
#define RILPOP_HEADER

#include "poptype.h"

namespace ibd
{

class RILpop : public pop_base
{
public:
	RILpop(alleletype p1, alleletype p2, int generation);
	double p_trans(const Genotype& g1, const Genotype& g2, double r) const;
	std::vector<double> model(const Genotype& g) const;
	Genome make_ind(const std::vector<double>& chr_length) const;
	std::vector<Genome> make_pop(const std::vector<double>& chr_length, int nind) const;
	void MapQTL(std::ostream& outp, const ObsGeno& g) const;
	ObsGeno MapQTL(std::istream& inp) const;
private:
	const Genotype AA,AB,BB;
	int gen;
};

}

#endif