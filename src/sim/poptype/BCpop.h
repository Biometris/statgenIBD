#ifndef BCPOP_HEADER
#define BCPOP_HEADER

#include "poptype.h"

namespace ibd
{

class BCpop : public pop_base
{
public:
	BCpop(alleletype p1, alleletype p2) : AA(Genotype(p1,p1)), AB(Genotype(p1,p2)),
					pop_base(ObsGeno(Genotype(p1,p1),Genotype(p1,p2)),1,"BC1") {}
	double p_trans(const Genotype& g1, const Genotype& g2, double r) const;
	std::vector<double> model(const Genotype& g) const;
	Genome make_ind(const std::vector<double>& chr_length) const;
	std::vector<Genome> make_pop(const std::vector<double>& chr_length, int nind) const;
	void MapQTL(std::ostream& outp, const ObsGeno& g) const;
	ObsGeno MapQTL(std::istream& inp) const;
private:
	const Genotype AA,AB;
};

}

#endif