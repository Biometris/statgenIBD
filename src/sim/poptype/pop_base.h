// toevoegen van nieuwe populatietype NEW
// A. maak een nieuwe class NEWpop, afgeleid van pop_base
// B. pas constructor van PopulationType aan.

#ifndef POP_BASE_HEADER
#define POP_BASE_HEADER

#include <vector>
#include "util_genetics.h"
#include "Genotype.h"
#include "ObsGeno.h"
#include "Loc.h"
#include "Genome.h"

namespace ibd
{

class pop_base
{
public:
	pop_base(const ObsGeno& u, int p, const std::string& name_pop) 
		: U(u), npar(p), p_name(name_pop) {}
	virtual ~pop_base() { }
	virtual double p_trans(const Genotype& g1, const Genotype& g2, double r) const = 0;
	virtual double p_trans_QTL(const Genotype& g1, const Genotype& g2) const
		{ return (g1==g2) ? 1.0 : 0.0; } 
	virtual std::vector<double> model(const Genotype& g) const = 0;
	virtual Genome make_ind(const std::vector<double>& chr_length) const = 0;
	virtual std::vector<Genome> make_pop(const std::vector<double>& chr_length, 
		int nind) const = 0;
	virtual void MapQTL(std::ostream& outp, const ObsGeno& g) const = 0;
	virtual ObsGeno MapQTL(std::istream& inp) const = 0;

	std::string name() const { return p_name; }
	int nparam() const { return npar; }  
	ObsGeno unknown_genotype() const { return U; }
protected:
	const ObsGeno U;
	const int npar;
	const std::string p_name;
};

}


#endif