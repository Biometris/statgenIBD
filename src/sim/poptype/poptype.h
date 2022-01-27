// toevoegen van nieuwe populatietype NEW
// A. maak een nieuwe class NEWpop, afgeleid van abstracte classe pop_base
// B. pas constructor van PopulationType aan.

#ifndef POPTYPE_HEADER
#define POPTYPE_HEADER

#include <vector>
#include "util_genetics.h"
#include "Genotype.h"
#include "ObsGeno.h"
#include "Loc.h"
#include "Genome.h"
#include "pop_base.h"
#include "matvec.h"

namespace ibd
{

class PopulationType
{
public:
	PopulationType() : h_popt(0) {}
	PopulationType(const std::string& popname);

	double p_trans(const Genotype& g1, const Genotype& g2, double r) const 
		{ return h_popt->p_trans(g1,g2,r); }
	double p_trans_QTL(const Genotype& g_parent, const Genotype& g_prog) const
		{  return h_popt->p_trans_QTL(g_parent, g_prog); }
	std::vector<double> model(const Genotype& g) const { return h_popt->model(g); }
	Genome make_ind(const std::vector<double>& chr_length) const
		{ return h_popt->make_ind(chr_length);}
	void MapQTL(std::ostream& outp, const ObsGeno& g) const { h_popt->MapQTL(outp,g); }
	ObsGeno MapQTL(std::istream& inp) const { return h_popt->MapQTL(inp); }
	std::string name() const { return h_popt->name(); }
	int nparam() const { return h_popt->nparam(); } 
	ObsGeno unknown_genotype() const { return h_popt->unknown_genotype(); }
	//pop_base * get_ptr() const { return h_popt.get(); }

private:
	SmartPtr<pop_base> h_popt;
};

matrix<double> transition_matrix(const ObsGeno& g1, 
								 const ObsGeno& g2, 
								 double r, 
								 const PopulationType& poptype);


int calcprob_ind(matrix<Genotype>& geno_augm,
				 std::vector<double>& P,
				 const std::vector<double>& r,
				 const std::vector<ObsGeno>& obs_geno,
				 const PopulationType& poptype,
				 double threshold = 0.0);

void calcprob_pop(std::vector<int>& N_ext,
				  matrix<Genotype>& geno_augm,
				  std::vector<double>& P,
				  const LinkageMap& Markermap,
				  const matrix<ObsGeno>& obs_geno,
				  const PopulationType& poptype, 
				  double threshold = 0.0);

std::vector<double> calc_P_stepII(const matrix<Genotype> geno_augm, 
							 const std::vector<double>& P,
							 const LinkageMap& markermap,
							 const PopulationType& popt,
							 const Locus& QTLpos);

}


#endif