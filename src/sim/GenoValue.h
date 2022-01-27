#ifndef GENOVALUE_SIMQTL_HEADER
#define GENOVALUE_SIMQTL_HEADER

#include "Loc.h"
#include "matvec.h"
#include "Genotype.h"
#include "Genome.h"

class Phi
{
public:
	Phi(const std::map<Locus,std::vector<double> >& QTLs, 
		const std::map<std::string,std::string>& inbfnds, 
		const ibd::matrix<double>& epi_addxadd, double mu);
	double operator()(const ibd::Genome& genome) const;
	int add_indicator(const std::vector<ibd::Genotype>& g, int q) const;
	int dom_indicator(const std::vector<ibd::Genotype>& g, int q) const;
	int NumberOfQTL() const { return nQTL; } 
	std::vector<std::string> getQTLs(const ibd::Genome& genome) const;
	std::vector<int> getORG(const ibd::Genome& genome, const ibd::Parent& par) const;
private:
	double mu_;
	int nQTL,nFND;
	std::vector<double> add,dom;	
	ibd::matrix<double> epi_axa;
	ibd::matrix<bool> pos_allele;     // ncol: Nqtl, nrow: Nfnd
	LinkageMap QTLmap;        
};

#endif
