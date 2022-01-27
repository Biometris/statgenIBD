#ifndef MARKERTYPE_SIMQTL_HEADER
#define MARKERTYPE_SIMQTL_HEADER

#include "util_genetics.h"  

#include <string>
#include "simQTL.h"
#include "Genotype.h"

class MarkerType
{
public:
	MarkerType(int Nfnd, int Nalleles, double perc_missing);
	std::string operator()(ibd::Genotype g) const; 
	std::string First(ibd::Genotype g) const;
	std::string Second(ibd::Genotype g) const;
private:
	std::vector<std::string> table;
	double percentage_missing;
};

std::vector<MarkerType> 
generate_markertypes(int nloc,int nfnd, int nalleles, double percentage_missing);

#endif