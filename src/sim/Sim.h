#ifndef SIM_SIMQTL_HEADER
#define SIM_SIMQTL_HEADER

#include "util_genetics.h"
#include "SimPop.h"
#include <string>

typedef std::pair<IndProp,ibd::Genome> SimInd;
typedef std::vector<SimInd> SimPop;

SimPop sim_population(const PopProp& pop, const std::map<std::string,ibd::Genome>& fnd, const std::string& pop_name);

std::map<std::string,ibd::Genome>
sim_pedigree(const std::vector<IndProp>& pedigree, const std::map<std::string,ibd::Genome>& founders);

void sim_SS_NS_pedigree(const std::vector<std::string>& fnd_names,
						std::string filename, int Ngen, int N);


SimPop sim_multiple_populations(const std::vector<PopProp>& pops, const std::map<std::string,ibd::Genome>& genome_par);

#endif
