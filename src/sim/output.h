#ifndef OUTPUT_SIMQTL_HEADER
#define OUTPUT_SIMQTL_HEADER

#include "util_genetics.h"

#include <map>
#include <string>
#include "input.h"
#include "MarkerType.h"
#include "GenoValue.h"

// void make_eval_file(const std::string& eval_filename,
// 					std::vector<double>& chr_length,
// 					double dist_eval_pos);

void make_ped_file(const SimPop& sim_pop, const std::string& filename);

void make_loc_file(const SimPop& sim_pop,
				   const LinkageMap& markermap,
				   const std::vector<MarkerType>& markertype,
			       const std::string& filename);

void make_qua_file(const SimPop& sim_pop, const Phi& phi, double sigma,
				   const std::string& filename);

void print(std::ostream& outp, const std::map<Locus,std::vector<double> >& QTLs);

void print_genstat_table(const std::vector<PopProp>& pops,
						 const std::map<std::string,std::string>& inbfnds,
						 const std::string filename);

void make_part_ldf_file(int N,std::vector<double>& chr_length,
						double dist_eval_pos, std::string filename);

void make_flapjack_map_file(std::string name, const LinkageMap& marker_map);


#endif
