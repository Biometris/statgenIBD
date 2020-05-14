/*!
\file 
\brief  Definition of Isome aux. functions
\author Martin Boer, Biometris
\date   2006-2010
*/
#ifndef INDPROP_HEADER
#define INDPROP_HEADER

#include <set>
#include <vector>
#include <string>

#include "util_genetics.h"

std::vector<int> get_ndx_set(const std::vector<IndProp>& pop, 
							 const std::set<std::string>& setID);

std::set<std::string> get_ind_coa_file(const std::vector<IndProp>& pop, 
									   const std::string& filename_coa);

std::set<std::string> get_parents_families(const std::vector<IndProp>& pop);
std::set<std::string> get_ind_within_fam(const std::vector<IndProp>& pop);

bool single_family(const std::vector<IndProp>& pop);



#endif


