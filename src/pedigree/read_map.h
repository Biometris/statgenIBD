/*!
\file 
\brief  Read markermap and calculate recombinations.
\author Martin Boer, Biometris
\date   2006-2010
*/

#ifndef READ_MAP_PEDIGREE_HEADER
#define READ_MAP_PEDIGREE_HEADER

#include <iostream>
#include <vector>
#include <string>
#include <map>

// library file
#include "Loc.h"
#include "Args.h"

// output: LinkageMap with only the markers which are in the loc file
LinkageMap reduce_markermap(const LinkageMap& markermap, const std::vector<std::string>& markers);		

// if distance between markers is less than eps=1.0e-6 then distance is 
// assumed to be delta = 0.001
LinkageMap adjust_markermap(const LinkageMap& markermap);

LinkageMap read_map_file(const std::string& mapfile);
//int NrChromosomes(const LinkageMap& markermap);

void print_marker_warnings(const LinkageMap& markermap, const std::vector<std::string>& markers);

std::vector<double> make_rec_map_inbred_lines(const LinkageMap& markermap);

std::vector<bool> get_selected_chr(const mbl::Args& Argu, 
								   const LinkageMap& markermap);

LinkageMap read_eval_pos_file(const std::string& filename);

LinkageMap select_chr(const LinkageMap& markermap, int sel_chr);

#endif



