/*!
\file 
\brief  general output of pedigree program
\author Martin Boer, Biometris
\date   2006-2010
*/

#ifndef OUTPUT_PEDIGREE_HEADER
#define OUTPUT_PEDIGREE_HEADER

#include <iostream>
#include <vector>
#include "Args.h"
#include "Loc.h"
#include "pedigree.h"
#include "markerscore.h"
#include "analysis_ped.h"

// Three types of output:
// BIN: binary output (-bin option)
// TAB: tab delimited file 
// FLX: input file for FlexQTL file (-flx option)
enum OutputType {BIN,TAB,FLX};

void outp_header(std::ostream& outp, 
					 const LinkageMap& markermap,
					 const std::vector<IndProp>& pop,
					 const std::set<std::string>& sel_ind,
					 OutputType type);

void outp_fam(std::ostream& outp,
					const Pedigree& ped_fam, 
					const mbl::matrix3D<double>& IBD,
					const LinkageMap& markermap,
					OutputType type);

void outp_ped_ftp(std::ostream& outp,
						const std::vector<int>& ndx, 
						const IBDped& IBD,
						const Pedigree& ped,
						const LinkageMap& markermap,
						OutputType type);

void outp_ped_coa(std::ostream& outp,
					const std::vector<int>& ndx, 
					const IBDped& IBD,
					 const LinkageMap& markermap,
					OutputType type);

OutputType open_output(std::ofstream& outp, 
					   const mbl::Args& Argu,
					   const std::string& output, 
					   bool coa_parents_fam);

#endif

