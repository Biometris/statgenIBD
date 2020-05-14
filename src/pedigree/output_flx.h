/*!
\file 
\brief  Output for FlexQTL
\author Martin Boer, Biometris
\date   2006-2010
*/
#ifndef OUTPUT_FLX_PEDIGREE_HEADER
#define OUTPUT_FLX_PEDIGREE_HEADER

#include <iostream>
#include <vector>
#include "Loc.h"
#include "pedigree.h"
#include "markerscore.h"
#include "analysis_ped.h"

void outp_flx_header(std::ostream& outp, 
					 const LinkageMap& markermap,
					 const std::vector<IndProp>& pop,
					 const std::set<std::string>& sel_ind);

void outp_flx_ped_coa(std::ostream& outp,
						const std::vector<int>& ndx, 
						const IBDped& IBD,
						const LinkageMap& markermap);

#endif

