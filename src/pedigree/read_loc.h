/*!
\file 
\brief  Read marker scores
\author Martin Boer, Biometris
\date   2006-2010
*/
#ifndef READ_LOC_PEDIGREE_HEADER
#define READ_LOC_PEDIGREE_HEADER

#include <iostream>
#include <vector>
#include <string>
#include <map>

// library files
#include "matvec.h"
#include "Loc.h"

#include "markerscore.h"
#include "IndProp.h"

int read_loc_file(std::vector<std::string>& ID, 
				  std::vector<std::string>& marker_names,
				  mbl::matrix<score>& geno_obs,
				  const std::string& locfile);

int convert_loc_file();


#endif



