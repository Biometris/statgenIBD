// /*!
// \file
// \brief  Output in tab-delimited format
// \author Martin Boer, Biometris
// \date   2006-2010
// */
// #ifndef OUTPUT_TAB_PEDIGREE_HEADER
// #define OUTPUT_TAB_PEDIGREE_HEADER
//
// #include <iostream>
// #include <vector>
// #include "Loc.h"
// #include "pedigree.h"
// #include "markerscore.h"
// #include "analysis_ped.h"
//
// void outp_tab_header(std::ostream& outp,
// 					 const LinkageMap& markermap,
// 					 const std::vector<IndProp>& pop,
// 				     const std::set<std::string>& sel_ind);
//
// void outp_tab_fam(std::ostream& outp,
// 			      const Pedigree& ped_fam,
// 				  const mbl::matrix3D<double>& Z,
// 				  const LinkageMap& markermap);
//
// void outp_tab_ped_ftp(std::ostream& outp,
// 						const std::vector<int>& ndx,
// 						const IBDped& IBD,
// 						const Pedigree& ped,
// 						const LinkageMap& markermap);
//
// void outp_tab_ped_coa(std::ostream& outp,
// 						const std::vector<int>& ndx,
// 						const IBDped& IBD,
// 						const LinkageMap& markermap);
//
// #endif
//
