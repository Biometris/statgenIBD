// /*!
// \file
// \brief  IBD-analysis of pedigree
// \author Martin Boer, Biometris
// \date   2006-2010
//
// Analysis of complete MEPS pedigree, i.e. defined pedigree on top of the families
// , i.e. starting with the founders (F) on the top, and the parents of the families (P),
// see e.g. powerpoint presentations of Pioneer-WUR meeting Nov 2006 for more details.
//
// In this module we split the original pedigree in subpopulations, using the generalized
// tabular method. At the moment we analyse the IBD-transition on marker positions. The method can be
// further generalized by defining a functions object (functor) of the form:
//
// \code
// 	mbl::matrix<double> FndToPar(mbl::Locus& QTLpos)
// \endcode
// */
// #ifndef ANALYSIS_PEDIGREE_HEADER
// #define ANALYSIS_PEDIGREE_HEADER
//
// #include <set>
// #include <vector>
// #include <string>
//
// #include "matvec.h"
// #include "pedigree.h"
// #include "markerscore.h"
//
// class IBDped
// {
// public:
// 	IBDped() {}
//     IBDped(const Pedigree& ped,
// 			  const mbl::matrix3D<double>& IBD_founders,
// 			  const mbl::matrix<int>& geno,
// 			  const std::vector<bool>& ancestor,
// 			  const std::vector<double>& r,
// 			  double eps);
// 	mbl::matrix<double> operator()(int m) const;
// 	int Nloc() const { return nloc; }
// private:
// 	int nloc;
// 	int nfnd;
// 	mbl::matrix3D<double> IBD_fnd;
// 	std::vector< std::vector<int> > anc;
// 	std::vector< mbl::matrix<double> > prob;
// };
//
//
// mbl::matrix<double> calc_coa(const Pedigree& ped);
//
// void pre_analysis_tabular_method(const Pedigree& ped,
// 								 const std::vector<bool>& ancestor);
//
// // int analysis_pedigree(IBDped& IBD,				// output
// // 					  Pedigree& ped,		    // output
// // 					  std::vector<int>& ndx,    // output
// // 					  const mbl::matrix3D<double>& IBD_fnd,
// // 					  const std::vector<IndProp>& pop,
// // 					  const mbl::matrix<score>& geno,
// // 					  const std::vector<double>& r,
// // 					  const std::set<std::string>& sel_ind,
// // 					  double fraction);
//
// /*
// void generalized_tabular_method_coa(mbl::matrix3D<double>& X,
// 							   const Pedigree& ped,
// 							   const mbl::matrix<int>& geno,
// 							   const std::vector<bool>& ancestor,
// 							   const std::vector<double>& r,
// 							   double eps);
//
// */
//
// #endif
//
//
//
//
