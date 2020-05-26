/*!
 \file
 \brief  IBD-analysis of families
 \author Martin Boer, Biometris
 \date   2006-2010
 */
#ifndef ANALYSIS_FAMILIES_HEADER
#define ANALYSIS_FAMILIES_HEADER

#include <string>
#include <vector>
#include <map>

#include "matvec.h"
#include "popt.h"
#include "markerscore.h"
#include "OrdGeno.h"
#include "Loc.h"

//! calculation of IBD probabilities in families at QTL positions
/*!
 Function Object for calculation of IBD probabilities at putative QTL positions
 */
class IBD_fam
{
public:
  IBD_fam(const mbl::matrix<mbl::OrdGeno>& P,
          const std::vector<score>& Offspring,
          const LinkageMap& MarkerMap,
          const std::string& poptype);
  ~IBD_fam() { delete popt; }
  std::map<mbl::OrdGeno,double> operator()(const Locus& QTLpos) const;
  std::vector<double> check_scores(const std::vector<mbl::OrdGeno>& geno, const score& sc_off) const;
private:
  pop_base *popt;
  int len_inh; // length of inheritance vector
  std::vector<mbl::OrdGeno> par;
  LinkageMap markermap;
  mbl::matrix<double> l_cond, r_cond;
};

#endif


