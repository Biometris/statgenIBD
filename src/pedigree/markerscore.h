/*!
 \file
 \brief  Definition of marker scores.
 \author Martin Boer, Biometris
 \date   2006-2010
 */
#ifndef MARKERSCORES_HEADER
#define MARKERSCORES_HEADER

#include <string>
#include <iostream>

#include "matvec.h"
#include "OrdGeno.h"

//! Marker score
/*!
 Score for a marker. It is assumed that scores are pairs of non-negative integers.
 */
class score : public std::pair<int,int>
{
public:
  score() {}
  score(int a, int b);
  bool homozygous() const;
  std::string print_string() const;
};

std::ostream& operator<<(std::ostream&, const score& sc);

const int U_haplo_sc = -1;
const score Uscore = score(U_haplo_sc,U_haplo_sc);

score read_score(std::istream& s,
                 char delimit);
bool check_score(const mbl::OrdGeno& g,
                 const score& sc);

std::map<score,int> ndx_score(int nparents);

int read_flapjackfile(std::vector<std::string>& geno,
                      std::vector<std::string>& markers,
                      mbl::matrix<score>& scores,
                      const std::string filename);


#endif

