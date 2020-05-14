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

#include "pedigree.h"
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

//std::ostream& operator<<(std::ostream&, const mbl::OrdGeno& g);

const int U_haplo_sc = -1;
const score Uscore = score(U_haplo_sc,U_haplo_sc);

score read_score(std::istream& s,char delimit);
bool check_score(const mbl::OrdGeno& g, const score& sc);

mbl::matrix<int> gen_haplo_score_ped(std::vector<IndProp>& sel_pop,
								  const mbl::matrix<score>& obsgeno, 
								  const std::vector<IndProp>& pop);

// haplo scores for families + parents families
mbl::matrix<int> gen_haplo_score_fam(std::vector<IndProp>& sel_pop,
									 const mbl::matrix<score>& obsgeno, 
									 const std::vector<IndProp>& pop, 
									 const std::set<std::string>& Pfam);

mbl::matrix<score> select_families(std::vector<IndProp>& sel_pop,
									const mbl::matrix<score>& obsgeno, 
									const std::vector<IndProp>& pop, 
									const std::set<std::string>& Pfam);

std::vector<bool> detect_high_density_ind(const mbl::matrix<int>& geno, int min);
std::map<score,int> ndx_score(int nparents);

int test_read();

std::map<score,int> ndx_score(int nparents);
int tst_ndx_score();

int read_flapjackfile(std::vector<std::string>& geno, std::vector<std::string>& markers,
    mbl::matrix<score>& scores, const std::string filename);


#endif

