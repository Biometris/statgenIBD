/*!
 \file
 \brief  Analysis of biparental, 3-way and 4-way crosses
 \author Martin Boer, Biometris
 \date   2006-2010
 */
#ifndef CROSSES_HEADER
#define CROSSES_HEADER

#include <vector>
#include <RcppArmadillo.h>

#include "matvec.h"
#include "util_genetics.h"
#include "markerscore.h"

int count_parents(const std::vector<IndProp>& pop);

arma::cube analysis_cross(const std::vector<IndProp>& pop,
                          const ibd::matrix<score>& geno,
                          const LinkageMap& markermap,
                          const LinkageMap& eval_pos);

#endif
