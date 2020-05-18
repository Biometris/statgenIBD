/*!
 \file
 \brief  Analysis of biparental, 3-way and 4-way crosses
 \author Martin Boer, Biometris
 \date   2006-2010
 */
#ifndef CROSSES_HEADER
#define CROSSES_HEADER

#include <vector>
#include <map>

#include "Args.h"

#include "matvec.h"
#include "popt.h"
#include "pedigree.h"
#include "markerscore.h"
#include "IndProp.h"

int count_parents(const std::vector<IndProp>& pop);

mbl::matrix3D<double> analysis_cross(const std::vector<IndProp>& pop,
                                     const mbl::matrix<score>& geno,
                                     const LinkageMap& markermap,
                                     const LinkageMap& eval_pos,
                                     //const std::string& outp,
                                     const mbl::Args& Argu);

#endif
