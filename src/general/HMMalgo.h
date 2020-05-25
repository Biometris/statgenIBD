/*!
 \file
 \brief  Hidden Markov Models algorithms
 \author Martin Boer, Biometris
 \date   1998-2007
 */
#ifndef HMM_ALGORITHM
#define HMM_ALGORITHM

#include <iostream>
#include <numeric>
#include <iomanip>
#include <vector>
#include <cmath>
#include "matvec.h"
#include "TransMatSym2D.h"

namespace mbl
{

std::vector<double> forward_equation(const std::vector<double>& p_prev,
                                     const TransMatSym2D& T,
                                     const std::vector<double>& q);

std::vector<double> backward_equation(const std::vector<double>& p_next,
                                      const TransMatSym2D& T,
                                      const std::vector<double>& q);

matrix<double> calc_prob_left(const std::vector<double>& pi0,
                              const matrix<double>& q,
                              const std::vector<TransMatSym2D>& T);

matrix<double> calc_prob_right(const matrix<double>& q,
                               const std::vector<TransMatSym2D>& T);

matrix<double> calc_prob(const std::vector<double>& pi0,
                         const matrix<double>& q,
                         const std::vector<TransMatSym2D>& T);

}

#endif
