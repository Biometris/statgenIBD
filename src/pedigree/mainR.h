/*!
 \file
 \brief  header file with some defined constants
 \author Martin Boer, Biometris
 \date   2006-2011
 */
#ifndef HEADER_MAIN_PEDIGREE
#define HEADER_MAIN_PEDIGREE

#include <string>

// #include "util.h"
#include "matvec.h"
#include "Loc.h"

const std::string version = "2.75";
const std::string date    = "april 13, 2020";

//! indicates whether scores are given in parentheses
const bool scores_in_parentheses = true;


int main_forR(mbl::matrix3D<double>& Z,
              std::vector<std::string>& parents,
              std::vector<std::string>& offspring,
              LinkageMap& positions,
              const std::string& poptype,
              const std::string& locfile,
              const std::string& mapfile,
             // const std::string& output,
              const std::string& eval_pos,
              const double& max_step_size);

#endif

