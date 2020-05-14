/***********************************************************
*														   *
*   Martin Boer											   *
*   Biometris Wageningen-UR & NWO, grant no. 925-01-001    *
*   July 12, 2002										   *
*														   *
*	original C-code from Numerical Recipes in C			   *
*														   *
************************************************************/

#ifndef LINALG_HEADER
#define LINALG_HEADER

//#pragma warning(disable:4786)   // disable C4786 warning

#include <vector>
#include <complex>
#include "matvec.h"

namespace mbl
{

bool LUdecomposition(matrix<double>&, std::vector<int>&, int&);
void LUbacksubst(matrix<double>&, std::vector<int>&, std::vector<double>& );
bool Solve(const matrix<double>&, std::vector<double>&, const std::vector<double>&);
matrix<double> Inverse(const matrix<double>& A);
double Determinant(const matrix<double>& A);

bool LUdecomposition(matrix< std::complex<double> >&, std::vector<int>&, int&);
void LUbacksubst(matrix< std::complex<double> >&, std::vector<int>&, std::vector< std::complex<double> >& );
bool Solve(const matrix< std::complex<double> >&,
		   std::vector< std::complex<double> >&, const std::vector< std::complex<double> >&);
matrix<std::complex<double> > Inverse(const matrix< std::complex<double> >& A);
std::complex<double> Determinant(const matrix<std::complex<double> >& A);

}

#endif
