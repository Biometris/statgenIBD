#ifndef EIGEN_HEADER
#define EIGEN_HEADER

#include <complex>
#include <cmath>
#include <limits>
#include <algorithm>
#include <iostream>
#include <iomanip>
//#include "util.h"
#include "matvec.h"

namespace mbl
{
	// calculation of eigenvalues for a SYMMETRIC matrix A
	std::vector<double> eigenvals(matrix<double>& A);
    void Jacobi(matrix<double> a, std::vector<double>& d, matrix<double>& v);
	matrix<double> square_root(const matrix<double>& M);

	//typedef std::vector< std::complex<double> > VectorComplex;
	//typedef matrix<double> MatrixDouble;

	// calculation of eigenvalues for a NON-SYMMETRIC matrix A
	std::vector< std::complex<double> > eigenvalues(matrix<double> A);

	// test functions
	void test_calc_eigenvals_sym();
	void test_calc_eigenvals_nonsym();

	//void eigen_vals_vects(mbl::matrix<double>& U, std::vector<double>& lambda, const mbl::matrix<double>& B);

	// real eigenvals and eigenvalues
	class EigenReal
	{
	public:
		EigenReal() {}
		EigenReal(double value,const std::vector<double>& vec) : eigenval(value),eigenvec(vec) {}
		double EigenValue() const { return eigenval; }
		std::vector<double> Eigenvector() const { return eigenvec; }
	private:
		friend bool compare_EigenReal(const EigenReal& , const EigenReal&);
		double eigenval;
		std::vector<double> eigenvec;
	};

	inline bool compare_EigenReal(const EigenReal& x, const EigenReal& y) { return (x.eigenval > y.eigenval); }

	std::vector<EigenReal> CalcEigen(const mbl::matrix<double>& SymmetricMatrix);

}

#endif

