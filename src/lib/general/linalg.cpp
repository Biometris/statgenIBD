/***********************************************************
*														   *
*   Martin Boer											   *
*   Biometris Wageningen-UR & NWO, grant no. 925-01-001    *
*   July 12, 2002										   *
*														   *
************************************************************/

// #include <cmath>
// #include <algorithm>
// #include "linalg.h"
// #include "mblexcept.h"
//
// using namespace std;
//
// namespace aux
// {
// using namespace mbl;
//
// /*
// double ** conversion_matrix(const matrix<double>& A)
// {
// 	int i,j;
// 	int nr = A.NrRows();
// 	int nc = A.NrCols();
// 	double **ptr = new double *[nr];
// 	for (i=0;i<nr;i++)
// 		*(ptr+i) = new double[nc];
// 	for (i=0;i<nr;i++)
// 		for (j=0;j<nr;j++)
// 			ptr[i][j] = A[i][j]; // kan nog wat sneller!
// 	return ptr;
// }
//
// void freematrix(double **ptr, int nr)
// {
// 	for (int i=0;i<nr;i++)
// 	    delete[] *(ptr+i);
// 	delete[] ptr;
// }
// */
//
// // bool LUdecomposition(double **A, int n, vector<int>& index, int& sign)
// // {
// // 	int imax,i,j,k;
// // 	double big, dum, sum;
// // 	vector<double> vv(n);
// // 	sign = 1;
// //     for (i=0; i<n; i++)
// // 	{
// // 		big = 0.0;
// // 		for (j=0; j<n; j++)
// // 			if (fabs(A[i][j]) > big) big = fabs(A[i][j]);
// // 		if (big == 0)
// // 			return false;
// // 		vv[i] = 1.0/big;
// // 	}
// // 	for (j=0; j<n; j++)
// // 	{
// // 		for (i=0;i<j;i++)
// // 		{
// // 			sum = A[i][j];
// // 			for (k=0;k<i;k++)
// // 				sum -= A[i][k]*A[k][j];
// // 			A[i][j] = sum;
// // 		}
// // 		big = 0.0;
// // 		for (i=j;i<n;i++)
// // 		{
// // 			sum = A[i][j];
// // 			for (k=0;k<j;k++)
// // 				sum -= A[i][k]*A[k][j];
// // 			A[i][j] = sum;
// // 			dum = vv[i]*fabs(sum);
// // 			if (dum >= big)
// // 			{
// // 				big = dum;
// // 				imax = i;
// // 			}
// // 		}
// // 		if (j != imax)
// // 		{
// // 			swap(A[imax],A[j]);
// // 			sign *= -1;
// // 			vv[imax] = vv[j];
// // 		}
// // 		index[j] = imax;
// // 		if (A[j][j] == 0)
// // 			A[j][j] = 1.0e-20; // return false
// // 		for (i=j+1; i<n; i++)
// // 			A[i][j] /= A[j][j];
// // 	}
// // 	return true;
// // }
//
// // void LUbacksubst(double **LU, int n, vector<int>& index, vector<double>& b)
// // {
// // 	int j,ip,ii,i;
// // 	double sum;
// // 	ii = -1;
// // 	for (i=0;i<n;i++)
// // 	{
// // 		ip = index[i];
// // 		sum = b[ip];
// // 		b[ip] = b[i];
// // 		if (ii != -1)
// // 			for (j=ii;j<i;j++)
// // 				sum -= LU[i][j]*b[j];
// // 		else if (sum != 0)
// // 			ii = i;
// // 		b[i] = sum;
// // 	}
// // 	for (i=n-1;i>=0;i--)
// // 	{
// // 		sum = b[i];
// // 		for (j=i+1;j<n;j++)
// // 			sum -= LU[i][j]*b[j];
// // 		b[i] = sum/LU[i][i];
// // 	}
// // }
//
// }
// //
// // bool mbl::LUdecomposition(matrix<double>& A, vector<int>& index, int& sign)
// // {
// // 	int imax,i,j,k;
// // 	int n = A.NrCols();
// // 	double big, dum, sum;
// // 	vector<double> vv(n);
// // 	sign = 1;
// //     for (i=0; i<n; i++)
// // 	{
// // 		big = 0.0;
// // 		for (j=0; j<n; j++)
// // 			if (fabs(A[i][j]) > big) big = fabs(A[i][j]);
// // 		if (big == 0)
// // 			return false;
// // 		vv[i] = 1.0/big;
// // 	}
// // 	for (j=0; j<n; j++)
// // 	{
// // 		for (i=0;i<j;i++)
// // 		{
// // 			sum = A[i][j];
// // 			for (k=0;k<i;k++)
// // 				sum -= A[i][k]*A[k][j];
// // 			A[i][j] = sum;
// // 		}
// // 		big = 0.0;
// // 		for (i=j;i<n;i++)
// // 		{
// // 			sum = A[i][j];
// // 			for (k=0;k<j;k++)
// // 				sum -= A[i][k]*A[k][j];
// // 			A[i][j] = sum;
// // 			dum = vv[i]*fabs(sum);
// // 			if (dum >= big)
// // 			{
// // 				big = dum;
// // 				imax = i;
// // 			}
// // 		}
// // 		if (j != imax)
// // 		{
// // 			A[imax].swap(A[j]);
// // 			//for (k=0;k<n;k++) // 13 september 2000: verwisselen van rijen sneller ?!
// // 			//	swap(A[imax][k],A[j][k]);
// // 			sign *= -1;
// // 			vv[imax] = vv[j];
// // 		}
// // 		index[j] = imax;
// // 		if (A[j][j] == 0)
// // 			A[j][j] = 1.0e-20; // return false
// // 		for (i=j+1; i<n; i++)
// // 			A[i][j] /= A[j][j];
// // 	}
// // 	return true;
// // }
// //
// // void mbl::LUbacksubst(matrix<double>& LU, vector<int>& index, vector<double>& b)
// // {
// // 	int j,ip,ii,i;
// // 	double sum;
// // 	int n = LU.NrCols();
// // 	ii = -1;
// // 	for (i=0;i<n;i++)
// // 	{
// // 		ip = index[i];
// // 		sum = b[ip];
// // 		b[ip] = b[i];
// // 		if (ii != -1)
// // 			for (j=ii;j<i;j++)
// // 				sum -= LU[i][j]*b[j];
// // 		else if (sum != 0)
// // 			ii = i;
// // 		b[i] = sum;
// // 	}
// // 	for (i=n-1;i>=0;i--)
// // 	{
// // 		sum = b[i];
// // 		for (j=i+1;j<n;j++)
// // 			sum -= LU[i][j]*b[j];
// // 		b[i] = sum/LU[i][i];
// // 	}
// // }
//
// // bool mbl::Solve(const matrix<double>& A, vector<double>& x, const vector<double>& b)
// // {
// // 	int sign,n = A.NrCols();
// // 	vector<int> index(n);
// // 	matrix<double> LU = A;
// //     //double **LU = conversion_matrix(A);
// // 	x = b;
// // 	if (!LUdecomposition(LU,index,sign))
// // 		return false;
// // 	LUbacksubst(LU,index,x);
// //
// // 	//if (!LUdecomposition(LU, n, index, sign))
// // 	//	return false;
// // 	//LUbacksubst(LU, n, index, x);
// //
// // 	//freematrix(LU,n);
// // 	return true;
// // }
//
//
// // mbl::matrix<double> mbl::Inverse(const matrix<double>& A)
// // {
// // 	int sign;
// // 	int n = A.NrCols();
// // 	vector<int> index(n);
// // 	matrix<double> Ainv(n,n);
// // 	double **LU = aux::conversion_matrix(A);
// //     //matrix<double> LU = A;
// // 	if (!aux::LUdecomposition(LU, n, index, sign))
// // 		throw mblib_error("Inverse of singular matrix");
// // 	for (int j=0;j<n;j++)
// // 	{
// // 		vector<double> col(n,0.0);
// // 		col[j] = 1.0;
// // 		aux::LUbacksubst(LU, n, index, col);
// // 		for (int i=0;i<n;i++)
// // 			Ainv[i][j] = col[i];
// // 	}
// // 	aux::freematrix(LU,n);
// // 	return Ainv;
// // }
//
// double mbl::Determinant(const matrix<double>& A)
// {
// 	int sign;
// 	int n = A.NrCols();
// 	double result;
// 	vector<int> index(n);
// 	double **LU = aux::conversion_matrix(A);
// 	//matrix<double> LU = A;
// 	if (!aux::LUdecomposition(LU, n, index,sign))
// 		return 0.0;
// 	result = sign;
// 	for (int j=0;j<n;j++)
// 		result *= LU[j][j];
// 	aux::freematrix(LU,n);
// 	return result;
// }
//
// // inline double abs_complex(complex<double> x) { return abs(x); }
// //
// // bool mbl::LUdecomposition(matrix< complex<double> >& A, vector<int>& index, int& sign)
// // {
// // 	int imax,i,j,k;
// // 	int n = A.NrCols();
// // 	double big, dum;
// // 	complex<double> sum;
// // 	vector<double> vv(n);
// // 	sign = 1;
// //     for (i=0; i<n; i++)
// // 	{
// // 		big = 0.0;
// // 		for (j=0; j<n; j++)
// // 			if (abs_complex(A[i][j]) > big) big = abs_complex(A[i][j]);
// // 		if (big == 0)
// // 			return false;
// // 		vv[i] = 1.0/big;
// // 	}
// // 	for (j=0; j<n; j++)
// // 	{
// // 		for (i=0;i<j;i++)
// // 		{
// // 			sum = A[i][j];
// // 			for (k=0;k<i;k++)
// // 				sum -= A[i][k]*A[k][j];
// // 			A[i][j] = sum;
// // 		}
// // 		big = 0.0;
// // 		for (i=j;i<n;i++)
// // 		{
// // 			sum = A[i][j];
// // 			for (k=0;k<j;k++)
// // 				sum -= A[i][k]*A[k][j];
// // 			A[i][j] = sum;
// // 			dum = vv[i]*abs_complex(sum);
// // 			if (dum >= big)
// // 			{
// // 				big = dum;
// // 				imax = i;
// // 			}
// // 		}
// // 		if (j != imax)
// // 		{
// // 			for (k=0;k<n;k++) // 13 september 2000: verwisselen van rijen sneller ?!
// // 				std::swap(A[imax][k],A[j][k]);
// // 			sign *= -1;
// // 			vv[imax] = vv[j];
// // 		}
// // 		index[j] = imax;
// // 		if (A[j][j] == 0.0)
// // 			A[j][j] = 1.0e-20; // return false
// // 		for (i=j+1; i<n; i++)
// // 			A[i][j] /= A[j][j];
// // 	}
// // 	return true;
// // }
// //
// // void mbl::LUbacksubst(matrix< complex<double> >& LU, vector<int>& index, vector< complex<double> >& b)
// // {
// // 	int j,ip,ii,i;
// // 	complex<double> sum;
// // 	int n = LU.NrCols();
// // 	ii = -1;
// // 	for (i=0;i<n;i++)
// // 	{
// // 		ip = index[i];
// // 		sum = b[ip];
// // 		b[ip] = b[i];
// // 		if (ii != -1)
// // 			for (j=ii;j<i;j++)
// // 				sum -= LU[i][j]*b[j];
// // 		else if (sum != 0.0)
// // 			ii = i;
// // 		b[i] = sum;
// // 	}
// // 	for (i=n-1;i>=0;i--)
// // 	{
// // 		sum = b[i];
// // 		for (j=i+1;j<n;j++)
// // 			sum -= LU[i][j]*b[j];
// // 		b[i] = sum/LU[i][i];
// // 	}
// // }
// //
// // bool mbl::Solve(const matrix< complex<double> >& A, vector< complex<double> >& x, const vector< complex<double> >& b)
// // {
// // 	int sign;
// // 	int n = A.NrCols();
// // 	vector<int> index(n);
// // 	matrix< complex<double> > LU = A;
// // 	x = b;
// // 	if (!LUdecomposition(LU,index,sign))
// // 		return false;
// // 	LUbacksubst(LU,index,x);
// // 	return true;
// // }
// //
// //
// // mbl::matrix< complex<double> > mbl::Inverse(const matrix< complex<double> >& A)
// // {
// // 	int sign;
// // 	int n = A.NrCols();
// // 	vector<int> index(n);
// // 	matrix< complex<double> > Ainv(n,n);
// //     matrix< complex<double> > LU = A;
// // 	if (!LUdecomposition(LU, index, sign))
// // 		throw mblib_error("Inverse of singular matrix");
// // 	for (int j=0;j<n;j++)
// // 	{
// // 		vector< complex<double> > col(n,0.0);
// // 		col[j] = 1.0;
// // 		LUbacksubst(LU, index, col);
// // 		for (int i=0;i<n;i++)
// // 			Ainv[i][j] = col[i];
// // 	}
// // 	return Ainv;
// // }
// //
// // complex<double> mbl::Determinant(const matrix< complex<double> >& A)
// // {
// // 	int sign;
// // 	int n = A.NrCols();
// // 	complex<double> result;
// // 	vector<int> index(n);
// // 	matrix< complex<double> > LU = A;
// // 	if (!LUdecomposition(LU, index,sign))
// // 		return 0.0;
// // 	result = sign;
// // 	for (int j=0;j<n;j++)
// // 		result *= LU[j][j];
// // 	return result;
// // }


