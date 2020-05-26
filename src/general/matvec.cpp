#include "matvec.h"

using namespace std;

mbl::matrix<double> mbl::identity_matrix(int n)
{
	matrix<double> I(n,n,0.0);
	for (int i=0;i<n;i++)
		I[i][i] = 1.0;
	return I;
}
