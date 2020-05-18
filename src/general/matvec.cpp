#include "misc.h"
#include "matvec.h"

using namespace std;

mbl::matrix<double> mbl::identity_matrix(int n)
{
	matrix<double> I(n,n,0.0);
	for (int i=0;i<n;i++)
		I[i][i] = 1.0;
	return I;
}

void mbl::write_bin(ostream& outp, const matrix<double>& x)
{
	const int ncol = x.NrCols();
	const int nrow = x.NrRows();
	for (int r=0;r<nrow;r++)
		for (int c=0;c<ncol;c++)
			write_bin(outp,x[r][c]);
}
