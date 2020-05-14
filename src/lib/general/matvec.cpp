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

/*
void mbl::resize_matrix_example()
{
	matrix<int> G(2,3,5);          // 2 x 3 matrix with initial values 5
	cout << setw(3) << G << endl;  
	G.resize(3,4,1);	           // resize to 3x4 matrix, new elements value 1 
	cout << setw(3) << G << endl;
	G.resize(4,5,2);	           // resize to 4x5 matrix, new elements value 2
	cout << setw(3) << G << endl;
	G.resize(2,7);				   // resize to 2x7 matrix.
	cout << setw(3) << G << endl;
}
*/
