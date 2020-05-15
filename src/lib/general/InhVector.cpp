#include <iostream>
#include <iomanip>
#include <vector>

#include "misc.h"
#include "convert.h"
#include "InhVector.h"

using namespace std;
using namespace mbl;

InhVector::InhVector(int len,unsigned int init) : x(init), max(mbl::pow2(len)),n(len)
{
	if (n > 15)
		throw mblib_error("Inheritance vector too long: " + stringify(n));
}


bool InhVector::next_indic()
{
	bool bit = (x & 01);
	x >>= 1;
	max >>= 1;
	n--;
	return bit;
}

void InhVector::print(std::ostream& outp) const
{
	unsigned int z = 1;
	for (unsigned int i=0;i<n;i++)
	{
		bool r = (x & z);
		outp << r;
		z <<= 1;
	}
}

ostream& mbl::operator<<(ostream& outp, const InhVector& h)
{
	int nr_whitespaces = outp.width() - h.length();
	outp.width(0);
	if (nr_whitespaces > 0)
		outp << string(nr_whitespaces,' ');
	h.print(outp);
	return outp;
}

// int mbl::testInhVector()
// {
// 	const int n = 4;
// 	const unsigned int start = 3;
// 	int N = pow2(n);
// 	vector<int> c(N,-1);
// 	int cnt=0;
// 	for (InhVector y(n,start);!y.end();y++)
// 	{
// 		c[y] = cnt++;
// 		unsigned int u = y;
// 		cout << setw(5) << u << setw(8) << y;
// 		InhVector z = y;
// 		cout << "   ";
// 		for (int i=0;i<n;i++)
// 			cout << setw(2) << z.next_indic();
// 		unsigned int u2 = z;
// 		cout << setw(8) << u2 << "  z:" << z << ":" << endl;
// 	}
//
// 	cout << endl;
// 	for (int i=0;i<N;i++)
// 		cout << setw(3) << i << setw(6) << c[i] << endl;
//
// 	return 0;
// }

