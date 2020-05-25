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