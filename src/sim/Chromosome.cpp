#include <math.h>
#include "Chromosome.h"
#include "util_genetics.h"

using namespace std;

void ibd::Chromosome::print(ostream& outp) const
{	
	int nr_segments = segments.size();
	for (int i=0;i<nr_segments;i++)
	{
		const int nr = segments[i].allele;
		outp << nr << ": ";
		double left = (i==0) ? 0.0 : segments[i-1].right;
		double right = segments[i].right;
		outp << left << "-" << right << ";"; 
    }
	outp << endl;
}

ibd::alleletype ibd::Chromosome::GetGenotype(double cM) const
{
	const_iterator iter;
	
	for (iter = segments.begin(); iter != segments.end(); ++iter)
		if (cM <= iter->right)
			return iter->allele;
	
	throw ibd_error("Chromosome::GetGenotype range error");
}

bool ibd::eq(const Chromosome& chr1, const Chromosome& chr2)
{
	if (chr1.GetNrSegments() != chr2.GetNrSegments()) 
		return false;

	Chromosome::const_iterator iter1, iter2;
	iter1 = chr1.segments.begin();
	iter2 = chr2.segments.begin();
	for (; iter1 != chr1.segments.end(); ++iter1, ++iter2)
	{
		if ((iter1->allele != iter2->allele) || fabs(iter1->right - iter2->right) > 1.0e-10)
			return false;
	}
	return true;
}
	
ibd::oBinFile& ibd::operator<<(oBinFile& os, const Chromosome& x) 
{ os << x.segments; return os; }

ibd::iBinFile& ibd::operator>>(iBinFile& is, Chromosome& x)       
{ is >> x.segments; return is; }
