#include "ChromosomePair.h"
#include "util_genetics.h"
#include "Random.h"

using namespace std;

bool ibd::eq(const ChromosomePair& x, const ChromosomePair& y)
{
	if ((x.chr_a == y.chr_a) && (x.chr_b == y.chr_b)) return true;
	if ((x.chr_a == y.chr_b) && (x.chr_b == y.chr_a)) return true;
	return false;
}

void ibd::ChromosomePair::print(ostream& outp) const
{
	chr_a.print(outp);
	chr_b.print(outp);
}

ibd::Chromosome ibd::ChromosomePair::Meiosis() const
{
	const double lambda = 0.01;
	const double length = chr_a.Length();
	Chromosome result;
	Chromosome::const_iterator iter1, iter2, iter1_end, iter2_end;
	iter1 = chr_a.segments.begin();
	iter2 = chr_b.segments.begin();
	iter1_end = chr_a.segments.end();
	iter2_end = chr_b.segments.end();

	if (randuniform() < 0.5)
	{
		swap(iter1,iter2);
		swap(iter1_end,iter2_end);
	}
	double crossover_pos = randexp(lambda);
	while (crossover_pos < length)
	{
		while (iter1->right < crossover_pos)
			result.add_segment(*iter1++);	
		while (iter2->right < crossover_pos)
			iter2++;
		if (iter1->allele != iter2->allele)
			result.add_segment(Segment(iter1->allele,crossover_pos));
		
		crossover_pos += randexp(lambda); 
		swap(iter1,iter2);
		swap(iter1_end,iter2_end);
	}

	for (;iter1 != iter1_end; ++iter1)
		result.add_segment(*iter1);
    return result;
}

ibd::oBinFile& ibd::operator<<(oBinFile& os, const ChromosomePair& x) 
{ os << x.chr_a << x.chr_b; return os;}

ibd::iBinFile& ibd::operator>>(iBinFile& is, ChromosomePair& x)       
{ is >> x.chr_a >> x.chr_b; return is; }



