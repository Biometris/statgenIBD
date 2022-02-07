#include <Rcpp.h>
#include <algorithm>
#include "Genome.h"
#include "util_genetics.h"

using namespace std;

namespace 
{

// 19 oktober 2000: function object (Stroustrup p. 514),
// wordt gebruikt in Genome::GetGenotype.
class Position
{
public:
	Position(double p) : pos(p) {}
	bool operator()(const ibd::Segment& x) { return (pos <= x.right); }
private:
	double pos;
};

}


ibd::Genome::Genome(const vector<double>& chr_length, alleletype a, alleletype b) 
{
	vector<double>::const_iterator iter;
	for (iter=chr_length.begin(); iter!=chr_length.end(); ++iter)
		chr_pair.push_back(ChromosomePair(*iter,a,b));
}

ibd::Genome::Genome(int nr_chr, double length, alleletype a, alleletype b)
{
	for (int i=0; i<nr_chr; ++i)
		chr_pair.push_back(ChromosomePair(length,a,b));
}

ibd::Genome::Genome(const vector<double>& chr_length, const Genotype& g) 
{
	alleletype a = g.First();
	alleletype b = g.Second();
	vector<double>::const_iterator iter;
	for (iter=chr_length.begin(); iter!=chr_length.end(); ++iter)
		chr_pair.push_back(ChromosomePair(*iter,a,b));
}

ibd::Genome::Genome(int nr_chr, double length, const Genotype& g)
{
	alleletype a = g.First();
	alleletype b = g.Second();
	for (int i=0; i<nr_chr; ++i)
		chr_pair.push_back(ChromosomePair(length,a,b));
}

bool ibd::Genome::IsHomozygote() const
{
	const_iterator iter = chr_pair.begin();
	for (; iter != chr_pair.end(); ++iter)
		if (!iter->IsHomozygote())
			return false;
	return true;
}

ibd::Genotype ibd::Genome::GetGenotype(int chr_nr, double cM) const
{
	Rcpp::Rcout << chr_nr << std::endl;
	if (chr_nr < 0 || (unsigned)chr_nr >= chr_pair.size())
		throw ibd_error("Genome::GetGenotype");
	return chr_pair[chr_nr].GetGenotype(cM);
}

ibd::Genotype ibd::Genome::GetGenotype(const Locus& loc) const
{
	int chr_nr = std::stoi(loc.GetChr());
	if (chr_nr < 0 || (unsigned)chr_nr >= chr_pair.size())
		throw ibd_error("Genome::GetGenotype");
	double cM = loc.GetPosition();
	return chr_pair[chr_nr].GetGenotype(cM);
}

void ibd::Genome::print(ostream& outp) const
{
	int nr_chrom = chr_pair.size();
	for (int i=0;i<nr_chrom;i++)
	{
		cout << "Chromosome " << i << endl;
		chr_pair[i].print(outp);
		outp << endl;
	}
}

ostream& ibd::operator<<(ostream& outp, const Genome& g)
{
	g.print(outp);
	return outp;
}

bool ibd::eq(const Genome& a, const Genome& b)
{
	if (a.GetNrChromosomePairs() != b.GetNrChromosomePairs() )
		return false;

	Genome::const_iterator iter1, iter2, iter1_end;
	iter1 = a.chr_pair.begin();
	iter2 = b.chr_pair.begin();
	iter1_end = a.chr_pair.end();
	for (; iter1 != iter1_end; ++iter1, ++iter2)
		if (*iter1 != *iter2)
			return false;
	return true;
}

ibd::Genome ibd::Cross(const Genome& a, const Genome& b)
{
	Genome c;
	unsigned int nr_chrom = a.chr_pair.size();
	if (nr_chrom != b.chr_pair.size())
		throw ibd_error("Cross: unequal number of chromosomes");
	c.chr_pair.resize(nr_chrom);
	Genome::const_iterator iter_a, iter_a_end, iter_b;
	Genome::iterator iter_c;
	iter_a = a.chr_pair.begin();
	iter_a_end = a.chr_pair.end();
	iter_b = b.chr_pair.begin();
	iter_c = c.chr_pair.begin();
	for (;iter_a != iter_a_end; ++iter_a, ++iter_b, ++iter_c)
	{
		iter_c->chr_a = iter_a->Meiosis(); 
		iter_c->chr_b = iter_b->Meiosis(); 
		if (fabs(iter_c->chr_a.Length()-iter_c->chr_b.Length()) > 1.0e-10)
			throw ibd_error("Cross: unequal lengths of chromosomes");
	}
	return c;
}


ibd::Genome ibd::DoubledHaploid(const Genome& a)
{
	Genome c;
	int nr_chrom = a.chr_pair.size();
	c.chr_pair.resize(nr_chrom);
	Genome::const_iterator iter_a, iter_a_end;
	Genome::iterator iter_c;
	iter_a = a.chr_pair.begin();
	iter_a_end = a.chr_pair.end();
	iter_c = c.chr_pair.begin();
	for (;iter_a != iter_a_end; ++iter_a, ++iter_c)
		iter_c->chr_a = iter_c->chr_b = iter_a->Meiosis();
	return c;
}

int ibd::Genome::GetNrSegments() const
{
	int nr_chr = chr_pair.size();
	int sum = 0;
	for (int i=0;i<nr_chr;i++)
		sum += chr_pair[i].GetNrSegments();
	return sum;
}

vector<ibd::Genotype> ibd::Genome::GetGenotype(const LinkageMap& linkage_map) const
{
	LinkageMap::const_iterator iter1;
	Chromosome::const_iterator iter2_a, iter2_b,a_end,b_end;
	vector<Genotype>::iterator iter3;
	vector<Genotype> result(linkage_map.size());
	int nchr = chr_pair.size();
	int prev_chr_nr = -1, chr_nr;
	iter3 = result.begin();
	for (iter1 = linkage_map.begin(); iter1 != linkage_map.end(); iter1++)
	{
		chr_nr = std::stoi(iter1->GetChr());
		if (prev_chr_nr != chr_nr)
		{
			Rcpp::Rcout << chr_nr << std::endl;
			if (chr_nr < 0 || chr_nr >= nchr)
				throw ibd_error("Genome::GetGenotype ");
			const ChromosomePair& cur = chr_pair[chr_nr];
			iter2_a = cur.chr_a.segments.begin();
			iter2_b = cur.chr_b.segments.begin();
			a_end = cur.chr_a.segments.end();
			b_end = cur.chr_b.segments.end();
		}
		Position pos(iter1->GetPosition());
	    iter2_a = find_if(iter2_a,a_end,pos);
		iter2_b = find_if(iter2_b,b_end,pos);
		
		if (iter2_a == a_end || iter2_b == b_end)
			throw ibd_error("Genome::GetGenotype");
		*iter3++ = Genotype(iter2_a->allele, iter2_b->allele);
		prev_chr_nr = chr_nr;
	}
	return result;
}


vector<ibd::alleletype> ibd::Genome::GetGamete(const LinkageMap& linkage_map, Parent parent) const
{
	LinkageMap::const_iterator iter1;
	Chromosome::const_iterator iter2, end;
	vector<alleletype> result(linkage_map.size());
	vector<alleletype>::iterator iter3 = result.begin();
	int nchr = chr_pair.size();
	int prev_chr_nr = -1, chr_nr;
	for (iter1 = linkage_map.begin(); iter1 != linkage_map.end(); iter1++)
	{
		chr_nr = std::stoi(iter1->GetChr());
		if (prev_chr_nr != chr_nr)
		{
			if (chr_nr < 0 || chr_nr >= nchr)
				throw ibd_error("Genome::GetGamete ");
			const ChromosomePair& cur = chr_pair[chr_nr];
			const Chromosome& chr = (parent == Father) ? cur.chr_a : cur.chr_b;
			iter2 = chr.segments.begin();
			end   = chr.segments.end();
		}
		Position pos(iter1->GetPosition());
	    iter2 = find_if(iter2,end,pos);
		if (iter2 == end)
			throw ibd_error("Genome::GetGamete");
		*iter3++ = iter2->allele;
		prev_chr_nr = chr_nr;
	}
	return result;
}

std::multimap<ibd::alleletype,double> ibd::Genome::GetLengthFragments(int chr_nr, Parent parent) const
{
	const ChromosomePair& cur = chr_pair[chr_nr];
	const Chromosome& chr = (parent == Father) ? cur.chr_a : cur.chr_b;
	double left = 0.0;
	std::multimap<ibd::alleletype,double> fragments;
	for (Chromosome::const_iterator it = chr.segments.begin(); it != chr.segments.end(); it++)
	{
		fragments.insert(pair<ibd::alleletype,double>(it->allele, it->right - left));
		left = it->right;
	}
	return fragments;
}

std::vector<ibd::Segment> ibd::Genome::GetSegments(int chr_nr, Parent parent) const
{
	const ChromosomePair& cur = chr_pair[chr_nr];
	const Chromosome& chr = (parent == Father) ? cur.chr_a : cur.chr_b;
	return chr.GetSegments();	
}

ibd::Genome ibd::Selfing(const Genome& x, int Ngen)
{
	Genome cur_ind = x;
	for (int gen=0;gen<Ngen;gen++)
		cur_ind = cur_ind*cur_ind;
	return cur_ind;
}

ibd::oBinFile& ibd::operator<<(oBinFile& os, const Genome& x) { os << x.chr_pair; return os; }
ibd::iBinFile& ibd::operator>>(iBinFile& is, Genome& x)       { is >> x.chr_pair; return is; } 


