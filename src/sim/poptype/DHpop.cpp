#include "DHpop.h"
#include "util_genetics.h"

using namespace std;
using namespace ibd;

double DHpop::p_trans(const Genotype& g1, const Genotype& g2, double r) const
{
	double s = 1.0 - r;
	if (g1 == AA)
    {
		if (g2 == AA) return s;
		if (g2 == BB) return r;
	}
	if (g1 == BB)
	{
		if (g2 == AA) return r;
		if (g2 == BB) return s;
	}
	throw ibd_error("DHpop::p_trans : unknown genotype g1 or g2"); 
	return 0.0;
}

vector<double> DHpop::model(const Genotype& g) const
{
	vector<double> X(npar); 
	if (g == AA) 
		X[0] = -1.0;
	else if (g == BB)
		X[0] =  1.0;
	else
		throw ibd_error("DHpop::model : unknown genotype g");
	return X;
}

Genome DHpop::make_ind(const vector<double>& chr_length) const
{
	Genome parA(chr_length,AA);
	Genome parB(chr_length,BB);
	Genome F1 = parA*parB;
	return DoubledHaploid(F1);
}


vector<Genome> DHpop::make_pop(const vector<double>& chr_length, int nind) const
{
	Genome parA(chr_length,AA);
	Genome parB(chr_length,BB);
	Genome F1 = parA*parB;
	vector<Genome> pop(nind);
	for (int ind=0;ind<nind;ind++)
		pop[ind] = DoubledHaploid(F1);
	return pop;
}

void DHpop::MapQTL(ostream& outp, const ObsGeno& g) const
{
	if (g == AA)
		outp << "a";
	else if (g == BB)
		outp << "b";
	else if (g == U)
		outp << "u";
	else
		throw ibd_error("unknown genotype in function DHpop::MapQTL()");
}

ObsGeno DHpop::MapQTL(istream& inp) const
{
	char c;
	inp >> eatcomment >> c;
	c = ::tolower(c);
	if (c == 'a')
		return AA;
	else if (c == 'b')
		return BB;
	else if (c == 'u' || c == '.' || c == '-')
		return U;
	else 
		throw ibd_error("DHpop::MapQTL(istream& )");
}
