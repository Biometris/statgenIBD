#include "BCpop.h"
#include "util_genetics.h"

using namespace std;
using namespace ibd;

double BCpop::p_trans(const Genotype& g1, const Genotype& g2, double r) const
{
	double s = 1.0 - r;
	if (g1 == AA)
    {
		if (g2 == AA) return s;
		if (g2 == AB) return r;
	}
	if (g1 == AB)
	{
		if (g2 == AA) return r;
		if (g2 == AB) return s;
	}
	throw ibd_error("BCpop::p_trans : unknown genotype g1 or g2"); 
	return 0.0;
}

vector<double> BCpop::model(const Genotype& g) const
{
	vector<double> X(npar); 
	if (g == AA) 
		X[0] = -0.5;
	else if (g == AB)
		X[0] =  0.5;
	else
		throw ibd_error("BCpop::model : unknown genotype g");
	return X;
}

Genome BCpop::make_ind(const vector<double>& chr_length) const
{
	Genome F1(chr_length,AB);
	Genome parA(chr_length,AA);
	return F1*parA;
}


vector<Genome> BCpop::make_pop(const vector<double>& chr_length, int nind) const
{
	Genome F1(chr_length,AB);
	Genome parA(chr_length,AA);
	vector<Genome> pop(nind);
	for (int ind=0;ind<nind;ind++)
		pop[ind] = F1*parA;
	return pop;
}

void BCpop::MapQTL(ostream& outp, const ObsGeno& g) const
{
	if (g == AA)
		outp << "a";
	else if (g == AB)
		outp << "h";
	else if (g == U)
		outp << "u";
	else
		throw ibd_error("unknown genotype in function BCpop::MapQTL()");
}

ObsGeno BCpop::MapQTL(istream& inp) const
{
	char c;
	inp >> eatcomment >> c;
	c = ::tolower(c);
	if (c == 'a')
		return AA;
	else if (c == 'h')
		return AB;
	else if (c == 'u' || c == '.' || c == '-')
		return U;
	else 
		throw ibd_error("BCpop::MapQTL(istream& )");
}
