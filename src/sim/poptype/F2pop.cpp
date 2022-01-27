#include "F2pop.h"
#include "util_genetics.h"

using namespace std;
using namespace ibd;

double F2pop::p_trans(const Genotype& g1, const Genotype& g2, double r) const
{
	double s = 1.0 - r;
	if (g1 == AA)
	{
		if (g2 == AA) return   s*s;
		if (g2 == AB) return 2*r*s;
		if (g2 == BB) return   r*r;
	}
	if (g1 == AB)
	{
		if (g2 == AA) return   r*s;
		if (g2 == AB) return r*r + s*s;
		if (g2 == BB) return   r*s;
	}
	if (g1 == BB)
	{
		if (g2 == AA) return   r*r;
		if (g2 == AB) return 2*r*s;
		if (g2 == BB) return   s*s;
	}
	throw ibd_error("F2pop::p_trans : unknown genotype g1 or g2");
	return 0.0;
}


vector<double> F2pop::model(const Genotype& g) const
{
	vector<double> X(npar,0.0); // X[0] : add, X[1]: dom 
	if (g == AA) 
		X[0] = -1.0;
	else if (g == BB)
		X[0] =  1.0;
	else if (g == AB)
		X[1] =  1.0;
	else
		throw ibd_error("F2pop::model : unknown genotype g");
	return X;
}

Genome F2pop::make_ind(const vector<double>& chr_length) const
{
	Genome F1(chr_length,AB);
	return F1*F1;
}


vector<Genome> F2pop::make_pop(const vector<double>& chr_length, int nind) const
{
	Genome F1(chr_length,AB);
	vector<Genome> pop(nind);
	for (int ind=0;ind<nind;ind++)
		pop[ind] = F1*F1;
	return pop;
}


void F2pop::MapQTL(ostream& outp, const ObsGeno& g) const
{
	ObsGeno C(AB,BB); // genotype B*
	ObsGeno D(AA,AB); // genotype A*
	if (g == AA)
		outp << "a";
	else if (g == AB)
		outp << "h";
	else if (g == BB)
		outp << "b";
	else if (g == C)
		outp << "c";
	else if (g == D)
		outp << "d";
	else if (g == U)
		outp << "u";
	else
		throw ibd_error("unknown genotype in function BCpop::MapQTL()");
}


ObsGeno F2pop::MapQTL(istream& inp) const
{
	ObsGeno C(AB,BB); // genotype B*
	ObsGeno D(AA,AB); // genotype A*
	char c;
	inp >> eatcomment >> c;
	c = ::tolower(c);
	if (c == 'a')
		return AA;
	else if (c == 'h')
		return AB;
	else if (c == 'b')
		return BB;
	else if (c == 'c')
		return C;
	else if (c == 'd')
		return D;
	else if (c == 'u' || c == '.' || c == '-')
		return U;
	else 
		throw ibd_error("F2pop::MapQTL(istream& )");
}


