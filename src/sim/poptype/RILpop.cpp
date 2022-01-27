#include "RILpop.h"
#include "util_genetics.h"

using namespace std;
using namespace ibd;

namespace 
{

string make_RIL_name(int gen) 
{ 
	string name = "RI";
	char buffer[100]; 
	itoa(gen,buffer,10); 
	name += buffer;
	return name;
}

};

RILpop::RILpop(alleletype p1, alleletype p2, int generation)						
     : AA(Genotype(p1,p1)), 
	   AB(Genotype(p1,p2)), 
	   BB(Genotype(p2,p2)),
       pop_base(ObsGeno(Genotype(p1,p1),
					   Genotype(p1,p2),
					   Genotype(p2,p2)), 1,make_RIL_name(generation)), gen(generation) {}


double RILpop::p_trans(const Genotype& g1, const Genotype& g2, double r) const
{
	double p, rs, pow_0_5_gen, factor;

	rs= r*(1.0-r);
	pow_0_5_gen= pow(0.5,gen);

	if (g1 == AA)
	{
		if (g2 == AA)
		{
            factor = 0.5/(1.0+2.0*r);
                 p =   factor - pow_0_5_gen
                     - factor * pow(0.5-r,gen)
                     + 0.25   * pow(0.5-rs,gen-1);
        }
		else if (g2 == BB)
		{
            factor = 1.0+2.0*r;
                 p =     r/factor - pow_0_5_gen
                    + 0.5/factor * pow(0.5-r,gen)
                    + 0.25       * pow(0.5-rs,gen-1);
        }
		else if (g2 == AB)
		{
			     p =                pow_0_5_gen
                    - 0.5        * pow(0.5-rs,gen-1);
        }
		else throw ibd_error("Unknown genotype g2 ");
		p /= 0.5 - pow_0_5_gen;
		return p;
	}
    if (g1 == BB)
	{
        if (g2 == AA)
		{
			factor= 1.0+2.0*r;
                 p=     r/factor - pow_0_5_gen
                    + 0.5/factor * pow(0.5-r,gen)
                    + 0.25       * pow(0.5-rs,gen-1);
        } 
		else if (g2 == BB)
		{
            factor= 0.5/(1.0+2.0*r);
                 p=   factor - pow_0_5_gen
                    - factor * pow(0.5-r,gen)
                    + 0.25   * pow(0.5-rs,gen-1);
		}
		else if (g2 == AB)
        {
                 p=            pow_0_5_gen
                    - 0.5    * pow(0.5-rs,gen-1);
        }
		else throw ibd_error("Unknown genotype g2 ");
		p/= 0.5 - pow_0_5_gen;
		return p;
	}
	if (g1 == AB)
	{
		if (g2 == AA || g2 == BB)
			p =   0.5 - 0.5 * pow(1.0-2.0*rs,gen-1);
		else if (g2 == AB)
            p =               pow(1.0-2.0*rs,gen-1);
		else throw ibd_error("Unknown genotype g2 ");
		return p;
	}
	throw ibd_error("Unknown genotype g1");
	return 0.0;
} 




vector<double> RILpop::model(const Genotype& g) const
{
	vector<double> X(npar,0.0); // no dominance !
	if (g == AA) 
		X[0] = -1.0;
	else if (g == BB)
		X[0] =  1.0;
	else if (g == AB)
		X[0] =  0.0;
	else
		throw ibd_error("RILpop::model : unknown genotype g");
	return X;
}

Genome RILpop::make_ind(const vector<double>& chr_length) const
{
	Genome F1(chr_length,AB); // F1 population
	Genome F = F1;
	for (int i=1;i<gen;i++)
		F = F*F;
	return F;
}

vector<Genome> RILpop::make_pop(const vector<double>& chr_length, int nind) const
{
	Genome F1(chr_length,AB); // F1 population
	vector<Genome> pop(nind);
	for (int ind=0;ind<nind;ind++)
	{
		Genome F = F1;
		for (int i=1;i<gen;i++)
			F = F*F;
		pop[ind] = F;
	}
	return pop;
}

void RILpop::MapQTL(ostream& outp, const ObsGeno& g) const
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
		throw ibd_error("unknown genotype in function RILpop::MapQTL()");
}


ObsGeno RILpop::MapQTL(istream& inp) const
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
		throw ibd_error("RILpop::MapQTL(istream& )");
}


