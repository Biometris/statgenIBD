#include "GenoValue.h"
#include "Chromosome.h"

using namespace ibd;
using namespace std;

Phi::Phi(const map<Locus,vector<double> >& QTLs, const map<string,string>& inbfnds,
		 const matrix<double>& epi_addxadd, double mu) :  mu_(mu), epi_axa(epi_addxadd)
{
	nQTL = QTLs.size();
	nFND = inbfnds.size();
	typedef map<Locus,vector<double> >::const_iterator IterQTLs;
	for (IterQTLs iter=QTLs.begin();iter!=QTLs.end();++iter)
	{
		QTLmap.push_back(iter->first);
		vector<double> par = iter->second;
		add.push_back(par[0]);
		dom.push_back(par[1]);
	}
	pos_allele = matrix<bool>(nFND,nQTL);
	int k=0;
	typedef map<string,string>::const_iterator IterFnds;
	for (IterFnds iter=inbfnds.begin();iter!=inbfnds.end();++iter)
	{
		string str = iter->second;
		for (int q=0;q<nQTL;q++)
		{
			pos_allele[k][q] = (str[q] == '+');
		}
		k++;
	}
}
		
int Phi::add_indicator(const vector<Genotype>& geno, int q) const
{
	Genotype g = geno[q];
	bool g1 = pos_allele[g.First()][q];
	bool g2 = pos_allele[g.Second()][q];
	if ( g1 &&  g2) return  1;  // both alleles positive
	if (!g1 && !g2) return -1;  // both alleles negative
	return 0;					// heterozygous
}

int Phi::dom_indicator(const vector<Genotype>& geno, int q) const
{
	Genotype g = geno[q];
	bool g1 = pos_allele[g.First()][q];
	bool g2 = pos_allele[g.Second()][q];
	if ( g1 &&  g2) return 0;  // both alleles positive
	if (!g1 && !g2) return 0;  // both alleles negative
	return 1;				   // heterozygous
}

double Phi::operator()(const Genome& genome) const 
{ 
	double result = mu_;
	vector<Genotype> geno = genome.GetGenotype(QTLmap);
	for (int q=0;q<nQTL;q++)
	{
		int x = add_indicator(geno,q); // add: (--) -> -1; (-+) -> 0; (++) -> 1
		int y = dom_indicator(geno,q); // dom: (--) ->  0; (-+) -> 1; (++) -> 0
		result += x*add[q] + y*dom[q];
	}
	// two-way epistatic interactions (additive by additive)
	for (int q1=0;q1<nQTL;q1++)
	{
		int x1 = add_indicator(geno,q1);
		for (int q2=0;q2<q1;q2++)
		{
			int x2 = add_indicator(geno,q2);
			result += x1*x2*epi_axa[q1][q2];
		}
	}
	return result; 
}


std::vector<string> Phi::getQTLs(const Genome& genome) const
{
	vector<string> result(nQTL);
	vector<Genotype> geno = genome.GetGenotype(QTLmap);
	for (int q=0;q<nQTL;q++)
	{
		string score;
		Genotype g = geno[q];
		bool g1 = pos_allele[g.First()][q];
		bool g2 = pos_allele[g.Second()][q];
		if (g1 && g2)		// both positive
			score = "A";
		else if (!g1 && !g2) 
			score = "B";    // both negative
		else 
			score = "H";	// heterozygous
		result[q] = score;
	}
	return result;
}

std::vector<int> Phi::getORG(const Genome& genome,const Parent& par) const
{
	vector<alleletype> alleles = genome.GetGamete(QTLmap,par);
	vector<int> result(nQTL);
	for (int q=0;q<nQTL;q++)
		result[q] = (int) alleles[q];
	return result;
}