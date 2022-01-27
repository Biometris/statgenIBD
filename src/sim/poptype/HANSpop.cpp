#include "HANSpop.h"
#include "util_genetics.h"

#include <fstream>

using namespace std;
using namespace ibd;

poly ibd::operator*(const poly& f, const poly& g)
{
	poly res;
	for (poly::const_iterator it_f = f.equation.begin(); it_f != f.equation.end(); ++it_f)
	{
		int d_f = it_f->first.first;
		int k_f = it_f->first.second;
		long int c_f = it_f->second;

		for (poly::const_iterator it_g = g.equation.begin(); it_g != g.equation.end(); ++it_g)
		{
			int d_g = it_g->first.first;
			int k_g = it_g->first.second;
			long int c_g = it_g->second;
			//if (d_f + d_g < 2)
			res += poly(c_f*c_g, k_f + k_g, d_f + d_g);
		}
	}
	return res;
}

// generates the possible gametes with probabilities, given the genotype of individual ind.
map<haplo2,poly> ibd::meiosis(geno2& ind)
{
	map<haplo2,poly> res;
	haplo2 pH = ind.first;
	haplo2 mH = ind.second;

	haplo2 c1, c2;
	c1.first  = pH.first;
	c1.second = mH.second;

	c2.first  = mH.first;
	c2.second = pH.second;
  
	poly half(1,-1,0);
	poly min_halfr(-1,-1,1);
	poly halfr(1,-1,1);
	
	poly no_crossover(1,-1,0);
	no_crossover += poly(-1,-1,1);
	poly crossover(1,-1,1);

	res[pH] += no_crossover;
	res[mH] += no_crossover;
	res[c1] += crossover;
	res[c2] += crossover;

	return res;
}


map<geno2,poly> ibd::cross(geno2& parA, geno2& parB)
{
	map<geno2,poly> res;

    map<haplo2,poly> gamA = meiosis(parA);
    map<haplo2,poly> gamB = meiosis(parB);

	map<haplo2,poly>::const_iterator itA,itB;
	for (itA = gamA.begin(); itA != gamA.end(); itA++)
		for (itB = gamB.begin(); itB != gamB.end(); itB++)
			res[geno2(itA->first,itB->first)] += itA->second*itB->second; 

	return res;
}


map<brother_sister,poly> ibd::brother_sister_mating(map<brother_sister,poly>& par)
{
	map<brother_sister,poly> res;
	map<brother_sister,poly>::const_iterator it;
	for (it = par.begin(); it != par.end(); it++)
	{
		geno2 father = it->first.first;
		geno2 mother = it->first.second;
		map<geno2,poly> progeny = cross(father,mother);
		map<geno2,poly>::const_iterator it_b, it_s;
		for (it_b = progeny.begin(); it_b != progeny.end(); ++it_b)
		{
			for (it_s = progeny.begin(); it_s != progeny.end(); ++it_s)
			{
				poly prob = it_b->second*it_s->second*it->second;
				res[brother_sister(it_b->first,it_s->first)] += prob;
			}
		}
	}
	return res;
}

map<geno2,poly> ibd::brother_sister_mating_2(map<brother_sister,poly>& par)
{
	map<geno2,poly> res;
	map<brother_sister,poly>::const_iterator it;
	for (it = par.begin(); it != par.end(); it++)
	{
		geno2 father = it->first.first;
		geno2 mother = it->first.second;
		map<geno2,poly> progeny = cross(father,mother);
		map<geno2,poly>::const_iterator it2;
		for (it2 = progeny.begin(); it2 != progeny.end(); ++it2)
		{
			poly prob = it2->second*it->second;
			res[it2->first] += prob;
		}
	}
	return res;
}

map<geno2,poly> ibd::self_fert(map<geno2,poly>& par)
{
	map<geno2,poly> res;
	map<geno2,poly>::const_iterator it;
	for (it = par.begin(); it != par.end(); it++)
	{
		geno2 parent = it->first;
		map<geno2,poly> progeny = cross(parent,parent);
		map<geno2,poly>::const_iterator it2;
		for (it2 = progeny.begin(); it2 != progeny.end(); ++it2)
		{
			poly prob = it2->second*it->second;
			res[it2->first] += prob;
		}
	}
	return res;
}

map<geno2,poly> ibd::intermating(map<geno2,poly>& par)
{
	map<geno2,poly> res;
	map<geno2,poly>::const_iterator it1,it2;
	for (it1 = par.begin(); it1 != par.end(); it1++)
	{
		geno2 par1 = it1->first;
		for (it2 = par.begin(); it2 != par.end(); it2++)
		{
			geno2 par2 = it2->first;
			map<geno2,poly> progeny = cross(par1,par2);
			map<geno2,poly>::const_iterator it3;
			for (it3 = progeny.begin(); it3 != progeny.end(); ++it3)
			{
				poly prob = it3->second*it1->second*it2->second;
				res[it3->first] += prob;
			}
		}
	}
	return res;
}



// generates the possible gametes with probabilities, given the genotype of individual ind.
map<haplo2,double> ibd::meiosis(geno2& ind, double r)
{
	map<haplo2,double> res;
	haplo2 pH = ind.first;
	haplo2 mH = ind.second;

	haplo2 c1, c2;
	c1.first  = pH.first;
	c1.second = mH.second;

	c2.first  = mH.first;
	c2.second = pH.second;

	res[pH] += 0.5*(1.0-r);
	res[mH] += 0.5*(1.0-r);
	res[c1] += 0.5*r;
	res[c2] += 0.5*r;

	return res;
}

// generates the possible genotypes with probabilities, given the genotypes of the parents
map<geno2,double> ibd::cross(geno2& parA, geno2& parB, double r)
{
	map<geno2,double> res;

    map<haplo2,double> gamA = ibd::meiosis(parA,r);
    map<haplo2,double> gamB = ibd::meiosis(parB,r);

	map<haplo2,double>::const_iterator itA,itB;
	for (itA = gamA.begin(); itA != gamA.end(); itA++)
		for (itB = gamB.begin(); itB != gamB.end(); itB++)
			res[geno2(itA->first,itB->first)] += itA->second*itB->second; 

	return res;
}

map<brother_sister,double> ibd::brother_sister_mating(map<brother_sister,double>& par, 
													  double r)
{
	map<brother_sister,double> res;
	map<brother_sister,double>::const_iterator it;
	for (it = par.begin(); it != par.end(); it++)
	{
		geno2 father = it->first.first;
		geno2 mother = it->first.second;
		map<geno2,double> progeny = cross(father,mother,r);
		map<geno2,double>::const_iterator it_b, it_s;
		for (it_b = progeny.begin(); it_b != progeny.end(); ++it_b)
		{
			for (it_s = progeny.begin(); it_s != progeny.end(); ++it_s)
			{
				double prob = it_b->second*it_s->second*it->second;
				res[brother_sister(it_b->first,it_s->first)] += prob;
			}
		}
	}
	return res;
}

map<geno2,double> ibd::brother_sister_mating_2(map<brother_sister,double>& par, 
													  double r)
{
	map<geno2,double> res;
	map<brother_sister,double>::const_iterator it;
	for (it = par.begin(); it != par.end(); it++)
	{
		geno2 father = it->first.first;
		geno2 mother = it->first.second;
		map<geno2,double> progeny = cross(father,mother,r);
		map<geno2,double>::const_iterator it2;
		for (it2 = progeny.begin(); it2 != progeny.end(); ++it2)
		{
			double prob = it2->second*it->second;
			res[it2->first] += prob;
		}
	}
	return res;
}

map<geno2,double> ibd::self_fert(map<geno2,double>& par, double r)
{
	map<geno2,double> res;
	map<geno2,double>::const_iterator it;
	for (it = par.begin(); it != par.end(); it++)
	{
		geno2 parent = it->first;
		map<geno2,double> progeny = cross(parent,parent,r);
		map<geno2,double>::const_iterator it2;
		for (it2 = progeny.begin(); it2 != progeny.end(); ++it2)
		{
			double prob = it2->second*it->second;
			res[it2->first] += prob;
		}
	}
	return res;
}

map<geno2,double> ibd::intermating(map<geno2,double>& par, double r)
{
	map<geno2,double> res;
	map<geno2,double>::const_iterator it1,it2;
	for (it1 = par.begin(); it1 != par.end(); it1++)
	{
		geno2 par1 = it1->first;
		for (it2 = par.begin(); it2 != par.end(); it2++)
		{
			geno2 par2 = it2->first;
			map<geno2,double> progeny = cross(par1,par2,r);
			map<geno2,double>::const_iterator it3;
			for (it3 = progeny.begin(); it3 != progeny.end(); ++it3)
			{
				double prob = it3->second*it1->second*it2->second;
				res[it3->first] += prob;
			}
		}
	}
	return res;
}



HANSpop::HANSpop(alleletype p1, alleletype p2)						
     : AA(Genotype(p1,p1)), 
	   AB(Genotype(p1,p2)), 
	   BB(Genotype(p2,p2)),
       pop_base(ObsGeno(Genotype(p1,p1),
					   Genotype(p1,p2),
					   Genotype(p2,p2)), 1,"HANS") {}


double HANSpop::p_trans(const Genotype& g1, const Genotype& g2, double r) const
{
	throw ibd_error("Not implemented yet");
	return 0.0;
} 

vector<double> HANSpop::model(const Genotype& g) const
{
	throw ibd_error("Not implemented yet");
	vector<double> X(npar,0.0); // no dominance !
	return X;
}


Genome HANSpop::make_ind(const vector<double>& chr_length) const
{
	Genome F1(chr_length,AB); // F1 population

	// 32 F1
	// 16 G2
	//  8 G3
	//  4 G4
	//  2 G5
	//  1 G6 
	int N = 32;
	vector<Genome> G(N,F1);
	while (N > 1)
	{
		N /= 2;
		vector<Genome> G_next(N);
		for (int i=0;i<N;i++)
			G_next[i] = G[2*i]*G[2*i+1];

		G = G_next;		
	}
	Genome F = G[0];
	for (int i=0;i<8;i++)
		F = F*F;
	return F;
}

vector<Genome> HANSpop::make_pop(const vector<double>& chr_length, int nind) const
{
	throw ibd_error("Not implemented yet");	

	vector<Genome> pop(nind);
	return pop;
}

void HANSpop::MapQTL(ostream& outp, const ObsGeno& g) const
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
		throw ibd_error("unknown genotype in function HANSpop::MapQTL()");
}


ObsGeno HANSpop::MapQTL(istream& inp) const
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
		throw ibd_error("HANSpop::MapQTL(istream& )");
}


int ibd::test_HANS()
{
	cout.setf(ios::fixed, ios::floatfield);

	haplo2 hap1('A','A');
	haplo2 hap2('B','B');
	geno2 F1(hap1,hap2);

	map<geno2,poly> prob = cross(F1,F1);

	int N1 = 2;
	int N2 = 8;
	
	for (int i=0;i<N1;i++)
	{
		cout << "N1: " << i << endl;	
		prob = intermating(prob);
	}
	for (int j=0;j<N2;j++)
	{
		cout << "N2: " << j << endl;
		prob = self_fert(prob);
	}
	poly sum;
	for (map<geno2,poly>::const_iterator it = prob.begin(); it!=prob.end(); ++it)
	{
		cout << it->first.first.first 
			 << it->first.first.second << " "
			 << it->first.second.first  
			 << it->first.second.second << "   ";
		const poly& F = it->second;
		cout << F(0.1) << "   ";
		F.print();
		cout << endl;
		sum += F;
	}

	cout << endl << "sum : ";
	sum.print();
	cout << endl;

	return 0;

	ofstream outp("recomb3.dat");
	outp.setf(ios::fixed, ios::floatfield);
	outp << setprecision(8) << endl;

	for (double r = 0.0; r <= 0.0001; r += 0.00001)
	{
		map<geno2,double> prob2 = cross(F1,F1,r);
		for (int i=0;i<N1;i++)
			prob2 = intermating(prob2,r);
		for (int j=0;j<N2;j++)
			prob2 = self_fert(prob2,r);
	
		map<geno2,double>::const_iterator it2;
		double R = 0.0;
		for (it2 = prob2.begin(); it2 != prob2.end(); it2++)
		{
			if (it2->first.first.first != it2->first.first.second)
				R += it2->second;
		}

		map<geno2,poly>::const_iterator it3;
		poly sum;
		for (it3 = prob.begin(); it3 != prob.end(); it3++)
			if (it3->first.first.first != it3->first.first.second)
				sum += it3->second;
		sum.simplify();

		cout << setw(12) << r << setw(12) << R << setw(12) << sum(r) << endl;
		outp << setw(12) << r << setw(12) << R << setw(12) << sum(r) << endl;
	}
	return 0;
}

