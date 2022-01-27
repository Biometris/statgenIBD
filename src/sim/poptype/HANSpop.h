#ifndef HANSPOP_HEADER
#define HANSPOP_HEADER

#include "poptype.h"

namespace ibd
{

class HANSpop : public pop_base
{
public:
	HANSpop(alleletype p1, alleletype p2);
	double p_trans(const Genotype& g1, const Genotype& g2, double r) const;
	std::vector<double> model(const Genotype& g) const;
	Genome make_ind(const std::vector<double>& chr_length) const;
	std::vector<Genome> make_pop(const std::vector<double>& chr_length, int nind) const;
	void MapQTL(std::ostream& outp, const ObsGeno& g) const;
	ObsGeno MapQTL(std::istream& inp) const;
private:
	const Genotype AA,AB,BB;
};


typedef std::pair<char,char> haplo2;
typedef std::pair<haplo2,haplo2> geno2;
typedef std::pair<geno2,geno2> brother_sister;

std::map<haplo2,double> meiosis(geno2& ind, double r);
std::map<geno2,double> cross(geno2& parA, geno2& parB, double r);
std::map<brother_sister,double> brother_sister_mating(std::map<brother_sister,double>& par, 
													  double r);
std::map<geno2,double> brother_sister_mating_2(std::map<brother_sister,double>& par,
											   double r);

std::map<geno2,double> self_fert(std::map<geno2,double>& par, double r);

class poly
{
public:
	typedef std::map< std::pair<int,int>, long int>::const_iterator const_iterator;
	typedef std::map< std::pair<int,int>, long int>::iterator iterator;
	
	poly() {}
	poly(long int coeff, int n, int d) { equation[std::pair<int,int>(d,n)] = coeff; }

	void operator+=(const poly& F)
	{
		for (const_iterator it = F.equation.begin(); it != F.equation.end(); ++it)
		{
			iterator it2 = equation.find(it->first);
			if (it2 == equation.end() || it2->second + it->second != 0)
				equation[it->first] += it->second;
			else 
				equation.erase(it2);
		}
		simplify();
	}
	void simplify()
	{
		for (;;)
		{
			std::map< std::pair<int,int>, long int> new_equation;
			const_iterator it;
			for (it = equation.begin(); it != equation.end(); ++it)
			{
				int d = it->first.first;
				int k = it->first.second;
				long int c = it->second;
				while (c % 2 == 0)
				{
					c /= 2;
					k++;
				}
				new_equation[std::pair<int,int>(d,k)] += c;
			}
			if (new_equation == equation) return;
			equation.clear();
			for (it = new_equation.begin(); it != new_equation.end(); ++it)
				if (it->second != 0)
					equation[it->first] = it->second;
		}
	}

	void print(int degree = 1000) const
	{
		for (const_iterator it = equation.begin(); it != equation.end(); ++it)
		{
			int d = it->first.first;
			int k = it->first.second;
			long int c = it->second;
			if (d > degree) break;
			std::cout << c << " * 2^(" << k << ") * r^" << d;
			const_iterator it2 = it;
			it2++;
			if (it2 != equation.end())
				std::cout << " +";
		}
	}
	double operator() (double r) const 
	{
		double sum = 0.0;
		for (const_iterator it = equation.begin(); it != equation.end(); ++it)
		{
			int d = it->first.first;
			int k = it->first.second;
			long int c = it->second;
			sum += c*pow(2.0,k)*pow(r,d);
		}
		return sum;
	}
private:
	std::map< std::pair<int,int>, long int> equation;
	friend poly operator*(const poly& f, const poly& g);
};

poly operator*(const poly& f, const poly& g);


std::map<haplo2,poly> meiosis(geno2& ind);
std::map<geno2,poly> cross(geno2& parA, geno2& parB);

std::map<brother_sister,poly> brother_sister_mating(std::map<brother_sister,poly>& par);
std::map<geno2,poly> brother_sister_mating_2(std::map<brother_sister,poly>& par);
std::map<geno2,poly> self_fert(std::map<geno2,poly>& par);

std::map<geno2,poly> intermating(std::map<geno2,poly>& par);

std::map<geno2,double> intermating(std::map<geno2,double>& par, double r);

int test_HANS();


}

#endif