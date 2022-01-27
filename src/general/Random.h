#ifndef RANDOMPROB_HEADER
#define RANDOMPROB_HEADER

#include <cmath>
#include <vector>
#include "misc.h"

namespace ibd
{
	
const double TWO_PI			 = 6.283185307179586;

void   randstart();
double randnormal();
double randuniform();
double randexp(double lambda);
double randexp(double lambda, double max); // truncated exponential distribution, CtB
int	   rand_poisson(double lambda);
int    randint(int a, int b);
void set_randomseed(long int);
long int get_randomseed();

double random_chi_sqr(int df);

// inline functions:
inline double randnormal(double mu, double sigma) 
{ return mu + sigma*randnormal();}

inline double randuniform(double a, double b)
{ return a + (b-a)*randuniform(); }

// returns a random integer x,
// P(X = x) = 1/(b-a+1) for a <= x <= b;
inline int randint(int a, int b)
{ return (int) (a + floor((b-a+1)*randuniform())); }

// exponential distribution: f(x)=lambda*exp(-lambda*x) 
inline double randexp (double lambda)
{ return -log(randuniform())/lambda; }

// return a random integer x, 0 <= x < N
inline int randint(int N)
{ return randint(0,N-1); }

std::vector<int> random_order_index(int N);

/*
// see Gelman p. 478,479 (first edition)
class rand_multinormal : public ibd::void_function< std::vector<double> >
{
public:
	rand_multinormal(const std::vector<double>& mu, const matrix<double>& Sigma)
		: _mu(mu), _A(CholeskyDecomp(Sigma)) {} 
	void set(const std::vector<double>& mu, const matrix<double>& Sigma)
		{_mu = mu; _A = CholeskyDecomp(Sigma);} 
	std::vector<double> operator()() const;
private:
	std::vector<double> _mu;
	matrix<double> _A;
};

matrix<double> wishart(const matrix<double>& PriorW, int df); // Joao Paulo

double rand_sc_inv_chi_sqr(double nu, double s_sqr); // CtB

// use this function instead of the STL function test_random_shuffle(),
// see also the test function test_random_permutation()
template <class RandomAccessIterator>
void random_permutation(RandomAccessIterator first, RandomAccessIterator last)
{	
	for (int n = last - first;n>1;n--)
		swap(*(first + randint(n)), *--last);
}

*/


template <class C>
void random_permutation(C& c) { random_permutation(c.begin(),c.end()); }



// alias method, see The Cross-Entropy Method (Rubinstein & Kroese), p 21

template<class T>
bool compare(const std::pair<T,double>& p1, const std::pair<T,double>& p2)
{
	return p1.second < p2.second;
}

template <class T>
struct ProbTable : public std::vector< std::pair<T,double> >
{
	void add(const T& key, double pval)
	{ push_back(std::pair<T,double>(key,pval)); }
};

/*
template <class T>
class DiscreteDistr
{
public:
	DiscreteDistr() {}
	DiscreteDistr(ProbTable<T> table)
	{
		typedef ProbTable<T>::iterator Iter;
		typedef ProbTable<T>::const_iterator ConstIter;
		
		if (table.size() < 2)
			throw MQMlib_error("DiscreteDistr<Type>: # elem < 2");
		M = table.size()-1;
		// make sum of elements in table equal to 1.
		double sum = 0.0;
		for (ConstIter it=table.begin(); it!=table.end(); it++)
			sum += it->second;
		for (Iter it=table.begin(); it!=table.end(); it++)
			it->second/=sum;

		// generate alpha, L and R
		for (int i=0;i<M;i++)
		{
			std::sort(table.begin(),table.end(),compare<T>);
			Iter it_min = table.begin();
			Iter it_max = table.end()-1;
			double gamma = M*it_min->second;
			alpha.push_back(gamma);
			it_max->second -= (1-gamma)/M;
			L.push_back(it_min->first);
			R.push_back(it_max->first);
			table.erase(it_min);
		}
	}

	T operator()() const 
	{
		// see Algorithm 1.7.3 Rubinstein & Kroese
		int u = ibd::randint(0,M-1); // M = N-1
		return (ibd::randuniform() < alpha[u]) ? L[u] : R[u];
	}

	double prob(const T& x) const 
	{
		double p = 0.0;
		for (int i=0;i<M;i++)
		{
			if (L[i] == x) 
				p += alpha[i]/M;
			else if (R[i] == x)
				p += (1-alpha[i])/M;
		}
		return p;
	}
private:
	std::vector<T> L,R;
	std::vector<double> alpha;
	int M;
};

template <typename T>
class Composition : public DiscreteDistr<T>
{
public:
	Composition(const ProbTable<T>& tab) : distr(tab) {}  
	typename T::result_type operator()() const 
	{ 
		T H = distr();
		return H(); 
	}
private:
    DiscreteDistr<T> distr;
};

void test_multinormal();
void test_scaled_inv_chi_sqr();
int test_random_permutation(int K, bool test_shuffle);
int test_AliasMethod();
int test_Composition();
*/


}


#endif