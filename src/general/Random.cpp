#include <time.h>
#include <float.h>
#include <stdlib.h>
#include <limits.h>
#include <algorithm>
#include <functional>
#include <fstream>
#include "misc.h"
#include "Random.h"

#define IM1     2147483563L
#define IM2     2147483399L
#define AM        (1.0/IM1)
#define IMM1        (IM1-1)
#define IA1          40014L
#define IA2          40692L
#define IQ1          53668L
#define IQ2          52774L
#define IR1          12211L
#define IR2           3791L
#define NTAB             32
#define NDIV  (1+IMM1/NTAB)
#define EPS     DBL_EPSILON
#define RNMX      (1.0-EPS)

using namespace std;
using namespace ibd;

namespace {

long int ran2_idum = 0;

}

// 12 maart 2002: get_randomseed en set_randomseed werken alleen
// als get_randomseed eerste getal van serie is (negatief getal)
void ibd::set_randomseed(long int k)
{  
	ran2_idum = k;
}

long int ibd::get_randomseed()
{
	return ran2_idum;
}

double ibd::randuniform()
{
	long *idum=&ran2_idum;
	long j, k;
	static long idum2=123456789L;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0) 
	{
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) 
		{
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[(int)j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[(int)j]-idum2;
	iv[(int)j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

void ibd::randstart()
{
	ran2_idum= -time(NULL);
} 

double ibd::randnormal()
{
	static char chifirst= 1;
	static double chisecond;
	double angle, radius, cosang, sinang;

	if (chifirst)
	{
		radius= sqrt(-2*log(randuniform()));
		angle= TWO_PI*randuniform();
		cosang= cos(angle);
		sinang= sqrt(1.0-cosang*cosang);
		if (angle>M_PI) sinang= -sinang;
		chisecond= radius*sinang;
		chifirst= 0;
		return  radius*cosang;
    }
	else
    {
		chifirst= 1;
		return  chisecond;
    }
} 

int ibd::rand_poisson(double lambda)
{
	double U = randuniform();
	double p = exp(-lambda);
	double sum = p;
	int i = 0;
	while (sum < U)
	{
		i++;
		p *= (lambda/i);
		sum += p;
	}
	return i;
}

double ibd::random_chi_sqr(int df)
{
	double sum = 0.0;
	for (int i=0;i<df;i++)
		sum += sqr(randnormal());
	return sum;
}

/*
vector<double> ibd::rand_multinormal::operator()() const
{
	const int dim = _mu.size();
	vector<double> z(dim);
	for (int i=0;i<dim;i++)
		z[i] = randnormal();
	return _mu + _A*z;
}

// Wishart function, Joao Paulo, juli 4 2003
ibd::matrix<double> ibd::wishart(const matrix<double>& PriorW, int df)
{
	int p = PriorW.NrCols();
	vector<double> muNull(p,0.0);				 // vector muNull
	matrix<double> Z(df,p);				         // matrix of df * p 	
	rand_multinormal  MultiNormal(muNull,PriorW); 
	
	for (int j = 0; j < df; j++) 
		Z[j] = MultiNormal();

	return transpose(Z)*Z;

}

// random number from a scaled inverse chi-square distribution 
// (defined as in Gelman et al 2003)
// last modified by CtB November 6 2003
// uses new edgamma definition using rgs_ for nu/2 <=1 and as91 for nu/2 >1
double ibd::rand_sc_inv_chi_sqr(double nu, double s_sqr)
{
	double X,r;
	do // only accept positive values for r
	{
		X = 2.0 * edgamma(nu/2.0, 1);  
		r = nu*(s_sqr+DBL_MIN)/X;
	}
	while (r < DBL_MIN); // min positive value

	return r;
}

// Cajo ter Braak
double ibd::randexp(double lambda, double max)
{
    if (fabs(lambda) < DBL_MIN)			// uniform if lambda == 0
		return randuniform(0.0, max);
	else
	{
		double p = randuniform();
		return  -log(p + (1.0-p)*exp(-lambda*max) )/lambda;  
	}
}

vector<int> ibd::random_order_index(int N)
{
	vector<int> index(N);
	for (int i=0;i<N;i++)
		index[i] = i;
	random_permutation(index);
	return index;
}
*/

