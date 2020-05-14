// Martin Boer, Biometris
#include "InhVector.h"
#include "analysis_sub.h"
#include "TransMatSym2D.h"
#include "HMMalgo.h"
#include "markerscore.h"

using namespace mbl;
using namespace std;

// n is the length of the inheritance vector (number of inheritance indicactors h_i = {0,1})
vector<TransMatSym2D> generate_transition_matrices(const vector<double>& vec_recomb, int n)
{
	int Nintervals = vec_recomb.size();
	vector<TransMatSym2D> T(Nintervals); // Transition matrices between the loci.	
	for (int k=0;k<Nintervals;k++)
	{
		double r = vec_recomb[k];
		// For the moment we assume that R is equal for all inheritance indicators, i.e. R[i] = r.
		// Further extensions are possible, e.g. a mixture of DH and RIL individuals, e.g.
		// R[0]   = r;          
		// R[1]   = 2*r/(1+2*r); (RIL)
		// R[2]   = 2*r/(1+2*r); (RIL)
		// ....
		// R[n-1] = r;
		vector<double> R(n,r);  
		T[k] = TransMatSym2D(R);
	}
	return T;
}

bool consistent_score(vector<int>& g, int i, int par)
{
	if (g[i] == g[par] || g[par] == U_haplo_sc) 
		return true;

	if (g[i] == U_haplo_sc)
	{
		g[i] = g[par];
		return true; 
	}
	return false;
}

bool possible_inh_vector(InhVector x, const matrix<int>& geno, int m, const Pedigree& ped)
{
	int N = geno.NrRows();
	vector<int> g(N);
	for (int i=0;i<N;i++)
	{
		g[i] = geno[i][m];
		if (!ped[i].IsFounder()) // if non-founder
		{
			int par = ped[i][x.next_indic()]; // get parent
			if (!consistent_score(g,i,par))
				return false;
		}
	}
	return true;
}

// n is number of non-founders.
matrix<double> generate_Q_matrix(const matrix<int>& geno, const Pedigree& ped, int n,
								 double eps)
{
	int M = geno.NrCols();
	unsigned int dim_state_space = pow2(n); 
	matrix<double> Q(M,dim_state_space);

	for (int m=0;m<M;m++)
	{
		vector<int> x(dim_state_space);
		for (InhVector u(n);!u.end();u++)
			x[u] = possible_inh_vector(u,geno,m,ped);

		int sum = accumulate(x.begin(),x.end(),0);
		if (sum > 0 && sum < (int) dim_state_space)
		{
			for (InhVector u(n);!u.end();u++)
			{
				Q[m][u] = (x[u] == 1) ? 1.0-eps/sum : eps/(dim_state_space-sum);
			}
		}
		else
			Q[m] = vector<double>(dim_state_space,1.0); 
	}
	return Q;
}

matrix<double> calc_inheritance_prob(const Pedigree& ped,
							         const matrix<int>& geno,
							         const vector<double>& r, double eps)
{
	const int N = geno.NrRows();
	int Nfnd = ped.Nfounders();
	int n = N - Nfnd;

	unsigned int dim_state_space = pow2(n);
	vector<double> pi0(dim_state_space,1.0/dim_state_space);
	vector<TransMatSym2D> T = generate_transition_matrices(r,n);
	matrix<double> Q = generate_Q_matrix(geno,ped,n,eps);
	return calc_prob(pi0,Q,T);
}

int ancestor_last_nonfounder(const Pedigree& ped, InhVector x)
{
	int N = ped.size();
	vector<int> Z(N);
	int y = -1; 
	vector<int> Y;
	int fnd_nr = 0;
	for (int i=0;i<N;i++)
	{
		if (ped[i].IsFounder())
		{
			Z[i] = fnd_nr++;
		}
		else
		{
			Z[i] = Z[ped[i][x.next_indic()]];
			y = Z[i];
		}
	}
	return y;
}

matrix<double> calc_IBD_prob_last_ind(const Pedigree& ped,
							   const matrix<int>& geno,
							   const vector<double>& r,double eps)
{
	matrix<double> P = calc_inheritance_prob(ped,geno,r,eps);
	int M = P.NrRows();
	int N = ped.size();
	int Nfnd = ped.Nfounders();
	int n = N-Nfnd; // number of nonfounders

	matrix<double> IBD(M,Nfnd,0.0);
	for (InhVector u(n);!u.end();u++)
	{
		int fndnr = ancestor_last_nonfounder(ped,u);
		for (int m=0;m<M;m++)
		{
			IBD[m][fndnr] += P[m][u];
		}
	}
	return IBD;
}


