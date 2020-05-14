// Martin Boer, Biometris
#include <map>

#include "misc.h"
#include "mblexcept.h"

#include "markerscore.h"
#include "analysis_ped.h"
#include "analysis_sub.h"

using namespace mbl;
using namespace std;

matrix<double> calc_coa(const Pedigree& ped)
{
	const int N = ped.size();
	matrix<double> X(N,N,0.0);
	for (int i=0;i<N;i++)
	{
		X[i][i] = 1.0;
		if (!(ped[i].IsFounder()))
		{
			int p1 = ped[i][0];
			int p2 = ped[i][1];
			for (int j=0;j<i;j++)
			{
				X[i][j] = X[j][i] = 0.5*X[p1][j] + 0.5*X[p2][j];
			}
		}
	}
	return X;
}

vector<bool> member_sub_pedigree(const Pedigree& ped,
								 const vector<bool>& hd, const int N)
{
	matrix<int> Z(N,N,0);
	for (int i=0;i<N;i++)
	{
		Z[i][i] = 1;
		if (!(ped[i].IsFounder() || hd[i]))
			Z[i] += Z[ped[i][0]] + Z[ped[i][1]];
	}
	vector<int> v = Z[ped[N][0]] + Z[ped[N][1]];

	vector<bool> result(N);
	for (int i=0;i<N;i++)
		result[i] = (v[i] != 0);
	return result;
}

vector<int> sub_pedigree(const Pedigree& ped, const vector<bool>& hd, const int N)
{
	vector<bool> member = member_sub_pedigree(ped,hd,N);
	vector<int> ndx;
	for (int i=0;i<N;i++)
		if (member[i])
			ndx.push_back(i);
	return ndx;
}

int cnt_non_founders_sub_ped(const Pedigree& ped, const vector<bool>& hd, const int N)
{
	vector<bool> member = member_sub_pedigree(ped,hd,N);
	int cnt=1;
	for (int i=0;i<N;i++)
		if (member[i] && !ped[i].IsFounder() && !hd[i])
			cnt++;
	return cnt;
}

void pre_analysis_tabular_method(const Pedigree& ped,
								const vector<bool>& ancestor)
{
	const int N = ped.size();
	map<int,int> hist;
	int max_sz = -1;
	for (int i=0;i<N;i++)
	{
		if (!ped[i].IsFounder())
		{
			int sz = cnt_non_founders_sub_ped(ped,ancestor,i);
			if (sz > max_sz)
				max_sz = sz;
			hist[sz]++;
		}
	}
	cout << "Histogram: " << endl;
	cout << setw(4) << "sz" << setw(8) << "number" << endl;
	for (map<int,int>::const_iterator it=hist.begin();it!=hist.end();it++)
		cout << setw(4) << it->first << setw(8) << it->second << endl;

	// cout << "Max dimension subpopulation " << "2^" << max_sz << endl;
	if (max_sz > 15)
		throw mblib_error("Pedigree too complex to analyse!");
}

// output: indices of the sub_pedigree
vector<int> make_selection(Pedigree& ped_sel,
						   matrix<int>& geno_sel,
						   const Pedigree& ped,
						   const matrix<int>& geno,
						   const vector<bool>& hd,
						   int indnr)
{
	vector<int> sel = sub_pedigree(ped,hd,indnr);
	map<int,int> ndx_sel = make_index(sel);
	int Nsel = sel.size();
	vector<ParentsInd> sel_parents(Nsel+1);
	geno_sel = matrix<int>(Nsel+1,0);

	vector<int> ndx_fnd;
	for (int i=0;i<Nsel;i++)
	{
		int k = sel[i];
		geno_sel[i] = geno[k];
		if (ped[k].IsFounder() || hd[k])
		{
			// default value of sel_parents[i] is (-1,-1), i.e. a founder;
			ndx_fnd.push_back(k);
      	}
		else
		{
			sel_parents[i] = ParentsInd(ndx_sel[ped[k][0]],ndx_sel[ped[k][1]]);
		}
	}

	geno_sel[Nsel] = geno[indnr];
	sel_parents[Nsel] = ParentsInd(ndx_sel[ped[indnr][0]],ndx_sel[ped[indnr][1]]);

	ped_sel = Pedigree(sel_parents);

	return ndx_fnd;
}

vector<int> calc_IBD_subpedigree(matrix<double>& P,
								 const Pedigree& ped,
								 const matrix<int>& geno,
								 const vector<bool>& ancestor,
								 const vector<double>& r,
								 int indnr,double eps)
{
	Pedigree ped_sel;
	matrix<int> geno_sel;
	vector<int> a = make_selection(ped_sel,geno_sel,ped,geno,ancestor,indnr);
	P = calc_IBD_prob_last_ind(ped_sel,geno_sel,r,eps);
	return a;
}

// generalized 'tabular' method for inbred lines
//   input: pop, vector of individuals
//		    geno,  matrix of genotypic scores
//			ancestor, use individual as ancestor
//						all true:  use for all the non-founders the parents
//				        all false: exact IBD-prob for last individual
//`         r, vector of recombinations
//          eps, allow for errors in genotypic scores.
//  output: See IBDped::operator()(int m)
IBDped::IBDped(const Pedigree& ped,
			  const mbl::matrix3D<double>& IBD_founders,
			  const mbl::matrix<int>& geno,
			  const std::vector<bool>& ancestor,
			  const std::vector<double>& r,
			  double eps) : nloc(geno.NrCols()), nfnd(IBD_founders.Dim2()),IBD_fnd(IBD_founders)
{
	const int N = ped.size();
	for (int i=nfnd;i<N;i++)
	{
		//cout << "analysis ind " << i << endl;
		vector<int> a;		// vector with indices for ancestors
		matrix<double> P;   // rows: Number of markers; cols number of ancestors
		//if (!ped[i].IsFounder()) // calc a and P (if founder, both a and P empty)
		a = calc_IBD_subpedigree(P,ped,geno,ancestor,r,i,eps);
		anc.push_back(a);
		prob.push_back(P);
	}
}

// calculates the NxN coefficient coancestry matrix for marker position m:
mbl::matrix<double> IBDped::operator()(int m) const
{
	const int N = anc.size() + nfnd; // change this code!!
	matrix<double> X(N,N);

	const matrix<double>& IBD_fnd_cur = IBD_fnd[m];
	for (int i=0;i<nfnd;i++)
		for (int j=0;j<nfnd;j++)
			X[i][j] = IBD_fnd_cur[i][j];

	for (int i=nfnd;i<N;i++)
	{
		int ii=i-nfnd;
		X[i][i] = 1.0;
		const vector<int>& a = anc[ii];
		const int N_ancestors = a.size();

		const matrix<double>& P = prob[ii];
		const vector<double>& lambda = P[m];
		// generalized tabular method:
		for (int j=0;j<i;j++)
		{
			double sum = 0.0;
			for (int k=0;k<N_ancestors;k++)
			{
				sum += lambda[k]*X[a[k]][j];
			}
			X[i][j] = X[j][i] = sum;
		}
	}
	return X;
}

int analysis_pedigree(IBDped& IBD,				    // output
					  Pedigree& ped,				// output
					  vector<int>& ndx,				// output
					  const matrix3D<double>& IBD_fnd,
					  const vector<IndProp>& pop,
				      const matrix<score>& geno,
					  const vector<double>& r,
					  const set<string>& sel_ind,
					  double fraction)
{
	cout << "start analysis pedigree ..." << endl;
	//const double eps = 0.0001;
	const double eps = 0.00;
	vector<IndProp> sel_pop;
	matrix<int> haplo = gen_haplo_score_ped(sel_pop,geno,pop);
	ped = Pedigree(sel_pop);

	const int M = haplo.NrCols();
	int min_scored_markers = int(fraction*M);
	vector<bool> ancestor = detect_high_density_ind(haplo, min_scored_markers);

	pre_analysis_tabular_method(ped,ancestor);

	IBD = IBDped(ped,IBD_fnd,haplo,ancestor,r,eps);

	ndx = get_ndx_set(sel_pop,sel_ind);

	return 0;
}

/*
double gen_formula_tab_method(const vector<double>& lambda,
	   						  const vector<int>& a,
							  const matrix<double>& Xm, int j)
{
	const int dim = a.size();
	double sum = 0.0;
	for (int k=0;k<dim;k++)
		sum += lambda[k]*Xm[a[k]][j];
	return sum;
}

void update_coa_matrix(matrix<double>& Xm,
					   const vector<double>& lambda,
					   const vector<int>& a,
					   int i)
{
	for (int j=0;j<i;j++)
	{
		Xm[i][j] = Xm[j][i] = gen_formula_tab_method(lambda,a,Xm,j);
	}
	Xm[i][i] = 1.0;
}

// generalized 'tabular' method for inbred lines
//   input: pop, vector of individuals
//		    geno,  matrix of genotypic scores
//			ancestor, use individual as ancestor
//						all true:  use for all the non-founders the parents
//				        all false: exact IBD-prob for last individual
//`         r, vector of recombinations
//          eps, allow for errors in genotypic scores.
//   ouput: 3D-dimensional matrix X, where X[m][i][j] is the (approximate)
//          probability ind. i is IBD to another individual j at marker position m.
void generalized_tabular_method_coa(matrix3D<double>& X,
	                               const Pedigree& ped,
							       const matrix<int>& geno,
								   const vector<bool>& ancestor,
								   const vector<double>& r,
								   double eps)
{
	int N = ped.size();
	int M = geno.NrCols();
	matrix<double> zero_mat(N,N,0.0);
	X = matrix3D<double>(M,0,0);
	for (int m=0;m<M;m++)
		X[m] = zero_mat;
	//X = matrix3D<double>(M,N,N,0.0);

	for (int i=0;i<N;i++)
	{
		if (ped[i].IsFounder())
		{
			for (int m=0;m<M;m++)
				X[m][i][i] = 1.0;
		}
		else
		{
			matrix<double> P;
			vector<int> a = calc_IBD_subpedigree(P,ped,geno,ancestor,r,i,eps);
			for (int m=0;m<M;m++)
			{
				update_coa_matrix(X[m],P[m],a,i);
			}
		}
	}

}
*/

