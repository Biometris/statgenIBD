#include <fstream>
#include <sstream>
#include "GenoProb.h"

using namespace std;
using namespace ibd;

bool ibd::EXACT_RIL = true;

matrix<double> ibd::calc_prob_left(const vector<ObsGeno>& obs_geno,
							  const PopulationType& popt,
							  const vector<double>& r)
{
	ObsGeno U = popt.unknown_genotype();
	int dimU = U.size();
	int nloc = obs_geno.size();

	map<Genotype,int> ndx;
	for (int q=0;q<dimU;q++)
		ndx[U[q]] = q;

	matrix<double> L(nloc,dimU,0.0);

	ObsGeno g = obs_geno[0]; 
	for (unsigned int k = 0; k < g.size(); k++)
		L[0][ndx[g[k]]] = popt.p_trans(U[0],g[k],0.5);
	make_conditional(L[0]);
	for (int loc=1; loc < nloc; loc++)
	{
		g = obs_geno[loc];
		for (unsigned int k = 0; k < g.size(); k++)
			for (int q = 0; q < dimU; q++)
				L[loc][ndx[g[k]]] += L[loc-1][q]*popt.p_trans(U[q],g[k],r[loc-1]);
		make_conditional(L[loc]);
	}
	return L;
}	

matrix<double> ibd::calc_prob_right(const vector<ObsGeno>& obs_geno,
							  const PopulationType& popt,
							  const vector<double>& r)
{
	ObsGeno U = popt.unknown_genotype();
	int dimU = U.size();
	int nloc = obs_geno.size();

	map<Genotype,int> ndx;
	for (int q=0;q<dimU;q++)
		ndx[U[q]] = q;

	matrix<double> R(nloc,dimU,0.0);

	ObsGeno	g = obs_geno[nloc-1]; 
	for (unsigned int k = 0; k < g.size(); k++)
		R[nloc-1][ndx[g[k]]] = popt.p_trans(g[k],U[0],0.5);
	make_conditional(R[nloc-1]);
	for (int loc=nloc-2; loc >= 0; loc--)
	{
		g = obs_geno[loc];
		for (unsigned int k = 0; k < g.size(); k++)
			for (int q = 0; q < dimU; q++)
				R[loc][ndx[g[k]]] += R[loc+1][q]*popt.p_trans(g[k],U[q],r[loc]);
		make_conditional(R[loc]);
	}
	return R;
}


GenoProb::GenoProb(const matrix<ObsGeno>& obs_geno,
			const PopulationType& popt,
			const LinkageMap& markermap)
			: markermap_(markermap), popt_(popt)
{
	U = popt.unknown_genotype();
	dimU = U.size();
	nind = obs_geno.NrRows();
	nloc = obs_geno.NrCols();
	Tr_l = matrix3D<double>(nind,nloc,dimU);
	Tr_r = matrix3D<double>(nind,nloc,dimU);
	
	vector<double> r = make_rec_map(markermap);

	// if poptype is a RIL, use algorithms from calcRILexact.*
	const int gen = IsRIpop(popt);
	const int N = gen-1;
	const int M = 2*N;
	for (int ind=0;ind<nind;ind++)
	{
		if (gen > 0 && EXACT_RIL)
		{
			vector<ObsInhVec> obs_inh_vect = obs_inh_vectors(obs_geno[ind],N);
			Tr_l[ind] = calc_prob_left(obs_inh_vect,M,r);
			Tr_r[ind] = calc_prob_right(obs_inh_vect,M,r);
		}
		else
		{
			Tr_l[ind] = calc_prob_left(obs_geno[ind],popt,r);
			Tr_r[ind] = calc_prob_right(obs_geno[ind],popt,r);
		}
	}
}

matrix<double> GenoProb::operator()(const Locus& QTLpos) const
{
	// if poptype is a RIL, use algorithms from calcRILexact.*
	const int gen = IsRIpop(popt_);
	if (gen > 0 && EXACT_RIL)
		return calc_QTL_prob(U,markermap_,Tr_l,Tr_r,nind,2*(gen-1),QTLpos);

	matrix<double> P(nind,dimU);
	const int left = pos_qtl(markermap_,QTLpos);
	const int right = left + 1;
	
	double r_left = recomb(markermap_[left],QTLpos);
	double r_right = recomb(QTLpos,markermap_[right]);

	matrix<double> T1 = transition_matrix(U,U,r_left,popt_);
	matrix<double> T2 = transition_matrix(U,U,r_right,popt_);

	for (int i=0;i<nind;i++)
	{
		const vector<double>& L = Tr_l[i][left];
		const vector<double>& R = Tr_r[i][right];
		vector<double>& prob_QTL = P[i];

		vector<double> left_prob(dimU,0.0);
		vector<double> right_prob(dimU,0.0);
		for (int q1=0;q1<dimU;q1++)
		{
			for (int q2=0;q2<dimU;q2++)
			{
				left_prob[q1]  += L[q2]*T1[q2][q1];
				right_prob[q1] += T2[q1][q2]*R[q2];
			}
		}
		for (int q=0;q<dimU;q++)
			prob_QTL[q] = left_prob[q]*right_prob[q];
		make_conditional(prob_QTL);
	}

	return P;
}

// check population : -1 if its not a RIL, otherwise generation of RIL 
int ibd::IsRIpop(const PopulationType& popt)
{
	const string name = popt.name();
	if (name.find("RI")==0)
	{
		int gen;
		istringstream str(name.c_str());
		str.ignore(2);
		str >> gen;
		return gen;
	}
	return -1;
}

/*
vector<double> GenoProb::operator()(int ind, int loc) const
{
	vector<double> prob_QTL(dimU);
	const vector<double>& L = Tr_l[ind][loc];
	const vector<double>& R = Tr_r[ind][loc];
	for (int q=0;q<dimU;q++)
		prob_QTL[q] = L[q]*R[q];
	make_conditional(prob_QTL);
	return prob_QTL;
}
*/
