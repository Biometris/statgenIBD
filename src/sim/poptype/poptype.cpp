#include "poptype.h"
#include "BCpop.h"
#include "F2pop.h"
#include "DHpop.h"
#include "RILpop.h"
#include "PIONpop.h"
#include "HANSpop.h" // added april 17, 2003
#include <strstream>
#include <numeric>

using namespace std;

namespace 
{
using namespace ibd;

int calcprob_aux(matrix<Genotype>& geno_augm, 
				 vector<double>& P,
				 double prob, 
				 int level, 
				 int prev,
				 int nloc,
				 vector<Genotype>& geno,
				 const vector<ObsGeno>& obs_geno,
				 const vector< matrix<double> >& P_intra,
				 double total_prob,
				 double threshold)
{
	double new_prob;
	if (level == nloc)
	{
		new_prob = prob*P_intra[level][prev][0];
		if (new_prob <= threshold*total_prob)
			return 0;
		geno_augm.push_back(geno);
		P.push_back(new_prob/total_prob);
		return 1;
	}
	int sum = 0;
	ObsGeno M = obs_geno[level]; // observation of current marker M 
	int nr = M.size();
	for (int k=0;k<nr;k++)
	{
		geno[level] = M[k];
		new_prob = prob*P_intra[level][prev][k];
		if (new_prob > threshold*total_prob)
			sum += calcprob_aux(geno_augm,P,new_prob,level+1,k,
							nloc,geno,obs_geno,P_intra,total_prob,threshold);
    }
	return sum;
}

bool no_missing_marker_info(matrix<Genotype>& geno_augm,
				  vector<double>& P, const vector<ObsGeno>& obs_geno)
{
	int nloc = obs_geno.size();
	vector<Genotype> geno(nloc);
	for (int i=0;i<nloc;i++)
	{
		const ObsGeno& M = obs_geno[i]; 
		if (M.size() > 1) // if marker information is uncomplete
			return false;
		else
			geno[i] = M[0];
	}
	// no markers with uncomplete information
	P.push_back(1.0);
	geno_augm.push_back(geno);
	return true;
};

};

ibd::PopulationType::PopulationType(const string& pop_name)
{
	h_popt = 0;
	if (pop_name == "BC1" || pop_name == "BC") 
		h_popt = new BCpop('A','B'); // BC1
	else if (pop_name == "F2") 
		h_popt = new F2pop('A','B'); // F2
	else if (pop_name[0] == 'D' && pop_name[1] == 'H')
		h_popt = new DHpop('A','B');
	else if (pop_name == "HANS")
		h_popt = new HANSpop('A','B'); // added april 17, 2003
	else if (pop_name[0] == 'R' && pop_name[1] == 'I') // RIX
	{
		int gen;
		istrstream str(pop_name.c_str());
		str.ignore(2);
		str >> gen;
		if (str.eof())
		   h_popt = new RILpop('A','B',gen);
	}
	else if (pop_name[0] == 'P' && pop_name[1] == 'I' &&
		     pop_name[2] == 'O' && pop_name[3] == 'N') // PIONx
	{
		int gen;
		istrstream str(pop_name.c_str());
		str.ignore(4);
		str >> gen;
		if (str.eof())
		   h_popt = new PIONpop('A','B',gen);
	}
	if (!h_popt)
		throw ibd_error("Unknown population type " + pop_name);
}


ibd::matrix<double> ibd::transition_matrix(const ObsGeno& g1, const ObsGeno& g2, 
									 double r, const PopulationType& poptype)
{
	const int nrows = g1.size();
	const int ncols = g2.size();
	matrix<double> result(nrows,ncols);
	for (int i=0;i<nrows;i++)
		for (int j=0;j<ncols;j++)
			result[i][j] = poptype.p_trans(g1[i],g2[j],r);
	return result;
}


int ibd::calcprob_ind(matrix<Genotype>& geno_augm,
				 vector<double>& P,
				 const vector<double>& r,
				 const vector<ObsGeno>& obs_geno,
				 const PopulationType& poptype,
				 double threshold)
{
	if (no_missing_marker_info(geno_augm,P,obs_geno))
		return 1;
	const int nloc = obs_geno.size();
	ObsGeno g_aux(obs_geno[0][0]); 
	vector< matrix<double> > P_intra(nloc+1); 
	matrix<double> A;
	P_intra[0] = transition_matrix(g_aux,obs_geno[0],0.5,poptype);
	A = P_intra[0];
	for (int i=0;i<nloc-1;i++)
	{
		P_intra[i+1] = transition_matrix(obs_geno[i],obs_geno[i+1],r[i],poptype);
		A = A*P_intra[i+1];
	}
	P_intra[nloc] = transition_matrix(obs_geno[nloc-1],g_aux,0.5,poptype);
	A = A*P_intra[nloc];
	double total_prob = A[0][0];

	vector<Genotype> geno(nloc);
	
	int nr = calcprob_aux(geno_augm,P,1.0,0,0,nloc,geno,obs_geno,P_intra,total_prob,threshold);
	vector<double>::iterator end   = P.end();
	vector<double>::iterator start = end - nr;
	double total_sum = accumulate(start,end,0.0);
	for (vector<double>::iterator iter=start; iter != end; ++iter)
		*iter /= total_sum;
	return nr;

}

void ibd::calcprob_pop(vector<int>& N_ext,
			  matrix<Genotype>& geno_augm,
			  vector<double>& P,
			  const LinkageMap& Markermap,
			  const matrix<ObsGeno>& obs_geno,
			  const PopulationType& poptype,
			  double threshold)
{
	const int nind = obs_geno.NrRows();
	geno_augm.clear();
	P.clear();
	vector<double> r = make_rec_map(Markermap);
	for (int ind=0;ind<nind;ind++)
		N_ext[ind] = calcprob_ind(geno_augm,P,r,obs_geno[ind],poptype,threshold);
}

vector<double> ibd::calc_P_stepII(const matrix<Genotype> geno_augm, 
							 const vector<double>& P,
							 const LinkageMap& markermap,
							 const PopulationType& popt,
							 const Locus& QTLpos)
{
	ObsGeno U = popt.unknown_genotype();
	int dimU = U.size();
	int nrow = P.size();
	vector<double> P2(dimU*nrow);
	const int left  = pos_qtl(markermap,QTLpos);
	const int right = left+1;
	Locus M_left  = markermap[left];
	Locus M_right = markermap[right];
	double r_left  = recomb(M_left,QTLpos);
	double r_right = recomb(QTLpos,M_right);

	vector<double> QTLprob(dimU);
	for (int row = 0; row < nrow; row++)
	{
		int k;
		double sum = 0.0;
		Genotype g_left  = geno_augm[row][left];
		Genotype g_right = geno_augm[row][right];
		for (k=0;k<dimU;k++)
		{
			Genotype QTL_prog = U[k];

			double p = 0.0;
			for (int k2=0;k2<dimU;k2++)
			{
				Genotype QTL_par = U[k2];
				p += popt.p_trans(g_left,QTL_par,r_left)*
					 popt.p_trans(QTL_par,g_right,r_right)*
					 popt.p_trans_QTL(QTL_par,QTL_prog);
			}
			QTLprob[k] = p;
			sum += p;
		}
		for (k=0;k<dimU;k++)
			P2[dimU*row + k] = P[row]*QTLprob[k]/sum;
	}
	return P2;
}

