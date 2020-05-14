// Martin Boer, Biometris
#include <set>
#include <fstream>
#include <sstream>

#include "misc.h"
#include "Loc.h"
#include "HMMalgo.h"
#include "convert.h"
#include "TransMatSym2D.h"

#include "popt.h"
#include "OrdGeno.h"
#include "InhVector.h"
#include "markerscore.h"
#include "analysis_fam.h"

using namespace mbl;
using namespace std;

IBD_fam::IBD_fam(const mbl::matrix<OrdGeno>& P,
			   const std::vector<score>& offspring,
			   const LinkageMap& MarkerMap,
			   const std::string& poptype)
				: markermap(MarkerMap)
{
	popt = init_pop(poptype);
	len_inh = popt->get_len();

	const unsigned int npar = P.NrRows();
	const unsigned int M = markermap.size();
	const unsigned int N = pow2(len_inh); // dimension of state space

	for (unsigned int i=0;i<npar;i++)
		par.push_back(OrdGeno(i,i));

	// init pi0
	vector<double> pi0(N,1.0/N);

	// init T
	vector<double> r = make_rec_map(markermap);

	int Nintervals = r.size();
	vector<TransMatSym2D> T(Nintervals); // Transition matrices between the loci.
	for (int k=0;k<Nintervals;k++)
		T[k] = TransMatSym2D(len_inh,r[k]);

	// init Q
	matrix<double> Q(M,N);
	for (unsigned int m=0;m<M;m++)
	{
		vector<OrdGeno> geno(npar);
		for (unsigned int par=0;par<npar;par++)
			geno[par] = P[par][m];
		Q[m] = check_scores(geno,offspring[m]);
	}

	// calculate left- and right- conditional probabilities.
	l_cond = calc_prob_left(pi0,Q,T);
	r_cond = calc_prob_right(Q,T);
}

std::vector<double> IBD_fam::check_scores(const std::vector<OrdGeno>& geno, const score& sc_off) const
{
	bool all_inconsistent = true;
	const unsigned int N = pow2(len_inh); // dimension of state space
	vector<double> q(N);
	for (InhVector u(len_inh);!u.end();u++)
	{
		OrdGeno g = popt->gen_off(geno,u);
		if (check_score(g,sc_off))
		{
			all_inconsistent = false;
			q[u] = 1.0;
		}
		else
			q[u] = 0.0;
	}
	if (all_inconsistent)
		return vector<double>(N,1.0);
	else
		return q;
}

std::map<OrdGeno,double> IBD_fam::operator()(const Locus& QTLpos) const
{
	const int left = pos_qtl(markermap,QTLpos);
	const int right = left+1;
	double r_left = recomb(markermap[left],QTLpos);
	double r_right = recomb(QTLpos,markermap[right]);
	const vector<double>& L = l_cond[left];
	const vector<double>& R = r_cond[right];
	TransMatSym2D T1(len_inh,r_left);
	TransMatSym2D T2(len_inh,r_right);
	vector<double> p = elem_prod(L*T1,T2*R);
	make_conditional(p);

	map<OrdGeno,double> IBDprob;
	for (InhVector u(len_inh);!u.end();u++)
	{
		OrdGeno g = popt->gen_off(par,u);
		IBDprob[g] += p[u];
	}
	return IBDprob;
}

IBD_fam make_IBD_biparental(const matrix<score>& sel_score, const LinkageMap& markermap,
						 int ind, int par1, int par2, string type)
{
	const int npar = 2; // biparental cross
	const int M = markermap.size();
	matrix<OrdGeno> par(npar,M);
	vector<score> offs(M);
	for (int m=0;m<M;m++)
	{
		score sc1 = sel_score[par1][m]; // par1
		score sc2 = sel_score[par2][m]; // par2
		if (sc1.homozygous() && sc2.homozygous())
		{
			// if there is no inheritance vector consistent with the scores
			// (e.g. P1: [1,1]	  P2: [2,2] -> offspring: (3,3)),
			// the constructor of IBD_fam will make the score non-informative
			// (see member function check_scores).
			par[0][m] = OrdGeno(sc1.first,sc1.first);
			par[1][m] = OrdGeno(sc2.first,sc2.first);
			offs[m] = sel_score[ind][m];
		}
		else
		{
			par[0][m] = OrdGeno(0,0);
			par[1][m] = OrdGeno(0,0);
			offs[m] = Uscore;
		}
	}
	return IBD_fam(par,offs,markermap,type);
}

matrix<double> calc_P_biparental(const matrix<score>& sel_score, const LinkageMap& markermap,
								 int ind, int par1, int par2, string type)
{
	IBD_fam IBD_ind = make_IBD_biparental(sel_score,markermap,ind,par1,par2,type);
	const int M = markermap.size();
	const int npar = 2; // biparental cross
	map<score,int> ndx = ndx_score(npar);
	matrix<double> P(M,ndx.size(),0.0);
	for (int m=0;m<M;m++)
	{
		map<OrdGeno,double> IBD = IBD_ind(markermap[m]);
		for (map<OrdGeno,double>::const_iterator it=IBD.begin();it!=IBD.end();it++)
		{
			OrdGeno g = it->first;
			score sc(g.first,g.second);
			int k = ndx[sc];
			P[m][k] += it->second;
		}
	}
	return P;
}

int analysis_families(matrix3D<double>& Z,         // output (Nind, M markers, Nfnd)
					  Pedigree& ped,               // output pedigree
					  const vector<IndProp>& pop,
					  const matrix<score>& geno,
					  const LinkageMap& markermap)
{
	cout << "start analysis families ..." << endl;
	set<string> Pfam = get_parents_families(pop);
	vector<IndProp> sel_pop;
	matrix<score> sel_score = select_families(sel_pop,geno,pop,Pfam);
	ped = Pedigree(sel_pop);
	const int N = ped.size();
	Z.clear();
	matrix<double> P;
	for (int i=0;i<N;i++)
	{
		ParentsInd ind = ped[i];
		if (!ind.IsFounder())
		{
			string type = sel_pop[i].GetType();
			P = calc_P_biparental(sel_score,markermap,i,ind[0],ind[1],type);
			Z.push_back(P);
		}
	}
	return 0;
}

/*
int test_analysis_F2_ind()
{
	cout.setf(ios::fixed, ios::floatfield);
	const double r1 = 0.20;
	const double r2 = 0.10;
	const double d1 = -50.0*log(1.0-2*r1);
	const double d2 = -50.0*log(1.0-2*r2);

	LinkageMap markermap;
	markermap.push_back(Locus(0,0.0,"M1"));
	markermap.push_back(Locus(0,d1,"M2"));
	markermap.push_back(Locus(0,d1+2*d2,"M3"));
	const int M = markermap.size();

	Locus QTLpos(0,d1+d2,"QTL");

	cout << setprecision(6) << endl;
	cout << "r M1 M2:  " << recomb(markermap[0],markermap[1]) << endl;
	cout << "r M2 M3:  " << recomb(markermap[1],markermap[2]) << endl;
	cout << "r M2 QTL: " << recomb(markermap[1],QTLpos) << endl;
	cout << "r QTL M3: " << recomb(QTLpos,markermap[2]) << endl << endl;

	//  "MapQTL"-score:       A     C     B
	istringstream inp_par1("(1,1) (1,1) (3,3)");
	istringstream inp_par2("(2,2) (2,2) (2,2)");
	//istringstream inp_offs("(2,2) (1,*) (3,2)");
	istringstream inp_offs("(1,1) (2,*) (2,2)");

	const int npar = 2; // biparental
	matrix<OrdGeno> par(npar,M);
	score sc;
	vector<score> offs(M);
	for (int m=0;m<M;m++)
	{
		sc = read_score(inp_par1);
		par[0][m].first = par[0][m].second = sc.first;
		sc = read_score(inp_par2);
		par[1][m].first = par[1][m].second = sc.first;
		offs[m] = read_score(inp_offs);
	}
	for (int i=0;i<npar;i++)
		cout << "par " << i << ": " << setw(7) << par[i] << endl;
	cout << "offs:  " << setw(7) << offs << endl;

	string poptype = "F2";
	cout << "poptype: "; cin >> poptype;
	IBD_fam tst(par,offs,markermap,poptype);

	cout << endl << endl;
	cout << "IBD: " << endl;
	map<OrdGeno,double> IBD = tst(QTLpos);
	for (map<OrdGeno,double>::const_iterator it = IBD.begin(); it!=IBD.end();it++)
		cout << setw(5) << it->first << setw(8) << setprecision(3) << it->second << endl;
	return 0;
}
*/

// This simple example shows that we can use the combination of the type
// OrdGeno and InheritanceVector for:
//  1) search for consistent marker scores
//  2) calculation of total IBD probabilities
int test_analysis_fam()
{
	int Ngen = 1;		  // # generations selfing
	int n = 2+2*Ngen + 1; // cross + Ngen selfing + DH
	int N = pow2(n);	  // dimension state space

	OrdGeno g1(0,1);
	OrdGeno g2(2,3);
	score sc(2,2);
	map<OrdGeno,int> distr;
	for (InhVector x(n); !x.end(); x++)
	{
		InhVector u = x; // make a copy of x
		bool b1 = u.next_indic();
		bool b2 = u.next_indic();
		OrdGeno g   = cross(g1,b1,g2,b2);
		OrdGeno gs  = selfing(g,u,Ngen);
		OrdGeno gDH = DH(gs,u);
		if (u.length() != 0)
			cout << "Warning: Length InheritanceVector u should be zero!" << endl;

		if (check_score(gDH,sc))
		{
			unsigned int v = x;
			cout << setw(8) << v << setw(8) << x << setw(10) << gDH << endl;
		}
		distr[gDH]++;
	}
	cout << endl << "List with IBDs (not using marker information): " << endl;
	for (map<OrdGeno,int>::const_iterator it=distr.begin();it!=distr.end();it++)
	{
		double IBD = it->second/(1.0*N);
		cout << setw(8) << it->first << setw(8) << IBD << endl;
	}
	return 0;
}

// march 3, 2007: this old implementation gives other results,
// probably because of the handling of missing values:
//    e.g. (AxB) -> C
//          * 1     2  --> old implementation concludes that 2 is not from B, thus from A
// So far, I don't use this information in the new implementation (i.e. I'm using a different
// analysis for families.

/*
#include "analysis_sub.h"

// output m x 2 matrix (2 founders)
matrix<double> analyse_ind_fam(int ind,const Pedigree& ped, matrix<int>& haplo,
							   const vector<double>& r, double eps)
{
	vector<ParentsInd> selpar(3);
	// 0 and 1 are defined as founders
	selpar[2] = ParentsInd(0,1);
	Pedigree ped_sel(selpar);

	matrix<int> haplo_sel(3,0);
	const int P1 = ped[ind][0];
	const int P2 = ped[ind][1];
	haplo_sel[0] = haplo[P1];
	haplo_sel[1] = haplo[P2];
	haplo_sel[2] = haplo[ind];
	cout << ind << endl;
	cout << setw(5) << haplo_sel[0] << endl;
	cout << setw(5) << haplo_sel[1] << endl;
	cout << setw(5) << haplo_sel[2] << endl;

	return calc_IBD_prob_last_ind(ped_sel,haplo_sel,r,eps);
}

int analysis_families(matrix3D<double>& Z,         // output (Nind, M markers, Nfnd)
					  Pedigree& ped,               // output pedigree
					  const vector<IndProp>& pop,
					  const matrix<score>& geno,
				      const LinkageMap& markermap)
{
	vector<double> r = make_rec_map(markermap);
	cout << "start analysis families ..." << endl;
	//const double eps = 0.0001;
	const double eps = 0.00;
	set<string> Pfam = get_parents_families(pop);

	vector<IndProp> sel_pop;
	matrix<int> haplo = gen_haplo_score_fam(sel_pop,geno,pop,Pfam);
	ped = Pedigree(sel_pop);

	const int N = ped.size();

	Z.clear();
	for (int i=0;i<N;i++)
	{
		if (!ped[i].IsFounder())
		{
			matrix<double> P = analyse_ind_fam(i,ped,haplo,r,eps);
			Z.push_back(P);
		}
	}
	return 0;
}
*/
