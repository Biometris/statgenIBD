#include <set>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "misc.h"
#include "HMMalgo.h"
#include "convert.h"
#include "TransMatSym2D.h"

#include "Loc.h"
#include "popt.h"
#include "OrdGeno.h"
#include "InhVector.h"
#include "markerscore.h"
#include "analysis_fam.h"
#include "crosses.h"
#include "mainR.h"

using namespace mbl;
using namespace std;

int count_parents(const vector<IndProp>& pop)
{
	int npar = 0;
	for (vector<IndProp>::const_iterator it=pop.begin();it!=pop.end();it++)
	{
		if (it->IsInbredParent())
			npar++;
	}
	return npar;
}

int count_progeny(const vector<IndProp>& pop)
{
	int nfam = 0;
	for (vector<IndProp>::const_iterator it=pop.begin();it!=pop.end();it++)
	{
		if (it->IsMemberFamily())
			nfam++;
	}
	return nfam;
}

IndProp find_first_progeny(const vector<IndProp>& pop)
{
	for (vector<IndProp>::const_iterator it=pop.begin();it!=pop.end();it++)
	{
		if (it->IsMemberFamily())
			return *it;
	}
	throw mblib_error("Cannot find progeny in functio find_first_progeny");
	return pop[0]; // dummy
}

string find_type(const vector<IndProp>& pop)
{
	return find_first_progeny(pop).GetType();
}

// print headers for cur position on the genome.
void print_cur_eval_pos(ostream& outp, const string pre,
						const LinkageMap& eval_pos, int m)
{
	Locus loc = eval_pos[m];
	outp << "# " << pre << ", evaluation point " << m+1
		 << ", chr " << loc.GetChr()
		 << ", " << loc.GetPosition() << "cM" << endl;
}

void print_progeny(ostream& outp, const vector<IndProp>& pop)
{
	const int nfam = count_progeny(pop);
	outp << "# Number of individuals: " << nfam << endl;
	outp << "# List with individuals: " << endl;
	outp << "indnr" << '\t' << "name" << endl;

	int cnt=1;
	for (vector<IndProp>::const_iterator it=pop.begin();it!=pop.end();it++)
	{
		if (it->IsMemberFamily())
			outp << cnt++ << '\t' << it->GetID() << '\n';
	}
}

void print_IBDs(ostream& outp, const matrix3D<double>& Z, const LinkageMap& eval_pos)
{
	int nfam = Z.Dim1();
	int M = Z.Dim2();
	for (int m=0;m<M;m++)
	{
		print_cur_eval_pos(outp,"IBD-probabilities ", eval_pos,m);
		for (int r=0;r<nfam;r++)
			write_tab_delimit_line(outp,Z[r][m]);
	}
}

void print_eval_pos(ostream& outp, const LinkageMap& eval_pos)
{
	const int M = eval_pos.size();
	outp << "# Number of markers or evaluation points: " << M << endl;
	outp << "posnr" << setw(5) << "chr" << setw(12) << "posmap" << '\n';
	for (int m=0;m<M;m++)
	{
		int chr = eval_pos[m].GetChr();
		double pos = eval_pos[m].GetPosition();
		outp << setw(5) << m+1 << setw(5) << chr << setw(12) << pos << '\n';
	}
}

void print_IBD_header(ostream& outp, int npar)
{
	outp << "#" << endl
		 << "# description of columns with IBD probabilities: " << endl;
	int c=1;
	for (int i=0;i<npar;i++)
	{
		outp << "#  column" << setw(3) << c++
			 << ": both alleles inherited from parent " << i+1 << endl;
	}
	if (npar == 2)
	{
		outp << "#  column" << setw(3) << c++ << ": heterozygous" << endl;
	}
	else if (npar == 3)
	{
		for (int i=0;i<2;i++)
		{
			outp << "#  column" << setw(3) << c++
			     << ": heterozygous, parents " << i+1 << " and 3" << endl;
		}
	}
	else if (npar == 4)
	{
		for (int i=0;i<2;i++)
		{
			for (int j=2;j<4;j++)
			{
				outp << "#  column" << setw(3) << c++
					 << ": heterozygous, parents " << i+1 << " and " << j+1 << endl;
			}
		}
	}
	outp << "# " << endl
		 << "# probability random allele of an individual derived from parent i, " << endl
         << "# where c[k] refers to the columns as defined above: " << endl
	     << "#" << endl;
	if (npar == 2)
	{
		outp << "#  P(parent=1) = c[1] + 0.5*c[3] " << endl
			 << "#  P(parent=2) = c[2] + 0.5*c[3] " << endl;
	}
	else if (npar == 3)
	{
		outp << "#  P(parent=1) = c[1] + 0.5*c[4] " << endl
			 << "#  P(parent=2) = c[2] + 0.5*c[5] " << endl
			 << "#  P(parent=3) = c[3] + 0.5*c[4] + 0.5*c[5] " << endl;
	}
	else if (npar == 4)
	{
		outp << "#  P(parent=1) = c[1] + 0.5*c[5] + 0.5*c[6]" << endl
			 << "#  P(parent=2) = c[2] + 0.5*c[7] + 0.5*c[8]" << endl
			 << "#  P(parent=3) = c[3] + 0.5*c[5] + 0.5*c[7]" << endl
			 << "#  P(parent=4) = c[4] + 0.5*c[6] + 0.5*c[8]" << endl;
	}
	outp << "# " << endl;
}

void print_header(ostream& outp,const vector<IndProp>& pop,
				  const vector<int>& ndx_par, const string& type)
{
	int npar = ndx_par.size();
	outp << "# Output of program pedigree.exe (v" << version << ")" << endl;
	outp << "#" << endl;
    switch (npar)
    {
		case 2:
			outp << "# analysis of biparental cross"  << endl;
			break;
		case 3:
            outp << "# analysis of three-way cross"  << endl;
            break;
		case 4:
            outp << "# analysis of four-way cross"  << endl;
            break;
    }
	outp << "# population type: " << type << endl;
	outp << "# parent list: " << endl;
	outp << "par" << '\t' << "name" << '\n';
	for (int i=0;i<npar;i++)
		outp << i+1 << '\t' << pop[ndx_par[i]].GetID() << '\n';
}

matrix<double> calc_P(const LinkageMap& eval_pos, int nparents, const IBD_fam& IBD_ind)
{
	const int M = eval_pos.size();
	map<score,int> ndx = ndx_score(nparents);
	matrix<double> P(M,ndx.size(),0.0);
	for (int m=0;m<M;m++)
	{
		map<OrdGeno,double> IBD = IBD_ind(eval_pos[m]);
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

matrix3D<double> calc_IBDs(const vector<IndProp>& pop,
						   const vector<int>& ndx_par,
						   matrix<score> geno,
						   const LinkageMap& markermap,
						   const LinkageMap& eval_pos,
						   const string& type)
{
	const int Nrow = geno.NrRows();
	const int npar = ndx_par.size();
	const int nloc = markermap.size();

	matrix<OrdGeno> par(npar,nloc);
	for (int m=0;m<nloc;m++)
	{
		bool parents_scores_OK = true;
		for (int i=0;i<npar;i++)
		{
			score sc = geno[ndx_par[i]][m];
			if (sc.homozygous())
				par[i][m] = OrdGeno(sc.first,sc.first);
			else
			{
				parents_scores_OK = false;
				par[i][m] = OrdGeno(0,0);
			}
		}
		if (!parents_scores_OK)
		{
			for (int i=0;i<Nrow;i++)
				geno[i][m] = Uscore;
		}
	}

	int r=0;
	matrix3D<double> Z;
	for (vector<IndProp>::const_iterator it=pop.begin();it!=pop.end();it++)
	{
		if (it->IsMemberFamily())
		{
			IBD_fam IBD_ind(par,geno[r],markermap,type);
			matrix<double> P = calc_P(eval_pos,npar,IBD_ind);
			Z.push_back(P);
		}
		r++;
	}
	return Z;
}

vector<int> get_ndx_par(const vector<IndProp>& pop)
{
	int npar = count_parents(pop);
	vector<int> ndx_par;
	IndProp ind = find_first_progeny(pop);
	string P1 = ind.GetP1();
	string P2 = ind.GetP2();
	int ndx_P1 = ndxID(pop,P1);
	int ndx_P2 = ndxID(pop,P2);
	IndProp par1 = pop[ndx_P1];
	IndProp par2 = pop[ndx_P2];
	if (npar == 2)
	{
		ndx_par.push_back(ndx_P1);
		ndx_par.push_back(ndx_P2);
	}
	else if (npar == 3)
	{
		if (par1.IsHybrid())
		{
			ndx_par.push_back(ndxID(pop,par1.GetP1()));
			ndx_par.push_back(ndxID(pop,par1.GetP2()));
			ndx_par.push_back(ndx_P2);
		}
		else
		{
			ndx_par.push_back(ndxID(pop,par2.GetP1()));
			ndx_par.push_back(ndxID(pop,par2.GetP2()));
			ndx_par.push_back(ndx_P1);
		}
	}
	else // npar == 4
	{
		ndx_par.push_back(ndxID(pop,par1.GetP1()));
		ndx_par.push_back(ndxID(pop,par1.GetP2()));
		ndx_par.push_back(ndxID(pop,par2.GetP1()));
		ndx_par.push_back(ndxID(pop,par2.GetP2()));
	}
	return ndx_par;
}

matrix3D<double> analysis_cross(const vector<IndProp>& pop,
					const matrix<score>& geno,
					const LinkageMap& markermap,
					const LinkageMap& eval_pos,
					//const string& filename,
					const Args& Argu)
{
	string tmp1;
	double tmp2;
	if (Argu.GetOption("bin") || Argu.GetOption("flx") ||
		Argu.GetOption("coa",tmp1) || Argu.GetOption("frac",tmp2))
			throw mblib_error("Analysis of family: one of the options cannot be used");
	int prec=3;
	Argu.GetOption("prec",prec);

	cout << "analysis of family ........" << endl;

	string type = find_type(pop);
	vector<int> ndx_par = get_ndx_par(pop);
	//const int npar = ndx_par.size();
	matrix3D<double> Z = calc_IBDs(pop,ndx_par,geno,markermap,eval_pos,type);

	//ofstream outp;
	//string outputfile = filename + ".txt";
	//OpenFile(outp,outputfile);
	//outp << setprecision(prec);
	//print_header(outp,pop,ndx_par,type);
	//print_progeny(outp,pop);
	//print_eval_pos(outp,eval_pos);
	//print_IBD_header(outp,npar);
	//print_IBDs(outp,Z,eval_pos);
	return Z;
}

