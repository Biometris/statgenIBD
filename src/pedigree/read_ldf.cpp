// Martin Boer, Biometris

#include "read_ldf.h"
#include "manipulators.h"
#include "convert.h"
#include "eigen.h"

using namespace mbl;
using namespace std;

// Jan 21 2009: Check the correct format!!
void read_ldf_file(vector<int>& fndname,
				   LinkageMap& markermap,
				   matrix3D<double>& Q, const string filename)
{
	ifstream inp;
	OpenFile(inp,filename);

	int nind,npos;
	string line,keyword;
	inp >> keyword >> nind;
	if (keyword != "nind")
		throw mblib_error("Error!");
	skip_rest_of_line(inp);
	skip_lines(inp,1);

	vector<int> nr;
	read_columns(nr,fndname,nind,inp);
	skip_lines(inp,1);

	inp >> keyword >> npos;
	if (keyword != "npos")
		throw mblib_error("Error!");
	skip_rest_of_line(inp);
	skip_lines(inp,1);

	vector<int> posnr,chr;
	vector<double> pos;
	read_columns(posnr,chr,pos,npos,inp);
	skip_lines(inp,2);

	Q = matrix3D<double>(npos,nind,nind);
	for (int m=0;m<npos;m++)
	{
		markermap.push_back(Locus(chr[m],pos[m],EVAL_POS));
		skip_lines(inp,1);
		inp >> Q[m];
	}
}

vector< matrix<double> > decomposition_ldf_file(const LinkageMap& evalpos, string filename,
											  double eps)
{
	vector<int> fndname;
	matrix3D<double> Q;
	LinkageMap markermap_ldf;
	read_ldf_file(fndname,markermap_ldf,Q,filename);

	map<Locus, matrix<double> > Q_all;
	for (unsigned int m=0;m<markermap_ldf.size();m++)
		Q_all[markermap_ldf[m]] = Q[m];

	vector< matrix<double> > U_all;
	int cnt = 0;

	const double num_lim_eps = numeric_limits<double>::epsilon();

	for (unsigned int m=0;m<evalpos.size();m++)
	{
		map<Locus,matrix<double> >::const_iterator it=Q_all.find(evalpos[m]);
		if (it == Q_all.end())
		{
			string error = "marker " + stringify(m) + " not found in ldf_file";
			throw mblib_error(error);
		}
		matrix<double> Q_cur = it->second;
		//cout << setw(5) << m << setw(5) << evalpos[m].GetChr()
		//	 << setw(12) << evalpos[m].GetPosition() << endl;
		int dim = Q_cur.NrRows();

		vector<EigenReal> eigen = CalcEigen(Q_cur + num_lim_eps*identity_matrix(dim));
		int ncol=0;
		for (vector<EigenReal>::const_iterator it = eigen.begin();it!=eigen.end();it++)
			if (it->EigenValue() > eps)
				ncol++;
		//cout << setw(4) << m << setw(5) << ncol << endl;
		matrix<double> Ut(ncol,dim);
		for (int k=0;k<ncol;k++)
			Ut[k] = sqrt(eigen[k].EigenValue())*eigen[k].Eigenvector();
		U_all.push_back(transpose(Ut));
		cnt += ncol;
	}
	cout << "ncols spectral decomposition: " << cnt << endl;
	return U_all;
}
