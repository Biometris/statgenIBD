// Martin Boer, Biometris
#include "output_flx.h"

using namespace std;
using namespace mbl;

void outp_flx_header(ostream& outp, 
					 const LinkageMap& markermap,
					 const vector<IndProp>& pop,
					 const set<string>& sel_ind)
{
	vector<int> ndx = get_ndx_set(pop,sel_ind);
	int Nind = ndx.size();
	outp << "nind " << Nind << endl;
    outp << "indnr    indname" << endl;
	for (int i=0;i<Nind;i++)
		outp << setw(5) << i+1 << setw(10) << pop[ndx[i]].GetID() << endl;
	outp << endl;

	int M = markermap.size();

	int n_eval_pos = 0;
	for (int m=0;m<M;m++)
		if (eval_pos(markermap[m]))
			n_eval_pos++;

	outp << "npos " << n_eval_pos << endl;
	outp << "posnr   chr   posmap" << endl;

	int k=1;
	for (int m=0;m<M;m++)
	{
		if (eval_pos(markermap[m]))
		{
			int chr = markermap[m].GetChr();
			double pos = markermap[m].GetPosition();
			outp << setw(4) << k++ << setw(5) << chr << setw(10) << pos << endl;
		}
	}
	outp << endl;
	outp << "inbred 1" << endl;
}

void outp_flx_ped_coa(ostream& outp,const vector<int>& ndx, 
					const IBDped& IBD,const LinkageMap& markermap)
{
	int M = IBD.Nloc();
	const int dim_ndx = ndx.size();
	outp << setprecision(2); 
	for (int m=0;m<M;m++)
	{
		if (eval_pos(markermap[m]))
		{
			// coefficients of coancestry for sel ind.
			outp << endl; // << endl;
			const matrix<double> IBD_m = IBD(m);
			for (int k=0;k<dim_ndx;k++)
			{
				for (int j=0;j<dim_ndx;j++)
				{
					outp << setw(5) << IBD_m[ndx[k]][ndx[j]];
				}
				outp << endl;
			}
		}
	}
}
