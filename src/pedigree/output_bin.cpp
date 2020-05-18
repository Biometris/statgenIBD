// // Martin Boer, Biometris
// #include "output_bin.h"
//
// using namespace std;
// using namespace mbl;
//
// void outp_bin_header(ostream& outp,
// 					 const LinkageMap& markermap,
// 					 const vector<IndProp>& pop)
// {
// 	int M = markermap.size();
// 	int Nfnd = Pedigree(pop).Nfounders();
// 	int Npar = get_parents_families(pop).size();
// 	int Nfam_ind = get_ind_within_fam(pop).size();
//
// 	write_bin(outp,M);		  // # markers
// 	write_bin(outp,Nfam_ind); // # ind within families (total number)
// 	write_bin(outp,Npar);	  // # Npar
// 	write_bin(outp,Nfnd);     // # fnd
//
// 	for (int m=0;m<M;m++)
// 	{
// 		double chr = markermap[m].GetChr();
// 		double pos = markermap[m].GetPosition();
//  		write_bin(outp,chr);
// 		write_bin(outp,pos);
// 	}
// }
//
//
// void outp_bin_fam(ostream& outp,const Pedigree& ped_fam, const matrix3D<double>& Z)
// {
// 	const int N = ped_fam.size();
// 	const int Nsel = Z.Dim1();
// 	const int M	   = Z.Dim2();
//
// 	for (int i=0;i<N;i++)
// 	{
// 		const ParentsInd& ind = ped_fam[i];
// 		if (!ind.IsFounder())
// 		{
// 			write_bin(outp, ind[0]+1); // ndx first parent
// 			write_bin(outp, ind[1]+1); // ndx second parent
// 		}
// 	}
//
// 	for (int m=0;m<M;m++)
// 	{
// 		for (int k=0;k<Nsel;k++)
// 		{
// 			write_bin(outp,Z[k][m][0]);
// 			//write_bin(outp,Z[k][m][1]);
// 		}
// 	}
// }
//
// void outp_bin_ped_ftp(ostream& outp,const vector<int>& ndx,
// 					    const IBDped& IBD,const Pedigree& ped)
//
// {
// 	int M = IBD.Nloc();
// 	vector<int> ndx_fnd = ped.GetNdxFnd();
// 	const int Nfnd = ndx_fnd.size();
// 	const int dim_ndx = ndx.size();
// 	for (int m=0;m<M;m++)
// 	{
// 		// prob. Founder to Parents
// 		const matrix<double> IBD_m = IBD(m);
// 		for (int k=0;k<dim_ndx;k++)
// 			for (int j=0;j<Nfnd;j++)
// 				write_bin(outp,IBD_m[ndx[k]][ndx_fnd[j]]);
// 	}
// }
//
// void outp_bin_ped_coa(ostream& outp,const vector<int>& ndx,
// 					const IBDped& IBD)
//
// {
// 	int M = IBD.Nloc();
// 	const int dim_ndx = ndx.size();
// 	for (int m=0;m<M;m++)
// 	{
// 		// coefficients of coancestry for sel ind.
// 		const matrix<double> IBD_m = IBD(m);
// 		for (int k=0;k<dim_ndx;k++)
// 			for (int j=0;j<dim_ndx;j++)
// 				write_bin(outp,IBD_m[ndx[k]][ndx[j]]);
// 	}
// }
