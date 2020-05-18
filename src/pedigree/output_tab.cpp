// // Martin Boer, Biometris
// #include "output_tab.h"
// #include "mainR.h"
// #include <algorithm>
//
// using namespace std;
// using namespace mbl;
//
// void outp_tab_header(ostream& outp,
// 					 const LinkageMap& markermap,
// 					 const vector<IndProp>& pop,
// 					 const set<string>& sel_ind)
// {
// 	outp << "# Output of program pedigree.exe (v" << version << ")" << endl;
// 	outp << "#" << endl;
// 	vector<string> fnd;
// 	int Nind_fam = 0;
// 	for (vector<IndProp>::const_iterator it=pop.begin();it!=pop.end();it++)
// 	{
// 		if (it->IsFounder())
// 			fnd.push_back(it->GetID());
// 		if (it->IsMemberFamily())
// 			Nind_fam++;
// 	}
// 	const int Nfnd = fnd.size();
// 	outp << "# Number of founders: " << Nfnd << endl;
// 	outp << "# Founder list" << endl;
// 	outp << "fndnr" << '\t' << "name" << endl;
// 	for (int i=0;i<Nfnd;i++)
// 		outp << i+1 << '\t' << fnd[i] << endl;
//
// 	vector<int> ndx = get_ndx_set(pop,sel_ind);
// 	int Nind = ndx.size();
// 	outp << "# Number of selected individuals (i.e. dimension of coa matrices): "
// 		 << Nind << endl;
// 	outp << "# Selected individuals list: " << endl;
// 	outp << "indnr" << '\t' << "name" << endl;
// 	for (int i=0;i<Nind;i++)
// 		outp << i+1 << '\t' << pop[ndx[i]].GetID() << endl;
// 	outp << "# Number of individuals in families: " << Nind_fam << endl;
// 	outp << "# List with individuals in families: " << endl;
// 	outp << "indnr" << '\t' << "name" << '\t' << "fam"
// 		 << '\t' << "type" << '\t' << "P1" << '\t' << "P2" << endl;
//
// 	int i=1;
// 	for (vector<IndProp>::const_iterator it=pop.begin();it!=pop.end();it++)
// 	{
// 		if (it->IsMemberFamily())
// 		{
// 			outp << i++ << '\t' << it->GetID()
// 						<< '\t' << it->GetFam()
// 						<< '\t' << it->GetType()
// 						<< '\t' << it->GetP1()
// 						<< '\t' << it->GetP2() << endl;
// 		}
// 	}
// 	int Neval = count_if(markermap.begin(),markermap.end(),eval_pos);
// 	outp << "# Number of markers or evaluation points: " << Neval << endl;
// 	outp << "# List with loci: " << endl;
// 	outp << "posnr" << setw(5) << "chr" << setw(12) << "posmap" << endl;
// 	int cnt=0;
// 	int M = markermap.size();
// 	for (int m=0;m<M;m++)
// 	{
// 		if (eval_pos(markermap[m]))
// 		{
// 			int chr = markermap[m].GetChr();
// 			double pos = markermap[m].GetPosition();
// 			outp << ++cnt << '\t' << chr << '\t' << pos << '\n';
// 			//outp << setw(5) << ++cnt << setw(5) << chr << setw(12) << pos << endl;
// 		}
// 	}
// }
//
// // print headers for cur position on the genome.
// void print_cur_pos(ostream& outp, const string pre, const Locus& loc, int cnt)
// {
// 	outp << "# " << pre << ", evaluation point " << cnt
// 		 << ", chr " << loc.GetChr()
// 		 << ", " << loc.GetPosition() << "cM" << endl;
// }
//
// void outp_tab_fam(ostream& outp,const Pedigree& ped_fam,
// 				  const matrix3D<double>& Z, const LinkageMap& markermap)
// {
// 	const int N = ped_fam.size();
// 	const int Nsel = Z.Dim1();
// 	const int M	   = Z.Dim2();
// 	outp << "# " << endl;
// 	outp << "# calculated IBD-probabilities for the families along the genome: " << endl;
// 	outp << "#   column 1: probability that both alleles are derived from first parent." << endl;
// 	outp << "#   column 2: probability that both alleles are derived from second parent." << endl;
// 	outp << "#   column 3: probability that individual is heterozygous." << endl;
// 	outp << "# " << endl;
// 	outp << "# parents of the families: "  << endl;
// 	outp << "indnr" << '\t' << "par1" << '\t' << "par2" << endl;
// 	int k=1;
// 	for (int i=0;i<N;i++)
// 	{
// 		const ParentsInd& ind = ped_fam[i];
// 		if (!ind.IsFounder())
// 		{
// 			outp << k++ << '\t' << ind[0] + 1 << '\t' << ind[1] + 1 << endl;
// 		}
// 	}
//
// 	int cnt=0;
// 	for (int m=0;m<M;m++)
// 	{
// 		if (eval_pos(markermap[m]))
// 		{
// 			print_cur_pos(outp,"probabilities families", markermap[m],++cnt);
// 			for (int k=0;k<Nsel;k++)
// 			{
// 				write_tab_delimit_line(outp,Z[k][m]);
// 			}
// 		}
// 	}
// }
//
// void outp_tab_ped_ftp(ostream& outp,const vector<int>& ndx,
// 					    const IBDped& IBD,
// 						const Pedigree& ped,
// 						const LinkageMap& markermap)
// {
// 	int M = IBD.Nloc();
// 	vector<int> ndx_fnd = ped.GetNdxFnd();
// 	const int Nfnd = ndx_fnd.size();
// 	const int dim_ndx = ndx.size();
// 	outp << "# " << endl;
// 	outp << "# Calculated founder to parents probabilities for positions " << endl
// 		 << "# along the genome" << endl;
// 	outp << "# " << endl;
//
// 	int cnt=0;
// 	for (int m=0;m<M;m++)
// 	{
// 		if (eval_pos(markermap[m]))
// 		{
// 			print_cur_pos(outp,"founder to parents", markermap[m],++cnt);
// 			matrix<double> IBD_m = IBD(m);
// 			for (int k=0;k<dim_ndx;k++)
// 			{
// 				for (int j=0;j<Nfnd;j++)
// 				{
// 					char delimiter = (j != Nfnd-1) ? '\t' : '\n';
// 					outp << IBD_m[ndx[k]][ndx_fnd[j]] << delimiter;
// 				}
// 			}
// 		}
// 	}
// }
//
// void outp_tab_ped_coa(ostream& outp,const vector<int>& ndx,
// 					const IBDped& IBD,
// 					const LinkageMap& markermap)
//
// {
// 	int M = IBD.Nloc();
// 	const int dim_ndx = ndx.size();
// 	outp << "# " << endl;
// 	outp << "# calculated coefficient of coancestry along the genome  " << endl;
// 	outp << "# " << endl;
//
// 	int cnt=1;
// 	for (int m=0;m<M;m++)
// 	{
// 		if (eval_pos(markermap[m]))
// 		{
// 			print_cur_pos(outp,"coefficients of coancestry",markermap[m],cnt++);
// 			matrix<double> IBD_m = IBD(m);
// 			for (int k=0;k<dim_ndx;k++)
// 			{
// 				for (int j=0;j<dim_ndx;j++)
// 				{
// 					char delimiter = (j != dim_ndx-1) ? '\t' : '\n';
// 					outp << IBD_m[ndx[k]][ndx[j]] << delimiter;
// 				}
// 			}
// 		}
// 	}
// }
