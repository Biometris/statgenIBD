#include <string>
#include <RcppArmadillo.h>

#include "mainR.h"
#include "matvec.h"
#include "Loc.h"
#include "read_map.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include "sim.h"
#include "input.h"
#include "output.h"
#include "popt.h"
#include "util_genetics.h"
#include "dir.h"
#include "poptype/mapqtl.h"

using namespace Rcpp;
using namespace std;
using namespace ibd;

void sim_pops(const vector<PopProp>& pops, map<string,Genome>& genome_inbred_parents,
			  const LinkageMap& markermap, const vector<MarkerType>& markertype,
			  const Phi& phi, double sigma, const LinkageMap& QTLmap, int nrep)
{
	Rcout << endl << "Start of simulations ..... " << endl;
	const int width = 5; // maximum number of simulations: 00001 - 99999
	MakeLabel make_rep("rep",width);
	for (int r=0;r<nrep;r++)
	{
		if (nrep > 1) ChangeDir(make_rep(r));
		//if (GLOBAL_FLAPJACK)
		//{
		//	const string pop_name = "star_design"; // 24 march, for the moment this name
		//	SimPop sim_pop = sim_multiple_populations(pops,genome_inbred_parents);
		//	make_ped_file(sim_pop,pop_name);
		//	make_loc_file(sim_pop,markermap,markertype,pop_name);
		//	make_qua_file(sim_pop,phi,sigma,pop_name);
		//}
		//else
		//{
			for (vector<PopProp>::const_iterator it=pops.begin();it!=pops.end();++it)
			{
				string pop_name = it->GetName();
				SimPop sim_pop = sim_population(*it,genome_inbred_parents,pop_name);
				make_ped_file(sim_pop,pop_name);
				make_loc_file(sim_pop,markermap,markertype,pop_name);
				make_qua_file(sim_pop,phi,sigma,pop_name);
			}
		//}
		if (nrep > 1) ChangeDir("..");
	}
	Rcout << "end of simulations........" << endl;
}

map<string,Genome> make_genome_inbred_founders(const map<string,string>& inbfnds,
										vector<double> chr_length, const Phi& phi)
{
	int k=0;
	Rcout << endl << "Founders: " << endl; 
	map<string,Genome> genome_founders;
	typedef map<string,string>::const_iterator IterFnds;
	for (IterFnds iter=inbfnds.begin();iter!=inbfnds.end();++iter)
	{
		Genome fnd = Genome(chr_length,k,k);
		genome_founders[iter->first] = fnd;
		Rcout << setw(2) << k << setw(8) << iter->first 
			 << setw(5) << iter->second << setw(12) << phi(fnd) << endl;
		k++;
	}
	return genome_founders;
}

void print_coa_file(const vector<string>& parname, const string& filename)
{
	ofstream outp;
	OpenFile(outp,filename);
	for (vector<string>::const_iterator it=parname.begin();it!=parname.end();it++)
		outp << *it << endl;
}

//' simQTL
//'
//' @export
// [[Rcpp::export]]
int simQTL(CharacterVector& inputfile,
		   Nullable<CharacterVector&> dir_name = R_NilValue,
           const int& nrep = 1,
           const bool& print_flexQTL = false)
{
	cout.setf(ios::fixed, ios::floatfield);
	
	//Argu.GetArgument("scriptfile",inputfile);
	//Argu.GetOption("r",nrep);
	string inputfile_str = Rcpp::as<string>(inputfile);
	
	if (dir_name.isNotNull()) 
	{
	  string dir_name_str = Rcpp::as<string>(dir_name);
	  ChangeDir(dir_name_str);
	}    	
	
	int nr_alleles;
	long int start_seed;
	double fr_miss = 0.0;
	double mu,sigma2_e,dist_eval_pos;
	LinkageMap markermap;
	vector<double> chr_length;
	string eval_filename,mapfile;
	Commands commands = read_input_file(inputfile_str);
	read(mu,commands,"mu");
	read(sigma2_e,commands,"var");
	read(eval_filename,dist_eval_pos,commands,"dist");
	read_seed(start_seed,commands);
	read_genome(chr_length,commands);
    read_marker(markermap,mapfile,chr_length,nr_alleles,fr_miss,commands);

	map<Locus, vector<double> > QTLs = read_QTLs(commands);
	map<string,string> inbfnds = read_inbfnd(commands,QTLs.size());

	vector<string> fnd_names;
	for (map<string,string>::const_iterator it=inbfnds.begin();it!=inbfnds.end();it++)
		fnd_names.push_back(it->first);
	read_makeped(commands,fnd_names);

	vector<PopProp> pops = read_pop(commands);
	matrix<double> A = read_epi(commands,QTLs);

	vector<IndProp> ped;
	string pedfile,locfile,outputdir;
	bool pedigree_defined = read(pedfile,locfile,outputdir,commands,"pedigree");
	if (pedigree_defined)
		ped = read_ped_file(pedfile);

	vector<string> parfam;
	bool inbpar_defined = multi_read(parfam,commands,"inbpar");
	if (inbpar_defined)
		print_coa_file(parfam,"coa.txt");

	const int M = markermap.size();
	const int Nfnd = inbfnds.size();
	const int nQTL = QTLs.size();

	// Test output
	Rcout << "Test output: " << endl << endl;
	Rcout << "seed:    " << start_seed << endl;
	Rcout << "dist:    " << eval_filename << "  " << dist_eval_pos << endl;
	Rcout << "nr_all:  " << nr_alleles << endl;
	Rcout << "fr_mis:  " << fr_miss << endl;
	Rcout << "nloc:    " << M << endl << endl;
	Rcout << "mu:      " << mu << endl;
	Rcout << "var:     " << sigma2_e << endl;
	Rcout << "genome:  " << setw(12) << chr_length << endl << endl;
	Rcout << "Epistatic interactions: " << endl << setw(12) << A << endl;
	Rcout << endl << "Populations: " << endl;
	for (vector<PopProp>::const_iterator it=pops.begin();it!=pops.end();it++)
		Rcout << it->GetName() << endl;
	Rcout << endl << "QTLs: (nQTL = " << nQTL << ")" << endl ;
	print(Rcout,QTLs);

	LinkageMap QTLmap;
	typedef map<Locus,vector<double> >::const_iterator IterQTLs;
	for (IterQTLs iter=QTLs.begin();iter!=QTLs.end();++iter)
		QTLmap.push_back(iter->first);

	Phi phi(QTLs,inbfnds,A,mu);
	double sigma = sqrt(sigma2_e);	
	vector<MarkerType> markertype = generate_markertypes(M,Nfnd,nr_alleles,fr_miss);
	map<string,Genome> genome_founders = make_genome_inbred_founders(inbfnds,chr_length,phi);

	map<string,Genome> genome_par_fam;
	if (pedigree_defined)
	{
		// 18 febr 2008: Here we assume only 1 simulation for pedigree->parents families.
		genome_par_fam = sim_pedigree(ped,genome_founders);
		const int N = genome_par_fam.size();
		SimPop vecSimInd(N);
		int k=0;
		for (vector<IndProp>::const_iterator it=ped.begin();it!=ped.end();it++)
		{
			string id = it->GetID();
			vecSimInd[k++] = SimInd(*it,genome_par_fam.find(id)->second);
		}
		string cur_dir = get_current_dir();
		ChangeDir(outputdir);
 		make_loc_file(vecSimInd,markermap,markertype,locfile);
		ChangeDir(cur_dir);
	}
	else
		genome_par_fam = genome_founders;

	string inputfile_hybr,outputfile_hybr;
	if (read(inputfile_hybr,outputfile_hybr,outputdir,commands,"hybridsfile"))
	{
		Rcout << "outputdir: " << outputdir << " ...." << endl;
		sim_hybridsfile(inputfile_hybr,outputdir,outputfile_hybr,ped,genome_par_fam,
						markertype,markermap,phi,sigma,print_flexQTL);
	}

	// ?! check this
	//make_hybrids(commands,genome_par_fam,phi,sigma);

	//print_genstat_table(pops,inbfnds,"table.dat");	

	// write output and simulations:
	if (GLOBAL_FLAPJACK)
		make_flapjack_map_file(mapfile,markermap);
	else
		make_map_file(mapfile,markermap);
	make_eval_file(eval_filename,chr_length,dist_eval_pos);
	
	sim_pops(pops,genome_par_fam,markermap,markertype,phi,sigma,QTLmap,nrep);
	return 0;
}

/*
int main(int argc, char *argv[])
{
	const string arg = "scriptfile";
	const string optional_arg = "";
	Options options;
	options.Add("H",			"help");
	options.Add("flx",			"output flexQTL files hybrids");
	options.Add("r <int>",		"number of repetitions");
	options.Add("d <string>",   "working directory");
	Rcout << "SimQTL, version " << version << " (" << date << ")" << endl
		 << "Martin Boer, Biometris; martin.boer@wur.nl" << endl << endl;	
	Args Argu(options,arg,optional_arg,argc,argv);
	if (Argu.GetOption("H") || !Argu.State())
	{
		Rcout << "format: SimQTL <scriptfile> [opt]"
			 << endl << endl << options << endl << endl;
		return 1;
	}
	return exception_handler(main_sim,Argu);
}
*/