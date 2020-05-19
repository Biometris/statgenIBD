// Martin Boer, Biometris

#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <Rcpp.h>

// library files
#include "Args.h"
#include "excepthnd.h"
#include "convert.h"
#include "read_map.h"
#include "read_ldf.h"
#include "crosses.h"
#include "analysis_fam.h"
#include "mainR.h"
#include "Loc.h"
#include "matvec.h"
#include "InhVector.h"

using namespace std;
using namespace mbl;
using namespace Rcpp;

vector<IndProp> make_ped_file(const string& poptype,
                              const vector<string>& ID)
{
  // only to check poptype has correct format:
  const pop_base *popt = init_pop(poptype);

  vector<IndProp> result;
  int npar;
  if (poptype.find("C4") == 0)
    npar = 4;
  else if (poptype.find("C3") == 0)
    npar = 3;
  else
    npar = 2;
  // header, parents
  string par1, par2;
  for (int i=0; i<npar;i++)
  {
    IndProp tmp = IndProp(ID[i],"*","INBPAR","0","0");
    result.push_back(tmp);
  }
  if (npar == 4) {
    IndProp tmp1 = IndProp("H1","*","HYBRID",ID[0],ID[1]);
    IndProp tmp2 = IndProp("H2","*","HYBRID",ID[2],ID[3]);
    result.push_back(tmp1);
    result.push_back(tmp2);
    par1 = "H1";
    par2 = "H2";
  } else if (npar == 3)
  {
    IndProp tmp = IndProp("H","*","HYBRID",ID[0],ID[1]);
    result.push_back(tmp);
    par1 = "H";
    par2 = ID[2];
  } else
  {
    par1 = ID[0];
    par2 = ID[1];
  }
  for (unsigned int i=npar;i<ID.size();i++)
  {
    IndProp tmp = IndProp(ID[i],"ID_FAM",poptype,par1,par2);
    result.push_back(tmp);
  }
  return result;
}

void marker_selection(LinkageMap& markermap,
                      LinkageMap& eval_pos,
                      int sel_chr,
                      bool analysis_fam,
                      bool pos_option)
{
  markermap = adjust_markermap(markermap);
  if (sel_chr != -1)
  {
    eval_pos = select_chr(eval_pos,sel_chr);
    markermap = select_chr(markermap,sel_chr);
  }
  const int M = eval_pos.size();
  if (M < 1)
    throw mblib_error("no evaluation points!");
  if (analysis_fam)
  {
    if (markermap.empty())
    {
      markermap.push_back(eval_pos.front());
      markermap.push_back(eval_pos.back());
    }
    else if (eval_pos.front() < markermap.front())
    {
      markermap.push_back(eval_pos.front());
      sort(markermap.begin(),markermap.end());
    }
    if (eval_pos.back() > markermap.back())
      markermap.push_back(eval_pos.back());
  }
  else
  {
    if (pos_option)
    {
      const int M = eval_pos.size();
      for (int m=0;m<M;m++)
        markermap.push_back(eval_pos[m]);
      sort(markermap.begin(),markermap.end());
    }
  }
  if (markermap.size() == 1)
  {
    Locus loc = markermap[0];
    markermap.push_back(Locus(loc.GetChr()+1,loc.GetPosition(),EXTR_POS));
  }
}

int linking_data(matrix<score>& geno,
                 const LinkageMap& markermap,
                 const vector<IndProp>& pop,
                 const matrix<score>& geno_locfile,
                 const vector<string>& ID,
                 const vector<string>& marker_names)
{
  // put all data in right format
  const int M = markermap.size();
  const int N = pop.size();
  map<string,int> ID_ndx = make_index(ID);

  map<string,int> marker_names_ndx = make_index(marker_names);

  geno = matrix<score>(N,M,Uscore);
  for (int i=0;i<N;i++)
  {
    map<string,int>::const_iterator it=ID_ndx.find(pop[i].GetID());
    if (it != ID_ndx.end())
    {
      if (pop[i].IsHybrid())
      {
        Rcout << "!Warning: Genotypic data for " << pop[i].GetID()
             << " will be ignored" << endl;
      }
      else
      {
        int ndx_ind = it->second;
        for (int m=0;m<M;m++)
        {
          string locname = markermap[m].GetName();
          map<string,int>::const_iterator it = marker_names_ndx.find(locname);
          if (it != marker_names_ndx.end())
          {
            int ndx_M = it->second;
            geno[i][m] = geno_locfile[ndx_ind][ndx_M];
          }
        }
      }
    }
  }
  return 0;
}

int count_scores(const vector<score>& geno)
{
  int M = geno.size();
  int cnt = 0;
  for (int m=0;m<M;m++)
  {
    if (geno[m] != Uscore)
      cnt++;
  }
  return cnt;
}

void print_diagnostics(const vector<IndProp>& pop,
                       const matrix<score>& geno,
                       const string& filename)
{
  const int N = pop.size();
  ofstream outp;
  OpenFile(outp,filename);
  for (int i=0;i<N;i++)
  {
    IndProp ind = pop[i];
    int cnt = count_scores(geno[i]);
    int cntP1, cntP2;
    if (ind.IsFounder())
    {
      cntP1 = 0;
      cntP2 = 0;
    }
    else
    {
      cntP1 = count_scores(geno[ndxID(pop,ind.GetP1())]);
      cntP2 = count_scores(geno[ndxID(pop,ind.GetP2())]);
    }

    outp << setw(4) << i
         << setw(10) << ind.GetID()
         << setw(10) << ind.GetP1()
         << setw(10) << ind.GetP2()
         << setw(10) << cnt
         << setw(10) << cntP1
         << setw(10) << cntP2 << endl;
  }
}

matrix3D<double> set_IBD_founders(const matrix3D<double>& A,
                                  const LinkageMap& markermap,
                                  double alpha)
{
  Rcout << "alpha:  " << alpha << endl;

  const int M = markermap.size();
  const int Nfnd = A.Dim2();
  matrix3D<double> IBD_fnd(M,Nfnd,Nfnd);
  matrix<double> I_Nfnd = identity_matrix(Nfnd);
  int k=0;
  for (int m=0;m<M;m++)
  {
    if (eval_pos(markermap[m]))	// if current position is evaluation point.
      IBD_fnd[m] = alpha*A[k++] + (1.0-alpha)*I_Nfnd;
    else
      IBD_fnd[m] = I_Nfnd;	// actually not used.
  }
  return IBD_fnd;
}

// main part of the program (reading data, calculating IBD-prob, save results)
// Argu contains the values of the command line arguments and options
int main_pedigreeR(matrix3D<double>& Z,
                   vector<string>& parents,
                   vector<string>& offspring,
                   LinkageMap& positions,
                   const Args& Argu)
{
  // get arguments (i.e. get the filenames for input)
  string poptype,locfile,mapfile; //,output;
  Argu.GetArgument("poptype",poptype);
  Argu.GetArgument("locfile",locfile);
  Argu.GetArgument("mapfile",mapfile);

  // read all the data
  matrix<score> geno;
  LinkageMap markermap = read_map_file(mapfile);
  bool analysis_family = true; // single_family(pop);

  int sel_chr = -1; // default: analyse all chromosomes
  Argu.GetOption("chr",sel_chr);

  matrix<score> geno_locfile;
  vector<string> ID, marker_names;
  Rcout << "reading data .............." << endl;
  read_flapjackfile(ID,marker_names,geno_locfile,locfile);

  vector<IndProp> pop = make_ped_file(poptype, ID);

  markermap = reduce_markermap(markermap,marker_names);
  print_marker_warnings(markermap,marker_names);
  LinkageMap eval_pos;
  string eval_pos_file;
  bool pos_option = Argu.GetOption("pos",eval_pos_file);

  string max_step_size_str;
  bool mss = Argu.GetOption("max_step_size",max_step_size_str);
  double max_step_size = std::stod(max_step_size_str);
  if (pos_option)
  {
    if (Argu.GetOption("flx") || Argu.GetOption("bin"))
      throw mblib_error("cannot use option '-pos' in flexqtl (-flx) or binary mode (-bin)");
    eval_pos = read_eval_pos_file(eval_pos_file);
  }
  else if (max_step_size > 0)
    eval_pos = generate_extended_map(markermap, max_step_size);
  else
    eval_pos = markermap;

  // Count number of inbred founders:
  int Nfnd = 0;
  const int Npop = pop.size();
  for (int i=0;i<Npop;i++)
  {
    if (pop[i].IsFounder())
      Nfnd++;
  }
  string IBD_founder_file;
  matrix3D<double> IBD_fnd;
  vector<int> fndname; // assumes here that fndnames are integers!!
  if (Argu.GetOption("rdldf",IBD_founder_file))
  {
    matrix3D<double> IBD_fnd_eval_pos;
    if (pos_option)
      throw mblib_error("cannot use both options '-pos' and '-rdldf'");
    pos_option = true;

    // here further checks ar needed for checking right order fndname!
    eval_pos.clear(); // !?
    read_ldf_file(fndname,eval_pos,IBD_fnd_eval_pos,IBD_founder_file);

    // only simple check for correct number of founders:
    const int fndname_size = fndname.size();
    if (Nfnd != fndname_size)
      throw mblib_error("wrong number of founders in '-rdldf' option");

    marker_selection(markermap,eval_pos,sel_chr,analysis_family,pos_option);

    double alpha = 1.0;
    Argu.GetOption("alpha",alpha);
    IBD_fnd = set_IBD_founders(IBD_fnd_eval_pos,markermap,alpha);
  }
  else
  {
    marker_selection(markermap,eval_pos,sel_chr,analysis_family,pos_option);
    int M = markermap.size();
    IBD_fnd = matrix3D<double>(M,Nfnd,Nfnd);
    for (int m=0;m<M;m++)
      IBD_fnd[m] = identity_matrix(Nfnd);
  }

  linking_data(geno,markermap,pop,geno_locfile,ID,marker_names);

  string diag_file;
  if (Argu.GetOption("diag",diag_file))
    print_diagnostics(pop,geno,diag_file);
  // start of analysis and write results to output:
  Z = analysis_cross(pop,geno,markermap,eval_pos,Argu);

  const int npar = count_parents(pop);
  for (int i=0;i<npar;i++)
    parents.push_back(ID[i]);
  for (unsigned int i=npar;i<ID.size();i++)
    offspring.push_back(ID[i]);

  positions = eval_pos;

  return 0;
}

// This function for making connection with R/Rcpp
//  1. It reads the command-line and saves the arguments in class Args
//  2. Exception handling (see proc. exception_handler in files excepthnd.*)
int main_forR(matrix3D<double>& Z,
              vector<string>& parents,
              vector<string>& offspring,
              LinkageMap& positions,
              const string& poptype,
              const string& locfile,
              const string& mapfile,
              const string& eval_pos,
              const double& max_step_size)
{
  string max_step_size_str = stringify(max_step_size);
  int argc;
  if (eval_pos.length() > 0)
    argc = 8;
  else
    argc = 6;
  const char *argv[argc];
  argv[0] = "main_forR";
  argv[1] = poptype.c_str();
  argv[2] = locfile.c_str();
  argv[3] = mapfile.c_str();
  argv[4] = "-max_step_size";
  argv[5] = max_step_size_str.c_str();
  if (eval_pos.length() > 0)
  {
    argv[6] = "-pos";
    argv[7] = eval_pos.c_str();
  }
  const string arg = "poptype locfile mapfile";
  const string optional_arg = "";

  Options options;
  options.Add("H"				  , "help");
  options.Add("bin"             , "save all results in binary file");
  options.Add("flx"			  , "generate LD-input file for FlexQTL");
  options.Add("chr <int>"       , "select chromosome number");
  options.Add("coa <string>"    , "coa matrices for individuals defined in file");
  options.Add("prec <int>"      , "precision in tab-delimited output files");
  options.Add("frac <double>"   , "threshold for ancestors (frac. scored loci)");
  options.Add("pos <string>"    , "IBD-prob. at positions defined in file");
  options.Add("max_step_size <string>"   , "maximum distance for in between marker observations");
  options.Add("diag <string>"   , "some diagnostics, extra output");
  options.Add("nocoa"           , "write no coa matrices to output");
  options.Add("rdldf <string>"  , "read IBD probabilities for founders");
  options.Add("alpha <double>"  , "weighting factor for rdldf option");
  options.Add("eigen <string>"  , "use eigen-decomposition of ldf-file");
  options.Add("epseig <double>" , "precision for spectral decomposition");
  options.Add("savecoa <string>", "save coeff. of coancestry in file");
  options.Add("binformat <int> ", "format in binary file");

  Rcout << "Pedigree, version " << version << " (" << date << ")" << endl
       << "Martin Boer, Biometris; martin.boer@wur.nl" << endl << endl;

  Args Argu(options,arg,optional_arg,argc,argv);
  if (Argu.GetOption("H") || !Argu.State())
  {
    Rcout << "format: Pedigree <pedfile> <locfile> <mapfile>  [options]"
         << endl << endl << options << endl << endl;
    return 1;
  }
  int result = main_pedigreeR(Z, parents, offspring, positions, Argu);

  return result;
}

