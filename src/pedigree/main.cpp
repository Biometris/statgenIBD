// Martin Boer, Biometris
// #define TEST_PEDIGREE
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>

// library files
#include "Args.h"
#include "excepthnd.h"
#include "convert.h"

#include "read_map.h"
#include "read_loc.h"
#include "read_ped.h"
#include "read_ldf.h"
#include "crosses.h"
#include "analysis_fam.h"
#include "analysis_ped.h"
#include "output.h"
#include "main.h"
#include "Loc.h"
#include "matvec.h"

#include "InhVector.h"

#include "eigen.h"


using namespace std;
using namespace mbl;

vector<IndProp> make_ped_file(const string& poptype, const vector<string>& ID)
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


void marker_selection(LinkageMap& markermap,LinkageMap& eval_pos,
                      int sel_chr, bool analysis_fam, bool pos_option)
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

void change_markernames(LinkageMap& markermap, bool pos_option)
{
  if (!pos_option)
  {
    const int M = markermap.size();
    for (int m=0;m<M;m++)
    {
      Locus loc = markermap[m];
      double pos = loc.GetPosition();
      int chr = loc.GetChr();
      string name = EVAL_POS + loc.GetName();
      markermap[m] = Locus(chr,pos,name);
    }
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
        cout << "!Warning: Genotypic data for " << pop[i].GetID()
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

void save_coa_pedigree_based(const Pedigree& ped,const vector<IndProp>& pop,
                             const vector<int>& ndx, const string filename)
{
  ofstream outp;
  OpenFile(outp,filename);
  matrix<double> coa = calc_coa(ped);
  outp << "# Output of program pedigree.exe (v" << version << ")" << endl;
  const int nind = ndx.size();
  outp << "nind = " << nind << endl;
  outp << setw(5) << "indnr" << setw(12) << "name" << endl;
  for (int i=0;i<nind;i++)
    outp << setw(5) << i+1 << setw(12) << pop[ndx[i]].GetID() << endl;
  outp << "# coefficients of coancestry matrix, pedigree based:" << endl;
  for (int i=0;i<nind;i++)
  {
    for (int j=0;j<nind;j++)
      outp << setw(12) << coa[ndx[i]][ndx[j]];
    outp << endl;
  }
}


// analysis of pedigree
void analysis_ped(const vector<IndProp>& pop,
                  const matrix<score>& geno,
                  LinkageMap markermap,
                  const matrix3D<double>& IBD_fnd,
                  const string& output,
                  const Args& Argu)
{
  string eval_pos_file;
  bool pos_option = (Argu.GetOption("pos",eval_pos_file) || Argu.GetOption("rdldf",eval_pos_file));
  //cout << "Pos_option " << pos_option << endl;
  change_markernames(markermap,pos_option);
  int M = count_if(markermap.begin(),markermap.end(),eval_pos);
  if (M == 0)
    throw mblib_error("Number of evaluation points empty!");
  //cout << "M: " << M  << endl;

  vector<double> r = make_rec_map(markermap);
  vector<double> r_inbred = make_rec_map_inbred_lines(markermap);

  // get selected individuals (default: parents of the families)
  set<string> sel_ind = get_parents_families(pop);
  string filename_coa;
  bool coa_parents_fam = !Argu.GetOption("coa",filename_coa);
  if (!coa_parents_fam)
  {
    sel_ind = get_ind_coa_file(pop,filename_coa);
  }

  // get option which defines the threshold for individuals
  // with relative high number of genotyped markers
  double frac = 0.80;
  Argu.GetOption("frac",frac);

  // Open output file and choose the type of output. There are three types:
  // 1. tab-delimited files (and comments starting with '#'): default option
  // 2. binary, used for communication with genstat: option -bin
  // 3. generate input file for flexqtl: option -flx
  ofstream outp;
  OutputType outp_type = open_output(outp,Argu,output,coa_parents_fam);

  // write 'header' information to output
  outp_header(outp,markermap,pop,sel_ind,outp_type);

  Pedigree ped;
  vector<int> ndx;
  matrix3D<double> IBD;

  if (outp_type != FLX && coa_parents_fam)
  {
    // analysis of families
    analysis_families(IBD,ped,pop,geno,markermap);
    // write probabilities to output
    outp_fam(outp,ped,IBD,markermap,outp_type);
  }

  // analysis of pedigree
  IBDped ibd_ped;
  analysis_pedigree(ibd_ped,ped,ndx,IBD_fnd,pop,geno,r_inbred,sel_ind,frac);

  string save_coa_file;
  if (Argu.GetOption("savecoa",save_coa_file))
  {
    save_coa_pedigree_based(ped,pop,ndx,save_coa_file);
  }

  cout << "writing output to file ...." << endl;

  // read and write ldf-file, using eigen-decomposition (april 7, 2009):
  string ldf_file;
  if (Argu.GetOption("eigen",ldf_file))
  {
    cout << "Reading ldf file and eigen decomposition ...." << endl;
    //cout << "markermap.size() " << markermap.size() << endl;

    LinkageMap eval_positions;
    for (unsigned int m=0;m<markermap.size();m++)
    {
      if (eval_pos(markermap[m]))
      {
        eval_positions.push_back(markermap[m]);
      }
    }

    int Neval_points = eval_positions.size();
    //cout << "evalsize: " << Neval_points << endl;

    vector<int> ndx_fnd = ped.GetNdxFnd();
    const int Nfnd = ndx_fnd.size();
    const int dim_ndx = ndx.size();
    //cout << "dim_ndx " << ndx.size() << endl;

    double eps = 1.0e-8;
    Argu.GetOption("epseig",eps);
    //cout << "eps: " << setprecision(20) << eps << endl;
    vector< matrix<double> > U = decomposition_ldf_file(eval_positions,ldf_file,eps);
    vector< matrix<double> > X;
    for (unsigned int m=0;m<markermap.size();m++)
    {
      if (eval_pos(markermap[m]))
      {
        matrix<double> X_cur(dim_ndx,Nfnd);
        // prob. Founder to Parents
        const matrix<double> IBD_m = ibd_ped(m);
        for (int k=0;k<dim_ndx;k++)
          for (int j=0;j<Nfnd;j++)
            X_cur[k][j] = IBD_m[ndx[k]][ndx_fnd[j]];
        X.push_back(X_cur);
      }
    }

    int binformat = 0;
    Argu.GetOption("binformat",binformat);

    string output;
    Argu.GetArgument("output",output);
    string outputfile = output + ".bin";
    ofstream outp;
    outp.open(outputfile.c_str(),ios::binary);
    double dbl_version = convertTo<double>(version);
    write_bin(outp,dbl_version);  // version of pedigree software
    write_bin(outp,binformat);    // format in binary file (see below)
    write_bin(outp,Neval_points); // nloc
    write_bin(outp,Nfnd);         // nfnd
    write_bin(outp,ndx.size());   // npar
    for (int m=0;m<Neval_points;m++) // mkchr
      write_bin(outp,eval_positions[m].GetChr());
    for (int m=0;m<Neval_points;m++) // mkpos
      write_bin(outp,eval_positions[m].GetPosition());
    for (int m=0;m<Neval_points;m++) // ncol
      write_bin(outp,U[m].NrCols());

    if (binformat == 0)
    {
      for (int m=0;m<Neval_points;m++) // write XU matrices
        write_bin(outp,X[m]*U[m]);
    }
    else if (binformat == 1)
    {
      for (int m=0;m<Neval_points;m++) // write U matrices
        write_bin(outp,U[m]);
      for (int m=0;m<Neval_points;m++) // write U matrices
        write_bin(outp,X[m]);
    }
    else
      throw mblib_error("binformat " + stringify(binformat) + " not defined!");
    return;
  }

  // write IBD-probabilities to file
  outp_ped_ftp(outp,ndx,ibd_ped,ped,markermap,outp_type);

  if (!Argu.GetOption("nocoa"))
    outp_ped_coa(outp,ndx,ibd_ped,markermap,outp_type);
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

void print_diagnostics(const vector<IndProp>& pop, const matrix<score>& geno, const string& filename)
{
  const int N = pop.size();
  // const int M = geno.NrCols();
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

matrix3D<double> set_IBD_founders(const matrix3D<double>& A,const LinkageMap& markermap, double alpha)
{
  cout << "alpha:  " << alpha << endl;

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
int main_pedigree(const Args& Argu)
{
  // get arguments (i.e. get the filenames for input and ouput)
  string poptype,locfile,mapfile,output;
  Argu.GetArgument("poptype",poptype);
  Argu.GetArgument("locfile",locfile);
  Argu.GetArgument("mapfile",mapfile);
  Argu.GetArgument("output",output);

  // read all the data
  matrix<score> geno;
  //vector<IndProp> pop = read_ped_file(pedfile);
  LinkageMap markermap = read_map_file(mapfile);
  bool analysis_family = true; // single_family(pop);

  int sel_chr = -1; // default: analyse all chromosomes
  Argu.GetOption("chr",sel_chr);

  matrix<score> geno_locfile;
  vector<string> ID, marker_names;
  cout << "reading data .............." << endl;
  read_flapjackfile(ID,marker_names,geno_locfile,locfile);

  vector<IndProp> pop = make_ped_file(poptype, ID);

  //read_loc_file(ID,marker_names,geno_locfile,locfile);
  markermap = reduce_markermap(markermap,marker_names);
  print_marker_warnings(markermap,marker_names);
  LinkageMap eval_pos;
  string eval_pos_file;
  bool pos_option = Argu.GetOption("pos",eval_pos_file);
  if (pos_option)
  {
    if (Argu.GetOption("flx") || Argu.GetOption("bin"))
      throw mblib_error("cannot use option '-pos' in flexqtl (-flx) or binary mode (-bin)");
    eval_pos = read_eval_pos_file(eval_pos_file);
  }
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

  //cout << "eval_pos.size(): " << eval_pos.size() << endl;
  //cout << "markermap.size(): " << markermap.size() << endl;

  //cout << "Nfnd: " << Nfnd << endl;

  string diag_file;
  if (Argu.GetOption("diag",diag_file))
    print_diagnostics(pop,geno,diag_file);

  matrix3D<double> Z;

  // start of analysis and write results to output:
  if (analysis_family)
    Z = analysis_cross(pop,geno,markermap,eval_pos,output,Argu);
  else
    analysis_ped(pop,geno,markermap,IBD_fnd,output,Argu);
  return 0;
}


// main part of the program (reading data, calculating IBD-prob, save results)
// Argu contains the values of the command line arguments and options
int main_pedigreeR(matrix3D<double>& Z,
                   vector<string>& parents,
                   vector<string>& offspring,
                   LinkageMap& positions,
                   const Args& Argu)
{
  // get arguments (i.e. get the filenames for input and ouput)
  string poptype,locfile,mapfile,output;
  Argu.GetArgument("poptype",poptype);
  Argu.GetArgument("locfile",locfile);
  Argu.GetArgument("mapfile",mapfile);
  Argu.GetArgument("output",output);

  // read all the data
  matrix<score> geno;
  //vector<IndProp> pop = read_ped_file(pedfile);
  LinkageMap markermap = read_map_file(mapfile);
  bool analysis_family = true; // single_family(pop);

  int sel_chr = -1; // default: analyse all chromosomes
  Argu.GetOption("chr",sel_chr);

  matrix<score> geno_locfile;
  vector<string> ID, marker_names;
  cout << "reading data .............." << endl;
  read_flapjackfile(ID,marker_names,geno_locfile,locfile);

  vector<IndProp> pop = make_ped_file(poptype, ID);

  //read_loc_file(ID,marker_names,geno_locfile,locfile);
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

  //cout << "eval_pos.size(): " << eval_pos.size() << endl;
  //cout << "markermap.size(): " << markermap.size() << endl;

  //cout << "Nfnd: " << Nfnd << endl;

  string diag_file;
  if (Argu.GetOption("diag",diag_file))
    print_diagnostics(pop,geno,diag_file);
  // start of analysis and write results to output:
  if (analysis_family)
    Z = analysis_cross(pop,geno,markermap,eval_pos,output,Argu);
  else
    analysis_ped(pop,geno,markermap,IBD_fnd,output,Argu);

  const int npar = count_parents(pop);
  for (int i=0;i<npar;i++)
    parents.push_back(ID[i]);
  for (unsigned int i=npar;i<ID.size();i++)
    offspring.push_back(ID[i]);

  positions = eval_pos;

  return 0;
}


// The main program is a kind of wrapper around main_pedigree:
//  1. It reads the command-line and saves the arguments in class Args
//  2. Exception handling (see proc. exception_handler in files excepthnd.*)
int main(int argc, char const *argv[])
{
#ifdef TEST_PEDIGREE
  //return exception_handler(tst_ndx_score);
  //return exception_handler(convert_ped_file);
  //return exception_handler(testInhVector);
  //return exception_handler(test_read);
  //return exception_handler(test_analysis_fam);
  return exception_handler(test_analysis_F2_ind);
#else
  const string arg = "poptype locfile mapfile output";
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
  options.Add("max_step_size <double>"  , "Maximum distance for in between marker observations");
  options.Add("diag <string>"   , "some diagnostics, extra output");
  options.Add("nocoa"           , "write no coa matrices to output");
  options.Add("rdldf <string>"  , "read IBD probabilities for founders");
  options.Add("alpha <double>"  , "weighting factor for rdldf option");
  options.Add("eigen <string>"  , "use eigen-decomposition of ldf-file");
  options.Add("epseig <double>" , "precision for spectral decomposition");
  options.Add("savecoa <string>", "save coeff. of coancestry in file");
  options.Add("binformat <int> ", "format in binary file");

  cout << "Pedigree, version " << version << " (" << date << ")" << endl
       << "Martin Boer, Biometris; martin.boer@wur.nl" << endl << endl;

  Args Argu(options,arg,optional_arg,argc,argv);
  
  if (Argu.GetOption("H") || !Argu.State())
  {
    cout << "format: Pedigree <pedfile> <locfile> <mapfile> <output> [options]"
         << endl << endl << options << endl << endl;
    return 1;
  }
  return exception_handler(main_pedigree,Argu);
#endif
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
			  const string& output,
			  const string& eval_pos,
			  const double& max_step_size)
{
  stringstream max_step_size_ss;
  max_step_size_ss << max_step_size;
  string max_step_size_str = max_step_size_ss.str();
  int argc;
  if (eval_pos.length() > 0)
	  argc = 9;
  else
	  argc = 7;
  const char *argv[argc];
  argv[0] = "main_forR";
  argv[1] = poptype.c_str();
  argv[2] = locfile.c_str();
  argv[3] = mapfile.c_str();
  argv[4] = output.c_str();
  argv[5] = "-max_step_size";
  argv[6] = max_step_size_str.c_str();
  if (eval_pos.length() > 0)
  {
	  argv[7] = "-pos";
	  argv[8] = eval_pos.c_str();
  }
  const string arg = "poptype locfile mapfile output";
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
  

  cout << "Pedigree, version " << version << " (" << date << ")" << endl
       << "Martin Boer, Biometris; martin.boer@wur.nl" << endl << endl;

  //  for (int i = 0; i<argc;i++)
  //  cout << argv[i] << endl;

  Args Argu(options,arg,optional_arg,argc,argv);
  if (Argu.GetOption("H") || !Argu.State())
  {
    cout << "format: Pedigree <pedfile> <locfile> <mapfile> <output> [options]"
         << endl << endl << options << endl << endl;
    return 1;
  }
  int result = main_pedigreeR(Z, parents, offspring, positions, Argu);

  return result;
}


// auxiliary function, not used at the moment
void eigenvalue_diagnostics(const IBDped& ibd_ped, const LinkageMap& markermap,
                            const vector<int>& ndx, int Napprox)
{
  Stopwatch time;
  ofstream outp;
  OpenFile(outp,"dimension.txt");
  cout << "Mtotal: " << markermap.size() << endl;
  const int K=ndx.size();
  cout << "ndx.size() : " << K << endl;
  for (unsigned int m=0;m<markermap.size();m++)
  {
    if (eval_pos(markermap[m]))
    {
      matrix<double> IBD_m = ibd_ped(m);
      matrix<double> M(K,K);
      for (int i=0;i<K;i++)
        for (int j=0;j<K;j++)
          M[i][j] = IBD_m[ndx[i]][ndx[j]];

      vector<EigenReal> eigen = CalcEigen(M);
      int nNeg = 0;
      int nPos = 0;
      int nZero = 0;
      const double eps = 1.0e-3;
      for (unsigned int i=0;i<eigen.size();i++)
      {
        if (eigen[i].EigenValue() < -eps)
          nNeg++;
        else if (eigen[i].EigenValue() >  eps)
          nPos++;
        else
          nZero++;
      }
      double sum = 0.0;
      for (int i=0;i<Napprox;i++)
        sum += eigen[i].EigenValue();
      outp << setw(4) << markermap[m].GetChr() << setw(12) << markermap[m].GetPosition()
           << setw(5) << nNeg << setw(5) << nZero << setw(5) << nPos
           << setw(12) << eigen[K-1].EigenValue() << setw(12) << eigen[0].EigenValue()
           << setw(12) << sum << endl;
      cout << setw(4) << m << setw(5) << nNeg << setw(5) << nZero
           << setw(5) << nPos << setw(12) << sum << endl;
    }
  }
  cout << "Total time: " << time <<  endl;
}

// aux, function, not used at the moment
// mean IBD calculated along the genome (using the evaluation points)
// N: number of individuals
// M: number of evaluation points.
void mean_IBD(const IBDped& ibd_ped, const LinkageMap& markermap, int N, int M)
{
  ofstream outp;
  OpenFile(outp,"mean_IBD.txt");
  Stopwatch time;
  cout << "M: " << M << endl;
  cout << "N: " << N << endl;
  matrix<double> A(N,N,0.0);

  for (unsigned int m=0;m<markermap.size();m++)
  {
    if (eval_pos(markermap[m]))
    {
      cout << setw(3) << markermap[m].GetChr() << setw(8) << markermap[m].GetPosition() << endl;
      A += ibd_ped(m);
    }
  }
  cout << "Total time: " << time << endl;
  for (int i=0;i<N;i++)
    for (int j=0;j<N;j++)
      A[i][j] /= (1.0*M);

  outp << setw(5) << setprecision(2) << A << endl;
}



