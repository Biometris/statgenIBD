#include "input.h"
#include "GenoValue.h"
#include "MarkerType.h"
#include "output.h"

#include "Sim.h"
#include "dir.h"
#include "Loc.h"
#include "Random.h"
#include "convert.h"
#include "urn.h"
#include "poptype/mapqtl.h"

using namespace ibd;
using namespace std;

Commands read_input_file(const string& filename)
{
  DefCommands defined_commands;
  AddCommand(defined_commands,"seed");
  AddCommand(defined_commands,"genome");
  AddCommand(defined_commands,"markers");
  AddCommand(defined_commands,"dist");
  AddCommand(defined_commands,"mu");
  AddCommand(defined_commands,"var");

  AddOptionalCommand(defined_commands,"makeped");
  AddOptionalCommand(defined_commands,"pedigree");

  // option hybrids is for a random selection of hybrids
  // option hybridsfile reads a file with defined hybrids
  AddOptionalCommand(defined_commands,"hybrids");
  AddOptionalCommand(defined_commands,"hybridsfile");

  AddOptionalCommand(defined_commands,"inbfndfile");

  const int max = 100;
  AddMultiCommand(defined_commands,"inbfnd",0,max);
  AddMultiCommand(defined_commands,"pop",0,max);
  AddMultiCommand(defined_commands,"qtl",0,max);
  AddMultiCommand(defined_commands,"epi",0,max); // two-way add x add interactions.

  AddMultiCommand(defined_commands,"inbpar",0,max); // inbred parents
  return read_command_file(defined_commands,filename);
}

// map<Locus, vector<double> > read_QTLs(const Commands& commands)
// {
//   map<Locus, vector<double> > QTLs;
//   typedef Commands::const_iterator Iter;
//   pair<Iter,Iter> range = commands.equal_range("qtl");
//   for (Iter iter=range.first;iter!=range.second;++iter)
//   {
//     const line_command& lc = iter->second;
//     istringstream line_stream(lc.argument);
//     string chr;
//     double pos;
//     string name;
//     vector<double> par(2);
//     line_stream >> name >> chr >> pos >> par;
//     //Locus loc(chr-1,pos,name);
//     Locus loc(chr,pos,name);
//     QTLs[loc] = par;
//   }
//   return QTLs;
// }

map<Locus, vector<double> > read_QTLs(const Rcpp::DataFrame& QTLposdf)
{
  Rcpp::NumericVector pos = QTLposdf["pos"];
  Rcpp::CharacterVector chr = QTLposdf["chr"];
  Rcpp::CharacterVector name = QTLposdf["name"];
  Rcpp::NumericVector add = QTLposdf["add"];
  Rcpp::NumericVector dom = QTLposdf["dom"];
  map<Locus, vector<double> > QTLs;
  for (int m=0;m<QTLposdf.nrows();m++)
  {
    string chrM = Rcpp::as<std::string>(chr[m]);
    string nameM = Rcpp::as<std::string>(name[m]);
    double posM = pos[m];
    vector<double> parM = {add[m], dom[m]};
    Locus loc(chrM,posM,nameM);
    QTLs[loc] = parM;
  }
  return(QTLs);
}

matrix<double> read_epi(const Commands& commands, const map< Locus,vector<double> >& QTLs)
{
  int nQTL = QTLs.size();
  matrix<double> A(nQTL,nQTL,0.0);
  typedef Commands::const_iterator Iter;
  pair<Iter,Iter> range = commands.equal_range("epi");
  for (Iter iter=range.first;iter!=range.second;++iter)
  {
    const line_command& lc = iter->second;
    istringstream line_stream(lc.argument);
    string qtl1,qtl2;
    double eff;
    line_stream >> qtl1 >> qtl2 >> eff;
    int k=0;
    int qtl1_nr = -1;
    int qtl2_nr = -1;
    typedef map<Locus,vector<double> >::const_iterator Iter;
    for (Iter it = QTLs.begin();it != QTLs.end();++it)
    {
      string name = it->first.GetName();
      if (name == qtl1)
        qtl1_nr = k;
      if (name == qtl2)
        qtl2_nr = k;
      k++;
    }
    if (qtl1_nr < 0 || qtl2_nr < 0 || qtl1_nr == qtl2_nr)
      lc.error("error while reading epistatic interaction");
    else
      A[qtl1_nr][qtl2_nr] = A[qtl2_nr][qtl1_nr] = eff;
  }
  return A;
}

map<string,string> read_inbfnd(const Commands& commands,unsigned int nqtl)
{
  map<string,string> inbfnd;
  typedef Commands::const_iterator Iter;
  pair<Iter,Iter> range = commands.equal_range("inbfnd");
  for (Iter iter=range.first;iter!=range.second;++iter)
  {
    const line_command& lc = iter->second;
    istringstream line_stream(lc.argument);
    string fnd_name,fnd_effects;
    line_stream >> fnd_name >> fnd_effects;
    if (!line_stream.eof() || line_stream.bad())
      lc.error("format error");
    if (fnd_effects.length() != nqtl)
      lc.error("wrong number of qtls");
    for (unsigned int i=0;i<nqtl;i++)
    {
      char c = fnd_effects[i];
      if (c != '+' && c != '-')
        lc.error("error in sign of qtl effects");
    }
    inbfnd[fnd_name] = fnd_effects;
  }

  string filename;
  if (!read(filename,commands,"inbfndfile"))
    return inbfnd;

  string line;
  ifstream inp;
  OpenFile(inp,filename);
  while (getline(inp,line))
  {
    if (line.empty()) continue;
    istringstream line_stream(line);
    string fnd_name,fnd_effects;
    line_stream >> fnd_name >> fnd_effects;

    if (fnd_effects.length() != nqtl)
      throw ibd_error("wrong number of qtls in file " + filename);
    for (unsigned int i=0;i<nqtl;i++)
    {
      char c = fnd_effects[i];
      if (c != '+' && c != '-')
        throw ibd_error("error in file " + filename);
    }
    inbfnd[fnd_name] = fnd_effects;
  }

  return inbfnd;
}


map<string,string> read_inbfnd(const Rcpp::DataFrame& inbfnddf,
                               unsigned int nqtl)
{
  map<string,string> inbfnd;
  Rcpp::CharacterVector name = inbfnddf["name"];
  Rcpp::DataFrame fnd_effects = inbfnddf;
  if (fnd_effects.size() != (nqtl + 1))
    throw ibd_error("wrong number of qtls");
  for (int m=0;m<fnd_effects.nrows();m++)
  {
    Rcpp::Rcout<< m << endl;

    string nameM = Rcpp::as<std::string>(name[m]);
    string fnd_effectsM = "";
    Rcpp::CharacterVector fnd_effectsI;
    for (unsigned int i=1;i<=nqtl;i++)
    {
      Rcpp::Rcout<< i << endl;

      fnd_effectsI = fnd_effects[i];

      Rcpp::Rcout << fnd_effectsI << endl;

      string s = Rcpp::as<std::string>(fnd_effectsI[m]);

      Rcpp::Rcout << s << endl;

      // if (s != "+" && s != "-")
      //   throw ibd_error("error in sign of qtl effects");
      fnd_effectsM = fnd_effectsM + s;
    }
    Rcpp::Rcout << fnd_effectsM << endl;

    inbfnd[nameM] = fnd_effectsM;
  }
  return inbfnd;
}






vector<PopProp> read_pop(const Commands& commands)
{
  vector<PopProp> pops;
  typedef Commands::const_iterator Iter;
  pair<Iter,Iter> range = commands.equal_range("pop");
  for (Iter iter=range.first;iter!=range.second;++iter)
  {
    const line_command& lc = iter->second;
    istringstream line_stream(lc.argument);
    string name,popt;
    line_stream >> name >> popt;
    int nind,npar;
    if (popt.find("C4") == 0)
      npar = 4;
    else if (popt.find("C3") == 0)
      npar = 3;
    else
      npar = 2;
    vector<string> P(npar);
    line_stream >> P >> nind;
    PopProp pop(nind,name,popt,P);
    if (!line_stream.eof() || line_stream.bad())
      lc.error("format error");
    pops.push_back(pop);
  }
  return pops;
}

// void read_seed(long int& start_seed, const Commands& commands)
// {
//   const line_command& lc = GetCommand(commands,"seed");
//   istringstream line_stream(lc.argument);
//   line_stream >> start_seed;
//   if (!line_stream.eof() || line_stream.bad())
//     lc.error("format error");
//
//   if (start_seed == 0)
//   {
//     randstart();
//     start_seed = get_randomseed();
//   }
//   else if (start_seed < 0)
//     set_randomseed(start_seed);
//   else
//     lc.error("Start value of random generator ran2 must be negative!");
// }

// void read_genome(vector<double>& chr_length, const Commands& commands)
// {
//   const line_command& lc = GetCommand(commands,"genome");
//   istringstream f(lc.argument);
//
//   char c;
//   f >> c;
//   if (c == '(')
//   {
//     for (;;)
//     {
//       double length;
//       f >> length >> c;
//       chr_length.push_back(length);
//       if (c != ',')
//       {
//         if (c != ')')
//           lc.error("Error while reading chromosome lengths");
//         else
//           break;
//       }
//     }
//   }
//   else
//   {
//     int nchr;
//     double length;
//     f.putback(c);
//     f >> nchr >> length;
//     chr_length = vector<double>(nchr,length);
//   }
//
//   eatcomment(f); // ?!
//   if (!f.eof() || f.bad())
//     lc.error("format error");
// }

// void read_number_markers(vector<int>& nr_markers_per_chr, istream& f)
// {
//   const int nchr = nr_markers_per_chr.size();
//   char c;
//   f >> c;
//   if (c == '(')
//   {
//     for (int i=0;i<nchr-1;i++)
//     {
//       f >> nr_markers_per_chr[i] >> c;
//       if (c != ',')
//         f.clear(ios::badbit);
//     }
//     f >> nr_markers_per_chr[nchr-1] >> c;
//     if (c != ')')
//       f.clear(ios::badbit);
//   }
//   else
//   {
//     int nr_loc;
//     f.putback(c);
//     f >> nr_loc;
//     for (int i=0;i<nchr;i++)
//       nr_markers_per_chr[i] = nr_loc;
//   }
//   eatcomment(f);
// }

void read_marker(LinkageMap& Markermap,
                 const string& filename,
                 vector<double> chr_length,
                 vector<double> nloc_chr,
                 const int& nr_alleles) //,
                 //double& fr_miss,
                 //const Commands& commands)
{
  //const line_command& lc = GetCommand(commands,"markers");
  //istringstream f(lc.argument);

  const int nchr = chr_length.size();
  //vector<int> nloc_chr(nchr);
  //f >> filename;
  //read_number_markers(nloc_chr,f);
  //if (f.bad())
//    lc.error("format error");
  //f >> nr_alleles;
  //if (f.bad())
//    lc.error("format error");
  //if (f.eof())
  //		fr_miss = 0.0;
  //else
  //		f >> fr_miss;

  //if (!f.eof() || f.bad())
//    lc.error("format error");

  for (int chr=0;chr<nchr;chr++)
  {
    double length = chr_length[chr];
    int nloc_cur_chr = nloc_chr[chr];
    if (nloc_cur_chr == 0)
      continue;
    vector<double> marker_pos(nloc_cur_chr);
    if (nloc_cur_chr == 1)
      marker_pos[0] = 0.5*length;
    else
    {
      marker_pos[0] = 0.0;
      double dist = length/(nloc_cur_chr-1);
      for (int m=1;m<nloc_cur_chr-1;m++)
        marker_pos[m] = dist*m;
      marker_pos[nloc_cur_chr-1] = length;
    }
    for (int m=0;m<nloc_cur_chr;m++)
    {
      string name = "M" + stringify(chr+1) + "_" + stringify(m+1);
      Markermap.push_back(Locus(stringify(chr),marker_pos[m],name));
    }
  }
}

vector<IndProp> read_pedigree(const Commands& commands)
{
  vector<IndProp> result;
  string filename="not defined";
  read(filename,commands,"pedigree");
  if (filename == "not defined")
    return result;
  return read_ped_file(filename);
}

vector<string> make_labels(const string& pre, int N, int width)
{
  MakeLabel label(pre,width);
  vector<string> result;
  for (int i=0;i<N;i++)
    result.push_back(label(i));
  return result;
}

void read_makeped(const Commands& commands, const vector<string>& fnd_names)
{
  int N, Ngen;
  string filename, coafile;

  if (!read(filename,N,Ngen,coafile,commands,"makeped"))
    return;
  sim_SS_NS_pedigree(fnd_names, filename, Ngen, N);

  vector<string> grp_A = make_labels("A" + stringify(Ngen) + "_",N,2);
  vector<string> grp_B = make_labels("B" + stringify(Ngen) + "_",N,2);

  ofstream outp1,outp2,outp3;
  OpenFile(outp1,coafile+".txt");
  OpenFile(outp2,coafile+"1.txt");
  OpenFile(outp3,coafile+"2.txt");

  for (vector<string>::const_iterator it=grp_A.begin();it!=grp_A.end();it++)
  {
    outp1 << *it << endl;
    outp2 << *it << endl;
  }
  for (vector<string>::const_iterator it=grp_B.begin();it!=grp_B.end();it++)
  {
    outp1 << *it << endl;
    outp3 << *it << endl;
  }

}


Genome get(const map<string,Genome>& pop, string ID)
{
  map<string,Genome>::const_iterator it = pop.find(ID);
  if (it == pop.end())
    throw ibd_error("cannot find ID " + ID + " in map");
  return it->second;
}

vector<string> get_strings(string filename)
{
  vector<string> result;

  ifstream inp;
  OpenFile(inp,filename);
  while (inp)
  {
    string tmp;
    inp >> tmp;
    if (tmp.empty())
      break;
    result.push_back(tmp);
  }
  return result;
}


void make_hybrids(const Commands& commands,
                  const map<string,Genome>& simpop,
                  const Phi& phi, double sigma)
{
  int N;
  string file1,file2,filename;
  if (!read(file1,file2,N,filename,commands,"hybrids"))
    return;

  vector<string> grpA = get_strings(file1);
  vector<string> grpB = get_strings(file2);
  cout << "Hybrids: " << endl
       << " size A: " << grpA.size() << endl
       << " size B: " << grpB.size() << endl << endl;

  vector<int> nr;
  const int Ntot = grpA.size()*grpB.size();
  for (int i=0;i<Ntot;i++)
    nr.push_back(i);
  Urn<int> urn(nr);
  vector<int> indicator(Ntot,0);
  for (int i=0;i<N;i++)
    indicator[urn.random_draw()] = 1;

  ofstream outp;
  OpenFile(outp,filename);

  outp <<  "b" << setw(12)  << "ID" << setw(8)  << "P1!" << setw(8)  << "P2!"
       << setw(12) << "pheno" << setw(12) << "geno" << endl;

  int k=0;
  MakeLabel hybr_label("HYBR",4);
  ///typedef pair<string,string> str_pr;
  //vector<str_pr> par_comb; // all parent combinations
  for (vector<string>::const_iterator itA=grpA.begin();itA!=grpA.end();itA++)
  {
    for (vector<string>::const_iterator itB=grpB.begin();itB!=grpB.end();itB++)
    {
      Genome g1 = get(simpop,*itA);
      Genome g2 = get(simpop,*itB);
      Genome g = g1*g2;

      double geno_val = phi(g);
      double error = randnormal(0.0,sigma);
      double pheno_val = geno_val + error;

      outp << indicator[k] << setw(12)  << hybr_label(k)
           << setw(8)  << *itA << setw(8)  << *itB
           << setw(12) << pheno_val  << setw(12) << geno_val << endl;

      k++;
    }
  }
}


// void make_hybrids_old(const Commands& commands,
//                       const map<string,Genome>& simpop,
//                       const Phi& phi, double sigma)
// {
//   int N;
//   string file1,file2,filename;
//   if (!read(file1,file2,N,filename,commands,"hybrids"))
//     return;
//
//   vector<string> grpA = get_strings(file1);
//   vector<string> grpB = get_strings(file2);
//   cout << "Hybrids: " << endl
//        << " size A: " << grpA.size() << endl
//        << " size B: " << grpB.size() << endl << endl;
//
//   ofstream outp;
//   OpenFile(outp,filename);
//
//   typedef pair<string,string> str_pr;
//   vector<str_pr> par_comb; // all parent combinations
//   for (vector<string>::const_iterator itA=grpA.begin();itA!=grpA.end();itA++)
//     for (vector<string>::const_iterator itB=grpB.begin();itB!=grpB.end();itB++)
//       par_comb.push_back(make_pair(*itA,*itB));
//   Urn<str_pr> urn(par_comb);
//
//   outp << setw(8)  << "ID" << setw(8)  << "P1!" << setw(8)  << "P2!"
//        << setw(12) << "pheno" << setw(12) << "geno" << endl;
//
//   MakeLabel hybr_label("HYBR",4);
//   for (int i=0;i<N;i++)
//   {
//     str_pr h = urn.random_draw();
//     Genome g1 = get(simpop,h.first);
//     Genome g2 = get(simpop,h.second);
//     Genome g = g1*g2;
//
//     double geno_val = phi(g);
//     double error = randnormal(0.0,sigma);
//     double pheno_val = geno_val + error;
//
//     outp << setw(8)  << hybr_label(i)
//          << setw(8)  << h.first << setw(8)  << h.second
//          << setw(12) << pheno_val  << setw(12) << geno_val << endl;
//   }
// }


// april 11: use the following two functions for generating hybrids:

vector<IndProp> read_hybrids(const Commands& commands)
{
  int N;
  string file1,file2,filename;
  if (!read(file1,file2,N,filename,commands,"hybrids"))
    return vector<IndProp>();

  vector<string> grpA = get_strings(file1);
  vector<string> grpB = get_strings(file2);
  cout << "Hybrids: " << endl
       << " size A: " << grpA.size() << endl
       << " size B: " << grpB.size() << endl << endl;

  typedef pair<string,string> str_pr;
  vector<str_pr> par_comb; // all parent combinations
  for (vector<string>::const_iterator itA=grpA.begin();itA!=grpA.end();itA++)
    for (vector<string>::const_iterator itB=grpB.begin();itB!=grpB.end();itB++)
      par_comb.push_back(make_pair(*itA,*itB));
  Urn<str_pr> urn(par_comb);

  vector<IndProp> result(N);
  MakeLabel hybr_label("HYBR",4);
  for (int i=0;i<N;i++)
  {
    str_pr h = urn.random_draw();
    result[i] = IndProp(hybr_label(i),"*","*",h.first,h.second);
  }
  return result;
}


void make_hybrids(const vector<IndProp>& hybr,const map<string,Genome>& simpop,
                  string filename, const Phi& phi, double sigma)
{
  ofstream outp;
  OpenFile(outp,filename);

  outp << setw(8)  << "ID" << setw(8)  << "P1!" << setw(8)  << "P2!"
       << setw(12) << "pheno" << setw(12) << "geno" << endl;

  MakeLabel hybr_label("HYBR",4);
  const int N = hybr.size();
  for (int i=0;i<N;i++)
  {
    IndProp ind = hybr[i];
    Genome g1 = get(simpop,ind.GetP1());
    Genome g2 = get(simpop,ind.GetP2());
    Genome g = g1*g2;

    double geno_val = phi(g);
    double error = randnormal(0.0,sigma);
    double pheno_val = geno_val + error;

    outp << setw(8)  << ind.GetID()
         << setw(8)  << ind.GetP1() << setw(8)  << ind.GetP2()
         << setw(12) << pheno_val  << setw(12) << geno_val << endl;
  }
}

void print_FlexQTL_ind(ofstream& outp, const IndProp& indprop, const Genome& genome,
                       const vector<MarkerType>& markertype,
                       const LinkageMap& markermap, double phenotype, double prob_SS_P1, int sparse_map)
{
  const int M = markermap.size();
  const string ID = indprop.GetID();
  const string P1 = indprop.GetP1();
  const string P2 = indprop.GetP2();

  if (ID.find("G") == 0) // assuming here that starting with "G" are extra simulated generations
    return;

  if (P1.find("G")==0 && P2.find("G")==0)  // intermediate founder:
  {
    outp << setw(4) << "1000" << setw(12) << ID
         << setw(12) << "0" << setw(12) << "-1000"
         << setw(5) << "1" << setw(5) << "*" << setw(12) << "*";
  }
  else if (indprop.IsRIL())
  {
    string F1_ind = "F1" + ID;
    // print F1 individual
    outp << setw(4) << "999" << setw(12) << F1_ind
         << setw(12) << indprop.GetP1() << setw(12) << indprop.GetP2()
         << setw(5) << "0" << setw(5) << "*" << setw(12) << "*";
    for (int m=0;m<M;m++)
      if (m % sparse_map == 0)
        outp << setw(5) << "-" << setw(2) << "-";
      outp << endl;
      // print offspring F1_ind
      outp << setw(4) << "1" << setw(12) << ID
           << setw(12) << F1_ind << setw(12) << "-777"
           << setw(5) << "1" << setw(5) << "*" << setw(12) << "*";
  }
  else // hybrid
  {
    int P1_SS = (prob_SS_P1 > 0.5) ? 1 : 0;

    outp << setw(4) << "100" << setw(12) << ID
         << setw(12) << indprop.GetP1() << setw(12) << indprop.GetP2()
         << setw(5) << "1" << setw(5) << P1_SS << setw(12) << phenotype;
  }

  vector<Genotype> geno = genome.GetGenotype(markermap);
  for (int m=0;m<M;m++)
  {
    if (m % sparse_map == 0)
    {
      Genotype g = geno[m];
      outp << setw(5) << markertype[m].First(g) << setw(2) << markertype[m].Second(g);
    }
  }
  outp << endl;
}

void sim_hybridsfile(const string inputfile,
                     const string outputdir,
                     const string outputfile,const vector<IndProp>& ped,
                     const map<string,Genome>& simpop,
                     const vector<MarkerType>& markertype,
                     const LinkageMap& markermap,
                     const Phi& phi, double sigma, bool FlexQTL)
{
  const int sparse_map = 1; // for flexQTL, marker 0,5,10 etc.
  const int M = markermap.size();

  ifstream inp;
  ofstream outp,outp2;

  OpenFile(inp,inputfile);

  ChangeDir(outputdir);
  OpenFile(outp,outputfile);

  if (FlexQTL)
  {
    LinkageMap flexqtl_map;
    for (int m=0;m<M;m++)
    {
      if (m % sparse_map == 0)
        flexqtl_map.push_back(markermap[m]);
    }
    make_map_file("flexqtl",flexqtl_map);
    OpenFile(outp2,"flexqtl.dat");
  }

  skip_lines(inp,1); // header
  int line_nr = 0;
  string line;

  const int nQTL = phi.NumberOfQTL();

  outp << setw(12)  << "ID" << setw(12)  << "TYPE"
       << setw(12) << "P1!" << setw(12)  << "P2!"
       << setw(12) << "pheno" << setw(12) << "geno";

  for (int q=0;q<nQTL;q++)
  {
    string qtl = "Q[" + stringify(q+1) + "]!" ;
    outp << setw(8) << qtl;
  }
  for (int q=0;q<nQTL;q++)
  {
    string str = "F[" + stringify(q+1) + "]!";
    outp << setw(8) << str;
  }
  for (int q=0;q<nQTL;q++)
  {
    string str = "M[" + stringify(q+1) + "]!";
    outp << setw(8) << str;
  }
  outp << endl;

  for (vector<IndProp>::const_iterator it=ped.begin();it!=ped.end();it++)
  {
    Genome g=get(simpop,it->GetID());
    double geno_val = phi(g);
    outp << setw(12) << it->GetID()
         << setw(12) << it->GetType()
         << setw(12) << it->GetP1()
         << setw(12) << it->GetP2()
         << setw(12) << "*"
         << setw(12) << geno_val
         << setw(8) << phi.getQTLs(g)
         << setw(8) << phi.getORG(g,Father)
         << setw(8) << phi.getORG(g,Mother) << endl;

    if (FlexQTL)
      print_FlexQTL_ind(outp2,*it,g,markertype,markermap,-10e10,0.0,sparse_map);
  }

  vector<double> hybrids_error, hybrids_pheno, hybrids_geno;

  ofstream outp3;
  OpenFile(outp3,"sim_hybr.qua");
  outp3 << setw(4) << "nr" << setw(12) << "ID"
        << setw(12) << "P1" << setw(12) << "P2"
        << setw(12) << "probSS_P1" << setw(12) << "probSS_P2"
        << setw(12) << "y" << endl;

  int hybr_nr = 1;
  while (getline(inp,line))
  {
    line_nr++;
    if (line.empty()) continue;
    istringstream line_stream(line);
    string ID,P1,P2;
    line_stream >> ID >> P1 >> P2;

    Genome g1 = get(simpop,P1);
    Genome g2 = get(simpop,P2);
    Genome g = g1*g2;

    double geno_val = phi(g);
    double error = randnormal(0.0,sigma);
    double pheno_val = geno_val + error;

    outp << setw(12) << ID << setw(12) << "HYBRID"
         << setw(12) << P1 << setw(12) << P2
         << setw(12) << pheno_val  << setw(12) << geno_val
         << setw(8) << phi.getQTLs(g)
         << setw(8) << phi.getORG(g,Father)
         << setw(8) << phi.getORG(g,Mother) << endl;

    double allele_QTL_P1 = phi.getORG(g,Father)[0];
    double allele_QTL_P2 = phi.getORG(g,Mother)[0];

    double probSS_P1,probSS_P2;
    if (allele_QTL_P1 > allele_QTL_P2)
    {
      probSS_P1 = 1.0; probSS_P2 = 0.0;
    }
    else
    {
      probSS_P1 = 0.0; probSS_P2 = 1.0;
    }

    outp3 << setw(4) << hybr_nr
          << setw(12) << ID << setw(12) << P1 << setw(12) << P2
          << setw(12) << probSS_P1 << setw(12) << probSS_P2
          << setw(12) << pheno_val << endl;

    IndProp prop = IndProp(ID,"*","HYBRID",P1,P2);

    if (FlexQTL)
      print_FlexQTL_ind(outp2,prop,g,markertype,markermap,pheno_val,probSS_P1, sparse_map);

    hybr_nr++;

    hybrids_error.push_back(error);
    hybrids_pheno.push_back(pheno_val);
    hybrids_geno.push_back(geno_val);
  }

  double sigma2_g = variance(hybrids_geno);
  double sigma2_p = variance(hybrids_pheno);

  cout << endl;
  cout << "sim sigma^2_y:    " << variance(hybrids_pheno) << endl;
  cout << "sim sigma^2_g:    " << variance(hybrids_geno) << endl;
  cout << "sim sigma^2_e:    " << variance(hybrids_error) << endl;
  cout << "sim heritability: " << sigma2_g/sigma2_p << endl << endl;
}
