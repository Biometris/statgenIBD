#include <float.h> // for DBL_MIN
#include "mapqtl.h"
#include "util_genetics.h"
#include "misc.h"
#include "convert.h"

using namespace std;

namespace {
using namespace ibd;

string MapQTLvar(istream& inp)
{
  char c;
  string x = "xxxx";
  eatcomment(inp);
  for (int i=0;i<4;i++)
    inp >> x[i];
  tolower(x);
  inp >> eatcomment >> c;
  if (c != '=' || !inp)
    throw ibd_error("MapQTLvar" + x + " " + c);
  return x;
}


void read_loc_header(string& name, string& popt, int& nind, int& nloc, istream& inp)
{
  popt = "";
  name = "";
  nind = -1;
  nloc = -1;
  for (int i=0;i<4;i++)
  {
    string var = MapQTLvar(inp);
    if (var == "popt")
      inp >> eatcomment >> popt;
    else if (var == "name")
      inp >> eatcomment >> name;
    else if (var == "nind")
      inp >> eatcomment >> nind;
    else if (var == "nloc")
      inp >> eatcomment >> nloc;
    else
      throw ibd_error("unknown variable " + var);
  }
  if (nind <=0 || nloc <= 0 || popt == "" || name == "")
    throw ibd_error("read_loc_header");
}

}

string ibd::read_qua_header(vector<string>& column_name, int& nind, int& ntrt, istream& inp)
{
  string miss = "";
  nind = -1;
  ntrt = -1;
  for (int i=0;i<3;i++)
  {
    string var = MapQTLvar(inp);
    if (var == "nind")
      inp >> eatcomment >> nind;
    else if (var == "ntrt")
      inp >> eatcomment >> ntrt;
    else if (var == "miss")
      inp >> eatcomment >> miss;
    else
      throw ibd_error("unknown variable " + var);
  }
  string str;
  for (int j=0;j<ntrt;j++)
  {
    inp >> eatcomment >> str;
    column_name.push_back(str);
  }
  if (nind <=0 || ntrt <= 0 || miss == "")
    throw ibd_error("read_qua_header");

  return miss;
}

void ibd::make_map_file(string name, const LinkageMap& marker_map)
{
  string cur_chr = "-1";
  string filename = name;
  ofstream outp(filename.c_str());
  if (!outp)
    throw ibd_error("Cannot open file " + filename);
  outp.setf(ios::fixed, ios::floatfield);
  outp.precision(4);
  int nloc = marker_map.size();
  for (int m=0; m < nloc; m++)
  {
    const Locus& loc = marker_map[m];
    if (cur_chr != loc.GetChr())
    {
      cur_chr = loc.GetChr();
      // 20 december: cur_chr --> cur_chr + 1
      outp << "chrom " << stringify(std::stoi(cur_chr) + 1) << endl;
    }
    outp << loc.GetName() << '\t' << loc.GetPosition() << endl;
  }
}

void ibd::make_pheno_file(string name,
                          const vector<double>& y,
                          const vector<double>& w,
                          const vector<int>& env_factor)
{
  string filename = name + "_pheno.txt";
  ofstream outp(filename.c_str());
  if (!outp)
    throw ibd_error("Cannot open file " + filename);
  outp.setf(ios::fixed, ios::floatfield);
  int nind = y.size();
  outp << "ntrt = 4" << endl;
  outp << "nind = "  << nind << endl;
  outp << "miss = *" << endl << endl;
  outp << "nr" << endl;
  outp << "env_factor" << endl;
  outp << "phenotype" << endl ;
  outp << "weight" << endl << endl;
  for (int ind=0; ind < nind; ind++)
    outp << setw(4) << ind
         << setw(4) << env_factor[ind]
         << setw(12) << y[ind]
         << setw(12) << w[ind] << endl;
}

void ibd::make_pheno_file(string name,
                          const vector<double>& y,
                          const string& trait_name)
{
  string filename = name + "_pheno.txt";
  ofstream outp(filename.c_str());
  if (!outp)
    throw ibd_error("Cannot open file " + filename);
  outp.setf(ios::fixed, ios::floatfield);
  int nind = y.size();
  outp << "ntrt = 2" << endl;
  outp << "nind = "  << nind << endl;
  outp << "miss = *" << endl << endl;
  outp << "ID" << endl;
  outp << trait_name << endl ;
  for (int ind=0; ind < nind; ind++)
    outp << setw(4) << ind
         << setw(12) << y[ind] << endl;
}


void ibd::make_loc_file(string name, const matrix<ObsGeno>& geno,
                        const PopulationType& poptype,
                        const LinkageMap& markermap)
{
  int nind = geno.NrRows();
  int nloc = geno.NrCols();
  string filename = name + ".loc";
  ofstream outp(filename.c_str());
  if (!outp)
    throw ibd_error("Cannot open file " + filename);

  outp << "name = " << name << endl;
  outp << "popt = " << poptype.name() << endl;
  outp << "nloc = " << nloc << endl;
  outp << "nind = " << nind << endl;
  for (int m = 0; m < nloc; m++)
  {
    outp << endl << endl << markermap[m].GetName();
    for (int ind = 0; ind < nind; ind++)
    {
      if (ind % 50 == 0) outp << endl;
      if (ind %  5 == 0) outp << " ";
      poptype.MapQTL(outp, geno[ind][m]);
    }
  }
}

void ibd::make_MapQTL_files(string name,
                            const vector<double>& pheno,
                            const vector<double>& weight,
                            const vector<int>& env_factor,
                            const matrix<ObsGeno>& geno,
                            const LinkageMap& markermap,
                            const PopulationType& poptype)
{
  make_map_file(name,markermap);
  make_pheno_file(name,pheno,weight,env_factor);
  make_loc_file(name,geno,poptype,markermap);
}

void ibd::make_MapQTL_files(string name,
                            const vector<double>& pheno,
                            const string& trait_name,
                            const matrix<ObsGeno>& geno,
                            const LinkageMap& markermap,
                            const PopulationType& poptype)
{
  make_map_file(name,markermap);
  make_pheno_file(name,pheno,trait_name);
  make_loc_file(name,geno,poptype,markermap);
}


string ibd::read_pheno_file(vector<string>& column_name,
                            matrix<string>& qua_info,
                            const string filename)
{
  ifstream inp(filename.c_str());
  if (!inp)
    throw ibd_error("Cannot read file" + filename);
  string str;
  int nind, ntrt;
  string miss = read_qua_header(column_name,nind,ntrt,inp);
  matrix<string> result(nind,ntrt);
  for (int i=0;i<nind;i++)
  {
    for (int j=0;j<ntrt;j++)
    {
      inp >> eatcomment >> str;
      result[i][j] = str;
    }
  }
  qua_info = result;
  return miss;
}


void ibd::read_map_file(LinkageMap&  markermap,
                        vector<string>& chr_name,
                        const string filename)
{
  ifstream inp(filename.c_str());
  if (!inp)
    throw ibd_error("Cannot read file" + filename);
  string str1, str2;
  int chr_nr = -1;
  while (inp)
  {
    double pos;
    inp >> eatcomment >> str1;
    tolower(str1);
    if (str1 == "group" || str1 == "chrom") // new chromosome
    {
      chr_nr++;
      inp >> str2;
      chr_name.push_back(str1 + str2);
    }
    else
    {
      inp >> pos;
      Locus loc(stringify(chr_nr),pos,str1);
      markermap.push_back(loc);
    }
    eatcomment(inp);
  }
}

PopulationType ibd::read_loc_file(string& experiment_name,
                                  map<string, vector<ObsGeno> >& marker_obs,
                                  string filename)
{
  ifstream inp(filename.c_str());
  if (!inp)
    throw ibd_error("Cannot open file" + filename);
  string pop_name;
  int nloc, nind;
  read_loc_header(experiment_name, pop_name, nind, nloc, inp);

  PopulationType poptype(pop_name);
  string marker_name;
  vector<ObsGeno> obs(nind);
  for (int m=0;m<nloc;m++)
  {
    inp >> eatcomment >> marker_name;
    tolower(marker_name);
    for (int i=0;i<nind;i++)
      obs[i] = poptype.MapQTL(inp);
    marker_obs[marker_name] = obs;
  }
  return poptype;
}

PopulationType ibd::read_loc_file(string& experiment_name,
                                  map<string, vector<ObsGeno> >& marker_obs,
                                  int& nind,
                                  string filename)
{
  ifstream inp(filename.c_str());
  if (!inp)
    throw ibd_error("Cannot open file" + filename);
  string pop_name;
  int nloc;
  read_loc_header(experiment_name, pop_name, nind, nloc, inp);

  PopulationType poptype(pop_name);
  string marker_name;
  vector<ObsGeno> obs(nind);
  for (int m=0;m<nloc;m++)
  {
    inp >> eatcomment >> marker_name;
    tolower(marker_name);
    for (int i=0;i<nind;i++)
      obs[i] = poptype.MapQTL(inp);
    marker_obs[marker_name] = obs;
  }
  return poptype;
}


PopulationType ibd::read_MapQTLfiles(string& name_quant_trait,
                                     vector<double>& pheno,
                                     vector<double>& weight,
                                     vector<int>& env_factor,
                                     matrix<ObsGeno>& geno,
                                     LinkageMap& markermap,
                                     const string& locfile,
                                     const string& quafile,
                                     const string& mapfile)
{
  pheno.clear();
  weight.clear();
  env_factor.clear();
  geno.clear();
  markermap.clear();

  string name;
  LinkageMap complete_markermap;
  vector<string> chr_name;
  vector<string> column_name;
  matrix<string> qua_info;
  map<string, vector<ObsGeno> > marker_obs;

  PopulationType poptype = read_loc_file(name,marker_obs,locfile);
  string miss = read_pheno_file(column_name,qua_info,quafile);
  read_map_file(complete_markermap,chr_name,mapfile);
  string qtrait_name;
  cout << "pheno file" << endl;
  for (unsigned int c=0;c<column_name.size();c++)
    cout << setw(3) << c << setw(20) << column_name[c] << endl;
  int col_nr_pheno, col_nr_weight, col_nr_env;
  cout << "Column number of quantitative trait: ";
  cin >> col_nr_pheno;
  name_quant_trait = column_name[col_nr_pheno];

  cout << "Column number of weighting factor: (-1 : no weight): ";
  cin >> col_nr_weight;

  cout << "Column number of env. factor: (-1 : no env. factor): ";
  cin >> col_nr_env;

  int nind = qua_info.NrRows();
  vector<bool> miss_ind(nind);

  for (int i=0;i<nind;i++)
  {
    if (qua_info[i][col_nr_pheno] == miss)
      miss_ind[i] = true;
    else
    {
      miss_ind[i] = false;
      pheno.push_back(strtod(qua_info[i][col_nr_pheno].c_str(), NULL));
      double w = (col_nr_weight >= 0) ? strtod(qua_info[i][col_nr_weight].c_str(), NULL) : 1.0;
      weight.push_back(w);
      string& str = qua_info[i][col_nr_env];
      int env = (col_nr_env >= 0) ? atoi(str.c_str()) : 0;
      env_factor.push_back(env);
    }
  }
  cout << endl;
  for (unsigned int m=0;m<complete_markermap.size();m++)
  {
    Locus& loc = complete_markermap[m];
    map<string, vector<ObsGeno> >::iterator iter = marker_obs.find(loc.GetName());
    if (iter != marker_obs.end())
    {
      markermap.push_back(loc);
      vector<ObsGeno> v = iter->second;
      vector<ObsGeno> v_sel;
      for (int i=0;i<nind;i++)
      {
        if (miss_ind[i] == false)
          v_sel.push_back(v[i]);
      }
      geno.push_back(v_sel);
    }
  }
  geno = transpose(geno);
  return poptype;
}

/*
 PopulationType ibd::read_MapQTLfiles(string& name_quant_trait,
 vector<double>& pheno,
 vector<double>& weight,
 vector<int>& env_factor,
 matrix<ObsGeno>& geno,
 LinkageMap& markermap,
 const string& locfile,
 const string& quafile,
 const string& mapfile)
 {
 pheno.clear();
 weight.clear();
 env_factor.clear();
 geno.clear();
 markermap.clear();

 string name;
 LinkageMap complete_markermap;
 vector<string> chr_name;
 vector<string> column_name;
 matrix<string> qua_info;
 map<string, vector<ObsGeno> > marker_obs;

 PopulationType poptype = read_loc_file(name,marker_obs,locfile);
 string miss = read_qua_file(column_name,qua_info,quafile);
 read_map_file(complete_markermap,chr_name,mapfile);
 string qtrait_name;
 cout << "qua file" << endl;
 for (int c=0;c<column_name.size();c++)
 cout << setw(3) << c << setw(20) << column_name[c] << endl;
 int col_nr_pheno, col_nr_weight, col_nr_env;
 cout << "Column number of quantitative trait: ";
 cin >> col_nr_pheno;
 name_quant_trait = column_name[col_nr_pheno];

 cout << "Column number of weighting factor: (-1 : no weight): ";
 cin >> col_nr_weight;

 cout << "Column number of env. factor: (-1 : no env. factor): ";
 cin >> col_nr_env;

 int nind = qua_info.NrRows();
 vector<bool> miss_ind(nind);

 for (int i=0;i<nind;i++)
 {
 if (qua_info[i][col_nr_pheno] == miss)
 miss_ind[i] = true;
 else
 {
 miss_ind[i] = false;
 pheno.push_back(strtod(qua_info[i][col_nr_pheno]));
 double w = (col_nr_weight >= 0) ? strtod(qua_info[i][col_nr_weight]) : 1.0;
 weight.push_back(w);
 string& str = qua_info[i][col_nr_env];
 int env = (col_nr_env >= 0) ? atoi(str.c_str()) : 0;
 env_factor.push_back(env);
 }
 }
 cout << endl;
 for (int m=0;m<complete_markermap.size();m++)
 {
 Locus& loc = complete_markermap[m];
 map<string, vector<ObsGeno> >::iterator iter = marker_obs.find(loc.GetName());
 if (iter != marker_obs.end())
 {
 markermap.push_back(loc);
 vector<ObsGeno> v = iter->second;
 vector<ObsGeno> v_sel;
 for (int i=0;i<nind;i++)
 {
 if (miss_ind[i] == false)
 v_sel.push_back(v[i]);
 }
 geno.push_back(v_sel);
 }
 }
 geno = transpose(geno);
 return poptype;
 }
 */


PopulationType ibd::read_MapQTLfiles(string& name_quant_trait,
                                     vector<double>& pheno,
                                     vector<double>& weight,
                                     vector<int>& env_factor,
                                     matrix<ObsGeno>& geno,
                                     LinkageMap& markermap,
                                     const string& locfile,
                                     const string& quafile,
                                     const string& mapfile,
                                     int col_nr_pheno,
                                     int col_nr_weight,
                                     int col_nr_env)
{
  pheno.clear();
  weight.clear();
  env_factor.clear();
  geno.clear();
  markermap.clear();

  string name;
  LinkageMap complete_markermap;
  vector<string> chr_name;
  vector<string> column_name;
  matrix<string> qua_info;
  map<string, vector<ObsGeno> > marker_obs;

  PopulationType poptype = read_loc_file(name,marker_obs,locfile);
  string miss = read_pheno_file(column_name,qua_info,quafile);
  read_map_file(complete_markermap,chr_name,mapfile);
  int nind = qua_info.NrRows();
  vector<bool> miss_ind(nind);
  name_quant_trait = column_name[col_nr_pheno];

  for (int i=0;i<nind;i++)
  {
    if (qua_info[i][col_nr_pheno] == miss)
      miss_ind[i] = true;
    else
    {
      miss_ind[i] = false;
      pheno.push_back(strtod(qua_info[i][col_nr_pheno].c_str(), NULL));
      double w = (col_nr_weight >= 0) ? strtod(qua_info[i][col_nr_weight].c_str(), NULL) : 1.0;
      weight.push_back(w);
      string& str = qua_info[i][col_nr_env];
      int env = (col_nr_env >= 0) ? atoi(str.c_str()) : 0;
      env_factor.push_back(env);
    }
  }
  cout << endl;
  for (unsigned int m=0;m<complete_markermap.size();m++)
  {
    Locus& loc = complete_markermap[m];
    map<string, vector<ObsGeno> >::iterator iter = marker_obs.find(loc.GetName());
    if (iter != marker_obs.end())
    {
      markermap.push_back(loc);
      vector<ObsGeno> v = iter->second;
      vector<ObsGeno> v_sel;
      for (int i=0;i<nind;i++)
      {
        if (miss_ind[i] == false)
          v_sel.push_back(v[i]);
      }
      geno.push_back(v_sel);
    }
  }
  geno = transpose(geno);
  return poptype;
}


PopulationType ibd::read_MapQTLfiles(vector<double>& pheno,
                                     matrix<ObsGeno>& geno,
                                     LinkageMap& markermap,
                                     const string& locfile,
                                     const string& quafile,
                                     const string& mapfile,
                                     const string& trait_name)
{
  pheno.clear();
  geno.clear();
  markermap.clear();

  string name;
  LinkageMap complete_markermap;
  vector<string> chr_name;
  vector<string> column_name;
  matrix<string> qua_info;
  map<string, vector<ObsGeno> > marker_obs;

  PopulationType poptype = read_loc_file(name,marker_obs,locfile);
  string miss = read_pheno_file(column_name,qua_info,quafile);
  read_map_file(complete_markermap,chr_name,mapfile);
  int nind = qua_info.NrRows();
  vector<bool> miss_ind(nind);

  vector<string>::const_iterator it;
  it = find(column_name.begin(), column_name.end(), trait_name);
  if (it == column_name.end()) throw ibd_error("Cannot find trait!");
  int col_nr_pheno = it - column_name.begin();

  for (int i=0;i<nind;i++)
  {
    if (qua_info[i][col_nr_pheno] == miss)
      miss_ind[i] = true;
    else
    {
      miss_ind[i] = false;
      pheno.push_back(strtod(qua_info[i][col_nr_pheno].c_str(), NULL));
    }
  }
  for (unsigned int m=0;m<complete_markermap.size();m++)
  {
    Locus& loc = complete_markermap[m];
    map<string, vector<ObsGeno> >::iterator iter = marker_obs.find(loc.GetName());
    if (iter != marker_obs.end())
    {
      markermap.push_back(loc);
      vector<ObsGeno> v = iter->second;
      vector<ObsGeno> v_sel;
      for (int i=0;i<nind;i++)
      {
        if (miss_ind[i] == false)
          v_sel.push_back(v[i]);
      }
      geno.push_back(v_sel);
    }
  }
  geno = transpose(geno);
  return poptype;
}


// nog lezen van commentaar toevoegen
LinkageMap ibd::read_cov_file(const string& filename, LinkageMap& markermap)
{
  ifstream inp(filename.c_str());
  if (!inp) return markermap;
  if (MapQTLvar(inp) != "ncof")
    throw ibd_error("command ncof not found");
  unsigned int ncof;
  inp >> ncof;

  LinkageMap selected_map;
  string marker;
  while (inp)
  {
    inp >> marker;
    if (marker.empty()) break;
    tolower(marker);

    LinkageMap::const_iterator it = markermap.begin();
    for (it = markermap.begin(); it != markermap.end(); ++it)
      if (it->GetName() == marker)
        break;
      if (it != markermap.end())
        selected_map.push_back(*it);
      else
        throw ibd_error(marker + " not in map");
  }

  if (selected_map.size() != ncof)
    throw ibd_error("ncof != number of cofactors");

  sort(selected_map.begin(), selected_map.end());

  return selected_map;
}


bool ibd::check_obsgeno(const matrix<ObsGeno>& geno, const LinkageMap& markermap)
{
  bool error = false;
  const int nind = geno.NrRows();
  const int nloc = geno.NrCols();
  for (int m=0;m<nloc-1;m++)
  {
    if (recomb(markermap[m],markermap[m+1]) < DBL_MIN)
    {
      for (int i=0;i<nind;i++)
      {
        ObsGeno g = intersection(geno[i][m],geno[i][m+1]);
        if (g.empty())
        {
          if (error == false)
          {
            error = true;
            cout.setf(ios::left,ios::adjustfield);
            cout << "Error: The following markers have conflicting   " << endl
                 << "scores (same position on the markermap, different scores): "
                 << endl << endl
                 << setw(20) << "marker 1"
                 << setw(20) << "marker 2"
                 << setw(15) << "individual" << endl << endl;
          }
          cout << setw(20) << markermap[m].GetName()
               << setw(20) << markermap[m+1].GetName()
               << setw(15) << i+1 << endl;
        }
      }
    }
  }
  cout << endl;
  return !error;
}

void ibd::remove_inconsistent_scores(matrix<ObsGeno>& geno, const LinkageMap& markermap,
                                     const ObsGeno& U)
{
  const int nind = geno.NrRows();
  const int nloc = geno.NrCols();
  for (int m1=0;m1<nloc-1;m1++)
  {
    for (int m2=m1+1;m2<nloc;m2++)
    {
      if (recomb(markermap[m1],markermap[m2]) > DBL_MIN) break;
      for (int i=0;i<nind;i++)
      {
        ObsGeno g = intersection(geno[i][m1],geno[i][m2]);
        if (g.empty())
          geno[i][m1] = geno[i][m2] = U;
      }
    }
  }
}



bool ibd::get_geno_and_markermap(matrix<ObsGeno>& geno, LinkageMap& markermap,
                                 const map<string, vector<ObsGeno> > marker_obs,
                                 const LinkageMap& complete_markermap, int nind)
{
  for (unsigned int m=0;m<complete_markermap.size();m++)
  {
    const Locus& loc = complete_markermap[m];
    map<string, vector<ObsGeno> >::const_iterator iter = marker_obs.find(loc.GetName());
    if (iter != marker_obs.end())
    {
      markermap.push_back(loc);
      vector<ObsGeno> v = iter->second;
      vector<ObsGeno> v_sel;
      for (int i=0;i<nind;i++)
        v_sel.push_back(v[i]);
      geno.push_back(v_sel);
    }
  }
  geno = transpose(geno);

  return check_obsgeno(geno,markermap);
}


pair<int,int> ibd::dimension_table(const string& filename)
{
  string line;
  ifstream inp;
  OpenFile(inp,filename);
  getline(inp,line);
  int ncols = count(line.begin(),line.end(),'\t') + 1;
  int nrows = 1;
  while (getline(inp,line))
  {
    if (line.empty()) break;
    nrows++;
    if (ncols != count(line.begin(),line.end(),'\t') + 1)
      throw ibd_error("error!!");
  }
  return pair<int,int>(nrows,ncols);
}

PopulationType ibd::read_loc_file_format2(map<string,vector<ObsGeno> > & marker_obs,
                                          int& nind, const string& filename)
{
  pair<int,int> dim = dimension_table(filename); // first check of file
  nind = dim.first - 1;
  int nloc = dim.second -1;

  ifstream inp;
  OpenFile(inp,filename);

  vector<string> ind_name(nind);
  vector<string> marker_name(nloc);
  matrix<ObsGeno> geno(nloc,nind);

  string pop_name;
  inp >> pop_name >> marker_name;
  PopulationType poptype(pop_name);

  for (int m=0;m<nloc;m++)
    tolower(marker_name[m]);

  for (int ind=0;ind<nind;ind++)
  {
    inp >> ind_name[ind];
    for (int m=0;m<nloc;m++)
      geno[m][ind] = poptype.MapQTL(inp);
  }
  for (int m=0;m<nloc;m++)
    marker_obs[marker_name[m]] = geno[m];
  return poptype;
}

void ibd::make_loc_file_format2(string name, const matrix<ObsGeno>& geno,
                                const PopulationType& poptype,
                                const LinkageMap& markermap)
{
  int nind = geno.NrRows();
  int nloc = geno.NrCols();
  string filename = name + ".loc";
  ofstream outp(filename.c_str());
  if (!outp)
    throw ibd_error("Cannot open file " + filename);
  outp << poptype.name();
  for (int m = 0; m < nloc; m++)
    outp << '\t' << markermap[m].GetName();
  outp << '\n';
  for (int ind=0;ind<nind;ind++)
  {
    outp << ind+1;
    for (int m=0;m<nloc;m++)
    {
      outp << '\t';
      poptype.MapQTL(outp,geno[ind][m]);
    }
    outp << '\n';
  }
}

void ibd::make_loc_file_format2(string name, const matrix<ObsGeno>& geno,
                                const vector<string>& ID,
                                const PopulationType& poptype,
                                const LinkageMap& markermap)
{
  int nind = geno.NrRows();
  int nloc = geno.NrCols();
  string filename = name + ".loc";
  ofstream outp(filename.c_str());
  if (!outp)
    throw ibd_error("Cannot open file " + filename);
  outp << poptype.name();
  for (int m = 0; m < nloc; m++)
    outp << '\t' << markermap[m].GetName();
  outp << '\n';
  for (int ind=0;ind<nind;ind++)
  {
    outp << ID[ind];
    for (int m=0;m<nloc;m++)
    {
      outp << '\t';
      poptype.MapQTL(outp,geno[ind][m]);
    }
    outp << '\n';
  }
}

void ibd::convert2new_format_loc_file(const vector<string>& ID,
                                      const string& old_locfile,
                                      const string& new_locfile)
{
  ofstream outp;
  ifstream inp;
  OpenFile(inp,old_locfile);
  OpenFile(outp,new_locfile);

  string pop_name,experiment_name;
  int nloc,nind;
  read_loc_header(experiment_name, pop_name, nind, nloc, inp);

  PopulationType poptype(pop_name);
  vector<string> marker_names(nloc);
  matrix<ObsGeno> obs(nind,nloc);
  for (int m=0;m<nloc;m++)
  {
    inp >> eatcomment >> marker_names[m];
    for (int i=0;i<nind;i++)
      obs[i][m] = poptype.MapQTL(inp);
  }

  outp << pop_name;
  for (int m = 0; m < nloc; m++)
    outp << '\t' << marker_names[m];
  outp << '\n';
  for (int ind=0;ind<nind;ind++)
  {
    outp << ID[ind];
    for (int m=0;m<nloc;m++)
    {
      outp << '\t';
      poptype.MapQTL(outp,obs[ind][m]);
    }
    outp << '\n';
  }
}



/*
 pop_base * read_MapQTLfiles(vector<double>& pheno,
 matrix<ObsGeno>& geno,
 LinkageMap& markermap,
 int& dimE,
 const string& locfile,
 const string& quafile,
 const string& mapfile)
 {
 int i,env;
 string name;
 vector<Locus> loci;
 vector<string> chr_name;
 vector<string> loc_name;
 vector<string> column_name;
 matrix<string> column;
 map<string, vector<ObsGeno> > marker_obs;

 pop_base *poptype = read_loc_file(name,marker_obs,locfile);
 string miss = read_qua_file(column_name,column,quafile);
 read_map_file(loci,loc_name,chr_name,mapfile);
 for (int k=0;k<loc_name.size();k++)
 cout << setw(3) << k << "   " <<  loc_name[k] << endl;
 string qtrait_name;
 cout << "qua file" << endl;
 for (int c=0;c<column_name.size();c++)
 cout << setw(3) << c+1 << setw(20) << column_name[c] << endl;
 cout << "Give number of traits: "; cin >> dimE;
 vector<int> ndx(dimE);
 for (env = 0; env<dimE;env++)
 {
 cout << "Give column nr: "; cin >> ndx[env];
 ndx[env] -= 1;
 }

 int nind = column.NrCols();
 vector<bool> miss_ind(nind);

 cout << "Missing individuals: ";
 for (i=0;i<nind;i++)
 {
 miss_ind[i] = false;
 for (env = 0; env<dimE;env++)
 if (column[ndx[env]][i] == miss)
 miss_ind[i] = true;
 if (miss_ind[i] == true)
 cout << i << " ";
 }
 cout << endl;
 for (env=0;env<dimE;env++)
 {
 for (i=0;i<nind;i++)
 {
 if (miss_ind[i] == false)
 {
 double val = strtod(column[ndx[env]][i]);
 pheno.push_back(val);
 }
 }
 }

 for (int m=0;m<loci.size();m++)
 {
 Locus loc = loci[m];
 string locname = loc_name[m];
 map<string, vector<ObsGeno> >::iterator iter;
 iter = marker_obs.find(locname);
 if (iter != marker_obs.end())
 {
 markermap.AddLocus(loc.GetChr(),loc.GetPosition());
 vector<ObsGeno> v = iter->second;
 vector<ObsGeno> v_sel;
 for (int i=0;i<nind;i++)
 {
 if (miss_ind[i] == false)
 v_sel.push_back(v[i]);
 }
 geno.push_back(v_sel);
 }
 }
 geno = transpose(geno);
 return poptype;
 }


 int test_qua()
 {
 vector<string> col_name;
 matrix<string> col;
 read_qua_file(col_name,col,"BCtest.qua");
 for (int j=0;j<col_name.size();j++)
 cout << col_name[j] << endl;
 for (int i=0;i<10;i++)
 cout << setw(3) << i << setw(12) << strtod(col[1][i]) << endl;
 return 0;
 }

 int test_loc()
 {
 string name;
 map<string, vector<ObsGeno> > marker_obs;
 pop_base * poptype = read_loc_file(name,marker_obs,"BCtest.loc");
 cout << name << endl;
 cout << poptype->name() << endl;
 vector<ObsGeno> M1 = marker_obs["M1"];
 for (int i=0;i<20;i++)
 {
 if (i%5 == 0) cout << endl;
 cout << setw(10) << M1[i];
 }
 return 0;
 }


 int test_map()
 {
 vector<string> loc_name;
 vector<string> chr_name;
 vector<Locus> markermap;
 read_map_file(markermap,loc_name,chr_name,"BCtest.loc");
 for (int i=0;i<3;i++)
 {
 Locus loc = markermap[i];
 int chr = loc.GetChr();
 double pos = loc.GetPosition();
 cout << setw(3)  << i
      << setw(10) << loc_name[i]
      << setw(10) << chr
      << setw(15) << chr_name[chr]
      << setw(12) << pos << endl;
 }
 return 0;
 }

 int test_reading_MapQTL()
 {
 vector<double> pheno;
 matrix<ObsGeno> geno;
 LinkageMap markermap;
 pop_base * poptype;
 poptype = read_MapQTLfiles(pheno,geno,markermap,"BCtest.loc","BCtest.qua","BCtest.map");
 cout << "Poptype: " << poptype->name() << endl;
 for (int i=0;i<5;i++)
 cout << setw(3) << setw(12) << pheno[i] << endl;
 cout << "Genotype ind 3, marker M1: " << geno[3][1] << endl;
 for (int j=0;j<15;j++)
 {
 Locus loc = markermap[j];
 cout << setw(5) << loc.GetChr() << setw(15) << loc.GetPosition() << endl;
 }
 return 0;
 }

 */
