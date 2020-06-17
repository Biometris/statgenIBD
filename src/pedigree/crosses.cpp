#include <Rcpp.h>

#include "Loc.h"
#include "OrdGeno.h"
#include "analysis_fam.h"
#include "crosses.h"

using namespace ibd;
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

IndProp find_first_progeny(const vector<IndProp>& pop)
{
  for (vector<IndProp>::const_iterator it=pop.begin();it!=pop.end();it++)
  {
    if (it->IsMemberFamily())
      return *it;
  }
  throw ibd_error("Cannot find progeny in function find_first_progeny");
  return pop[0]; // dummy
}

string find_type(const vector<IndProp>& pop)
{
  return find_first_progeny(pop).GetType();
}

matrix<double> calc_P(const LinkageMap& eval_pos,
                      int nparents,
                      const IBD_fam& IBD_ind)
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
  int ndx_P1 = ndxID(pop, P1);
  int ndx_P2 = ndxID(pop, P2);
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
      ndx_par.push_back(ndxID(pop, par1.GetP1()));
      ndx_par.push_back(ndxID(pop, par1.GetP2()));
      ndx_par.push_back(ndx_P2);
    }
    else
    {
      ndx_par.push_back(ndxID(pop, par2.GetP1()));
      ndx_par.push_back(ndxID(pop, par2.GetP2()));
      ndx_par.push_back(ndx_P1);
    }
  }
  else // npar == 4
  {
    ndx_par.push_back(ndxID(pop, par1.GetP1()));
    ndx_par.push_back(ndxID(pop, par1.GetP2()));
    ndx_par.push_back(ndxID(pop, par2.GetP1()));
    ndx_par.push_back(ndxID(pop, par2.GetP2()));
  }
  return ndx_par;
}

matrix3D<double> analysis_cross(const vector<IndProp>& pop,
                                const matrix<score>& geno,
                                const LinkageMap& markermap,
                                const LinkageMap& eval_pos)
{
  Rcpp::Rcout << "analysis of family ........" << endl;

  string type = find_type(pop);
  vector<int> ndx_par = get_ndx_par(pop);
  matrix3D<double> Z = calc_IBDs(pop, ndx_par, geno, markermap, eval_pos, type);
  return Z;
}

