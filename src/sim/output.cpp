#include "output.h"
#include "convert.h"
#include <Rcpp.h>

using namespace ibd;
using namespace std;

void make_ped_file(const SimPop& sim_pop, const string& filename)
{
  ofstream outp;
  OpenFile(outp, filename + ".ped");
  for (SimPop::const_iterator it=sim_pop.begin();it!=sim_pop.end();it++)
  {
    const IndProp& indprop = it->first;
    outp << setw(12) << indprop.GetID()
         << setw(12) << indprop.GetType()
         << setw(12) << indprop.GetP1()
         << setw(12) << indprop.GetP2() << endl;
  }
}

void make_loc_file(const SimPop& sim_pop,
                   const LinkageMap& markermap,
                   const vector<MarkerType>& markertype,
                   const string& filename)
{
  const int M = markermap.size();
  ofstream outp;
  OpenFile(outp,filename + ".loc");
  outp << "ID" << '\t';
  for (int m=0;m<M;m++)
  {
    char delimit = (m==M-1) ? '\n' : '\t';
    outp << markermap[m].GetName() << delimit;
  }
  for (SimPop::const_iterator it=sim_pop.begin();it!=sim_pop.end();it++)
  {
    if (!it->first.IsHybrid())
    {
      const Genome& genome = it->second;
      outp << it->first.GetID() << '\t';
      vector<Genotype> geno = genome.GetGenotype(markermap);
      for (int m=0;m<M;m++)
      {
        char delimit = (m==M-1) ? '\n' : '\t';
        Genotype g = geno[m];
        outp << markertype[m](g) << delimit;
      }
    }
  }
}

void make_pheno_file(const SimPop& sim_pop,
                     const Phi& phi,
                     double sigma,
                     const string& filename)
{
  ofstream outp;
  OpenFile(outp,filename + "_pheno.txt");
  outp << setw(12) << "ID"
       << setw(15) << "pheno"
       << setw(15) << "geno"
       << setw(15) << "error" << endl;

  for (SimPop::const_iterator it=sim_pop.begin();it!=sim_pop.end();it++)
  {
    if (it->first.IsMemberFamily()) // progeny
    {
      const Genome& genome = it->second;
      double geno_val = phi(genome);
      double error = R::rnorm(0.0,sigma);
      double pheno_val = geno_val + error;
      string ID = it->first.GetID();
      outp << setw(12) << ID
           << setw(15) << pheno_val
           << setw(15) << geno_val
           << setw(15) << error << endl;
    }
  }
}

// void make_eval_file(const string& eval_filename, vector<double>& chr_length,
// 					double dist_eval_pos)
// {
// 	ofstream outp;
// 	OpenFile(outp,eval_filename);
// 	const int Nchr = chr_length.size();
// 	for (int chr=0;chr<Nchr;chr++)
// 	{
// 		double len = chr_length[chr];
// 		for (double pos = 0.0; pos <=len; pos+=dist_eval_pos)
// 			outp << setw(3) << chr+1 << setw(12) << pos << endl;
// 	}
// }

void make_part_ldf_file(int N,vector<double>& chr_length,
                        double dist_eval_pos, string filename)
{
  LinkageMap eval_map;
  const int Nchr = chr_length.size();
  for (int chr=0;chr<Nchr;chr++)
  {
    double len = chr_length[chr];
    for (double pos = 0.0; pos <=len; pos+=dist_eval_pos)
      eval_map.push_back(Locus(stringify(chr+1),pos));
  }

  ofstream outp;
  OpenFile(outp,filename);
  const int nloc = eval_map.size();
  outp << "nind " << N << endl;
  outp << "indnr  indname " << endl;

  outp << "npos " << nloc << endl;
  outp << "posnr " << endl;
  for (int m=0;m<nloc;m++)
  {
    outp << setw(5) << m+1 << setw(6) << eval_map[m].GetChr()
         << setw(12) << eval_map[m].GetPosition() << endl;
  }

  outp << endl << "inbred 1" << endl;

  matrix<double> Id = identity_matrix(N);
  for (int m=0;m<nloc;m++)
  {
    outp << endl << setprecision(2) << setw(5) << Id;
  }

}


void make_flapjack_map_file(string name, const LinkageMap& marker_map)
{
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
    outp << loc.GetName() << '\t' << loc.GetChr() << '\t' << loc.GetPosition() << '\n';
  }
}



void print(ostream& outp, const map<Locus,vector<double> >& QTLs)
{
  typedef map<Locus,vector<double> >::const_iterator IterQTLs;
  for (IterQTLs iter=QTLs.begin();iter!=QTLs.end();++iter)
  {
    const Locus& loc = iter->first;
    const vector<double>& par = iter->second;
    outp << setw(3) << loc.GetChr() << setw(12) << loc.GetPosition()
         << setw(6) << loc.GetName() << setw(12) << par << endl;
  }
}

void print_genstat_table(const vector<PopProp>& pops, const map<string,string>& inbfnds,
                         const string filename)
{
  ofstream outp;
  OpenFile(outp,filename);

  outp << setw(5) << "crsnr" << setw(12) << "name" << endl;

  int crsnr = 1;
  for (vector<PopProp>::const_iterator it = pops.begin();it!=pops.end();it++)
  {
    outp << setw(5) << crsnr << setw(12) << it->GetName() << endl;
    crsnr++;
  }

  /*
   outp << setw(5) << "crsnr" << setw(12) << "name"
        << setw(8) << "type" << setw(6) << "nind"
        << setw(8) << "ndx[1]" << setw(8) << "ndx[2]"
        << setw(8) << "ndx[3]" << setw(8) << "ndx[4]" << endl;

   set<string> inbfnd_used; // only select founders used in crosses
   for (map<string,string>::const_iterator iter1=inbfnds.begin();iter1!=inbfnds.end();iter1++)
   {
   string fnd_name = iter1->first;
   for (vector<PopProp>::const_iterator iter2 = pops.begin();iter2!=pops.end();iter2++)
   {
   if (iter2->P[0] == fnd_name || iter2->P[1] == fnd_name ||
   iter2->P[2] == fnd_name || iter2->P[3] == fnd_name)
   inbfnd_used.insert(fnd_name);
   }
   }

   int crsnr = 1;
   for (vector<PopProp>::const_iterator it = pops.begin();it!=pops.end();it++)
   {
   int c1,c2,c3,c4;
   int k=1;
   for (set<string>::const_iterator it2=inbfnd_used.begin();it2!=inbfnd_used.end();it2++)
   {
   if (*it2 == it->P[0]) c1 = k;
   if (*it2 == it->P[1]) c2 = k;
   if (*it2 == it->P[2]) c3 = k;
   if (*it2 == it->P[3]) c4 = k;
   k++;
   }
   outp << setw(5) << crsnr << setw(12) << it->name
        << setw(8) << it->type << setw(6) << it->nind
        << setw(8) << c1 << setw(8) << c2
        << setw(8) << c3 << setw(8) << c4 << endl;
   crsnr++;
   }
   */
}
