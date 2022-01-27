#ifndef SIMUL_HEADER
#define SIMUL_HEADER

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include "poptype.h"
#include "Random.h"

namespace ibd 
{

template<class G>
double calc_genetic_var(const G& QTLmodel, const LinkageMap& QTLmap, 
						const PopulationType& poptype)
{
	const int nqtl = QTLmap.size();
	ObsGeno U = poptype.unknown_genotype();
	std::vector<ObsGeno> obs_geno(nqtl,U);

	std::vector<double> P;
	matrix<Genotype> geno_augm;
	std::vector<double> r = make_rec_map(QTLmap);
	int n = calcprob_ind(geno_augm,P,r,obs_geno,poptype,0.0);

	double E_x = 0.0;
	double E_sqr_x = 0.0;
	for (int i=0;i<n;i++)
	{
		double g_val = QTLmodel(geno_augm[i]);
		E_x += P[i]*g_val;
		E_sqr_x += P[i]*sqr(g_val);
	}
	return E_sqr_x - sqr(E_x);
}

template<class G>
double calc_genetic_mean(const G& QTLmodel, const LinkageMap& QTLmap, 
						 const PopulationType& poptype)
{
	const int nqtl = QTLmap.size();
	ObsGeno U = poptype.unknown_genotype();
	std::vector<ObsGeno> obs_geno(nqtl,U);

	std::vector<double> P;
	matrix<Genotype> geno_augm;
	std::vector<double> r = make_rec_map(QTLmap);
	int n = calcprob_ind(geno_augm,P,r,obs_geno,poptype,0.0);

	double E_x = 0.0;
	for (int i=0;i<n;i++)
	{
		double g_val = QTLmodel(geno_augm[i]);
		E_x += P[i]*g_val;
	}
	return E_x;
}


// Algemene classe voor omzetting van QTL genotype in genotypische waarde

class QTLModel
{
public:
	QTLModel(){}
	~QTLModel(){}
	virtual double operator()(const std::vector<Genotype>& g) const = 0;
	virtual double operator()(const std::vector<Genotype>& g, int env) const = 0;
};


class GenoValue : public QTLModel
{
public:
	GenoValue(){}
	GenoValue(PopulationType poptype, double mu, const std::vector<double>& main,
										const std::vector<double>& inter);
	double value(const std::vector<Genotype>& g) const;
	double operator()(const std::vector<Genotype>& g) const { return value(g); }
	double operator()(const std::vector<Genotype>& g, int env) const
	{ if (env == 0) return value(g);
	  ibd_error("Env factor != 0 in class Genovalue");
	  return 0.0;
	}
private:
	PopulationType poptype;
	std::vector<double> b;
};

class NO_QTL : public QTLModel
{
public:
	NO_QTL(double mu) :_mu(mu) {};
	double operator()(const std::vector<Genotype>& g) const { return 0.0; }
private:
	double _mu;
};

std::vector<ObsGeno> generate_markerdata(const Genome& genome, 
									const LinkageMap& MarkerMap,
									double fraction_missing,
									PopulationType poptype);


void sim(std::vector<double>& pheno,									    // phenotypic values
		 std::vector<double>& geno_val,										// genotypic values
		 matrix<ObsGeno>& geno,                                             // marker observations
		 const QTLModel& QTLmodel, const LinkageMap& QTLmap,			    // QTL info 
		 double sigma_env, const std::vector<double>& weight,				// env. var + weights 
		 const LinkageMap& markermap, double fr_miss,					    // marker info
		 std::vector<double> chr_length, int nind, PopulationType poptype);	// population info

// void read_chromosome_lengths(std::vector<double>& chr_length, std::ifstream &f);
void read_markermap(LinkageMap& markermap, std::ifstream& f);

void read_qtlmap(LinkageMap& QTLmap, std::vector<double>& main_effect, 
				 const std::string& filename,
				 PopulationType poptype);

PopulationType  read_param(int& nind,
				double& herit,
				double& fraction_missing,
				bool& random_start,
				int& nchr, double& length, int& nloc_chr,
				double& win_sz,
				double& eff_dim_epi,
				const std::string& filename);

std::vector<double> read_interactions(const std::string& filename, int nqtl, 
								 PopulationType poptype);

}
#endif
