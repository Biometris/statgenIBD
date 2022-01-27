#include "simul.h"
#include "util_genetics.h"
#include "Loc.h"
#include "convert.h"

using namespace std;

vector<ibd::ObsGeno> ibd::generate_markerdata(const Genome& genome, 
									const LinkageMap& MarkerMap,
									double fraction_missing,
									PopulationType poptype)
{
	const int nloc = MarkerMap.size();
	vector<Genotype> g = genome.GetGenotype(MarkerMap);
	vector<ObsGeno> result(nloc);
	ObsGeno U = poptype.unknown_genotype();
	for (int loc=0;loc<nloc;loc++)
	{
		if (randuniform() < fraction_missing)
			result[loc] = U;	
		else
			result[loc] = g[loc];
	}
	return result;
}


double ibd::GenoValue::value(const vector<Genotype>& g) const
{
	int nqtl = g.size();
	int npar = poptype.nparam();

	// compute genotypic value
	matrix<double> m(nqtl,npar);
	vector<double>::const_iterator iter = b.begin();
	double val = *iter++; 
	for (int q=0;q<nqtl;q++)
	{
		m[q] = poptype.model(g[q]);
		for (int j=0;j<npar;j++)
			val += m[q][j]*(*iter++);		
	}
	// 25 jan 2002, nog nagaan of dit goed gaat bij bijvoorbeeld F2,
	// zie ook functie read_interactions(). Wellicht is het beter om te werken met
	// een afgeleid type van vector<double>, waarbij we kunnen zetten:
	// vec.GetParam(int q1, int q2, int par_q1, int par_q2) o.i.d.
	// voor F2 doet de onderstaande code:
	// qtl 1	qtl 2	add/dom 1	 add/dom 2
	//	1		  0		  add		  add
	//	1		  0		  add		  dom
	//	1		  0		  dom		  add
	//  1		  0		  dom		  dom
	//	2		  0		  add		  add
	//	2		  0		  add		  dom
	//	2		  0		  dom		  add
	//  2		  0		  dom		  dom
	//	2		  1		  add		  add
	//	2		  1		  add		  dom
	//	2		  1		  dom		  add
	//  2		  1		  dom		  dom
	//  etc.
	for (int q1=0;q1<nqtl;q1++)
		for (int q2=0;q2<q1;q2++)
			for (int j1=0;j1<npar;j1++)
				for (int j2=0;j2<npar;j2++)
					val += m[q1][j1]*m[q2][j2]*(*iter++);

	if (iter != b.end())
		throw ibd_error("GenoValue::operator()");
	return val;
}


ibd::GenoValue::GenoValue(PopulationType pt, double mu, 
									const vector<double>& main,
									const vector<double>& inter) 
									: poptype(pt) 
{
	b.push_back(mu);
	vector<double>::const_iterator iter;
	for (iter = main.begin(); iter != main.end(); iter++)
		b.push_back(*iter);
	for (iter = inter.begin(); iter != inter.end(); iter++)
		b.push_back(*iter);
}


// The QTL positions in the input file must be given in ascending order!
void ibd::read_qtlmap(LinkageMap& QTLmap, vector<double>& main_effect, const string& filename,
				 PopulationType poptype)
{	
	ifstream f(filename.c_str());
	if (!f)
		throw ibd_error("Cannot open file " + filename);

	int chr;
	const string& str = "_";
	double pos, effect;

	int npar = poptype.nparam();

	f >> eatcomment;
	while (f)
	{
		f >> chr >> pos;
		QTLmap.push_back(Locus(stringify(chr),pos,str));
		for (int j=0;j<npar;j++)
		{
			f >> effect;
			main_effect.push_back(effect);
		}
		f >> skip_rest_of_line >> eatcomment; 
	}
	if (!f.eof())
		throw ibd_error("Error reading file " + filename);
}

vector<double> ibd::read_interactions(const string& filename, int nqtl, PopulationType poptype)
{
	ifstream f(filename.c_str());
	if (!f)
		throw ibd_error("Cannot open file " + filename);
	int npar = poptype.nparam();
	int npar_per_pair = npar*npar;
	int npairs = nqtl*(nqtl-1)/2;
	vector<double> result(npairs*npar_per_pair,0.0);
	int loc1,loc2;
	f >> eatcomment;
	while (f)
	{
		// let op, nog niet helemaal goed, zie commentaar bij GenoValue::operator()
		f >> loc1 >> loc2; 
		if (loc1 >= loc2)
			throw ibd_error("loc1 >= loc2 in read_interactions!");
		int k = (loc2*(loc2-1)/2 + loc1)*npar_per_pair;
		for (int j=0;j<npar_per_pair;j++)
			f >> result[k++];
	    f >> skip_rest_of_line >> eatcomment; 
	}
	if (!f.eof())
		throw ibd_error("Error reading file " + filename);

	return result;
}


ibd::PopulationType  ibd::read_param(int& nind,
				double& herit,
				double& fraction_missing, 
				bool& random_start,
				int& nchr,
				double& length,
				int& nloc_chr,
				double& win_sz,
				double& eff_dim_epi,
				const string& filename)
{
	ifstream f(filename.c_str());
	if (!f)
		throw ibd_error("Cannot open file " + filename);

	string str, pop_name;
	f >> eatcomment >> pop_name >> eatcomment;
	PopulationType popt(pop_name);
	f >> eatcomment;
	f >> nind >> eatcomment;
	f >> herit >> eatcomment;
	f >> fraction_missing >> eatcomment;
	f >> str >> eatcomment;
	if (str == "false")
		random_start = false;
	else if (str == "true")
		random_start = true;
	else
		throw ibd_error("Error reading file" + filename);
	f >> nchr >> eatcomment;
	f >> length >> eatcomment;
	f >> nloc_chr >> eatcomment;
	f >> win_sz >> eatcomment;
	f >> eff_dim_epi >> eatcomment;
	if (!f.eof())
		throw ibd_error("Error reading file " + filename);
	return popt;
}

/*
void ibd::read_chromosome_lengths(vector<double>& chr_length, ifstream &f)
{
	double l;
 	f >> eatcomment;
	while (f)
	{
		f >> l >> eatcomment;
		chr_length.push_back(l);
	}
	if (!f.eof())
		throw ibd_error("Error reading chromosome lenghts.");
}
*/

void ibd::read_markermap(LinkageMap& markermap, ifstream& f)
{
	skip_header(f);

	int chr;
	const string& str = "_";
	double pos;
	while (f >> chr)
	{
		f >> pos;
		markermap.push_back(Locus(stringify(chr),pos,str));		
	}
}

void ibd::sim(vector<double>& pheno,											// phenotypic values
		 vector<double>& geno_val,										// genotypic values
		 matrix<ObsGeno>& geno,                                         // marker observations
		 const QTLModel& QTLmodel, const LinkageMap& QTLmap,			// QTL info 
		 double sigma_env, const vector<double>& weight,				// env. var + weights 
		 const LinkageMap& markermap, double fr_miss,				    // marker info
		 vector<double> chr_length, int nind, PopulationType poptype)	// population info
{
	double pheno_value, geno_value, sigma;
	Genome genome;
	vector<Genotype> QTLgeno;
	vector<ObsGeno> locgeno;
	double sigma_env_sqr = sqr(sigma_env);
	for (int ind = 0; ind < nind; ind++)  // for each individual
	{
		// create genome (given a vector with chromosome lengths)
		genome = poptype.make_ind(chr_length);

		// find the QTL genotypes (given the positions of the QTL)
	    QTLgeno = genome.GetGenotype(QTLmap);

		// calculate the genotypic value (given the QTL genotypes)
		geno_value = QTLmodel(QTLgeno);

		// generate the phenotypic value
		sigma = sqrt(sigma_env_sqr/weight[ind]);
		pheno_value = geno_value + randnormal(0.0,sigma);

		// generate markerdata, given the markermap, fraction missing, and population type  
		locgeno = generate_markerdata(genome,markermap,fr_miss,poptype);

		// store the information of the new individual
		pheno.push_back(pheno_value);
		geno_val.push_back(geno_value);
		geno.push_back(locgeno);
	}		
}
