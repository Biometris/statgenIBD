#ifndef INPUT_SIMQTL_HEADER
#define INPUT_SIMQTL_HEADER

#include "util_genetics.h"
#include "Sim.h"
#include "GenoValue.h"
#include "MarkerType.h"
#include "Loc.h"
#include "CommandFile.h"

ibd::Commands read_input_file(const std::string& filename);
//void read_seed(long int& start_seed, const ibd::Commands& commands);
std::vector<PopProp> read_pop(const ibd::Commands& commands);
void read_number_markers(std::vector<int>& nr_markers_per_chr, std::istream& f);
void read_genome(std::vector<double>& chr_length, const ibd::Commands& commands);
std::map<Locus, std::vector<double> > read_QTLs(const ibd::Commands& commands);
std::map<std::string,std::string> read_inbfnd(const ibd::Commands& commands,unsigned int nqtl);

ibd::matrix<double> read_epi(const ibd::Commands& commands,
                             const std::map< Locus,std::vector<double> >& QTLs);

void read_marker(LinkageMap& Markermap,
                 const std::string& filename,
                 std::vector<double> chr_length,
                 std::vector<double> nloc_chr,
                 const int& nr_alleles); //,double& fr_miss,
//const ibd::Commands& commands);

std::vector<IndProp> read_pedigree(const ibd::Commands& commands);

void read_makeped(const ibd::Commands& commands,
                  const std::vector<std::string>& fnd_names);

void make_hybrids(const ibd::Commands& commands,
                  const std::map<std::string,ibd::Genome>& simpop,
                  const Phi& phi, double sigma);

void sim_hybridsfile(const std::string inputfile,
                     const std::string outputdir,
                     const std::string outputfile,
                     const std::vector<IndProp>& ped,
                     const std::map<std::string,ibd::Genome>& simpop,
                     const std::vector<MarkerType>& markertype,
                     const LinkageMap& markermap,
                     const Phi& phi, double sigma, bool FlexQTL) ;

#endif

