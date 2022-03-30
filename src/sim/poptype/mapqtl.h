#ifndef MAPQTL_HEADER
#define MAPQTL_HEADER

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <algorithm>
#include <limits>
#include "Loc.h"
#include "poptype.h"

namespace ibd
{

void make_map_file(std::string name, const LinkageMap& marker_map);
void make_pheno_file(std::string name, const std::vector<double>& y,
                     const std::vector<double>& w,
                     const std::vector<int>& env_factor);
void make_pheno_file(std::string name, const std::vector<double>& y,
                     const std::string& trait_name);
void make_loc_file(std::string name, const matrix<ObsGeno>& geno,
                   const PopulationType& poptype,
                   const LinkageMap& marker_map);

void make_MapQTL_files(std::string name,
                       const std::vector<double>& pheno,
                       const std::vector<double>& weight,
                       const std::vector<int>& env_factor,
                       const matrix<ObsGeno>& geno,
                       const LinkageMap& markermap,
                       const PopulationType& poptype);

void make_MapQTL_files(std::string name,
                       const std::vector<double>& pheno,
                       const std::string& trait_name,
                       const matrix<ObsGeno>& geno,
                       const LinkageMap& markermap,
                       const PopulationType& poptype);


std::string read_qua_header(std::vector<std::string>& column_name,
                            int& nind, int& ntrt, std::istream& inp);

std::string read_pheno_file(std::vector<std::string>& column_name,
                            matrix<std::string>& qua_info,
                            const std::string filename);
void read_map_file(LinkageMap&  markermap,
                   std::vector<std::string>& chr_name,
                   const std::string filename);
PopulationType read_loc_file(std::string& experiment_name,
                             std::map<std::string, std::vector<ObsGeno> >& marker_obs,
                             std::string filename);

PopulationType read_loc_file(std::string& experiment_name,
                             std::map<std::string, std::vector<ObsGeno> >& marker_obs,
                             int& nind,
                             std::string filename);


PopulationType read_MapQTLfiles(std::string& name_quant_trait,
                                std::vector<double>& pheno,
                                std::vector<double>& weight,
                                std::vector<int>& env_factor,
                                matrix<ObsGeno>& geno,
                                LinkageMap& markermap,
                                const std::string& locfile,
                                const std::string& quafile,
                                const std::string& mapfile);

PopulationType read_MapQTLfiles(std::string& name_quant_trait,
                                std::vector<double>& pheno,
                                std::vector<double>& weight,
                                std::vector<int>& env_factor,
                                matrix<ObsGeno>& geno,
                                LinkageMap& markermap,
                                const std::string& locfile,
                                const std::string& quafile,
                                const std::string& mapfile,
                                int col_pheno,
                                int col_weight,
                                int col_env);

PopulationType read_MapQTLfiles(std::vector<double>& pheno,
                                matrix<ObsGeno>& geno,
                                LinkageMap& markermap,
                                const std::string& locfile,
                                const std::string& quafile,
                                const std::string& mapfile,
                                const std::string& trait_name);


LinkageMap read_cov_file(const std::string& filename, LinkageMap& markermap);

bool check_obsgeno(const matrix<ObsGeno>& geno, const LinkageMap& markermap);

void remove_inconsistent_scores(matrix<ObsGeno>& geno, const LinkageMap& markermap,
                                const ObsGeno& U);

bool get_geno_and_markermap(matrix<ObsGeno>& geno, LinkageMap& markermap,
                            const std::map<std::string, std::vector<ObsGeno> > marker_obs,
                            const LinkageMap& complete_markermap, int nind);



std::pair<int,int> dimension_table(const std::string& filename);

PopulationType read_loc_file_format2(std::map<std::string,std::vector<ObsGeno> > & marker_obs,
                                     int& nind, const std::string& filename);

void make_loc_file_format2(std::string name, const matrix<ObsGeno>& geno,
                           const PopulationType& poptype,
                           const LinkageMap& markermap);

void make_loc_file_format2(std::string name, const matrix<ObsGeno>& geno,
                           const std::vector<std::string>& ID,
                           const PopulationType& poptype,
                           const LinkageMap& markermap);

void convert2new_format_loc_file(const std::vector<std::string>& ID,
                                 const std::string& old_locfile,
                                 const std::string& new_locfile);


}

/*
 pop_base * read_MapQTLfiles(vector<double>& pheno,
 matrix<ObsGeno>& geno,
 LinkageMap& markermap,
 int& dimE,
 const string& locfile,
 const string& quafile,
 const string& mapfile);


 // test functies voor lezen MapQTL files:
 int test_qua();
 int test_loc();
 int test_map();
 int test_reading_MapQTL();
 */

#endif
