#include <string>
#include <RcppArmadillo.h>

#include "mainR.h"
#include "matvec.h"
#include "Loc.h"
#include "read_map.h"

#include "popt.h"
#include "util_genetics.h"

using namespace Rcpp;
using namespace std;
using namespace ibd;

//' Calculate IBD probabilities
//'
//' Calculate IBD probabilities for different types of populations.
//'
//' IBD probabilities can be calculated for many different types of populations.
//' In the following table all supported populations are listed. Note that the
//' value of x in the population types is variable, with its maximum value
//' depicted in the last column.
//'
//' | __Population type__ | __Cross__ | __Description__ | __max. x__ |
//' | ------ | ----------------- | -------------------------------------- | --- |
//' | DH | biparental | doubled haploid population | |
//' | Fx | biparental | Fx population (F1, followed by x-1 generations of selfing) | 8 |
//' | FxDH | biparental | Fx, followed by DH generation | 8 |
//' | BCx | biparental | backcross, second parent is recurrent parent | 9 |
//' | BCxDH | biparental | BCx, followed by DH generation | 9 |
//' | BC1Sx | biparental | BC1, followed by x generations of selfing | 7 |
//' | BC1SxDH | biparental | BC1, followed by x generations of selfing and DH | 6 |
//' | C3 | three-way | three way cross: (AxB) x C |  |
//' | C3DH | three-way | C3, followed by DH generation |  |
//' | C3Sx | three-way | C3, followed by x generation of selfing | 7 |
//' | C3SxDH | three-way | C3, followed by x generation of selfing and DH generation | 6 |
//' | C4 | four-way | four-way cross: (AxB) x (C x D)	| |
//' | C4DH | four-way | C4, followed by DH generation |  |
//' | C4Sx | four-way | C4, followed by x generation of selfing | 6 |
//' | C4SxDH | four-way | C4, followed by x generation of selfing and DH generation | 6 |
//'
//' @param poptype A character string indicating the type of population. One of
//' DH, Fx, FxDH, BCx, BCxDH, BC1Sx, BC1SxDH, C3, C3DH, C3Sx, C3SxDH, C4, C4DH,
//' C4Sx, C4SxDH (see Details)
//' @param locfile A character string indicating the location of the file with
//' genotypic information for the population. The file should be in
//' tab-delimited format with a header containing marker names.
//' @param mapfile A character string indicating the location of the map file
//' for the population. The file should be in tab-delimited format. It should
//' consist of exactly three columns, marker, chromosome and position. There
//' should be no header.
//' @param evalposfile An optional character string indicating a tab separated
//' .txt file containing combinations of chromosomes and positions to which the
//' calculations should be limited.
//' @param evaldist An optional numerical value indicating the maximum
//' distance for in between marker evaluation positions.
//'
//' @return An object of class \code{calcIBD}, a \code{list} with four elements,
//' \describe{
//' \item{map}{a \code{data.frame} with chromosome and position of the markers.}
//' \item{markers}{a 3-dimensional \code{array} of IBD probabilities with
//' markers, genotypes and  parents as array dimensions.}
//' \item{poptype}{the population type}
//' \item{multiCross}{a boolean indicating if multiple crosses have been
//' combined in the \code{calcIBD} object}
//' }
//'
//' @examples
//' ## Compute IBD probabilities for Steptoe Morex.
//' SxMIBD <- calcIBD(poptype = "DH",
//'                   locfile = system.file("extdata", "SxM_geno.txt",
//'                                         package = "statgenIBD"),
//'                   mapfile = system.file("extdata", "SxM_map.txt",
//'                                         package = "statgenIBD"))
//'
//' ## Check summary.
//' summary(SxMIBD)
//'
//' ## Compute IBD probabilities for Steptoe Morex.
//' ## Add extra evaluation positions so positions are at most 5 apart.
//' SxMIBD_Ext <- calcIBD(poptype = "DH",
//'                       locfile = system.file("extdata", "SxM_geno.txt",
//'                                             package = "statgenIBD"),
//'                       mapfile = system.file("extdata", "SxM_map.txt",
//'                                             package = "statgenIBD"),
//'                       evaldist = 5)
//'
//' ## Check summary.
//' summary(SxMIBD_Ext)
//'
//' @export
// [[Rcpp::export]]
List calcIBD(CharacterVector& poptype,
             CharacterVector& locfile,
             CharacterVector& mapfile,
             Nullable<CharacterVector&> evalposfile = R_NilValue,
             Nullable<NumericVector&> evaldist = R_NilValue)
{
  string _poptype = Rcpp::as<std::string>(poptype);
  // only to check poptype has correct format:
  const pop_base *popt = init_pop(_poptype);

  int x;
  bool isDH = _poptype.find("DH") != std::string::npos;
  bool isBC = match(x, _poptype, "BCx");
  LinkageMap positions;
  arma::cube prob;
  vector<string> parents, offspring;
  double max_step_size = -1;
  if (evaldist.isNotNull())
    max_step_size = Rcpp::as<double>(evaldist);
  std::string _evalposfile;
  if (evalposfile.isNotNull())
    _evalposfile = Rcpp::as<std::string>(evalposfile);
  else
    _evalposfile = "";
  try
  {
    main_pedigreeR(prob, parents, offspring, positions,
                   _poptype,
                   Rcpp::as<std::string>(locfile),
                   Rcpp::as<std::string>(mapfile),
                   _evalposfile,
                   max_step_size);
  }
  catch (ibd_error& e)
  {
    forward_exception_to_r(e);
  }
  catch (std::exception& e)
  {
    forward_exception_to_r(e);
  }
  catch (...)
  {
    ::Rf_error("c++ exception (unknown reason)");
  }
  const int npar = parents.size();
  const int M = positions.size();
  unsigned int nSlices = prob.n_slices;
  // Remove slices with only zero
  for (arma::uword i = 0; i < nSlices; i++) {
    if (all(vectorise(prob.slice(nSlices - i - 1)) == 0)) {
      prob.shed_slice(nSlices - i - 1);
    }
  }
  // Construct vector of names for parents.
  CharacterVector parentNames (0);
  for (int i = 0; i < npar; i++)
  {
    if (!(isBC && i == 0))
    {
      parentNames.push_back("p" + parents[i]);
    }
  }
  if (!isDH)
  {
    if (npar == 2)
    {
      parentNames.push_back("p" + parents[0] + parents[1]);
    } else if (npar == 3)
    {
      parentNames.push_back("p" + parents[0] + parents[1]);
      parentNames.push_back("p" + parents[0] + parents[2]);
    } else if (npar == 4) {
      parentNames.push_back("p" + parents[0] + parents[2]);
      parentNames.push_back("p" + parents[0] + parents[3]);
      parentNames.push_back("p" + parents[1] + parents[2]);
      parentNames.push_back("p" + parents[1] + parents[3]);
    }
  }
  // Construct map file from positions.
  CharacterVector posNames = CharacterVector(M);
  IntegerVector chr = IntegerVector(M);
  NumericVector pos = NumericVector(M);
  for (int m = 0; m < M; m++) {
    posNames(m) = positions[m].GetName();
    chr(m) = positions[m].GetChr();
    pos(m) = positions[m].GetPosition();
  }
  DataFrame map = DataFrame::create(Named("chr") = chr, Named("pos") = pos);
  map.attr("row.names") = posNames;
  // Reshape prob to 3D array and add names to dimensions.
  NumericVector P = wrap(prob);
  P.attr("dimnames") = List::create(posNames, offspring, parentNames);
  // Create result list: map + markers.
  List res = List::create(Named("map") = map, Named("markers") = P,
                          Named("poptype") = poptype,
                          Named("multiCross") = false);
  res.attr("class") = "calcIBD";
  return res;
}
