#include <string>
#include <Rcpp.h>

#include "mainR.h"
#include "matvec.h"
#include "Loc.h"
#include "read_map.h"

using namespace Rcpp;
using namespace std;
using namespace mbl;

//' Calculate IBD probabilities
//'
//' A long description
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
//' @param locfile A character string indicating the .loc file for the
//' population. The file should be in Flapjack format.
//' @param mapfile A character string indicating the .map file for the
//' population. The file should be in Flapjack format.
//' @param evalposfile An optional character string indicating a tab separated
//' .txt file containing combinations of chromosomes and positions to which the
//' calculations should be limited.
//' @param evaldist An optional numerical value indicating the maximum
//' distance for in between marker evaluation positions.
//'
//' @return A data.frame with IBD probabilities.
//'
//' @export
// [[Rcpp::export]]
List calcIBD(CharacterVector& poptype,
             CharacterVector& locfile,
             CharacterVector& mapfile,
             Nullable<CharacterVector&> evalposfile = R_NilValue,
             Nullable<NumericVector&> evaldist = R_NilValue)
{
  LinkageMap positions;
  matrix3D<double> prob;
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
                   Rcpp::as<std::string>(poptype),
                   Rcpp::as<std::string>(locfile),
                   Rcpp::as<std::string>(mapfile),
                   _evalposfile,
                   max_step_size);
  }
  catch (mblib_error& e)
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
  const int npop = offspring.size();
  const int M = positions.size();
  const int ncol_prob = prob.Dim3();
  NumericMatrix P(M * npop, ncol_prob);
  int k = 0;
  for (int i = 0; i < npop; i++)
  {
    for (int m = 0; m < M; m++)
    {
      for (int j = 0; j < ncol_prob; j++) {
        P(k, j) = prob[i][m][j];
      }
      k++;
    }
  }
  // Construct vector of names for parents.
  CharacterVector parentNames;
  if (npar == 2) {
    parentNames = CharacterVector::create("p" + parents[0], "p" + parents[1],
                                          "p" + parents[0] + parents[1]);
  } else if (npar==3) {
    parentNames = CharacterVector::create("p" + parents[0], "p" + parents[1],
                                          "p" + parents[2],
                                                       "p" + parents[0] + parents[1],

                                                                                 "p" + parents[0] + parents[2]);

  } else if (npar == 4) {
    parentNames = CharacterVector::create("p" + parents[0], "p" + parents[1],
                                          "p" + parents[2], "p" + parents[3],
                                                                         "p" + parents[0] + parents[2],
                                                                                                   "p" + parents[0] + parents[3],
                                                                                                                             "p" + parents[1] + parents[2],
                                                                                                                                                       "p" + parents[1] + parents[3]);
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
  // Reshape P to 3D array and add names to dimensions.
  P.attr("dim") = IntegerVector {M, npop, ncol_prob};
  P.attr("dimnames") = List::create(posNames, offspring, parentNames);
  // Create result list: map + markers.
  List res = List::create(Named("map") = map, Named("markers") = P);
  res.attr("class") = "calcIBD";
  return res;
}
