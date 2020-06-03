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
//' @param poptype A character string indicating the type of population. One of
//' ....
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
  // Create result list: map + markers.
  P.attr("dim") = IntegerVector {M, npop, ncol_prob};
  P.attr("dimnames") = List::create(posNames, offspring, parentNames);
  List res = List::create(Named("map") = map, Named("markers") = P);
  res.attr("class") = "statgenIBD";
  return res;
}
