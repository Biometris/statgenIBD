// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// calcIBD
List calcIBD(CharacterVector& poptype, CharacterVector& locfile, CharacterVector& mapfile, Nullable<DataFrame&> evalpos, Nullable<NumericVector&> evaldist, const bool& verbose);
RcppExport SEXP _statgenIBD_calcIBD(SEXP poptypeSEXP, SEXP locfileSEXP, SEXP mapfileSEXP, SEXP evalposSEXP, SEXP evaldistSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector& >::type poptype(poptypeSEXP);
    Rcpp::traits::input_parameter< CharacterVector& >::type locfile(locfileSEXP);
    Rcpp::traits::input_parameter< CharacterVector& >::type mapfile(mapfileSEXP);
    Rcpp::traits::input_parameter< Nullable<DataFrame&> >::type evalpos(evalposSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector&> >::type evaldist(evaldistSEXP);
    Rcpp::traits::input_parameter< const bool& >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(calcIBD(poptype, locfile, mapfile, evalpos, evaldist, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_statgenIBD_calcIBD", (DL_FUNC) &_statgenIBD_calcIBD, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_statgenIBD(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
