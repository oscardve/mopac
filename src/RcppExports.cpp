// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// readFASTQC
void readFASTQC(const std::string library, const std::string id, const int spacer_start, const int spacer_length, const bool rev, const bool comp, const StringVector fastqIN);
RcppExport SEXP _MoPAC_readFASTQC(SEXP librarySEXP, SEXP idSEXP, SEXP spacer_startSEXP, SEXP spacer_lengthSEXP, SEXP revSEXP, SEXP compSEXP, SEXP fastqINSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type library(librarySEXP);
    Rcpp::traits::input_parameter< const std::string >::type id(idSEXP);
    Rcpp::traits::input_parameter< const int >::type spacer_start(spacer_startSEXP);
    Rcpp::traits::input_parameter< const int >::type spacer_length(spacer_lengthSEXP);
    Rcpp::traits::input_parameter< const bool >::type rev(revSEXP);
    Rcpp::traits::input_parameter< const bool >::type comp(compSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type fastqIN(fastqINSEXP);
    readFASTQC(library, id, spacer_start, spacer_length, rev, comp, fastqIN);
    return R_NilValue;
END_RCPP
}
// readFASTQC_shiny
void readFASTQC_shiny(const std::string library, const std::string id, const int spacer_start, const int spacer_length, const bool rev, const bool comp, const StringVector fastqIN, Rcpp::Function shinyF);
RcppExport SEXP _MoPAC_readFASTQC_shiny(SEXP librarySEXP, SEXP idSEXP, SEXP spacer_startSEXP, SEXP spacer_lengthSEXP, SEXP revSEXP, SEXP compSEXP, SEXP fastqINSEXP, SEXP shinyFSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type library(librarySEXP);
    Rcpp::traits::input_parameter< const std::string >::type id(idSEXP);
    Rcpp::traits::input_parameter< const int >::type spacer_start(spacer_startSEXP);
    Rcpp::traits::input_parameter< const int >::type spacer_length(spacer_lengthSEXP);
    Rcpp::traits::input_parameter< const bool >::type rev(revSEXP);
    Rcpp::traits::input_parameter< const bool >::type comp(compSEXP);
    Rcpp::traits::input_parameter< const StringVector >::type fastqIN(fastqINSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type shinyF(shinyFSEXP);
    readFASTQC_shiny(library, id, spacer_start, spacer_length, rev, comp, fastqIN, shinyF);
    return R_NilValue;
END_RCPP
}
// get_quantile
NumericMatrix get_quantile(NumericMatrix M);
RcppExport SEXP _MoPAC_get_quantile(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(get_quantile(M));
    return rcpp_result_gen;
END_RCPP
}
// getPValues
DataFrame getPValues(StringVector genes, StringVector sgrnas, NumericVector values, double alpha);
RcppExport SEXP _MoPAC_getPValues(SEXP genesSEXP, SEXP sgrnasSEXP, SEXP valuesSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< StringVector >::type genes(genesSEXP);
    Rcpp::traits::input_parameter< StringVector >::type sgrnas(sgrnasSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type values(valuesSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(getPValues(genes, sgrnas, values, alpha));
    return rcpp_result_gen;
END_RCPP
}
// RRA_1tail
DataFrame RRA_1tail(StringVector genes, StringVector sgrnas, NumericVector values);
RcppExport SEXP _MoPAC_RRA_1tail(SEXP genesSEXP, SEXP sgrnasSEXP, SEXP valuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< StringVector >::type genes(genesSEXP);
    Rcpp::traits::input_parameter< StringVector >::type sgrnas(sgrnasSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type values(valuesSEXP);
    rcpp_result_gen = Rcpp::wrap(RRA_1tail(genes, sgrnas, values));
    return rcpp_result_gen;
END_RCPP
}
// get_sorted
NumericMatrix get_sorted(NumericMatrix M);
RcppExport SEXP _MoPAC_get_sorted(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(get_sorted(M));
    return rcpp_result_gen;
END_RCPP
}
// get_sorted1
NumericMatrix get_sorted1(NumericVector V, int size);
RcppExport SEXP _MoPAC_get_sorted1(SEXP VSEXP, SEXP sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type V(VSEXP);
    Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(get_sorted1(V, size));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MoPAC_readFASTQC", (DL_FUNC) &_MoPAC_readFASTQC, 7},
    {"_MoPAC_readFASTQC_shiny", (DL_FUNC) &_MoPAC_readFASTQC_shiny, 8},
    {"_MoPAC_get_quantile", (DL_FUNC) &_MoPAC_get_quantile, 1},
    {"_MoPAC_getPValues", (DL_FUNC) &_MoPAC_getPValues, 4},
    {"_MoPAC_RRA_1tail", (DL_FUNC) &_MoPAC_RRA_1tail, 3},
    {"_MoPAC_get_sorted", (DL_FUNC) &_MoPAC_get_sorted, 1},
    {"_MoPAC_get_sorted1", (DL_FUNC) &_MoPAC_get_sorted1, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_MoPAC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
