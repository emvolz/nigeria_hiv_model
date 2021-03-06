// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// f_t
double f_t(double t, List ns_parmlist);
RcppExport SEXP gfmodel11_f_t(SEXP tSEXP, SEXP ns_parmlistSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< List >::type ns_parmlist(ns_parmlistSEXP);
    __result = Rcpp::wrap(f_t(t, ns_parmlist));
    return __result;
END_RCPP
}
// f_t_msm
double f_t_msm(double t, List ns_parmlist);
RcppExport SEXP gfmodel11_f_t_msm(SEXP tSEXP, SEXP ns_parmlistSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< List >::type ns_parmlist(ns_parmlistSEXP);
    __result = Rcpp::wrap(f_t_msm(t, ns_parmlist));
    return __result;
END_RCPP
}
// source_size
double source_size(double t, List parms);
RcppExport SEXP gfmodel11_source_size(SEXP tSEXP, SEXP parmsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< List >::type parms(parmsSEXP);
    __result = Rcpp::wrap(source_size(t, parms));
    return __result;
END_RCPP
}
// F_from_state
NumericMatrix F_from_state(double t, NumericVector y, List parms);
RcppExport SEXP gfmodel11_F_from_state(SEXP tSEXP, SEXP ySEXP, SEXP parmsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< List >::type parms(parmsSEXP);
    __result = Rcpp::wrap(F_from_state(t, y, parms));
    return __result;
END_RCPP
}
// G_from_state
NumericMatrix G_from_state(double t, NumericVector y, List parms);
RcppExport SEXP gfmodel11_G_from_state(SEXP tSEXP, SEXP ySEXP, SEXP parmsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< List >::type parms(parmsSEXP);
    __result = Rcpp::wrap(G_from_state(t, y, parms));
    return __result;
END_RCPP
}
// dydt
List dydt(double t, NumericVector y, List parms);
RcppExport SEXP gfmodel11_dydt(SEXP tSEXP, SEXP ySEXP, SEXP parmsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< List >::type parms(parmsSEXP);
    __result = Rcpp::wrap(dydt(t, y, parms));
    return __result;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP gfmodel11_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    __result = Rcpp::wrap(rcpp_hello_world());
    return __result;
END_RCPP
}
