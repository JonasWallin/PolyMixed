// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// WoodburySpecial
Eigen::MatrixXd WoodburySpecial(Eigen::VectorXd D, Eigen::MatrixXd U);
RcppExport SEXP _PolyMixed_WoodburySpecial(SEXP DSEXP, SEXP USEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type D(DSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type U(USEXP);
    rcpp_result_gen = Rcpp::wrap(WoodburySpecial(D, U));
    return rcpp_result_gen;
END_RCPP
}
// CholUpdate
Eigen::MatrixXd CholUpdate(Eigen::MatrixXd R, Eigen::MatrixXd U);
RcppExport SEXP _PolyMixed_CholUpdate(SEXP RSEXP, SEXP USEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type R(RSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type U(USEXP);
    rcpp_result_gen = Rcpp::wrap(CholUpdate(R, U));
    return rcpp_result_gen;
END_RCPP
}
// getTValues
Eigen::VectorXd getTValues(Eigen::VectorXi& index, Eigen::MatrixXd& Xu, Eigen::MatrixXd& uX, Eigen::MatrixXd& Sigma_U, Eigen::VectorXd& Y);
RcppExport SEXP _PolyMixed_getTValues(SEXP indexSEXP, SEXP XuSEXP, SEXP uXSEXP, SEXP Sigma_USEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXi& >::type index(indexSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type Xu(XuSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type uX(uXSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type Sigma_U(Sigma_USEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(getTValues(index, Xu, uX, Sigma_U, Y));
    return rcpp_result_gen;
END_RCPP
}
// getTValuesDiagS
Eigen::VectorXd getTValuesDiagS(Eigen::VectorXi& index, Eigen::MatrixXd& Xu, Eigen::MatrixXd& uX, Eigen::VectorXd& Sigma_U, Eigen::VectorXd& Y);
RcppExport SEXP _PolyMixed_getTValuesDiagS(SEXP indexSEXP, SEXP XuSEXP, SEXP uXSEXP, SEXP Sigma_USEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXi& >::type index(indexSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type Xu(XuSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type uX(uXSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type Sigma_U(Sigma_USEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(getTValuesDiagS(index, Xu, uX, Sigma_U, Y));
    return rcpp_result_gen;
END_RCPP
}
// crossSelect
void crossSelect(Eigen::Map<Eigen::MatrixXd> X1, Eigen::Map<Eigen::MatrixXd> X2, IntegerVector cross);
RcppExport SEXP _PolyMixed_crossSelect(SEXP X1SEXP, SEXP X2SEXP, SEXP crossSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cross(crossSEXP);
    crossSelect(X1, X2, cross);
    return R_NilValue;
END_RCPP
}
// mixing_population
void mixing_population(Eigen::Map<Eigen::MatrixXd> X1, Eigen::Map<Eigen::MatrixXd> X2, Eigen::Map<Eigen::MatrixXi> Z1, Eigen::Map<Eigen::MatrixXi> Z2, int n_gen, double lambda_cross);
RcppExport SEXP _PolyMixed_mixing_population(SEXP X1SEXP, SEXP X2SEXP, SEXP Z1SEXP, SEXP Z2SEXP, SEXP n_genSEXP, SEXP lambda_crossSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXi> >::type Z1(Z1SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXi> >::type Z2(Z2SEXP);
    Rcpp::traits::input_parameter< int >::type n_gen(n_genSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_cross(lambda_crossSEXP);
    mixing_population(X1, X2, Z1, Z2, n_gen, lambda_cross);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_PolyMixed_WoodburySpecial", (DL_FUNC) &_PolyMixed_WoodburySpecial, 2},
    {"_PolyMixed_CholUpdate", (DL_FUNC) &_PolyMixed_CholUpdate, 2},
    {"_PolyMixed_getTValues", (DL_FUNC) &_PolyMixed_getTValues, 5},
    {"_PolyMixed_getTValuesDiagS", (DL_FUNC) &_PolyMixed_getTValuesDiagS, 5},
    {"_PolyMixed_crossSelect", (DL_FUNC) &_PolyMixed_crossSelect, 3},
    {"_PolyMixed_mixing_population", (DL_FUNC) &_PolyMixed_mixing_population, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_PolyMixed(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
