// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// rcpp_dat
arma::mat rcpp_dat(const Rcpp::List& cfg);
RcppExport SEXP _orvacsim_rcpp_dat(SEXP cfgSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type cfg(cfgSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_dat(cfg));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_immu
Rcpp::List rcpp_immu(const arma::mat& d, const Rcpp::List& cfg, const int look);
RcppExport SEXP _orvacsim_rcpp_immu(SEXP dSEXP, SEXP cfgSEXP, SEXP lookSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type cfg(cfgSEXP);
    Rcpp::traits::input_parameter< const int >::type look(lookSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_immu(d, cfg, look));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_n_obs
int rcpp_n_obs(const arma::mat& d, const int look, const Rcpp::NumericVector looks, const Rcpp::NumericVector months, const double info_delay);
RcppExport SEXP _orvacsim_rcpp_n_obs(SEXP dSEXP, SEXP lookSEXP, SEXP looksSEXP, SEXP monthsSEXP, SEXP info_delaySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const int >::type look(lookSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type looks(looksSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type months(monthsSEXP);
    Rcpp::traits::input_parameter< const double >::type info_delay(info_delaySEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_n_obs(d, look, looks, months, info_delay));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_lnsero
Rcpp::List rcpp_lnsero(const arma::mat& d, const int nobs);
RcppExport SEXP _orvacsim_rcpp_lnsero(SEXP dSEXP, SEXP nobsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const int >::type nobs(nobsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_lnsero(d, nobs));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_immu_interim_post
arma::mat rcpp_immu_interim_post(const arma::mat& d, const int nobs, const int post_draw, const Rcpp::List& lnsero);
RcppExport SEXP _orvacsim_rcpp_immu_interim_post(SEXP dSEXP, SEXP nobsSEXP, SEXP post_drawSEXP, SEXP lnseroSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const int >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< const int >::type post_draw(post_drawSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type lnsero(lnseroSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_immu_interim_post(d, nobs, post_draw, lnsero));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_immu_interim_ppos
Rcpp::List rcpp_immu_interim_ppos(const arma::mat& d, const arma::mat& m, const int nobs, const int nimpute, const int post_draw, const Rcpp::List& lnsero, const Rcpp::List& cfg);
RcppExport SEXP _orvacsim_rcpp_immu_interim_ppos(SEXP dSEXP, SEXP mSEXP, SEXP nobsSEXP, SEXP nimputeSEXP, SEXP post_drawSEXP, SEXP lnseroSEXP, SEXP cfgSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const int >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< const int >::type nimpute(nimputeSEXP);
    Rcpp::traits::input_parameter< const int >::type post_draw(post_drawSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type lnsero(lnseroSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type cfg(cfgSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_immu_interim_ppos(d, m, nobs, nimpute, post_draw, lnsero, cfg));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_censoring
arma::mat rcpp_censoring(const arma::mat d, const int look, const int trtstatus, const int iend, const float curmonth, const float surveillancemonths);
RcppExport SEXP _orvacsim_rcpp_censoring(SEXP dSEXP, SEXP lookSEXP, SEXP trtstatusSEXP, SEXP iendSEXP, SEXP curmonthSEXP, SEXP surveillancemonthsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type d(dSEXP);
    Rcpp::traits::input_parameter< const int >::type look(lookSEXP);
    Rcpp::traits::input_parameter< const int >::type trtstatus(trtstatusSEXP);
    Rcpp::traits::input_parameter< const int >::type iend(iendSEXP);
    Rcpp::traits::input_parameter< const float >::type curmonth(curmonthSEXP);
    Rcpp::traits::input_parameter< const float >::type surveillancemonths(surveillancemonthsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_censoring(d, look, trtstatus, iend, curmonth, surveillancemonths));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_orvacsim_rcpp_dat", (DL_FUNC) &_orvacsim_rcpp_dat, 1},
    {"_orvacsim_rcpp_immu", (DL_FUNC) &_orvacsim_rcpp_immu, 3},
    {"_orvacsim_rcpp_n_obs", (DL_FUNC) &_orvacsim_rcpp_n_obs, 5},
    {"_orvacsim_rcpp_lnsero", (DL_FUNC) &_orvacsim_rcpp_lnsero, 2},
    {"_orvacsim_rcpp_immu_interim_post", (DL_FUNC) &_orvacsim_rcpp_immu_interim_post, 4},
    {"_orvacsim_rcpp_immu_interim_ppos", (DL_FUNC) &_orvacsim_rcpp_immu_interim_ppos, 7},
    {"_orvacsim_rcpp_censoring", (DL_FUNC) &_orvacsim_rcpp_censoring, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_orvacsim(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
