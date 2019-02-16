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
// rcpp_dat_small
void rcpp_dat_small(arma::mat& d, const Rcpp::List& cfg, const int look, const double l0, const double l1);
RcppExport SEXP _orvacsim_rcpp_dat_small(SEXP dSEXP, SEXP cfgSEXP, SEXP lookSEXP, SEXP l0SEXP, SEXP l1SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type cfg(cfgSEXP);
    Rcpp::traits::input_parameter< const int >::type look(lookSEXP);
    Rcpp::traits::input_parameter< const double >::type l0(l0SEXP);
    Rcpp::traits::input_parameter< const double >::type l1(l1SEXP);
    rcpp_dat_small(d, cfg, look, l0, l1);
    return R_NilValue;
END_RCPP
}
// rcpp_clin
Rcpp::List rcpp_clin(arma::mat& d, const Rcpp::List& cfg, const int look);
RcppExport SEXP _orvacsim_rcpp_clin(SEXP dSEXP, SEXP cfgSEXP, SEXP lookSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type cfg(cfgSEXP);
    Rcpp::traits::input_parameter< const int >::type look(lookSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_clin(d, cfg, look));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_clin_interim_ppos
Rcpp::List rcpp_clin_interim_ppos(arma::mat& d_new, const arma::mat& m, const int nimpute, const int look, const Rcpp::List& cfg);
RcppExport SEXP _orvacsim_rcpp_clin_interim_ppos(SEXP d_newSEXP, SEXP mSEXP, SEXP nimputeSEXP, SEXP lookSEXP, SEXP cfgSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type d_new(d_newSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const int >::type nimpute(nimputeSEXP);
    Rcpp::traits::input_parameter< const int >::type look(lookSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type cfg(cfgSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_clin_interim_ppos(d_new, m, nimpute, look, cfg));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_cens
Rcpp::List rcpp_cens(const arma::mat& d_new, const arma::vec& visits, const int i, const int look, const Rcpp::List& cfg);
RcppExport SEXP _orvacsim_rcpp_cens(SEXP d_newSEXP, SEXP visitsSEXP, SEXP iSEXP, SEXP lookSEXP, SEXP cfgSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type d_new(d_newSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type visits(visitsSEXP);
    Rcpp::traits::input_parameter< const int >::type i(iSEXP);
    Rcpp::traits::input_parameter< const int >::type look(lookSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type cfg(cfgSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_cens(d_new, visits, i, look, cfg));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_cens_interim
Rcpp::List rcpp_cens_interim(const arma::mat& d_new, const arma::vec& visits, const int i, const int look, const Rcpp::List& cfg);
RcppExport SEXP _orvacsim_rcpp_cens_interim(SEXP d_newSEXP, SEXP visitsSEXP, SEXP iSEXP, SEXP lookSEXP, SEXP cfgSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type d_new(d_newSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type visits(visitsSEXP);
    Rcpp::traits::input_parameter< const int >::type i(iSEXP);
    Rcpp::traits::input_parameter< const int >::type look(lookSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type cfg(cfgSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_cens_interim(d_new, visits, i, look, cfg));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_cens_final
Rcpp::List rcpp_cens_final(const arma::mat& d_new, const arma::vec& visits, const int i, const int look, const Rcpp::List& cfg);
RcppExport SEXP _orvacsim_rcpp_cens_final(SEXP d_newSEXP, SEXP visitsSEXP, SEXP iSEXP, SEXP lookSEXP, SEXP cfgSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type d_new(d_newSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type visits(visitsSEXP);
    Rcpp::traits::input_parameter< const int >::type i(iSEXP);
    Rcpp::traits::input_parameter< const int >::type look(lookSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type cfg(cfgSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_cens_final(d_new, visits, i, look, cfg));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_visits
arma::vec rcpp_visits(const arma::mat& d_new, const int i, const int look, const Rcpp::List& cfg);
RcppExport SEXP _orvacsim_rcpp_visits(SEXP d_newSEXP, SEXP iSEXP, SEXP lookSEXP, SEXP cfgSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type d_new(d_newSEXP);
    Rcpp::traits::input_parameter< const int >::type i(iSEXP);
    Rcpp::traits::input_parameter< const int >::type look(lookSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type cfg(cfgSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_visits(d_new, i, look, cfg));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_clin_set_obst
Rcpp::List rcpp_clin_set_obst(arma::mat& d, const Rcpp::List& cfg, const int look);
RcppExport SEXP _orvacsim_rcpp_clin_set_obst(SEXP dSEXP, SEXP cfgSEXP, SEXP lookSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type cfg(cfgSEXP);
    Rcpp::traits::input_parameter< const int >::type look(lookSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_clin_set_obst(d, cfg, look));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_clin_interim_post
void rcpp_clin_interim_post(arma::mat& m, const int n_uncen_0, const double tot_obst_0, const int n_uncen_1, const double tot_obst_1, const int post_draw, const Rcpp::List& cfg);
RcppExport SEXP _orvacsim_rcpp_clin_interim_post(SEXP mSEXP, SEXP n_uncen_0SEXP, SEXP tot_obst_0SEXP, SEXP n_uncen_1SEXP, SEXP tot_obst_1SEXP, SEXP post_drawSEXP, SEXP cfgSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const int >::type n_uncen_0(n_uncen_0SEXP);
    Rcpp::traits::input_parameter< const double >::type tot_obst_0(tot_obst_0SEXP);
    Rcpp::traits::input_parameter< const int >::type n_uncen_1(n_uncen_1SEXP);
    Rcpp::traits::input_parameter< const double >::type tot_obst_1(tot_obst_1SEXP);
    Rcpp::traits::input_parameter< const int >::type post_draw(post_drawSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type cfg(cfgSEXP);
    rcpp_clin_interim_post(m, n_uncen_0, tot_obst_0, n_uncen_1, tot_obst_1, post_draw, cfg);
    return R_NilValue;
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
// atest
void atest();
RcppExport SEXP _orvacsim_atest() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    atest();
    return R_NilValue;
END_RCPP
}
// rcpp_gamma
arma::vec rcpp_gamma(const int n, const double a, const double b);
RcppExport SEXP _orvacsim_rcpp_gamma(SEXP nSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_gamma(n, a, b));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_test_1
void rcpp_test_1(arma::mat& d);
RcppExport SEXP _orvacsim_rcpp_test_1(SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type d(dSEXP);
    rcpp_test_1(d);
    return R_NilValue;
END_RCPP
}
// rcpp_test_sub_1
void rcpp_test_sub_1(arma::mat& d);
RcppExport SEXP _orvacsim_rcpp_test_sub_1(SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type d(dSEXP);
    rcpp_test_sub_1(d);
    return R_NilValue;
END_RCPP
}
// rcpp_test_2
arma::mat rcpp_test_2(const arma::mat& d);
RcppExport SEXP _orvacsim_rcpp_test_2(SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_test_2(d));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_test_sub_2
arma::mat rcpp_test_sub_2(arma::mat& d);
RcppExport SEXP _orvacsim_rcpp_test_sub_2(SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_test_sub_2(d));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_orvacsim_rcpp_dat", (DL_FUNC) &_orvacsim_rcpp_dat, 1},
    {"_orvacsim_rcpp_dat_small", (DL_FUNC) &_orvacsim_rcpp_dat_small, 5},
    {"_orvacsim_rcpp_clin", (DL_FUNC) &_orvacsim_rcpp_clin, 3},
    {"_orvacsim_rcpp_clin_interim_ppos", (DL_FUNC) &_orvacsim_rcpp_clin_interim_ppos, 5},
    {"_orvacsim_rcpp_cens", (DL_FUNC) &_orvacsim_rcpp_cens, 5},
    {"_orvacsim_rcpp_cens_interim", (DL_FUNC) &_orvacsim_rcpp_cens_interim, 5},
    {"_orvacsim_rcpp_cens_final", (DL_FUNC) &_orvacsim_rcpp_cens_final, 5},
    {"_orvacsim_rcpp_visits", (DL_FUNC) &_orvacsim_rcpp_visits, 4},
    {"_orvacsim_rcpp_clin_set_obst", (DL_FUNC) &_orvacsim_rcpp_clin_set_obst, 3},
    {"_orvacsim_rcpp_clin_interim_post", (DL_FUNC) &_orvacsim_rcpp_clin_interim_post, 7},
    {"_orvacsim_rcpp_immu", (DL_FUNC) &_orvacsim_rcpp_immu, 3},
    {"_orvacsim_rcpp_n_obs", (DL_FUNC) &_orvacsim_rcpp_n_obs, 5},
    {"_orvacsim_rcpp_lnsero", (DL_FUNC) &_orvacsim_rcpp_lnsero, 2},
    {"_orvacsim_rcpp_immu_interim_post", (DL_FUNC) &_orvacsim_rcpp_immu_interim_post, 4},
    {"_orvacsim_rcpp_immu_interim_ppos", (DL_FUNC) &_orvacsim_rcpp_immu_interim_ppos, 7},
    {"_orvacsim_atest", (DL_FUNC) &_orvacsim_atest, 0},
    {"_orvacsim_rcpp_gamma", (DL_FUNC) &_orvacsim_rcpp_gamma, 3},
    {"_orvacsim_rcpp_test_1", (DL_FUNC) &_orvacsim_rcpp_test_1, 1},
    {"_orvacsim_rcpp_test_sub_1", (DL_FUNC) &_orvacsim_rcpp_test_sub_1, 1},
    {"_orvacsim_rcpp_test_2", (DL_FUNC) &_orvacsim_rcpp_test_2, 1},
    {"_orvacsim_rcpp_test_sub_2", (DL_FUNC) &_orvacsim_rcpp_test_sub_2, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_orvacsim(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
