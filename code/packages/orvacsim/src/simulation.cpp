
#include <RcppDist.h>
// [[Rcpp::depends(RcppDist)]]

#include <cmath>
#include <algorithm>


//#include <mcmc.hpp>

// column indices
#define COL_ID            0
#define COL_TRT           1
#define COL_ACCRT         2
#define COL_AGE           3

#define COL_SEROT2        4
#define COL_SEROT3        5
#define COL_PROBT3        6

#define COL_EVTT          7

#define COL_FU1           8
#define COL_FU2           9

//#define COL_CURAGE        10
//#define COL_CENT          11
#define COL_CEN           10
#define COL_OBST          11


#define NCOL              12


#define COL_THETA0        0
#define COL_THETA1        1
#define COL_DELTA         2


#define COL_LAMB0         0
#define COL_LAMB1         1
#define COL_RATIO         2



#define _DEBUG 1

#if _DEBUG
#define DBG( os, msg )                             \
(os) << "DBG: " << __FILE__ << "(" << __LINE__ << ") "\
     << msg << std::endl
#else
#define DBG( os, msg )
#endif

#define _INFO  0

#if _INFO
#define INFO( os, msg )                                \
   (os) << "INFO: " << __FILE__ << "(" << __LINE__ << ") "\
        << msg << std::endl
#else
#define INFO( os, msg )
#endif

// function prototypes
arma::mat rcpp_dat(const Rcpp::List& cfg);
void rcpp_dat_small(const arma::mat& d,
                         const Rcpp::List& cfg,
                         const int look,
                         const double l0,
                         const double l1);
Rcpp::List rcpp_clin(arma::mat& d, const Rcpp::List& cfg,
                     const int look);
Rcpp::List rcpp_cens(const arma::mat& d_new,
                     const arma::vec& visits,
                     const int i,
                     const int look,
                     const Rcpp::List& cfg);
Rcpp::List rcpp_cens_interim(const arma::mat& d_new,
                             const arma::vec& visits,
                             const int i,
                             const int look,
                             const Rcpp::List& cfg);
Rcpp::List rcpp_cens_final(const arma::mat& d_new,
                           const arma::vec& visits,
                           const int i,
                           const int look,
                           const Rcpp::List& cfg);
arma::vec rcpp_visits(const arma::mat& d_new,
                      const int i,
                      const int look,
                      const Rcpp::List& cfg);
Rcpp::List rcpp_clin_set_obst(arma::mat& d,
                              const Rcpp::List& cfg,
                              const int look);
void rcpp_clin_interim_post(arma::mat& m,
                            const int n_uncen_0,
                            const double tot_obst_0,
                            const int n_uncen_1,
                            const double tot_obst_1,
                            const int post_draw,
                            const Rcpp::List& cfg);
Rcpp::List rcpp_clin_interim_ppos(arma::mat& d_new,
                                  const arma::mat& m,
                                  const int nimpute,
                                  const int look,
                                  const Rcpp::List& cfg);
Rcpp::List rcpp_cens_interim_alt(const arma::mat& d_new,
                                 const int i,
                                 const int look,
                                 const Rcpp::List& cfg);

Rcpp::List rcpp_immu(const arma::mat& d, const Rcpp::List& cfg, const int look);
int rcpp_n_obs(const arma::mat& d,
               const int look,
               const Rcpp::NumericVector looks,
               const Rcpp::NumericVector months,
               const double info_delay);
Rcpp::List rcpp_lnsero(const arma::mat& d,
                       const int nobs);
void rcpp_immu_interim_post(const arma::mat& d,
                            arma::mat& m,
                             const int nobs,
                             const int post_draw,
                             const Rcpp::List& lnsero);
Rcpp::List rcpp_immu_interim_ppos(const arma::mat& d,
                             const arma::mat& m,
                             const int nobs,
                             const int nimpute,
                             const int post_draw,
                             const Rcpp::List& lnsero,
                             const Rcpp::List& cfg);

Rcpp::List rcpp_logrank(const arma::mat& d,
                        const int look,
                        const Rcpp::List& cfg);
void rcpp_outer(const arma::vec& z,
                const arma::vec& t,
                arma::mat& out);
arma::vec rcpp_gamma(const int n, const double a, const double b);
void rcpp_test_1(arma::mat& d);
void rcpp_test_sub_1(arma::mat& d);
arma::mat rcpp_test_2(const arma::mat& d) ;
arma::mat rcpp_test_sub_2(arma::mat& d);
// end function prototypes






// data generation



// [[Rcpp::export]]
arma::mat rcpp_dat(const Rcpp::List& cfg) {

  int n = cfg["nstop"];
  arma::mat d = arma::zeros(n, NCOL);
  double tpp = (double)cfg["months_per_person"];

  for(int i = 0; i < n; i++){

    d(i, COL_ID) = i+1;
    d(i, COL_TRT) = ((i-1)%2 == 0) ? 0 : 1;
    // simultaneous accrual of each next ctl/trt pair
    d(i, COL_ACCRT) = (i%2 == 0) ? ((i+1)*tpp)+tpp : (i+1)*tpp;

    d(i, COL_AGE) = r_truncnorm(cfg["age_months_mean"], cfg["age_months_sd"],
      cfg["age_months_lwr"], cfg["age_months_upr"]);

    d(i, COL_SEROT2) = R::rbinom(1, cfg["baselineprobsero"]);
    d(i, COL_SEROT3) = d(i, COL_SEROT2);
    d(i, COL_PROBT3) = d(i, COL_TRT) * (double)cfg["deltaserot3"];

    if(d(i, COL_SEROT2) == 0 && d(i, COL_TRT) == 1){
      d(i, COL_SEROT3) = R::rbinom(1, d(i, COL_PROBT3));
    }


    // tte - the paramaterisation of rexp uses SCALE NOTE RATE!!!!!!!!!!!
    // event time is the time from randomisation (not birth) at which first
    // medical presentation occurs
    if(d(i, COL_TRT) == 0){
      d(i, COL_EVTT) = R::rexp(1/(double)cfg["b0tte"])  ;
    } else {
      double beta = (double)cfg["b0tte"] + (double)cfg["b1tte"];
      d(i, COL_EVTT) = R::rexp(1/beta)  ;
    }

    // fu 1 and 2 times from time of accrual
    // fu 1 is between 14 and 21 days from accrual
    // fu 2 is between 28 and 55 days from accrual
    d(i, COL_FU1) = R::runif((double)cfg["fu1_lwr"], (double)cfg["fu1_upr"]);
    d(i, COL_FU2) = R::runif((double)cfg["fu2_lwr"], (double)cfg["fu2_upr"]);
  }

  d.col(COL_CEN) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["nstop"]));
  d.col(COL_OBST) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["nstop"]));

  return d;
}


// [[Rcpp::export]]
void rcpp_dat_small(arma::mat& d,
                    const Rcpp::List& cfg,
                    const int look,
                    const double l0,
                    const double l1) {

  Rcpp::NumericVector looks = cfg["looks"];
  Rcpp::NumericVector months = cfg["interimmnths"];
  int mylook = look - 1;
  int idxStart = looks[mylook];
  int n = cfg["nstop"];

  double tpp = (double)cfg["months_per_person"];

  //DBG(Rcpp::Rcout, "starting dat_small at " << idxStart );

  for(int i = idxStart; i < n; i++){

    //DBG(Rcpp::Rcout, "i " << i << "d(i, COL_EVTT) was " << d(i, COL_EVTT) );

    d(i, COL_AGE) = r_truncnorm(cfg["age_months_mean"], cfg["age_months_sd"],
      cfg["age_months_lwr"], cfg["age_months_upr"]);

    // tte - the paramaterisation of rexp uses SCALE NOTE RATE!!!!!!!!!!!
    if(d(i, COL_TRT) == 0){
      d(i, COL_EVTT) = R::rexp(1/l0)  ;
    } else {
      d(i, COL_EVTT) = R::rexp(1/l1)  ;
    }

    //DBG(Rcpp::Rcout, "i " << i << "d(i, COL_EVTT) now " << d(i, COL_EVTT) );

    // fu 1 and 2 times from time of accrual
    // fu 1 is between 14 and 21 days from accrual
    // fu 2 is between 28 and 55 days from accrual
    d(i, COL_FU1) = R::runif((double)cfg["fu1_lwr"], (double)cfg["fu1_upr"]);
    d(i, COL_FU2) = R::runif((double)cfg["fu2_lwr"], (double)cfg["fu2_upr"]);

  }

  d.col(COL_CEN) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["nstop"]));
  d.col(COL_OBST) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["nstop"]));

  return;
}






// clinical endpoint





// [[Rcpp::export]]
Rcpp::List rcpp_clin(arma::mat& d, const Rcpp::List& cfg,
                     const int look){

  int post_draw = (int)cfg["post_draw"];
  int n_uncen_0 = 0;
  int n_uncen_1 = 0;
  double tot_obst_0 = 0;
  double tot_obst_1 = 0;

  int mylook = look - 1;
  double ppos_max = 0;
  double fudge = 0.0001;

  Rcpp::NumericVector looks = cfg["looks"];
  Rcpp::NumericVector months = cfg["interimmnths"];
  Rcpp::List cens;
  Rcpp::List ret;
  Rcpp::List lsuffstat;

  arma::vec visits;
  arma::mat m = arma::zeros(post_draw , 3);

  d.col(COL_CEN) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["nstop"]));
  d.col(COL_OBST) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["nstop"]));


  // updates d(COL_CEN) and d(COL_OBST)
  lsuffstat = rcpp_clin_set_obst(d, cfg, look);

  n_uncen_0 = (double)lsuffstat["n_uncen_0"];
  tot_obst_0 = (double)lsuffstat["tot_obst_0"];
  n_uncen_1 = (double)lsuffstat["n_uncen_1"];
  tot_obst_1 = (double)lsuffstat["tot_obst_1"];


  // compute posterior
  // note the parameterisation of the gamma distribution is not the same as R
  // if the conjugate prior is gamma(k, q) then the posterior is:
  //
  // ~ gamma(k + sum(uncensored obs), q / (1 + q * total obs time))
  //
  // se ibrahim bayesian surv analysis and
  // https://cdn2.hubspot.net/hubfs/310840/VWO_SmartStats_technical_whitepaper.pdf

  // updates m
  rcpp_clin_interim_post(m, n_uncen_0, tot_obst_0,
                         n_uncen_1, tot_obst_1,
                         post_draw, cfg);

  arma::uvec tmp = arma::find(m.col(COL_RATIO) > 1);
  double ppos_n =  (double)tmp.n_elem / (double)post_draw;
  double mean_ratio =  arma::mean(m.col(COL_RATIO));

  // use posterior to do pp at final analysis
  // this involves simulating multiple datasets conditional on the
  // posterior and therefore requires multiple to the censoring funct

  arma::mat d_new = arma::mat(d);
  int nimpute = max(looks) - looks[mylook];

  // DBG(Rcpp::Rcout, " need to impute " << nimpute << " to get to end of trial");
  Rcpp::List lppos = rcpp_clin_interim_ppos(d_new, m, nimpute, look, cfg);

  if(_DEBUG == 1){
    ret = Rcpp::List::create(Rcpp::Named("ppn") = ppos_n,
                             Rcpp::Named("mean_ratio") = mean_ratio,
                             Rcpp::Named("ppmax") = (double)lppos["ppos"],
                             Rcpp::Named("ppmax_mean_ratio") = (double)lppos["mean_ratio"],
                             Rcpp::Named("n_uncen_0") = n_uncen_0,
                             Rcpp::Named("tot_obst_0") = tot_obst_0,
                             Rcpp::Named("n_uncen_1") = n_uncen_1,
                             Rcpp::Named("tot_obst_1") = tot_obst_1);

  } else {
    ret = Rcpp::List::create(Rcpp::Named("ppn") = ppos_n,
                             Rcpp::Named("mean_ratio") = mean_ratio,
                             Rcpp::Named("ppmax") = (double)lppos["ppos"],
                             Rcpp::Named("ppmax_mean_ratio") = (double)lppos["mean_ratio"],
                             Rcpp::Named("ppmax_sd_ratio") = (double)lppos["sd_ratio"]);
  }




  return ret;
}


// [[Rcpp::export]]
Rcpp::List rcpp_clin_interim_ppos(arma::mat& d_new,
                                  const arma::mat& m,
                                  const int nimpute,
                                  const int look,
                                  const Rcpp::List& cfg){


  int post_draw = (int)cfg["post_draw"];
  Rcpp::NumericVector looks = cfg["looks"];
  Rcpp::NumericVector months = cfg["interimmnths"];
  Rcpp::List lsuffstat;
  Rcpp::List llr;

  int win = 0;
  int mylook = look - 1;

  int n_uncen_0 = 0;
  double tot_obst_0 = 0;
  int n_uncen_1 = 0;
  double tot_obst_1 = 0;

  arma::vec postprob_ratio_gt1 = arma::zeros(post_draw);
  arma::vec mean_rat = arma::zeros(post_draw);
  arma::vec lrp = arma::zeros(post_draw);

  for(int i = 0; i < post_draw; i++){

    n_uncen_0 = 0;
    n_uncen_1 = 0;
    tot_obst_0 = 0;
    tot_obst_1 = 0;

    // updates records above current look in data using posterior values
    // changes evtt, fu1, fu2
    rcpp_dat_small(d_new, cfg, look, m[i, COL_LAMB0], m[i, COL_LAMB1]);

    lsuffstat = rcpp_clin_set_obst(d_new, cfg, looks.size());

    llr = rcpp_logrank(d_new, looks.size(), cfg);
    lrp(i) = (double)llr["pvalue"];

    n_uncen_0 = (double)lsuffstat["n_uncen_0"];
    tot_obst_0 = (double)lsuffstat["tot_obst_0"];
    n_uncen_1 = (double)lsuffstat["n_uncen_1"];
    tot_obst_1 = (double)lsuffstat["tot_obst_1"];

    // DBG(Rcpp::Rcout, "i " << i << " sufficient " << n_uncen_0
    //                       << "  " << n_uncen_1 << "  "
    //                       << tot_obst_0 << "  "  << tot_obst_1 );

    // compute posterior
    // see ibrahim bayesian surv analysis and
    // https://cdn2.hubspot.net/hubfs/310840/VWO_SmartStats_technical_whitepaper.pdf

    // updates m
    arma::mat m_new = arma::zeros(post_draw , 3);
    rcpp_clin_interim_post(m_new,
                           n_uncen_0, tot_obst_0,
                           n_uncen_1, tot_obst_1,
                           post_draw, cfg);

    // for(int j = 0; j < 10; j++){
    //  DBG(Rcpp::Rcout, "i " << i << " new posterior " << m_new(j, COL_RATIO));
    // }

    mean_rat(i) = arma::mean(m_new.col(COL_RATIO));

    // empirical posterior probability that ratio_lamb > 1
    arma::uvec tmp = arma::find(m_new.col(COL_RATIO) > 1);
    postprob_ratio_gt1(i) =  (double)tmp.n_elem / (double)post_draw;
    if(postprob_ratio_gt1(i) > (double)cfg["post_tte_thresh"]){
      win++;
    }
  }

  //DBG(Rcpp::Rcout, "pow " << std::pow((float)post_draw, 2.0) );

  double ppos = (double)win / (double)post_draw;

  Rcpp::List res ;
  if(_DEBUG == 1){

    res = Rcpp::List::create(Rcpp::Named("win") = win,
                                        Rcpp::Named("ppos") = ppos,
                                        Rcpp::Named("pp_ratio") = postprob_ratio_gt1,
                                        Rcpp::Named("mean_ratio") = arma::mean(mean_rat),
                                        Rcpp::Named("pvalue") = lrp);
  } else {

    res = Rcpp::List::create(Rcpp::Named("win") = win,
                             Rcpp::Named("ppos") = ppos,
                             Rcpp::Named("mean_ratio") = arma::mean(mean_rat),
                             Rcpp::Named("sd_ratio") = arma::stddev(mean_rat),
                             Rcpp::Named("pvalue") = lrp);
  }


  return res;

}


// [[Rcpp::export]]
Rcpp::List rcpp_cens(const arma::mat& d_new,
                     const arma::vec& visits,
                     const int i,
                     const int look,
                     const Rcpp::List& cfg) {

  Rcpp::List cens;

  Rcpp::NumericVector looks = cfg["looks"];
  Rcpp::NumericVector months = cfg["interimmnths"];

  int mylook = look - 1;
  int curmonth = months[mylook];

  double cen = 0;
  double obst = 0;

  double fudge = 0.0001;


  if(looks[mylook] != max(looks)){

    if((int)cfg["use_alt_censoring"] == 0){
      cens = rcpp_cens_interim(d_new, visits, i, look, cfg);
    } else {
      DBG(Rcpp::Rcout, "i " << i << " using alt censoring         : " << (int)cfg["use_alt_censoring"]);
      cens = rcpp_cens_interim_alt(d_new, i, look, cfg);
    }


  } else {
    cens = rcpp_cens_final(d_new, visits, i, look, cfg);
  }

  //DBG(Rcpp::Rcout, "i " << i << " cens indicator    : " << (double)cens["cen"]);
  //DBG(Rcpp::Rcout, "i " << i << " cens obst         : " << (double)cens["obst"]);

  return cens;
}


// [[Rcpp::export]]
Rcpp::List rcpp_cens_interim(const arma::mat& d_new,
                     const arma::vec& visits,
                     const int i,
                     const int look,
                     const Rcpp::List& cfg){

  Rcpp::List cens;

  Rcpp::NumericVector looks = cfg["looks"];
  Rcpp::NumericVector months = cfg["interimmnths"];

  int mylook = look - 1;
  int curmonth = months[mylook];

  double cen = NA_REAL;
  double obst = NA_REAL;

  double fudge = 0.0001;


  if(d_new(i, COL_ACCRT) <= months[mylook] + fudge){

    if(visits.n_elem == 0){

      if(i > looks[mylook]-1){
        cen = NA_REAL;
        obst = NA_REAL;
      } else {
        cen = 1;
        obst = months[mylook] - d_new(i, COL_ACCRT);
        obst = obst < 0 ? 0 : obst;
      }

      cens = Rcpp::List::create(Rcpp::Named("cen") = cen, Rcpp::Named("obst") = obst);
      //DBG(Rcpp::Rcout, "i " << i << " no visits     : " << obst);
      return cens;

    } else if (d_new(i, COL_ACCRT) <= (double)arma::max(visits)) {

      double age_at_this_analysis = d_new(i, COL_AGE) + months[mylook] - d_new(i, COL_ACCRT);
      double age_at_last_visit = d_new(i, COL_AGE) + (double)arma::max(visits) - d_new(i, COL_ACCRT);

      // event occurred prior to last visit and age is less than 36 at time of event
      if(d_new(i, COL_ACCRT) + d_new(i, COL_EVTT) <= (double)arma::max(visits) &&
         d_new(i, COL_AGE) + d_new(i, COL_EVTT) <= (double)cfg["max_age_fu_months"]){
        cen = 0;
        obst = d_new(i, COL_EVTT);
        //DBG(Rcpp::Rcout, "i " << i << " event     : " << obst);
        cens = Rcpp::List::create(Rcpp::Named("cen") = cen, Rcpp::Named("obst") = obst);
        return cens;
      }
      // event occurred prior to last visit and age at time of event was more than 36
      if(d_new(i, COL_ACCRT) + d_new(i, COL_EVTT) <= (double)arma::max(visits) &&
         d_new(i, COL_AGE) + d_new(i, COL_EVTT) > (double)cfg["max_age_fu_months"]){
        cen = 1;
        obst = (double)cfg["max_age_fu_months"];
        //DBG(Rcpp::Rcout, "i " << i << " cens 1     : " << obst);
        cens = Rcpp::List::create(Rcpp::Named("cen") = cen, Rcpp::Named("obst") = obst);
        return cens;
      }
      // event occurred after last visit  -- age is only relevant in the censoring time
      if(d_new(i, COL_ACCRT) + d_new(i, COL_EVTT) > (double)arma::max(visits)){
        cen = 1;

        if(d_new(i, COL_AGE) + months[mylook] - d_new(i, COL_ACCRT) > (double)cfg["max_age_fu_months"]){
          obst = (double)cfg["max_age_fu_months"];
          //DBG(Rcpp::Rcout, "i " << i << " cens 2a     : " << obst);
        } else {
          obst = months[mylook] - d_new(i, COL_ACCRT);
          obst = obst < 0 ? 0 : obst;
          //DBG(Rcpp::Rcout, "i " << i << " cens 2b     : " << obst);
        }

        cens = Rcpp::List::create(Rcpp::Named("cen") = cen, Rcpp::Named("obst") = obst);
        return cens;
      }

    } else { // accrual is before this look but is after last visit

      // but accrual happens before this analysis because we have one
      cen = 1;

      if(d_new(i, COL_AGE) + months[mylook] - d_new(i, COL_ACCRT) > (double)cfg["max_age_fu_months"]){
        obst = (double)cfg["max_age_fu_months"];
        //DBG(Rcpp::Rcout, "i " << i << " cens 3     : " << obst);
      } else {
        obst = months[mylook] - d_new(i, COL_ACCRT);
        //DBG(Rcpp::Rcout, "i " << i << " cens 4     : " << obst);
      }

      cens = Rcpp::List::create(Rcpp::Named("cen") = cen, Rcpp::Named("obst") = obst);
      return cens;

    }

  }

  cens = Rcpp::List::create(Rcpp::Named("cen") = cen, Rcpp::Named("obst") = obst);

  return cens;
}


// [[Rcpp::export]]
Rcpp::List rcpp_cens_final(const arma::mat& d_new,
                             const arma::vec& visits,
                             const int i,
                             const int look,
                             const Rcpp::List& cfg){
  Rcpp::List cens;

  Rcpp::NumericVector looks = cfg["looks"];
  Rcpp::NumericVector months = cfg["interimmnths"];

  int mylook = look - 1;
  int curmonth = months[mylook];

  double cen = NA_REAL;
  double obst = NA_REAL;

  double fudge = 0.0001;

  // should have definitely had some visits. if not there is something very wrong...

  // event occurred prior to last surveillance visit and age is less than 36 at time of event
  if(d_new(i, COL_EVTT) - d_new(i, COL_ACCRT) <= (double)arma::max(visits) &&
     d_new(i, COL_AGE) + d_new(i, COL_EVTT) <= (double)cfg["max_age_fu_months"]){

    cen = 0;
    obst = d_new(i, COL_EVTT);
    //DBG(Rcpp::Rcout, "i " << i << " event     : " << obst);
    cens = Rcpp::List::create(Rcpp::Named("cen") = cen, Rcpp::Named("obst") = obst);
    return cens;

  }
  // event occurred prior to last visit and age at time of event was more than 36
  if(d_new(i, COL_ACCRT) + d_new(i, COL_EVTT) <= (double)arma::max(visits) &&
     d_new(i, COL_AGE) + d_new(i, COL_EVTT) > (double)cfg["max_age_fu_months"]){

    cen = 1;
    obst = (double)cfg["max_age_fu_months"];
    //DBG(Rcpp::Rcout, "i " << i << " censor 1     : " << obst);
    cens = Rcpp::List::create(Rcpp::Named("cen") = cen, Rcpp::Named("obst") = obst);
    return cens;

  }
  // event occurred after the final follow up for each individual
  // -- age is only relevant in censoring time
  if(d_new(i, COL_ACCRT) + d_new(i, COL_EVTT) > (double)arma::max(visits)){
    cen = 1;

    if(d_new(i, COL_AGE) + (double)arma::max(visits) - d_new(i, COL_ACCRT) > (double)cfg["max_age_fu_months"]){
      obst = (double)cfg["max_age_fu_months"];
     // DBG(Rcpp::Rcout, "i " << i << " censor 2     : " << obst);
    } else {
      obst = d_new(i, COL_AGE) + (double)arma::max(visits) - d_new(i, COL_ACCRT);
      //DBG(Rcpp::Rcout, "i " << i << " censor 3    : " << obst);
    }

    cens = Rcpp::List::create(Rcpp::Named("cen") = cen, Rcpp::Named("obst") = obst);
    return cens;

  }

  return cens;
}


// [[Rcpp::export]]
arma::vec rcpp_visits(const arma::mat& d_new,
                      const int i,
                      const int look,
                      const Rcpp::List& cfg) {

  Rcpp::NumericVector looks = cfg["looks"];
  Rcpp::NumericVector months = cfg["interimmnths"];

  int mylook = look - 1;
  int nvisits = 0;

  double csum_visit_time = 0;
  double fu_36months = 0;

  // visits are the times (from the start of the trial) at which we review the
  // medical records data to deterime if a medical presentation has occurred.
  arma::vec visits = arma::vec();


  if(months[mylook] >= d_new(i, COL_ACCRT) + d_new(i, COL_FU1)  ||
     looks[mylook] == Rcpp::max(looks) ){
    nvisits++;
    visits.resize(nvisits);
    visits(nvisits - 1) = d_new(i, COL_ACCRT) + d_new(i, COL_FU1);
  }

  if(months[mylook] >= d_new(i, COL_ACCRT) + d_new(i, COL_FU2)  ||
     looks[mylook] == Rcpp::max(looks)){
    nvisits++;
    visits.resize(nvisits);
    visits(nvisits - 1) = d_new(i, COL_ACCRT) + d_new(i, COL_FU2);
  }

  csum_visit_time = R::runif((double)cfg["visit_lwr"], (double)cfg["visit_upr"]);

  while(d_new(i, COL_AGE) + csum_visit_time < (double)cfg["max_age_fu_months"] &&
        d_new(i, COL_ACCRT) + csum_visit_time < months[mylook]){

    nvisits++;
    visits.resize(nvisits);
    visits(nvisits - 1) = d_new(i, COL_ACCRT) + csum_visit_time;

    //DBG(Rcpp::Rcout, "i " << i << " nvisits                     " << nvisits );
    //DBG(Rcpp::Rcout, "i " << i << " COL_AGE                     " << d_new(i, COL_AGE)  );
    //DBG(Rcpp::Rcout, "i " << i << " COL_ACCRT                   " << d_new(i, COL_ACCRT)  );
    //DBG(Rcpp::Rcout, "i " << i << " csum_visit_time             " << csum_visit_time  );
    //DBG(Rcpp::Rcout, "i " << i << " age at visit                " << d_new(i, COL_AGE) + csum_visit_time );

    csum_visit_time += R::runif((float)cfg["visit_lwr"], (float)cfg["visit_upr"]);
  }

  //DBG(Rcpp::Rcout, "i " << i << " nvisits after loop          " << nvisits );
  // we are going to look at the records for an individual for the last time
  // when they are 36 months (+/- 4 weeks ~= 1 month).

  // if the last visit of those that have already been computed does not meet
  // the criteria for the last follow up then add another surveillance visit.
  if(visits.n_elem > 0 &&
     d_new(i, COL_AGE) + arma::max(visits) - d_new(i, COL_ACCRT) <
       (double)cfg["max_age_fu_months"]){

    // make a draw from 36 months +/- 4 weeks
    // if this is greater than the max visit then add otherwise the last follow up visit
    // is already generated
    fu_36months = R::runif((double)cfg["max_age_fu_months"] - 1,
                           (double)cfg["max_age_fu_months"] + 1);

    if(looks[mylook] == Rcpp::max(looks)){

      //DBG(Rcpp::Rcout, "i " << i << " adding fu_36months at max(looks) " << fu_36months );
      nvisits++;
      visits.resize(nvisits);
      visits(nvisits - 1) = d_new(i, COL_ACCRT) + fu_36months - d_new(i, COL_AGE);
    } else if(d_new(i, COL_ACCRT) + fu_36months - d_new(i, COL_AGE) <= months[mylook]){

      //DBG(Rcpp::Rcout, "i " << i << " adding fu_36months          " << fu_36months );
      nvisits++;
      visits.resize(nvisits);
      visits(nvisits - 1) = d_new(i, COL_ACCRT) + fu_36months - d_new(i, COL_AGE);
    }
  }

  // DBG(Rcpp::Rcout, "i " << i << " nvisits after fu36          " << nvisits );
  // DBG(Rcpp::Rcout, "i " << i << " COL_AGE                     " << d_new(i, COL_AGE)  );
  // DBG(Rcpp::Rcout, "i " << i << " COL_ACCRT                   " << d_new(i, COL_ACCRT)  );


  // DBG(Rcpp::Rcout, "i " << i << " has had " << visits.n_elem << " visits " );
  //
  // if(_DEBUG == 1){
  //   Rprintf("      time from start    time from accrual    age\n" );
  //
  //   for(int l = 0; l < visits.n_elem; l++){
  //     // DBG(Rcpp::Rcout, visits(l) << ",         "
  //     //                            << visits(l) - d_new(i, COL_ACCRT) <<  ",         "
  //     //                            << d_new(i, COL_AGE) + visits(l) - d_new(i, COL_ACCRT)  );
  //
  //     Rprintf("         %6.3f               %6.3f        %6.3f\n",
  //             visits(l),
  //             visits(l) - d_new(i, COL_ACCRT),
  //             d_new(i, COL_AGE) + visits(l) - d_new(i, COL_ACCRT));
  //
  //
  //   }
  // }


  // DBG(Rcpp::Rcout, std::endl );


  return visits;
}




// [[Rcpp::export]]
Rcpp::List rcpp_clin_set_obst(arma::mat& d, const Rcpp::List& cfg,
                              const int look){

  int mylook = look - 1;

  int n_uncen_0 = 0;
  int n_uncen_1 = 0;
  double tot_obst_0 = 0;
  double tot_obst_1 = 0;

  double fudge = 0.0001;

  Rcpp::NumericVector looks = cfg["looks"];
  Rcpp::NumericVector months = cfg["interimmnths"];
  Rcpp::List cens;
  Rcpp::List ret;

  arma::vec visits;

  for(int i = 0; i < d.n_rows; i++){

    if(d(i, COL_ACCRT) > months[mylook] + fudge && months[mylook] != max(months)){

      d(i, COL_CEN) = NA_REAL;
      d(i, COL_OBST) = NA_REAL;

      //DBG(Rcpp::Rcout, "i " << i << " accrual time " << d(i, COL_ACCRT)
      //                       << " occurs after current anlaysis month" << months[mylook]);
      continue;
    }

    // work out visits and censoring conditional on visit times
    visits = rcpp_visits(d, i, look, cfg);
    cens = rcpp_cens(d, visits, i, look, cfg);

    d(i, COL_CEN) = (double)cens["cen"];
    d(i, COL_OBST) = (double)cens["obst"];

    // this is an NA check in CPP
    // see https://stackoverflow.com/questions/570669/checking-if-a-double-or-float-is-nan-in-c
    if(d(i, COL_CEN) != d(i, COL_CEN) || d(i, COL_OBST) != d(i, COL_OBST)){
      d(i, COL_CEN) = NA_REAL;
      d(i, COL_OBST) = NA_REAL;
      continue;
    }

    // sufficient stats for comp posterior
    if(d(i, COL_CEN) == 0){
      if(d(i, COL_TRT) == 0) {
        n_uncen_0 += 1;
      } else {
        n_uncen_1 += 1;
      }
    }

    if(d(i, COL_TRT) == 0) {
      //DBG(Rcpp::Rcout, "i " << i << " tot_obst_0 " << tot_obst_0 << " adding " << d(i, COL_OBST));
      tot_obst_0 = tot_obst_0 + d(i, COL_OBST);
    } else {
      //DBG(Rcpp::Rcout, "i " << i << " tot_obst_1 " << tot_obst_1 << " adding " << d(i, COL_OBST));
      tot_obst_1 = tot_obst_1 + d(i, COL_OBST);
    }
  }

  ret = Rcpp::List::create(Rcpp::Named("n_uncen_0") = n_uncen_0,
                           Rcpp::Named("tot_obst_0") = tot_obst_0,
                           Rcpp::Named("n_uncen_1") = n_uncen_1,
                           Rcpp::Named("tot_obst_1") = tot_obst_1);

  return ret;
}









// [[Rcpp::export]]
void rcpp_clin_interim_post(arma::mat& m,
                                 const int n_uncen_0,
                                 const double tot_obst_0,
                                 const int n_uncen_1,
                                 const double tot_obst_1,
                                 const int post_draw,
                                 const Rcpp::List& cfg){

  double a = (double)cfg["prior_gamma_a"];
  double b = (double)cfg["prior_gamma_b"];

  for(int i = 0; i < post_draw; i++){

    // see VWO_SmartStats_technical_whitepaper.pdf page 22 formula 11.2
    //m(i, COL_LAMB0) = R::rgamma(a + n_uncen_0, b/(1 + b * tot_obst_0));
    //m(i, COL_LAMB1) = R::rgamma(a + n_uncen_1, b/(1 + b * tot_obst_1));

    m(i, COL_LAMB0) = R::rgamma(a + n_uncen_0, 1/(b + tot_obst_0));
    m(i, COL_LAMB1) = R::rgamma(a + n_uncen_1, 1/(b + tot_obst_1));


    m(i, COL_RATIO) = m(i, COL_LAMB0) / m(i, COL_LAMB1);
  }

  return;
}












// immunological endpoint


// [[Rcpp::export]]
Rcpp::List rcpp_immu(const arma::mat& d, const Rcpp::List& cfg, const int look){

  Rcpp::NumericVector looks = cfg["looks"];
  Rcpp::NumericVector months = cfg["interimmnths"];
  arma::mat m = arma::zeros((int)cfg["post_draw"] , 3);
  int mylook = look - 1;
  int nimpute1 = 0;
  int nimpute2 = 0;
  double ppos_n = 0;
  double ppos_max = 0;
  Rcpp::List lnsero;
  Rcpp::List pp1;
  Rcpp::List pp2;
  Rcpp::List ret;

  if(looks[mylook] <= (int)cfg["nmaxsero"]){

    DBG(Rcpp::Rcout, "immu interim analysis ");

    // how many records did we observe in total (assumes balance)
    int nobs = rcpp_n_obs(d, look, looks, months, (float)cfg["sero_info_delay"]);

    // how many successes in each arm?
    lnsero = rcpp_lnsero(d, nobs);

    // posterior at this interim
    arma::mat m = arma::zeros((int)cfg["post_draw"] , 3);
    rcpp_immu_interim_post(d, m, nobs, (int)cfg["post_draw"], lnsero);

    // therefore how many do we need to impute?
    nimpute1 = looks[mylook] - nobs;
    // predicted prob of success at interim
    pp1 = rcpp_immu_interim_ppos(d, m,
                                  nobs, nimpute1,
                                  (int)cfg["post_draw"],
                                  lnsero, cfg);

    // predicted prob of success at nmaxsero
    nimpute2 = (int)cfg["nmaxsero"] - nobs;
    pp2 = rcpp_immu_interim_ppos(d, m, nobs, nimpute2,
                                    (int)cfg["post_draw"],
                                            lnsero, cfg);

    if(_DEBUG == 1){
      ret = Rcpp::List::create(Rcpp::Named("nobs") = nobs,
                               Rcpp::Named("nimpute1") = nimpute1,
                               Rcpp::Named("nimpute2") = nimpute2,
                               Rcpp::Named("lnsero") = lnsero,
                               Rcpp::Named("posterior") = m,
                               Rcpp::Named("ppos_n") = (double)pp1["ppos"],
                               Rcpp::Named("ppos_max") = (double)pp2["ppos"]);
    } else {
      ret = Rcpp::List::create(Rcpp::Named("ppos_n") = (double)pp1["ppos"],
                               Rcpp::Named("ppos_max") = (double)pp2["ppos"]);
    }
  }

  return ret;
}


// [[Rcpp::export]]
int rcpp_n_obs(const arma::mat& d,
               const int look,
               const Rcpp::NumericVector looks,
               const Rcpp::NumericVector months,
               const double info_delay){

  // reset look to zero first element
  int mylook = look - 1;
  double obs_to_month = months[mylook] - info_delay;
  int nobs = 0;
  int flooraccrt = 0;
  float fudge = 0.0001;

  DBG(Rcpp::Rcout, "obs_to_month " << obs_to_month);
  DBG(Rcpp::Rcout, "info_delay " << info_delay);

  for(int i = 0; i < d.n_rows; i++){

    // have to fudge to work around inexact numeric representation :(
    flooraccrt = floor(d(i, COL_ACCRT) + info_delay);
    if(flooraccrt == months[mylook] && i%2 == 1){
      // we accrue ctl/trt pairs simultaneously.
      nobs = i + 1 ;
      DBG(Rcpp::Rcout, "(Equal to) ID " << d(i, COL_ID) << " ACCRT "
                                << d(i, COL_ACCRT) << " at i = "
                                << i << " nobs = " << nobs);
      break;
    }


    if(d(i, COL_ACCRT) + info_delay > months[mylook] + fudge  && i%2 == 0){
      // we accrue ctl/trt pairs simultaneously.
      nobs = i ;
      DBG(Rcpp::Rcout, "(Greater than) ID " << d(i, COL_ID) << " ACCRT "
                                               << d(i, COL_ACCRT) << " at i = "
                                               << i << " nobs = " << nobs);
      break;
    }

  }
  return nobs;
}


// [[Rcpp::export]]
Rcpp::List rcpp_lnsero(const arma::mat& d,
                       const int nobs){

  int n_sero_ctl = 0;
  int n_sero_trt = 0;

  // todo have now discovered this can probably be done with arma::find
  for(int i = 0; i < nobs; i++){

    if(d(i, COL_TRT) == 0){
      n_sero_ctl = n_sero_ctl + d(i, COL_SEROT3);
    } else {
      n_sero_trt = n_sero_trt + d(i, COL_SEROT3);
    }
  }

  //DBG(Rcpp::Rcout, "nobs " << nobs);
  //DBG(Rcpp::Rcout, "n_sero_ctl " << n_sero_ctl);
  //DBG(Rcpp::Rcout, "n_sero_trt " << n_sero_trt);

  Rcpp::List l = Rcpp::List::create(Rcpp::Named("n_sero_ctl") = n_sero_ctl,
                                    Rcpp::Named("n_sero_trt") = n_sero_trt);

  //Rcpp::List l = Rcpp::List::create(n_sero_ctl, n_sero_trt);

  //DBG(Rcpp::Rcout, "lnsero " << l);

  return l;
}


// [[Rcpp::export]]
void rcpp_immu_interim_post(const arma::mat& d,
                            arma::mat& m,
                            const int nobs,
                            const int post_draw,
                            const Rcpp::List& lnsero){

  for(int i = 0; i < post_draw; i++){
    m(i, COL_THETA0) = R::rbeta(1 + (int)lnsero["n_sero_ctl"], 1 + (nobs/2) - (int)lnsero["n_sero_ctl"]);
    m(i, COL_THETA1) = R::rbeta(1 + (int)lnsero["n_sero_trt"], 1 + (nobs/2) - (int)lnsero["n_sero_trt"]);
    m(i, COL_DELTA) = m(i, COL_THETA1) - m(i, COL_THETA0);
  }

  return;

}


// [[Rcpp::export]]
Rcpp::List rcpp_immu_interim_ppos(const arma::mat& d,
                                  const arma::mat& m,
                                  const int nobs,
                                  const int nimpute,
                                  const int post_draw,
                                  const Rcpp::List& lnsero,
                                  const Rcpp::List& cfg){

  int n_sero_ctl = 0;
  int n_sero_trt = 0;
  int win = 0;
  arma::vec t0 = arma::zeros(post_draw);
  arma::vec t1 = arma::zeros(post_draw);
  arma::vec delta1 = arma::zeros(post_draw);
  arma::vec postprobdelta_gt0 = arma::zeros(post_draw);
  int ntarget = nobs + nimpute;

  // create 1000 phony interims conditional on our current understanding
  // of theta0 and theta1.
  for(int i = 0; i < post_draw; i++){

    // This is a view of the total draws at a sample size of nobs + nimpute
    n_sero_ctl = lnsero["n_sero_ctl"] + R::rbinom((nimpute/2), m(i, COL_THETA0));
    n_sero_trt = lnsero["n_sero_trt"] + R::rbinom((nimpute/2), m(i, COL_THETA1));

    // update the posteriors
    for(int j = 0; j < post_draw; j++){

      t0(j) = R::rbeta(1 + n_sero_ctl, 1 + (ntarget/2) - n_sero_ctl);
      t1(j) = R::rbeta(1 + n_sero_trt, 1 + (ntarget/2) - n_sero_trt);

      delta1(j) = t1(j) - t0(j);
    }

    // empirical posterior probability that ratio_lamb > 1
    arma::uvec tmp = arma::find(delta1 > 0);
    postprobdelta_gt0(i) =  (double)tmp.n_elem / (double)post_draw;
    if(postprobdelta_gt0(i) > (double)cfg["post_sero_thresh"]){
      win++;
    }
  }

  //DBG(Rcpp::Rcout, "pow " << std::pow((float)post_draw, 2.0) );
  //DBG(Rcpp::Rcout, "mean difference mean_delta " << mean_delta );
  //DBG(Rcpp::Rcout, "mean   postprobdelta_gt0 " << (double)arma::mean(postprobdelta_gt0) );
  //DBG(Rcpp::Rcout, "median postprobdelta_gt0 " << (double)arma::median(postprobdelta_gt0) );
  //DBG(Rcpp::Rcout, "min    postprobdelta_gt0 " << (double)arma::min(postprobdelta_gt0) );
  //DBG(Rcpp::Rcout, "max    postprobdelta_gt0 " << (double)arma::max(postprobdelta_gt0) );
  //DBG(Rcpp::Rcout, "win " << (double)win );
  //DBG(Rcpp::Rcout, "post_draw " << (double)post_draw);
  //DBG(Rcpp::Rcout, "mean wins " << (double)win / (double)post_draw );

  Rcpp::List res = Rcpp::List::create(Rcpp::Named("ppos") = (double)win / (double)post_draw);

  return res;

}


// [[Rcpp::export]]
void rcpp_outer(const arma::vec& z,
                const arma::vec& t,
                arma::mat& out){

  for(int i = 0; i < z.n_elem; i++){
    for(int j = 0; j < t.n_elem; j++){
      out(i, j) = z(i) >= t(j) ? 1 : 0;
    }
  }
}


// [[Rcpp::export]]
Rcpp::List rcpp_logrank(const arma::mat& d,
                        const int look,
                        const Rcpp::List& cfg){

  int mylook = look - 1;
  Rcpp::NumericVector looks = cfg["looks"];
  Rcpp::NumericVector months = cfg["interimmnths"];

  arma::mat d_new = d.submat(0, COL_ID, looks[mylook]-1, COL_OBST);
  //DBG(Rcpp::Rcout, "  submatrix " << std::endl << d_new );

  // all observation time indices by group
  arma::uvec idxz1 = arma::find(d_new.col(COL_TRT) == 0);
  arma::uvec idxz2 = arma::find(d_new.col(COL_TRT) == 1);

  //DBG(Rcpp::Rcout, " idxz1 " <<  idxz1 );
  //DBG(Rcpp::Rcout, " idxz2 " <<  idxz2 );

  arma::vec tmp = d_new.col(COL_OBST);
  //DBG(Rcpp::Rcout, " tmp " <<  tmp );

  arma::vec z1 = tmp.elem(idxz1);
  arma::vec z2 = tmp.elem(idxz2);

  // uncensored observation time indices by group
  idxz1 = arma::find(d_new.col(COL_TRT) == 0 && d_new.col(COL_CEN) == 0);
  idxz2 = arma::find(d_new.col(COL_TRT) == 1 && d_new.col(COL_CEN) == 0);

  arma::vec t1 = tmp.elem(idxz1);
  arma::vec t2 = tmp.elem(idxz2);

  //DBG(Rcpp::Rcout, " z1 " << z1);
  //DBG(Rcpp::Rcout, " z2 " << z2);
  //DBG(Rcpp::Rcout, " t1 " << t1);
  //DBG(Rcpp::Rcout, " t2 " << t2);

  arma::mat z1t1 = arma::zeros(z1.n_elem, t1.n_elem);
  arma::mat z1t2 = arma::zeros(z1.n_elem, t2.n_elem);
  arma::mat z2t1 = arma::zeros(z2.n_elem, t1.n_elem);
  arma::mat z2t2 = arma::zeros(z2.n_elem, t2.n_elem);

  //DBG(Rcpp::Rcout, " rcpp_outer " );
  rcpp_outer(z1, t1, z1t1);
  rcpp_outer(z1, t2, z1t2);
  rcpp_outer(z2, t1, z2t1);
  rcpp_outer(z2, t2, z2t2);

  //DBG(Rcpp::Rcout, " cumsum " );
  arma::vec risk1t1 = arma::sum(z1t1.t(), 1);
  arma::vec risk1t2 = arma::sum(z1t2.t(), 1);
  arma::vec risk2t1 = arma::sum(z2t1.t(), 1);
  arma::vec risk2t2 = arma::sum(z2t2.t(), 1);

  // DBG(Rcpp::Rcout, " risk1t1 " << risk1t1);
  // DBG(Rcpp::Rcout, " risk1t2 " << risk1t2);
  // DBG(Rcpp::Rcout, " risk2t1 " << risk2t1);
  // DBG(Rcpp::Rcout, " risk2t2 " << risk2t2);

  double sum1 = 0;
  double sum2 = 0;
  double var1 = 0;
  double var2 = 0;

  for(int i=0; i < risk1t1.n_elem; i++){
    sum1 += risk2t1(i) / (risk1t1(i) + risk2t1(i));
    var1 += risk1t1(i) * risk2t1(i) / std::pow((risk1t1(i) + risk2t1(i)), 2.0);
  }
  for(int i=0; i < risk1t2.n_elem; i++){
    sum2 += risk1t2(i) / (risk1t2(i) + risk2t2(i));
    var2 += risk1t2(i) * risk2t2(i) / std::pow((risk1t2(i) + risk2t2(i)), 2.0);
  }

  // DBG(Rcpp::Rcout, " sum1 " << sum1 );
  // DBG(Rcpp::Rcout, " sum2 " << sum2 );
  // DBG(Rcpp::Rcout, " var1 " << var1 );
  // DBG(Rcpp::Rcout, " var2 " << var2 );
  //
  // DBG(Rcpp::Rcout, " (sum1 - sum2) " << (sum1 - sum2) );
  // DBG(Rcpp::Rcout, " std::pow(var1 + var2, 0.5) " << std::pow(var1 + var2, 0.5) );

  double logrank = (sum1 - sum2)/std::pow(var1 + var2, 0.5);
  double pvalue = R::pchisq(std::pow(logrank, 2), 1, 0, 0);

  //DBG(Rcpp::Rcout, " z1 " << std::endl << z1 );
  //DBG(Rcpp::Rcout, " z2 " << std::endl << z2 );

  //DBG(Rcpp::Rcout, " col of submatrix " << std::endl << d_new.col(COL_TRT)  );

  // z1 = d2$obst[d2$trt == 0],
  // delta1 = 1-d2$cen[d2$trt == 0],

  // z2 = d2$obst[d2$trt == 1],
  // delta2 = 1-d2$cen[d2$trt == 1])

  Rcpp::List res = Rcpp::List::create(Rcpp::Named("logrank") = logrank,
                                      Rcpp::Named("pvalue") = pvalue );

  return res;

}








// [[Rcpp::export]]
void atest(){


  arma::vec v = arma::zeros(10);
  arma::uvec gt1;

  for(int i = 0; i < 10; i++){
    v(i) = R::runif(-1, 1);
  }

  gt1 = arma::find(v > 0);

  DBG(Rcpp::Rcout, "v " << std::endl << v);
  DBG(Rcpp::Rcout, "gt1 " << std::endl << gt1);
  DBG(Rcpp::Rcout, "" );
  DBG(Rcpp::Rcout, "gt1 len " << gt1.n_elem);
  DBG(Rcpp::Rcout, "gt1 /10 " << (double) gt1.n_elem / (double)10);
}



// [[Rcpp::export]]
arma::vec rcpp_gamma(const int n, const double a, const double b) {

  arma::vec v = arma::zeros(n);

  for(int i = 0; i < n; i++){
    v(i) = R::rgamma(a, b);
  }

  return v;
}

// [[Rcpp::export]]
void rcpp_test_1(arma::mat& d) {

  for(int i = 1; i < d.n_rows; i++){
    d(i, 0) = 1;
  }
  for(int i = 1; i < d.n_cols; i++){
    d(0, i) = 1;
  }
  rcpp_test_sub_1(d);
  return;
}

// [[Rcpp::export]]
void rcpp_test_sub_1(arma::mat& d) {

  for(int i = 2; i < d.n_rows; i++){
    d(i, 1) = 2;
  }
  for(int i = 2; i < d.n_cols; i++){
    d(1, i) = 2;
  }
  return;
}

// [[Rcpp::export]]
arma::mat rcpp_test_2(const arma::mat& d) {

  arma::mat d_new = arma::mat(d);
  for(int i = 1; i < d_new.n_rows; i++){
    d_new(i, 0) = 1;
  }
  for(int i = 1; i < d_new.n_cols; i++){
    d_new(0, i) = 1;
  }
  d_new = rcpp_test_sub_2(d_new);
  return d_new;
}

// [[Rcpp::export]]
arma::mat rcpp_test_sub_2(arma::mat& d) {

  arma::mat d_new = arma::mat(d);
  for(int i = 2; i < d_new.n_rows; i++){
    d_new(i, 1) = 2;
  }
  for(int i = 2; i < d_new.n_cols; i++){
    d_new(1, i) = 2;
  }
  return d_new;
}






// [[Rcpp::export]]
Rcpp::List rcpp_cens_interim_alt(const arma::mat& d_new,
                             const int i,
                             const int look,
                             const Rcpp::List& cfg){

  Rcpp::List cens;

  Rcpp::NumericVector looks = cfg["looks"];
  Rcpp::NumericVector months = cfg["interimmnths"];

  int mylook = look - 1;
  int curmonth = months[mylook];

  double cen = NA_REAL;
  double obst = NA_REAL;

  double fudge = 0.0001;


  if(d_new(i, COL_ACCRT) <= months[mylook] + fudge){

    double age_at_this_analysis = d_new(i, COL_AGE) + months[mylook] - d_new(i, COL_ACCRT);

    // event occurred prior to last visit and age is less than 36 at time of event
    if(d_new(i, COL_AGE) + d_new(i, COL_EVTT) <= (double)cfg["max_age_fu_months"]){
      cen = 0;
      obst = d_new(i, COL_EVTT);
      DBG(Rcpp::Rcout, "i " << i << " event     : " << obst);
      cens = Rcpp::List::create(Rcpp::Named("cen") = cen, Rcpp::Named("obst") = obst);
      return cens;

    }
    // event occurred prior to last visit and age at time of event was more than 36
    if(d_new(i, COL_AGE) + d_new(i, COL_EVTT) > (double)cfg["max_age_fu_months"]){
      cen = 1;
      obst = (double)cfg["max_age_fu_months"];
      DBG(Rcpp::Rcout, "i " << i << " cens 1     : " << obst);
      cens = Rcpp::List::create(Rcpp::Named("cen") = cen, Rcpp::Named("obst") = obst);
      return cens;
    }
  } else{

    if(d_new(i, COL_AGE) + months[mylook] - d_new(i, COL_ACCRT) > (double)cfg["max_age_fu_months"]){
      obst = (double)cfg["max_age_fu_months"];
      DBG(Rcpp::Rcout, "i " << i << " cens 2a     : " << obst);
    } else {
      obst = d_new(i, COL_AGE) + months[mylook] - d_new(i, COL_ACCRT);
      DBG(Rcpp::Rcout, "i " << i << " cens 2b     : " << obst);
    }

    cens = Rcpp::List::create(Rcpp::Named("cen") = cen, Rcpp::Named("obst") = obst);
    return cens;
  }

  cens = Rcpp::List::create(Rcpp::Named("cen") = cen, Rcpp::Named("obst") = obst);

  return cens;
}

// if the event for this participant occurs after the last surveillance visit but
// before this visit, we are not aware that it has happened

// if(cen == 0 && obst == 0 &&
//    visits.n_elem > 0 &&
//    d_new(i, COL_EVTT) + d_new(i, COL_ACCRT) > (double)arma::max(visits)) {
//
//   double age_at_last_visit = d_new(i, COL_AGE) + (double)arma::max(visits) - d_new(i, COL_ACCRT);
//   double age_now = d_new(i, COL_AGE) + months[mylook] - d_new(i, COL_ACCRT);
//   DBG(Rcpp::Rcout, "i " << i << " event after max visit age_at_last_visit : " << age_at_last_visit);
//   cen = 1;
//   if(age_at_last_visit > (double)cfg["max_age_fu_months"] &&
//      (double)arma::max(visits) < months[mylook] ){
//     obst =  age_at_last_visit;
//   } else if (age_now < (double)cfg["max_age_fu_months"]) {
//     obst =  age_now;
//   } else {
//     obst = (double)cfg["max_age_fu_months"];
//   }
//   obst = obst < 0 ? 0 : obst;
//
//   DBG(Rcpp::Rcout, "i " << i << " event after max visit      : " << (double)arma::max(visits) << " after event at " << obst);
//   DBG(Rcpp::Rcout, "i " << i << " event after max visit accrt: " << d_new(i, COL_ACCRT) );
//   DBG(Rcpp::Rcout, "i " << i << " event after max visit obst : " << obst);
// }
// try {
//   if ((int)cens["obst"] == 0) {
//     DBG(Rcpp::Rcout, "i " << i << " ERROR missing censor     : " << obst);
//     //DBG(Rcpp::Rcout, "i " << i << " ERROR age_at_last_visit  : " << age_at_last_visit);
//     //DBG(Rcpp::Rcout, "i " << i << " ERROR age_now            : " << age_now);
//     DBG(Rcpp::Rcout, "i " << i << " ERROR d_new(i, COL_EVTT)     : " << d_new(i, COL_EVTT));
//     DBG(Rcpp::Rcout, "i " << i << " ERROR d_new(i, COL_ACCRT)     : " << d_new(i, COL_ACCRT));
//     DBG(Rcpp::Rcout, "i " << i << " ERROR d_new(i, COL_AGE)     : " << d_new(i, COL_AGE));
//
//     cen = NA_REAL;
//     obst =  NA_REAL;
//
//     throw std::range_error("missing censor rule");
//   }
//
//   cens = Rcpp::List::create(Rcpp::Named("cen") = cen,
//                             Rcpp::Named("obst") = obst);
//
//   return cens;
//
// } catch(std::exception &ex) {
//   forward_exception_to_r(ex);
// }



// // work out censoring to end of trial
// for(int j = 0; j < d.n_rows; j++){
//   //for(int i = 0; i < 200; i++){
//
//   // DBG(Rcpp::Rcout, "j " << j << " ppos starting analysis for month " << months[mylook] );
//
//   // work out visits and censoring conditional on visit times at study end
//   visits = rcpp_visits(d, j, looks.length(), cfg);
//   cens = rcpp_cens(d, visits, j, looks.length(), cfg);
//
//   d(j, COL_CEN) = (double)cens["cen"];
//   d(j, COL_OBST) = (double)cens["obst"];
//
//   // this is an NA check in CPP
//   // see https://stackoverflow.com/questions/570669/checking-if-a-double-or-float-is-nan-in-c
//   if(d(j, COL_OBST) != d(j, COL_OBST)){
//     // DBG(Rcpp::Rcout, "j " << j << " ppos tot_obst_0 " << d(j, COL_OBST) << " skipping to next iteration.");
//     continue;
//   }
//
//   // sufficient stats for comp posterior
//   if(d(j, COL_CEN) == 0){
//     if(d(j, COL_TRT) == 0) {
//       n_uncen_0 += 1;
//     } else {
//       n_uncen_1 += 1;
//     }
//   }
//
//   if(d(j, COL_TRT) == 0) {
//     tot_obst_0 = tot_obst_0 + d(j, COL_OBST);
//   } else {
//     tot_obst_1 = tot_obst_1 + d(j, COL_OBST);
//   }
// }
//
// DBG(Rcpp::Rcout, "i " << i << " ppos n_uncen_0 " << n_uncen_0);
// DBG(Rcpp::Rcout, "i " << i << " ppos n_uncen_1 " << n_uncen_1);
// DBG(Rcpp::Rcout, "i " << i << " ppos tot_obst_0 " << tot_obst_0);
// DBG(Rcpp::Rcout, "i " << i << " ppos tot_obst_1 " << tot_obst_1);
//
//
// // update the posteriors
// double a = (double)cfg["prior_gamma_a"];
// double b = (double)cfg["prior_gamma_b"];
//
// for(int k = 0; k < post_draw; k++){
//   t0(k) = R::rgamma(a + n_uncen_0, b / (1 + b * tot_obst_0));
//   t1(k) = R::rgamma(a + n_uncen_1, b / (1 + b * tot_obst_1));
//   ratio_lamb(k) = (double)t0(k) / (double)t1(k);
// }
