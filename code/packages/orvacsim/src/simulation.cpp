
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



#define _DEBUG  1

#if _DEBUG
#define DBG( os, msg )                             \
(os) << "DBG: " << __FILE__ << "(" << __LINE__ << ") "\
     << msg << std::endl
#else
#define DBG( os, msg )
#endif

#define _INFO  1

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
Rcpp::List rcpp_immu(const arma::mat& d, const Rcpp::List& cfg, const int look);
int rcpp_n_obs(const arma::mat& d,
               const int look,
               const Rcpp::NumericVector looks,
               const Rcpp::NumericVector months,
               const double info_delay);
Rcpp::List rcpp_lnsero(const arma::mat& d,
                       const int nobs);
arma::mat rcpp_immu_interim_post(const arma::mat& d,
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
Rcpp::List rcpp_clin(arma::mat& d, const Rcpp::List& cfg, const int look);
Rcpp::List rcpp_cens(const arma::mat& d_new,
                     const arma::vec& visits,
                     const int i,
                     const int look,
                     const Rcpp::List& cfg);
arma::vec rcpp_visits(const arma::mat& d_new,
                      const int i,
                      const int look,
                      const Rcpp::List& cfg);
arma::vec rcpp_gamma(const int n, const double a, const double b);
Rcpp::List rcpp_clin_interim_ppos(arma::mat& d,
                                  const arma::mat& m,
                                  const int nimpute,
                                  const int look,
                                  const Rcpp::List& cfg);
void rcpp_test_1(arma::mat& d);
void rcpp_test_sub_1(arma::mat& d);
arma::mat rcpp_test_2(const arma::mat& d) ;
arma::mat rcpp_test_sub_2(arma::mat& d);
// end function prototypes






// [[Rcpp::export]]
arma::vec rcpp_gamma(const int n, const double a, const double b) {

  arma::vec v = arma::zeros(n);

  for(int i = 0; i < n; i++){
    v(i) = R::rgamma(a, b);
  }

  return v;
}






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
    // tte - the paramaterisation of rexp uses SCALE NOTE RATE!!!!!!!!!!!
    // tte - the paramaterisation of rexp uses SCALE NOTE RATE!!!!!!!!!!!
    d(i, COL_EVTT) = R::rexp(1/((double)cfg["b0tte"] + d(i, COL_TRT) * (double)cfg["b1tte"]));

    // fu 1 and 2 times from time of accrual
    // fu 1 is between 14 and 21 days from accrual
    // fu 2 is between 28 and 55 days from accrual
    d(i, COL_FU1) = R::runif((double)cfg["fu1_lwr"], (double)cfg["fu1_upr"]);
    d(i, COL_FU2) = R::runif((double)cfg["fu2_lwr"], (double)cfg["fu2_upr"]);

    d(i, COL_CEN) = 0;
    d(i, COL_OBST) = 0;
  }

  return d;
}








// [[Rcpp::export]]
Rcpp::List rcpp_clin(arma::mat& d, const Rcpp::List& cfg, const int look){


  int post_draw = (int)cfg["post_draw"];

  Rcpp::NumericVector looks = cfg["looks"];
  Rcpp::NumericVector months = cfg["interimmnths"];

  arma::mat m = arma::zeros(post_draw , 3);

  int mylook = look - 1;
  int nimpute1 = 0;
  int nimpute2 = 0;
  double ppos_n = 0;
  double ppos_max = 0;
  double fudge = 0.0001;
  Rcpp::List lnsero;
  Rcpp::List cens;
  Rcpp::List pp1;
  Rcpp::List pp2;
  Rcpp::List ret = Rcpp::List::create(0);
  arma::vec visits;

  int n_uncen_0 = 0;
  double tot_obst_0 = 0;
  int n_uncen_1 = 0;
  double tot_obst_1 = 0;

  int i = 0;

  //for(i = 0; i < d.n_rows; i++){
  for(int i = 0; i < 200; i++){

    DBG(Rcpp::Rcout, "i " << i << " starting analysis for month " << months[mylook] );

    if(d(i, COL_ACCRT) > months[mylook] + fudge){
      d(i, COL_CEN) = NA_REAL;
      d(i, COL_OBST) = NA_REAL;

      DBG(Rcpp::Rcout, "i " << i << " accrual time " << d(i, COL_ACCRT)
                            << " occurs after current anlaysis month" << months[mylook]);
      continue;
    }

    // work out visits and censoring conditional on visit times
    visits = rcpp_visits(d, i, look, cfg);
    cens = rcpp_cens(d, visits, i, look, cfg);

    d(i, COL_CEN) = (double)cens["cen"];
    d(i, COL_OBST) = (double)cens["obst"];

    // this is an NA check in CPP
    // see https://stackoverflow.com/questions/570669/checking-if-a-double-or-float-is-nan-in-c
    if(d(i, COL_OBST) != d(i, COL_OBST)){
      DBG(Rcpp::Rcout, "i " << i << " tot_obst_0 " << d(i, COL_OBST) << " skipping to next iteration.");
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
      DBG(Rcpp::Rcout, "i " << i << " tot_obst_0 " << tot_obst_0 << " adding " << d(i, COL_OBST));
      tot_obst_0 = tot_obst_0 + d(i, COL_OBST);
    } else {
      DBG(Rcpp::Rcout, "i " << i << " tot_obst_1 " << tot_obst_1 << " adding " << d(i, COL_OBST));
      tot_obst_1 = tot_obst_1 + d(i, COL_OBST);
    }

  }


  // compute posterior
  // note the parameterisation of the gamma distribution is not the same as R
  // if the conjugate prior is gamma(k, q) then the posterior is:
  //
  // ~ gamma(k + sum(uncensored obs), q / (1 + q * total obs time))
  //
  // se ibrahim bayesian surv analysis and
  // https://cdn2.hubspot.net/hubfs/310840/VWO_SmartStats_technical_whitepaper.pdf

  double a = (double)cfg["prior_gamma_a"];
  double b = (double)cfg["prior_gamma_b"];

  for(int j = 0; j < post_draw; j++){

    m(j, COL_LAMB0) = R::rgamma(a + n_uncen_0, b / (1 + b * tot_obst_0));
    m(j, COL_LAMB1) = R::rgamma(a + n_uncen_1, b / (1 + b * tot_obst_1));
    m(j, COL_RATIO) = m(j, COL_LAMB0) / m(j, COL_LAMB1);

    DBG(Rcpp::Rcout, "i " << i << " post " << m.row(j) );
  }


  // use posterior to do pp at final analysis
  // this involves simulating multiple datasets conditional on the
  // posterior and therefore requires multiple to the censoring funct
  arma::mat d_new = arma::mat(d);
  int nimpute = max(looks) - looks[mylook];
  DBG(Rcpp::Rcout, " need to impute " << nimpute << " to get to end of trial");

  //double ppos_max = rcpp_clin_interim_ppos(d_new, m, nimpute, cfg);

  ret = Rcpp::List::create(Rcpp::Named("d") = d,
                           Rcpp::Named("posterior") = m,
                           Rcpp::Named("n_uncen_0") = n_uncen_0,
                           Rcpp::Named("tot_obst_0") = tot_obst_0,
                           Rcpp::Named("n_uncen_1") = n_uncen_1,
                           Rcpp::Named("tot_obst_1") = tot_obst_1);

  return ret;
}









// [[Rcpp::export]]
void rcpp_dat_small(arma::mat& d,
                    const Rcpp::List& cfg,
                    const int look,
                    const double l0,
                    const double l1) {

  Rcpp::NumericVector looks = cfg["looks"];
  Rcpp::NumericVector months = cfg["interimmnths"];

  int idxStart = looks[look];
  int n = cfg["nstop"];

  double tpp = (double)cfg["months_per_person"];

  for(int i = idxStart; i < n; i++){

    DBG(Rcpp::Rcout, "i " << i << "d(i, COL_EVTT) was " << d(i, COL_EVTT) );

    // tte - the paramaterisation of rexp uses SCALE NOTE RATE!!!!!!!!!!!
    if(d(i, COL_TRT) == 0){
      d(i, COL_EVTT) = R::rexp(1/l0);
    } else {
      d(i, COL_EVTT) = R::rexp(1/l1);
    }

    DBG(Rcpp::Rcout, "i " << i << "d(i, COL_EVTT) now " << d(i, COL_EVTT) );

    // fu 1 and 2 times from time of accrual
    // fu 1 is between 14 and 21 days from accrual
    // fu 2 is between 28 and 55 days from accrual
    d(i, COL_FU1) = R::runif((double)cfg["fu1_lwr"], (double)cfg["fu1_upr"]);
    d(i, COL_FU2) = R::runif((double)cfg["fu2_lwr"], (double)cfg["fu2_upr"]);

    d(i, COL_CEN) = 0;
    d(i, COL_OBST) = 0;
  }

  return;
}


// [[Rcpp::export]]
Rcpp::List rcpp_clin_interim_ppos(arma::mat& d,
                                  const arma::mat& m,
                                  const int nimpute,
                                  const int look,
                                  const Rcpp::List& cfg){


  Rcpp::NumericVector looks = cfg["looks"];

  int win = 0;
  int mylook = look - 1;
  int ntarget = looks[mylook] + nimpute;
  int post_draw = cfg["post_draw"];

  arma::vec t0 = arma::zeros(post_draw);
  arma::vec t1 = arma::zeros(post_draw);
  arma::vec delta_gt0 = arma::zeros(post_draw);
  arma::vec postprobdelta_gt0 = arma::zeros(post_draw);
  arma::mat delta1 = arma::zeros(post_draw, post_draw);


  for(int i = 0; i < post_draw; i++){

    // updates d
    rcpp_dat_small(d, cfg, look, m[i, COL_LAMB0], m[i, COL_LAMB1]);


    // // update the posteriors
    // for(int j = 0; j < post_draw; j++){
    //
    //   t0(j) = R::rbeta(1 + n_sero_ctl, 1 + (ntarget/2) - n_sero_ctl);
    //   t1(j) = R::rbeta(1 + n_sero_trt, 1 + (ntarget/2) - n_sero_trt);
    //
    //   delta1(j, i) = t1(j) - t0(j);
    //   mean_delta = mean_delta + delta1(j, i);
    //
    //   if(delta1(j, i) > 0){
    //     delta_gt0(j) = 1;
    //   }
    // }
    //
    // // empirical posterior probability that delta > 0
    // postprobdelta_gt0(i) = arma::mean(delta_gt0);
    // if(postprobdelta_gt0(i) > (float)cfg["post_sero_thresh"]){
    //   win++;
    // }
    //
    // //reset to zeros
    // delta_gt0 = arma::zeros(post_draw);
  }

  //DBG(Rcpp::Rcout, "pow " << std::pow((float)post_draw, 2.0) );

  // mean_delta = mean_delta / std::pow((float)post_draw, 2.0);

  //DBG(Rcpp::Rcout, "mean difference mean_delta " << mean_delta );
  //DBG(Rcpp::Rcout, "mean   postprobdelta_gt0 " << (double)arma::mean(postprobdelta_gt0) );
  //DBG(Rcpp::Rcout, "median postprobdelta_gt0 " << (double)arma::median(postprobdelta_gt0) );
  //DBG(Rcpp::Rcout, "min    postprobdelta_gt0 " << (double)arma::min(postprobdelta_gt0) );
  //DBG(Rcpp::Rcout, "max    postprobdelta_gt0 " << (double)arma::max(postprobdelta_gt0) );
  //DBG(Rcpp::Rcout, "win " << (double)win );
  //DBG(Rcpp::Rcout, "post_draw " << (double)post_draw);
  //DBG(Rcpp::Rcout, "mean wins " << (double)win / (double)post_draw );

  // Rcpp::List res = Rcpp::List::create(Rcpp::Named("ppos") = (double)win / (double)post_draw,
  // Rcpp::Named("delta") = mean_delta);

  Rcpp::List res ;

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

  // if there are no visits we cannot have seen
  if(cen == 0 && obst == 0 && visits.n_elem == 0) {

    DBG(Rcpp::Rcout, "i " << i << " no visits yet " );
    DBG(Rcpp::Rcout, "i " << i << " months[mylook] " << months[mylook] );
    if(d_new(i, COL_ACCRT) <= months[mylook]){
      cen = 1;
      // ensure non-negative (arises because of limits of precision representation)
      obst =  months[mylook] - d_new(i, COL_ACCRT) > 0 ? months[mylook] - d_new(i, COL_ACCRT) : 0;
    } else {
      cen = NA_REAL;
      obst =  NA_REAL;
    }
  }

  // if the event for this participant occurs after the last surveillance visit (excl
  // final age = 36 month) then it is not possible for us to have seen the event yet
  if(cen == 0 && obst == 0 &&
     d_new(i, COL_EVTT) + d_new(i, COL_ACCRT) > (double)arma::max(visits)) {

    if(visits.n_elem == 0){
      cen = 1;
      obst =  months[mylook] - d_new(i, COL_ACCRT);
    } else {
      cen = 1;
      obst = (double)arma::max(visits) - d_new(i, COL_ACCRT);
    }

    DBG(Rcpp::Rcout, "i " << i << " cen test 1 max obs: " << (double)arma::max(visits) );
    DBG(Rcpp::Rcout, "i " << i << " cen test 1 acc    : " << d_new(i, COL_ACCRT) );
    DBG(Rcpp::Rcout, "i " << i << " cen test 1 at     : " << obst);
  }

  // happy path
  if(cen == 0 && obst == 0 &&
     d_new(i, COL_EVTT) + d_new(i, COL_ACCRT) <= (double)arma::max(visits) &&
     d_new(i, COL_EVTT) + d_new(i, COL_AGE) <= (double)cfg["max_age_fu_months"]) {

    cen = 0;
    obst = d_new(i, COL_EVTT);

    DBG(Rcpp::Rcout, "i " << i << " happy path(1) max obs: " << (double)arma::max(visits) );
    DBG(Rcpp::Rcout, "i " << i << " happy path(1) acc    : " << d_new(i, COL_ACCRT) );
    DBG(Rcpp::Rcout, "i " << i << " happy path(1) at     : " << obst);

  } else if (d_new(i, COL_CEN) == 0 && d_new(i, COL_OBST) == 0 &&
    d_new(i, COL_EVTT) + d_new(i, COL_AGE) > (double)cfg["max_age_fu_months"]){

    // the event hasn't happened and the participant is now older than max fu age
    // or, more specifically, the last visit minus accrual which gives you the
    // max fu age observed
    if(visits.n_elem == 0){
      // ensure non-negative (arises because of limits of precision representation)
      // if they have enrolled with accrt exactly equal to the current look then do not
      // count them, i.e. set cen and obst to NA
      cen = months[mylook] - d_new(i, COL_ACCRT) > 0 ? 1 : NA_REAL;
      obst =  months[mylook] - d_new(i, COL_ACCRT) > 0 ? months[mylook] - d_new(i, COL_ACCRT) : NA_REAL;
    } else {
      cen = 1;
      obst = (double)arma::max(visits) - d_new(i, COL_ACCRT);
    }

    DBG(Rcpp::Rcout, "i " << i << " cen test 2     : ");

  }

  cens = Rcpp::List::create(Rcpp::Named("cen") = cen,
                            Rcpp::Named("obst") = obst);


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

  csum_visit_time = R::runif((float)cfg["visit_lwr"], (float)cfg["visit_upr"]);

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

  DBG(Rcpp::Rcout, "i " << i << " nvisits after loop          " << nvisits );

  if(visits.n_elem > 0 &&
     d_new(i, COL_AGE) + arma::max(visits) - d_new(i, COL_ACCRT) <
       (double)cfg["max_age_fu_months"] - 1){

    // make a draw from 36 months +/- 4 weeks
    // if this is greater than the max visit then add otherwise the last follow up visit
    // is already generated
    fu_36months = R::runif((double)cfg["max_age_fu_months"] - 1,
                           (double)cfg["max_age_fu_months"] + 1);


    if(d_new(i, COL_ACCRT) + fu_36months - d_new(i, COL_AGE) <= months[mylook] ||
       looks[mylook] == Rcpp::max(looks)){

      DBG(Rcpp::Rcout, "i " << i << " adding fu_36months          " << fu_36months );
      nvisits++;
      visits.resize(nvisits);
      visits(nvisits - 1) = d_new(i, COL_ACCRT) + fu_36months - d_new(i, COL_AGE);
    }
  }

  DBG(Rcpp::Rcout, "i " << i << " nvisits after fu36          " << nvisits );
  DBG(Rcpp::Rcout, "i " << i << " COL_AGE                     " << d_new(i, COL_AGE)  );
  DBG(Rcpp::Rcout, "i " << i << " COL_ACCRT                   " << d_new(i, COL_ACCRT)  );


  DBG(Rcpp::Rcout, "i " << i << " has had " << visits.n_elem << " visits " );

  if(_DEBUG == 1){
    Rprintf("      time from start    time from accrual    age\n" );

    for(int l = 0; l < visits.n_elem; l++){
      // DBG(Rcpp::Rcout, visits(l) << ",         "
      //                            << visits(l) - d_new(i, COL_ACCRT) <<  ",         "
      //                            << d_new(i, COL_AGE) + visits(l) - d_new(i, COL_ACCRT)  );

      Rprintf("         %6.3f               %6.3f        %6.3f\n",
              visits(l),
              visits(l) - d_new(i, COL_ACCRT),
              d_new(i, COL_AGE) + visits(l) - d_new(i, COL_ACCRT));


    }
  }


  DBG(Rcpp::Rcout, std::endl );


  return visits;
}






// [[Rcpp::export]]
arma::mat rcpp_clin_post(const arma::mat& d,
                         const int i,
                         const Rcpp::List& cfg,
                         const int look){


  arma::mat m = arma::zeros(1000 , 3);
  return m;


}








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
    m = rcpp_immu_interim_post(d, nobs, (int)cfg["post_draw"], lnsero);

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

    ret = Rcpp::List::create(Rcpp::Named("nobs") = nobs,
                             Rcpp::Named("nimpute1") = nimpute1,
                             Rcpp::Named("nimpute2") = nimpute2,
                             Rcpp::Named("lnsero") = lnsero,
                             Rcpp::Named("posterior") = m,
                             Rcpp::Named("ppn") = pp1,
                             Rcpp::Named("ppmax") = pp2);

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

  // the nobs - 1 is to adjust for the fact that we start at 0
  for(int i = 0; i < nobs - 1; i++){

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
arma::mat rcpp_immu_interim_post(const arma::mat& d,
                                 const int nobs,
                                 const int post_draw,
                                 const Rcpp::List& lnsero){

  arma::mat m = arma::zeros(post_draw, 3);

  for(int i = 0; i < post_draw; i++){
    m(i, COL_THETA0) = R::rbeta(1 + (int)lnsero["n_sero_ctl"], 1 + (nobs/2) - (int)lnsero["n_sero_ctl"]);
    m(i, COL_THETA1) = R::rbeta(1 + (int)lnsero["n_sero_trt"], 1 + (nobs/2) - (int)lnsero["n_sero_trt"]);
    m(i, COL_DELTA) = m(i, COL_THETA1) - m(i, COL_THETA0);
  }

  return m;

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
  arma::vec delta_gt0 = arma::zeros(post_draw);
  arma::vec postprobdelta_gt0 = arma::zeros(post_draw);
  arma::mat delta1 = arma::zeros(post_draw, post_draw);
  double mean_delta = 0;
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

      delta1(j, i) = t1(j) - t0(j);
      mean_delta = mean_delta + delta1(j, i);

      if(delta1(j, i) > 0){
        delta_gt0(j) = 1;
      }
    }

    // empirical posterior probability that delta > 0
    postprobdelta_gt0(i) = arma::mean(delta_gt0);
    if(postprobdelta_gt0(i) > (float)cfg["post_sero_thresh"]){
      win++;
    }

    //reset to zeros
    delta_gt0 = arma::zeros(post_draw);
  }

  //DBG(Rcpp::Rcout, "pow " << std::pow((float)post_draw, 2.0) );

  mean_delta = mean_delta / std::pow((float)post_draw, 2.0);

  //DBG(Rcpp::Rcout, "mean difference mean_delta " << mean_delta );
  //DBG(Rcpp::Rcout, "mean   postprobdelta_gt0 " << (double)arma::mean(postprobdelta_gt0) );
  //DBG(Rcpp::Rcout, "median postprobdelta_gt0 " << (double)arma::median(postprobdelta_gt0) );
  //DBG(Rcpp::Rcout, "min    postprobdelta_gt0 " << (double)arma::min(postprobdelta_gt0) );
  //DBG(Rcpp::Rcout, "max    postprobdelta_gt0 " << (double)arma::max(postprobdelta_gt0) );
  //DBG(Rcpp::Rcout, "win " << (double)win );
  //DBG(Rcpp::Rcout, "post_draw " << (double)post_draw);
  //DBG(Rcpp::Rcout, "mean wins " << (double)win / (double)post_draw );

  Rcpp::List res = Rcpp::List::create(Rcpp::Named("ppos") = (double)win / (double)post_draw,
                                      Rcpp::Named("delta") = mean_delta);

  return res;

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



