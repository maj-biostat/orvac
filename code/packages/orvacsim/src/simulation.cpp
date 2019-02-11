
#include <RcppDist.h>
// [[Rcpp::depends(RcppDist)]]

#include <cmath>
#include <algorithm>

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

#define COL_CURAGE        10
#define COL_CENT          11
#define COL_OBST          12
#define COL_CEN           13

#define NCOL              14


#define COL_THETA0        0
#define COL_THETA1        1
#define COL_DELTA         2

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
arma::mat rcpp_clin(const arma::mat& d, const Rcpp::List& cfg, const int look);
arma::mat rcpp_censoring(const arma::mat& d_new,
                         const int look,
                         const Rcpp::List& cfg);
// end function prototypes

















// [[Rcpp::export]]
arma::mat rcpp_dat(const Rcpp::List& cfg) {

  int n = cfg["nstop"];
  arma::mat d = arma::zeros(n, NCOL);
  float mu = 0;
  double tpp = (float)cfg["months_per_person"];

  for(int i = 0; i < n; i++){

    d(i, COL_ID) = i+1;
    d(i, COL_TRT) = ((i-1)%2 == 0) ? 0 : 1;
    // simultaneous accrual of each next ctl/trt pair
    d(i, COL_ACCRT) = (i%2 == 0) ? ((i+1)*tpp)+tpp : (i+1)*tpp;

    d(i, COL_AGE) = r_truncnorm(cfg["age_months_mean"], cfg["age_months_sd"],
      cfg["age_months_lwr"], cfg["age_months_upr"]);

    d(i, COL_SEROT2) = R::rbinom(1, cfg["baselineprobsero"]);
    d(i, COL_SEROT3) = d(i, COL_SEROT2);
    d(i, COL_PROBT3) = d(i, COL_TRT) * (float)cfg["deltaserot3"];

    if(d(i, COL_SEROT2) == 0 && d(i, COL_TRT) == 1){
      d(i, COL_SEROT3) = R::rbinom(1, d(i, COL_PROBT3));
    }

    // tte - the paramaterisation of rexp uses SCALE NOTE RATE!!!!!!!!!!!
    // tte - the paramaterisation of rexp uses SCALE NOTE RATE!!!!!!!!!!!
    // tte - the paramaterisation of rexp uses SCALE NOTE RATE!!!!!!!!!!!
    d(i, COL_EVTT) = R::rexp(1/((float)cfg["b0tte"] + d(i, COL_TRT) * (float)cfg["b1tte"]));

    // fu 1 and 2 times
    d(i, COL_FU1) = R::runif((float)cfg["fu1_lwr"], (float)cfg["fu1_upr"]);
    d(i, COL_FU2) = R::runif((float)cfg["fu2_lwr"], (float)cfg["fu2_upr"]);
  }

  return d;
}


// [[Rcpp::export]]
arma::mat rcpp_clin(const arma::mat& d, const Rcpp::List& cfg, const int look){


  arma::mat d_new = arma::mat(d);

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
  Rcpp::List ret = Rcpp::List::create(0);
  arma::mat surveill_res = arma::zeros(d.n_rows, 3);

  int nmaxvisit = std::floor((float)months[mylook]/(float)cfg["surveillance_mnths"]);
  arma::vec obs_mnths = arma::zeros(nmaxvisit);

  if(looks[mylook] < Rcpp::max(looks)){

    //DBG(Rcpp::Rcout, "clin interim analysis ");

    for(int i = 0; i < d_new.n_rows; i++){

      int cen = 0;
      double tte = 0;

      // if the current interim is prior to the enrollment of the individual then their
      // current age is set to NA
      d_new(i, COL_CURAGE) = months[mylook] - d_new(i, COL_ACCRT) + d_new(i, COL_AGE);
      d_new(i, COL_CURAGE) = d_new(i, COL_CURAGE) < 0 ? NA_REAL : d_new(i, COL_CURAGE);

      // We will observe participant x at these surveillance visits.
      // Note that these are adjusted by accrual time to give times from the start of the trial.
      obs_mnths(0) = d_new(i, COL_ACCRT) + d_new(i, COL_FU1);
      obs_mnths(1) = d_new(i, COL_ACCRT) + d_new(i, COL_FU2);

      for(int j = 2; j < nmaxvisit; j++){

        obs_mnths(j) = d_new(i, COL_ACCRT) + ((j-1) * (double)cfg["surveillance_mnths"]);
      }

      DBG(Rcpp::Rcout, "obs_mnths " << std::endl << obs_mnths);

      if(months[mylook] < (double)arma::min(obs_mnths)){
        d_new(i, COL_CEN) = 1;
        d_new(i, COL_OBST) = (double)arma::min(obs_mnths) - d_new(i, COL_ACCRT);
        continue;
      }



      // NEXT - FILL IN THE BLANKS FOR THE REMAINDER OF THIS CENSORING LOOP.


      // d_new(i, COL_CENT) = surveill_res(i, 1);
      // d_new(i, COL_OBST) = surveill_res(i, 2);
      // d_new(i, COL_CEN) = surveill_res(i, 3);

    }

    // ret = Rcpp::List::create(Rcpp::Named("nobs") = nobs,
    //                          Rcpp::Named("nimpute1") = nimpute1,
    //                          Rcpp::Named("nimpute2") = nimpute2,
    //                          Rcpp::Named("lnsero") = lnsero,
    //                          Rcpp::Named("posterior") = m,
    //                          Rcpp::Named("ppn") = pp1,
    //                          Rcpp::Named("ppmax") = pp2);

  }

  return d_new;
}


// [[Rcpp::export]]
arma::mat rcpp_censoring(const arma::mat& d_new,
                         const int look,
                         const Rcpp::List& cfg) {

  arma::mat res = arma::zeros(1, 3);
  arma::vec currentage = arma::zeros(d_new.n_rows);
  arma::vec obsmnths = arma::zeros(d_new.n_rows);

  Rcpp::NumericVector looks = cfg["looks"];
  Rcpp::NumericVector months = cfg["interimmnths"];

  int mylook = look - 1;
  int curmonth = months[mylook];
  // number of surveillance visits are per individual (excludes fu1 and fu2)
  // this is the maximum possile number of surveillance visits to date


  //DBG(Rcpp::Rcout, "curmonth " << curmonth);
  //DBG(Rcpp::Rcout, "nmaxvisit " << nmaxvisit);

  // int j = 0;
  // for(int i = 0; i < iend; i++){
  //
  //   Rcpp::Rcout << "test " << d_new(i, 0) << d_new(i, 1) << std::endl;
  //
  //   if(d_new(i, COL_TRT) == trtstatus){
  //
  //
  //     //d["current_age"] = d["age"] - d["accrt"] + curmonth;
  //
  //     currentage[j] = 3;
  //
  //
  //
  //   }
  //
  //
  //
  // }


  //d2 <-
  //
  //
  //NumericVector ageataccr = d["age"];
  //NumericVector accrt = d["accrt"];
  //NumericVector fu1 = d["fu1"];
  //NumericVector fu2 = d["fu2"];
  //NumericVector currentage = d["current_age"];
  //
  //
  //for(int i = istart-1; i < iend-1; i++){
  //
  //  //d["current_age"] = d["age"] - d["accrt"] + curmonth;
  //  currentage[i] = 3;
  //}


  //d["current_age"] = d["age"] - d["accrt"] + curmonth;

  //n_max_vis = std::floor(curmonth/surveillancemonths);

  //Rcpp::Rcout << "test" << d << std::endl;

  // allocate the output matrix
  //Rcpp::NumericMatrix output(x.nrow(), x.ncol());

  // SquareRoot functor (pass input and output matrixes)
  //SquareRoot squareRoot(x, output);

  // call parallelFor to do the work
  //parallelFor(0, x.length(), squareRoot);

  // return the output matrix

  //res = arma::join_rows(d, currentage);


  return res;
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














// mcens <- censoring(accrt = d[d$trt == trt_status,accrt][1:n_obs_grp],
//                    evtt = d[d$trt == trt_status,evtt][1:n_obs_grp],
//                    fu1 = d[d$trt == trt_status,fu1][1:n_obs_grp],
//                    fu2 = d[d$trt == trt_status,fu2][1:n_obs_grp],
//                    age = d$age_months[d$trt == trt_status][1:n_obs_grp],
//                    look = look)
//
//
//
//
// mcens <- censoring(accrt = d$accrt[d$trt == trt_status][1:(n_obs_grp + n_impute)],
//                    evtt = evtt_rep[[x]],
//                    fu1 = d[d$trt == trt_status,fu1][1:(n_obs_grp + n_impute)],
//                    fu2 = d[d$trt == trt_status,fu2][1:(n_obs_grp + n_impute)],
//                    age = d$age_months[d$trt == trt_status][1:(n_obs_grp + n_impute)],
//                    look = length(cfg$looks))

