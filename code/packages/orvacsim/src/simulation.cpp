
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
#define DBGVAR( os, msg )                             \
(os) << "DBG: " << __FILE__ << "(" << __LINE__ << ") "\
     << msg << std::endl
#else
#define DBGVAR( os, msg )
#endif




// [[Rcpp::export]]
arma::mat rcpp_gen_dat(const Rcpp::List cfg) {

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

    // tte
    d(i, COL_EVTT) = R::rexp((float)cfg["b0tte"] + d(i, COL_TRT) * (float)cfg["b1tte"]);

    // fu 1 and 2 times
    d(i, COL_FU1) = R::runif((float)cfg["fu1_lwr"], (float)cfg["fu1_upr"]);
    d(i, COL_FU2) = R::runif((float)cfg["fu2_lwr"], (float)cfg["fu2_upr"]);
  }

  return d;
}



// [[Rcpp::export]]
arma::mat rcpp_mod_immu(const arma::mat d, const Rcpp::List cfg, const int look){

  int post_draws = (int)cfg["post_draw"];
  arma::mat m = arma::zeros(post_draws, 3);

  Rcpp::NumericVector looks = cfg["looks"];
  Rcpp::NumericVector months = cfg["interimmnths"];

  int n_sero_ctl = 0;
  int n_sero_trt = 0;
  int nobs = 0;
  int n_per_grp = 0;

  // reset look to zero first element
  int mylook = look - 1;
  float obs_to_month = months[mylook] - (float)cfg["sero_info_delay"];


  DBGVAR(Rcpp::Rcout, "looks " << looks);

  if(looks[mylook] <= (int)cfg["nmaxsero"]){

    for(int i = 0; i < looks[mylook]; i++){

      if(d(i, COL_ACCRT) > obs_to_month){
        // no need to add 1 to this as we accrue
        // ctl/trt pairs simultaneously.
        nobs = i;
        n_per_grp = nobs/2;
        break;
      }

      if(d(i, COL_TRT) == 0){
        n_sero_ctl = n_sero_ctl + d(i, COL_SEROT3);
      } else {
        n_sero_trt = n_sero_trt + d(i, COL_SEROT3);
      }

    }

    for(int i = 0; i < post_draws; i++){

      m(i, COL_THETA0) = R::rbeta(1 + n_sero_ctl, 1 + n_per_grp - n_sero_ctl);
      m(i, COL_THETA1) = R::rbeta(1 + n_sero_trt, 1 + n_per_grp - n_sero_trt);
      m(i, COL_DELTA) = m(i, COL_THETA1) - m(i, COL_THETA0);
    }


  }


  return m;

}




// [[Rcpp::export]]
arma::mat rcpp_censoring(const arma::mat d,
                                   const int look,
                                   const int trtstatus,
                                   const int iend,
                                   const float curmonth,
                                   const float surveillancemonths) {

  arma::mat res = arma::zeros(d.n_rows, d.n_cols + 1);
  arma::vec currentage = arma::zeros(d.n_rows);
  arma::vec obsmnths = arma::zeros(d.n_rows);

  int nmaxvisit = std::floor(curmonth/surveillancemonths);

  int j = 0;
  for(int i = 0; i < iend; i++){

    Rcpp::Rcout << "test " << d(i, 0) << d(i, 1) << std::endl;

    if(d(i, COL_TRT) == trtstatus){


      //d["current_age"] = d["age"] - d["accrt"] + curmonth;

      currentage[j] = 3;



    }



  }


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

  res = arma::join_rows(d, currentage);


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
