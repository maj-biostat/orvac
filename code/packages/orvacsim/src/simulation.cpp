
#include <RcppDist.h>
// [[Rcpp::depends(RcppDist)]]

#include <cmath>
#include <algorithm>

// ese Makevars
// compiler flags
// https://stackoverflow.com/questions/42328346/changing-the-left-most-optimization-flag-during-compilation-of-code-from-rcpp

// file.path(R.home("etc"), "Makeconf")
// [1] "/usr/lib64/R/etc/Makeconf"


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
#define COL_CEN           10
#define COL_OBST          11
#define COL_REASON        12
#define COL_IMPUTE        13
#define COL_REFTIME       14
#define NCOL              15

#define COL_THETA0        0
#define COL_THETA1        1
#define COL_DELTA         2

#define COL_LAMB0         0
#define COL_LAMB1         1
#define COL_RATIO         2



#define _DEBUG 0

#if _DEBUG
#define DBG( os, msg )                             \
(os) << "DBG: " << __FILE__ << "(" << __LINE__ << ") "\
     << msg << std::endl
#else
#define DBG( os, msg )
#endif

#define _INFO  1

#if _INFO
#define INFO( os, i, msg )                                \
   (os) << "INFO: " << __FILE__ << "(" << __LINE__ << ") "\
        << " sim = " << i << " " << msg << std::endl
#else
#define INFO( os, i, msg )
#endif

// function prototypes

arma::mat rcpp_dat(const Rcpp::List& cfg);
void rcpp_dat_small(const arma::mat& d,
                         const Rcpp::List& cfg,
                         const int look,
                         const double l0,
                         const double l1);

Rcpp::List rcpp_clin(arma::mat& d, const Rcpp::List& cfg,
                     const int look, const int idxsim);
Rcpp::List rcpp_clin_med(arma::mat& d, const Rcpp::List& cfg,
                     const int look, const int idxsim);
Rcpp::List rcpp_clin_set_state(arma::mat& d, const int look,
                              const double ref_time,
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
                             const int look,
                             const int nobs,
                             const int nimpute,
                             const int post_draw,
                             const Rcpp::List& lnsero,
                             const Rcpp::List& cfg);
Rcpp::List rcpp_immu_ppos_test(const arma::mat& d,
                                  const arma::mat& m,
                                  const int look,
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
Rcpp::List rcpp_dotrial(const int idxsim, const Rcpp::List& cfg,
                        const bool rtn_trial_dat);

// end function prototypes




class Trial {
private:
  int stop_ven_samp = 0;
  int stop_immu_fut = 0;
  int stop_clin_fut = 0;
  int stop_clin_sup = 0;
  int inconclu = 0;
  int nmaxsero = 200;
  int nstartclin = 200;
  int immu_ss = 0;
  int clin_ss = 0;

  bool i_final_win = 0;
  bool c_final_win = 0;

public:
  Trial(Rcpp::List cfg)
  {
    nmaxsero = cfg["nmaxsero"];
    nstartclin = cfg["nstartclin"];
  }
  Trial(Rcpp::List cfg, int vstop, int ifut, int cfut, int csup, int inc)
  {
    nmaxsero = cfg["nmaxsero"];
    nstartclin = cfg["nstartclin"];
    stop_ven_samp = vstop;
    stop_immu_fut = ifut;
    stop_clin_fut = cfut;
    stop_clin_sup = csup;
    inconclu = inc;
  }
  int maxsero();
  int startclin_at_n();
  bool do_immu(int n_current);
  bool do_clin(int n_current);
  int is_v_samp_stopped(){return stop_ven_samp;}
  int is_immu_fut(){return stop_immu_fut;}
  int is_clin_fut(){return stop_clin_fut;}
  int is_clin_sup(){return stop_clin_sup;}
  int is_inconclusive(){return inconclu;}
  int getnmaxsero(){return nmaxsero;}
  int getnstartclin(){return nstartclin;}
  int get_immu_ss(){return immu_ss;}
  int get_clin_ss(){return clin_ss;}
  int immu_final(){return i_final_win;}
  int clin_final(){return c_final_win;}
  void immu_stopv();
  void immu_fut();
  void clin_fut();
  void clin_sup();
  void inconclusive();
  void immu_set_ss(int n);
  void clin_set_ss(int n);
  void immu_final_win(bool won);
  void clin_final_win(bool won);
  void immu_state(const int idxsim);
  void clin_state(const int idxsim);
};
int Trial::maxsero(){
  return nmaxsero;
}
int Trial::startclin_at_n(){
  return nstartclin;
}
bool Trial::do_immu(int n_current){
  if(stop_ven_samp == 1){
    return false;
  }
  if(n_current > nmaxsero){
    return false;
  }
  if(stop_immu_fut == 1){
    return false;
  }
  if(stop_clin_fut == 1){
    return false;
  }
  if(stop_clin_sup == 1){
    return false;
  }
  return true;
}
bool Trial::do_clin(int n_current){

  if(n_current < nstartclin){
    return false;
  }
  if(stop_clin_fut == 1){
    return false;
  }
  if(stop_clin_sup == 1){
    return false;
  }
  if(stop_immu_fut == 1){
    return false;
  }
  return true;
}
void Trial::immu_stopv(){stop_ven_samp = 1;}
void Trial::immu_fut(){
  stop_immu_fut = 1;
  return;
}
void Trial::clin_fut(){stop_clin_fut = 1;}
void Trial::clin_sup(){stop_clin_sup = 1;}
void Trial::inconclusive(){inconclu = 1;}
void Trial::immu_set_ss(int n){immu_ss = n;}
void Trial::clin_set_ss(int n){clin_ss = n;}
void Trial::immu_final_win(bool won){
  i_final_win = won;
}
void Trial::clin_final_win(bool won){
  c_final_win = won;
}
void Trial::immu_state(const int idxsim){
  INFO(Rcpp::Rcout, idxsim,  "immu ep state: intrm stop v samp " << stop_ven_samp <<
    " fut " << stop_immu_fut << " interm inconclu " << inconclu << " fin analy win " << i_final_win );
}
void Trial::clin_state(const int idxsim){
  INFO(Rcpp::Rcout, idxsim, "clin ep state: intrm sup " << stop_clin_sup <<
    " fut " << stop_clin_fut << " interm inconclu " << inconclu << " fin analy win " << c_final_win );
}


// dotrial loop





// [[Rcpp::export]]
Rcpp::List rcpp_dotrial(const int idxsim,
                        const Rcpp::List& cfg,
                        const bool rtn_trial_dat){

  INFO(Rcpp::Rcout, idxsim, "STARTED.");

  Rcpp::NumericVector looks = cfg["looks"];
  Rcpp::NumericVector months = cfg["interimmnths"];
  Rcpp::NumericVector post_tte_sup_thresh = cfg["post_tte_sup_thresh"];

  // used in assessing futility along with pp_tte_fut_thresh
  Rcpp::NumericVector post_tte_win_thresh = cfg["post_tte_win_thresh"];
  Rcpp::NumericVector post_sero_win_thresh = cfg["post_sero_win_thresh"];
  Rcpp::List m_immu_res;
  Rcpp::List m_clin_res;
  double current_sup;

  arma::mat d = rcpp_dat(cfg);
  int nobs = 0;

  //Trial t(cfg, vstop, ifut, cfut, csup, inc);
  Trial t(cfg);
  int look;
  int i = 0;
  for(i = 0; i < looks.length(); i++){
    // look is here because all the original methods were called from R with r indexing
    look = i + 1;

    // we may not have started analysing the clin ep yet, but
    // we still need to set ss here otherwise it would just be recorded as 0 and
    // we would therefore underestimate the avg
    t.clin_set_ss(looks[i]);

    if(t.do_immu(looks[i])){
      nobs = rcpp_n_obs(d, look, looks, months, (double)cfg["sero_info_delay"]);
      INFO(Rcpp::Rcout, idxsim, "doing immu, with " << looks[i]
            << " enrld and " << nobs << " test results."
            << " sup thresh (stop v samp) " << (double)cfg["pp_sero_sup_thresh"]
            << ", pp win thresh " << (double)post_sero_win_thresh[i]
            << ", fut thresh " << (double)cfg["pp_sero_fut_thresh"]);

      m_immu_res = rcpp_immu(d, cfg, look);

      if((double)m_immu_res["ppos_max"] < (double)cfg["pp_sero_fut_thresh"]){
         INFO(Rcpp::Rcout, idxsim, "immu futile - stopping now, ppos_max " << (double)m_immu_res["ppos_max"] << " with " << nobs << " test results.");
         t.immu_fut();
         t.immu_set_ss(nobs);
         break;
      }

      if ((double)m_immu_res["ppos_n"] > (double)cfg["pp_sero_sup_thresh"] && !t.is_immu_fut()){
        nobs = rcpp_n_obs(d, look, looks, months, (double)cfg["sero_info_delay"]);
        INFO(Rcpp::Rcout, idxsim, "immu sup - stopping v samp now, ppos_n " << (double)m_immu_res["ppos_n"] << " with " << nobs << " test results.");
        t.immu_stopv();
      }
      t.immu_set_ss(nobs);
    }


    if(t.do_clin(looks[i])){
      INFO(Rcpp::Rcout, idxsim, "doing clin, with " << looks[i]
            << " enrld and sup thresh " << (double)post_tte_sup_thresh[i]
            << ", pp win thresh " << (double)post_tte_win_thresh[i]
            << ", fut thresh " << (double)cfg["pp_tte_fut_thresh"]);

      m_clin_res = rcpp_clin(d, cfg, look, idxsim);

      // INFO(Rcpp::Rcout, idxsim, "blah " << (double)m_clin_res["ratio"] << " "
      // << (double)m_clin_res["lwr"] << " "
      // << (double)m_clin_res["upr"] );

      if((double)m_clin_res["ppmax"] < (double)cfg["pp_tte_fut_thresh"]){
        INFO(Rcpp::Rcout, idxsim, "clin futile - stopping now, ppmax " << (double)m_clin_res["ppmax"] << " fut thresh " << (double)cfg["pp_tte_fut_thresh"]);
        t.clin_fut();
        break;
      }

      if ((double)m_clin_res["ppn"] > (double)post_tte_sup_thresh[i]  && !t.is_clin_fut()){
        INFO(Rcpp::Rcout, idxsim, "clin sup - stopping now, ppn " << (double)m_clin_res["ppn"] << " sup thresh " << (double)post_tte_sup_thresh[i] );
        t.clin_sup();
        break;
      }
    }


    // if at last look set inconclusive
    if(i == looks.length()-1){
      t.inconclusive();
    }


  }



  // final analysis for sero
  //how many successes in each arm?
  Rcpp::List lnsero = rcpp_lnsero(d, (int)cfg["nmaxsero"]);
  // posterior at this interim
  arma::mat m = arma::zeros((int)cfg["post_draw"] , 3);
  rcpp_immu_interim_post(d, m, (int)cfg["nmaxsero"], (int)cfg["post_draw"], lnsero);
  arma::uvec tmp = arma::find(m.col(COL_DELTA) > 0);
  double post_prob_gt0 =  (double)tmp.n_elem / (double)cfg["post_draw"];
  double i_mym = arma::mean(m.col(COL_DELTA));
  double i_mysd = arma::stddev(m.col(COL_DELTA));
  double i_lwr = i_mym - 1.96 * i_mysd;
  double i_upr = i_mym + 1.96 * i_mysd;
  i_mym = round(i_mym * 1000) / 1000;
  i_lwr = round(i_lwr * 1000) / 1000;
  i_upr = round(i_upr * 1000) / 1000;
  INFO(Rcpp::Rcout, idxsim, "FINAL: immu postr: p0 " << arma::mean(m.col(COL_THETA0))
  << "  p1 " << arma::mean(m.col(COL_THETA1))
  << "  delta " << i_mym << " (" << i_lwr << ", " << i_upr
  << "). n delta gt0 " << tmp.n_elem  <<  " prob_gt0 " << post_prob_gt0);
  if(post_prob_gt0 > (double)cfg["post_final_thresh"]){
    t.immu_final_win(true);
  } else{
    t.immu_final_win(false);
  }
  t.immu_state(idxsim);






  // final analysis for tte
  d.col(COL_CEN) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["nstop"]));
  d.col(COL_OBST) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["nstop"]));
  // updates d(COL_CEN) and d(COL_OBST)
  Rcpp::List lss = rcpp_clin_set_state(d, looks.size(), max(looks) + 36, cfg);

  double n_evnt_0b = (double)lss["n_evnt_0"];
  double n_evnt_1b = (double)lss["n_evnt_1"];
  double tot_obst_0 = (double)lss["tot_obst_0"];
  double tot_obst_1 = (double)lss["tot_obst_1"];

  double a = (double)cfg["prior_gamma_a"];
  double b = (double)cfg["prior_gamma_b"];

  m = arma::zeros((int)cfg["post_draw"] , 3);

  for(int i = 0; i < (int)cfg["post_draw"]; i++){
    // compute the posterior based on the __observed__ data to the time of the interim
    // take single draw
    m(i, COL_LAMB0) = R::rgamma(a + n_evnt_0b, 1/(b + tot_obst_0));
    m(i, COL_LAMB1) = R::rgamma(a + n_evnt_1b, 1/(b + tot_obst_1));
    m(i, COL_RATIO) = m(i, COL_LAMB0) / m(i, COL_LAMB1);
  }

  tmp = arma::find(m.col(COL_RATIO) > 1);
  double post_prob_gt1 =  (double)tmp.n_elem / (double)cfg["post_draw"];

  double c_mym = arma::mean(m.col(COL_RATIO));
  double c_mysd = arma::stddev(m.col(COL_RATIO));
  double c_lwr = c_mym - 1.96 * c_mysd;
  double c_upr = c_mym + 1.96 * c_mysd;
  c_mym = round(c_mym * 1000) / 1000;
  c_lwr = round(c_lwr * 1000) / 1000;
  c_upr = round(c_upr * 1000) / 1000;

  INFO(Rcpp::Rcout, idxsim, "FINAL: clin postr: l0 " << arma::mean(m.col(COL_LAMB0))
    << "  l1 " << arma::mean(m.col(COL_LAMB1))
    << "  ratio " << c_mym << " (" << c_lwr << ", " << c_upr
    << "). n ratio gt1 " << tmp.n_elem  <<  "  prob_gt1 " << post_prob_gt1);

  if(post_prob_gt1 > (double)cfg["post_final_thresh"]){
    t.clin_final_win(true);
  } else{
    t.clin_final_win(false);
  }
  t.clin_state(idxsim);


  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("idxsim") = idxsim);
  ret["p0"] = (double)cfg["baselineprobsero"];
  ret["p1"] = (double)cfg["trtprobsero"];
  ret["m0"] = log(2)/(double)cfg["b0tte"];
  ret["m1"] = log(2)/((double)cfg["b0tte"] + (double)cfg["b1tte"]);
  ret["look"] = i < looks.length() ? looks[i] : max(looks);
  ret["ss_immu"] = t.get_immu_ss();
  ret["ss_clin"] = t.get_clin_ss();
  ret["stop_v_samp"] = t.is_v_samp_stopped();
  ret["stop_i_fut"] = t.is_immu_fut();
  ret["stop_c_fut"] = t.is_clin_fut();
  ret["stop_c_sup"] = t.is_clin_sup();
  ret["inconclu"] = t.is_inconclusive();
  ret["i_final"] = t.immu_final();
  ret["c_final"] = t.clin_final();
  ret["i_ppn"] = m_immu_res.length() > 0 ? (double)m_immu_res["ppos_n"] : NA_REAL;
  ret["i_ppmax"] = m_immu_res.length() > 0 ? (double)m_immu_res["ppos_max"] : NA_REAL;
  ret["c_ppn"] = m_clin_res.length() > 0 ? (double)m_clin_res["ppn"] : NA_REAL;
  ret["c_ppmax"] = m_clin_res.length() > 0 ? (double)m_clin_res["ppmax"] : NA_REAL;
  ret["i_mean"] = (double)i_mym;
  ret["i_lwr"] = (double)i_lwr;
  ret["i_upr"] = (double)i_upr;
  ret["c_mean"] = (double)c_mym;
  ret["c_lwr"] = (double)c_lwr;
  ret["c_upr"] = (double)c_upr;

  INFO(Rcpp::Rcout, idxsim, "FINISHED.");

  // if(rtn_trial_dat){
  //   ret["d"] = d;
  // }

  // Rcpp::List ret = Rcpp::List::create(Rcpp::Named("idxsim") = idxsim);

  return ret;
}






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

    // d(i, COL_AGE) = r_truncnorm(cfg["age_months_mean"], cfg["age_months_sd"],
    //   cfg["age_months_lwr"], cfg["age_months_upr"]);

    d(i, COL_AGE) = R::runif((double)cfg["age_months_lwr"], (double)cfg["age_months_upr"]);

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
    d(i, COL_FU1) = 0.575; //R::runif((double)cfg["fu1_lwr"], (double)cfg["fu1_upr"]);
    d(i, COL_FU2) = 1.36345; // R::runif((double)cfg["fu2_lwr"], (double)cfg["fu2_upr"]);

  }

  d.col(COL_CEN) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["nstop"]));
  d.col(COL_OBST) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["nstop"]));
  d.col(COL_REASON) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["nstop"]));
  d.col(COL_IMPUTE) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["nstop"]));
  d.col(COL_REFTIME) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["nstop"]));

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

  //DBG(Rcpp::Rcout, "starting dat_small at " << idxStart << " going up to " << n << " l0 " << l0 << " l1 " << l1);

  for(int i = idxStart; i < n; i++){

    d(i, COL_AGE) = R::runif((double)cfg["age_months_lwr"], (double)cfg["age_months_upr"]);

    // tte - the paramaterisation of rexp uses SCALE NOTE RATE!!!!!!!!!!!
    if(d(i, COL_TRT) == 0){
      d(i, COL_EVTT) = R::rexp(1/l0)  ;
    } else {
      d(i, COL_EVTT) = R::rexp(1/l1)  ;
    }
    //DBG(Rcpp::Rcout, "i " << i << "d(i, COL_EVTT) now " << d(i, COL_EVTT) );

    // fu 1 and 2 times from time of accrual
    // fu 1 is between 14 and 21 days from accrual
    // 365.25/12 = a, 14/a
    // fu 2 is between 28 and 55 days from accrual
    d(i, COL_FU1) = 0.575; //R::runif((double)cfg["fu1_lwr"], (double)cfg["fu1_upr"]);
    d(i, COL_FU2) = 1.36345; // R::runif((double)cfg["fu2_lwr"], (double)cfg["fu2_upr"]);
  }

  d.col(COL_CEN) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["nstop"]));
  d.col(COL_OBST) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["nstop"]));
  d.col(COL_REASON) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["nstop"]));
  d.col(COL_IMPUTE) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["nstop"]));
  d.col(COL_REFTIME) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["nstop"]));

  return;
}






// clinical endpoint





// [[Rcpp::export]]
Rcpp::List rcpp_clin(arma::mat& d, const Rcpp::List& cfg,
                     const int look, const int idxsim) {

  int post_draw = (int)cfg["post_draw"];
  int mylook = look - 1;
  double fudge = 0.0001;
  double fu_36 = 36;

  double ppos = 0;
  double ppos_int = 0;
  double ppos_max = 0;

  double a = (double)cfg["prior_gamma_a"];
  double b = (double)cfg["prior_gamma_b"];

  Rcpp::NumericVector looks = cfg["looks"];
  Rcpp::NumericVector months = cfg["interimmnths"];

  arma::mat d_new ;
  arma::mat m = arma::zeros(post_draw , 3);
  arma::mat m_pp_int = arma::zeros(post_draw , 3);
  arma::mat m_pp_max = arma::zeros(post_draw , 3);

  arma::uvec uimpute;
  arma::uvec ugt1;
  arma::vec ppos_int_ratio_gt1 = arma::zeros(post_draw);
  arma::vec ppos_max_ratio_gt1 = arma::zeros(post_draw);

  // compute suff stats (calls visits and censoring) for the current interim
  d.col(COL_CEN) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["nstop"]));
  d.col(COL_OBST) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["nstop"]));

  Rcpp::List lss_post = rcpp_clin_set_state(d, look, months[mylook], cfg);
  int n_evnt_0 = (double)lss_post["n_evnt_0"];
  int n_evnt_1 = (double)lss_post["n_evnt_1"];
  double tot_obst_0 = (double)lss_post["tot_obst_0"];
  double tot_obst_1 = (double)lss_post["tot_obst_1"];

  Rcpp::List lss_int;
  Rcpp::List lss_max;

  // keep a copy of the original state
  d_new = arma::mat(d);
  arma::mat d_orig = arma::zeros((int)cfg["nstop"], 6);
  d_orig.col(0) = arma::vec(d.col(COL_EVTT));
  d_orig.col(1) = arma::vec(d.col(COL_CEN));
  d_orig.col(2) = arma::vec(d.col(COL_OBST));
  d_orig.col(3) = arma::vec(d.col(COL_REASON));
  d_orig.col(4) = arma::vec(d.col(COL_IMPUTE));
  d_orig.col(5) = arma::vec(d.col(COL_REFTIME));

  // subjs that require imputation
  uimpute = arma::find(d.col(COL_IMPUTE) == 1);

  arma::mat d1 = arma::mat(d);
  arma::mat d2 ;
  arma::mat d3 ;

  // for i in postdraws do posterior predictive trials
  // 1. for the interim (if we are at less than 50 per qtr)
  // 2. for the max sample size
  for(int i = 0; i < post_draw; i++){

    // compute the posterior based on the __observed__ data to the time of the interim
    // take single draw
    m(i, COL_LAMB0) = R::rgamma(a + n_evnt_0, 1/(b + tot_obst_0));
    m(i, COL_LAMB1) = R::rgamma(a + n_evnt_1, 1/(b + tot_obst_1));
    m(i, COL_RATIO) = m(i, COL_LAMB0) / m(i, COL_LAMB1);

    // use memoryless prop of exponential and impute enrolled kids that have not
    // yet had event.
    for(int j = 0; j < (int)uimpute.n_elem; j++){
      int sub_idx = uimpute(j);
      if(d(sub_idx, COL_TRT) == 0){
        // this just assigns a new evtt time on which we can update state.
        d(sub_idx, COL_EVTT) = d(sub_idx, COL_OBST) + R::rexp(1/m(i, COL_LAMB0))  ;
      } else {
        d(sub_idx, COL_EVTT) = d(sub_idx, COL_OBST) + R::rexp(1/m(i, COL_LAMB1))  ;
      }
    }

    // update view of the sufficent stats using enrolled
    // kids that have all now been given an event time
    lss_int = rcpp_clin_set_state(d, look, months[mylook] + fu_36, cfg);
    d2 = arma::mat(d);

    for(int j = 0; j < post_draw; j++){
      m_pp_int(j, COL_LAMB0) = R::rgamma(a + (double)lss_int["n_evnt_0"], 1/(b + (double)lss_int["tot_obst_0"]));
      m_pp_int(j, COL_LAMB1) = R::rgamma(a + (double)lss_int["n_evnt_1"], 1/(b + (double)lss_int["tot_obst_1"]));
      m_pp_int(j, COL_RATIO) = m_pp_int(j, COL_LAMB0) / m_pp_int(j, COL_LAMB1);
    }
    // empirical posterior probability that ratio_lamb > 1
    ugt1 = arma::find(m_pp_int.col(COL_RATIO) > 1);
    ppos_int_ratio_gt1(i) =  (double)ugt1.n_elem / (double)post_draw;



    // impute the remaining kids
    //DBG(Rcpp::Rcout, "from " << looks[mylook] << " to " << max(looks));

    for(int k = looks[mylook]; k < max(looks); k++){
      if(d(k, COL_TRT) == 0){
        d(k, COL_EVTT) = R::rexp(1/m(i, COL_LAMB0))  ;
      } else {
        d(k, COL_EVTT) = R::rexp(1/m(i, COL_LAMB1))  ;
      }
    }
    //DBG(Rcpp::Rcout, "d last " << d(max(looks)-1, COL_OBST) );

    // set the state up to the max sample size at time of the final analysis
    lss_max = rcpp_clin_set_state(d, looks.length(), max(months) + fu_36, cfg);
    d3 = arma::mat(d);

    // what does the posterior at max sample size say?
    for(int j = 0; j < post_draw; j++){
      m_pp_max(j, COL_LAMB0) = R::rgamma(a + (double)lss_max["n_evnt_0"], 1/(b + (double)lss_max["tot_obst_0"]));
      m_pp_max(j, COL_LAMB1) = R::rgamma(a + (double)lss_max["n_evnt_1"], 1/(b + (double)lss_max["tot_obst_1"]));
      m_pp_max(j, COL_RATIO) = m_pp_max(j, COL_LAMB0) / m_pp_max(j, COL_LAMB1);
    }
    // empirical posterior probability that ratio_lamb > 1
    ugt1 = arma::find(m_pp_max.col(COL_RATIO) > 1);
    ppos_max_ratio_gt1(i) =  (double)ugt1.n_elem / (double)post_draw;

    // reset to original state ready for the next posterior draw
    d.col(COL_EVTT) = arma::vec(d_orig.col(0));
    d.col(COL_CEN) = arma::vec(d_orig.col(1));
    d.col(COL_OBST) = arma::vec(d_orig.col(2));
    d.col(COL_REASON) = arma::vec(d_orig.col(3));
    d.col(COL_IMPUTE) = arma::vec(d_orig.col(4));
    d.col(COL_REFTIME) = arma::vec(d_orig.col(5));

  }

  double ppn = arma::mean(ppos_int_ratio_gt1);
  double ppmax = arma::mean(ppos_max_ratio_gt1);

  INFO(Rcpp::Rcout, idxsim, "clin: ppn " << ppn << " ppmax " << ppmax);

  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("ppn") = ppn,
                                Rcpp::Named("ppmax") = ppmax,
                                Rcpp::Named("lss_post") = lss_post,
                                Rcpp::Named("lss_int") = lss_int,
                                Rcpp::Named("lss_max") = lss_max,
                                Rcpp::Named("uimpute") = uimpute,
                                Rcpp::Named("m") = m,
                                Rcpp::Named("m_pp_int") = m_pp_int,
                                Rcpp::Named("m_pp_max") = m_pp_max,
                                Rcpp::Named("d1") = d1,
                                Rcpp::Named("d2") = d2,
                                Rcpp::Named("d3") = d3,
                                Rcpp::Named("ppos_int_ratio_gt1") = ppos_int_ratio_gt1,
                                Rcpp::Named("ppos_max_ratio_gt1") = ppos_max_ratio_gt1);

  // Rcpp::List ret = Rcpp::List::create(Rcpp::Named("ppmax") = 0);

  return ret;
}




// [[Rcpp::export]]
Rcpp::List rcpp_clin_med(arma::mat& d, const Rcpp::List& cfg,
                     const int look, const int idxsim){

  int post_draw = (int)cfg["post_draw"];
  int mylook = look - 1;
  double fudge = 0.0001;
  double fu_36 = 36;

  double ppos = 0;
  double ppos_int = 0;
  double ppos_max = 0;

  double a = (double)cfg["prior_gamma_a"];
  double b = (double)cfg["prior_gamma_b"];

  Rcpp::NumericVector looks = cfg["looks"];
  Rcpp::NumericVector months = cfg["interimmnths"];

  arma::mat d_new ;
  arma::mat m = arma::zeros(post_draw , 3);
  arma::mat m_pp_int = arma::zeros(post_draw , 3);
  arma::mat m_pp_max = arma::zeros(post_draw , 3);

  arma::uvec uimpute;
  arma::uvec ugt1;
  arma::vec ppos_int_ratio_gt1 = arma::zeros(post_draw);
  arma::vec ppos_max_ratio_gt1 = arma::zeros(post_draw);

  // compute suff stats (calls visits and censoring) for the current interim
  d.col(COL_CEN) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["nstop"]));
  d.col(COL_OBST) = arma::vec(Rcpp::rep(NA_REAL, (int)cfg["nstop"]));

  Rcpp::List lss_post = rcpp_clin_set_state(d, look, months[mylook], cfg);
  int n_evnt_0 = (double)lss_post["n_evnt_0"];
  int n_evnt_1 = (double)lss_post["n_evnt_1"];
  double tot_obst_0 = (double)lss_post["tot_obst_0"];
  double tot_obst_1 = (double)lss_post["tot_obst_1"];
  double log2 = (double)0.6931472;

  Rcpp::List lss_int;
  Rcpp::List lss_max;

  Rcpp::List reslist1(post_draw);
  Rcpp::List reslist2(post_draw);

  // keep a copy of the original state
  d_new = arma::mat(d);
  arma::mat d_orig = arma::zeros((int)cfg["nstop"], 6);
  d_orig.col(0) = arma::vec(d.col(COL_EVTT));
  d_orig.col(1) = arma::vec(d.col(COL_CEN));
  d_orig.col(2) = arma::vec(d.col(COL_OBST));
  d_orig.col(3) = arma::vec(d.col(COL_REASON));
  d_orig.col(4) = arma::vec(d.col(COL_IMPUTE));
  d_orig.col(5) = arma::vec(d.col(COL_REFTIME));

  // subjs that require imputation
  uimpute = arma::find(d.col(COL_IMPUTE) == 1);

  arma::mat d1 = arma::mat(d);
  arma::mat d2 ;
  arma::mat d3 ;

  // for i in postdraws do posterior predictive trials
  // 1. for the interim (if we are at less than 50 per qtr)
  // 2. for the max sample size
  for(int i = 0; i < post_draw; i++){

    // compute the posterior based on the __observed__ data to the time of the interim
    // take single draw
    m(i, COL_LAMB0) = 1/R::rgamma(a + n_evnt_0, 1/(b + log2*tot_obst_0));
    m(i, COL_LAMB1) = 1/R::rgamma(a + n_evnt_1, 1/(b + log2*tot_obst_1));
    m(i, COL_RATIO) = m(i, COL_LAMB1) - m(i, COL_LAMB0);

    // use memoryless prop of exponential and impute enrolled kids that have not
    // yet had event.
    for(int j = 0; j < (int)uimpute.n_elem; j++){
      int sub_idx = uimpute(j);
      if(d(sub_idx, COL_TRT) == 0){
        // this just assigns a new evtt time on which we can update state.
        d(sub_idx, COL_EVTT) = d(sub_idx, COL_OBST) + R::rexp(1/(log2/m(i, COL_LAMB0)))  ;
      } else {
        d(sub_idx, COL_EVTT) = d(sub_idx, COL_OBST) + R::rexp(1/(log2/m(i, COL_LAMB1)))  ;
      }
    }

    // update view of the sufficent stats using enrolled
    // kids that have all now been given an event time
    lss_int = rcpp_clin_set_state(d, look, months[mylook] + fu_36, cfg);
    d2 = arma::mat(d);

    for(int j = 0; j < post_draw; j++){
      m_pp_int(j, COL_LAMB0) = 1/R::rgamma(a + (double)lss_int["n_evnt_0"],
               1/(b + log2*(double)lss_int["tot_obst_0"]));
      m_pp_int(j, COL_LAMB1) = 1/R::rgamma(a + (double)lss_int["n_evnt_1"],
               1/(b + log2*(double)lss_int["tot_obst_1"]));
      m_pp_int(j, COL_RATIO) = m_pp_int(j, COL_LAMB1) - m_pp_int(j, COL_LAMB0);
    }
    ugt1 = arma::find(m_pp_int.col(COL_RATIO) > 0);
    ppos_int_ratio_gt1(i) =  (double)ugt1.n_elem / (double)post_draw;

    reslist1[i] = lss_int;
    reslist2[i] = m_pp_int;

    // impute the remaining kids
    //DBG(Rcpp::Rcout, "from " << looks[mylook] << " to " << max(looks));

    for(int k = looks[mylook]; k < max(looks); k++){
      if(d(k, COL_TRT) == 0){
        d(k, COL_EVTT) = R::rexp(1/(log2/m(i, COL_LAMB0)))  ;
      } else {
        d(k, COL_EVTT) = R::rexp(1/(log2/m(i, COL_LAMB1)))  ;
      }
    }
    //DBG(Rcpp::Rcout, "d last " << d(max(looks)-1, COL_OBST) );

    // set the state up to the max sample size at time of the final analysis
    lss_max = rcpp_clin_set_state(d, looks.length(), max(months) + fu_36, cfg);
    d3 = arma::mat(d);

    // what does the posterior at max sample size say?
    for(int j = 0; j < post_draw; j++){
      m_pp_max(j, COL_LAMB0) = 1/R::rgamma(a + (double)lss_max["n_evnt_0"],
               1/(b + log2*(double)lss_max["tot_obst_0"]));
      m_pp_max(j, COL_LAMB1) = 1/R::rgamma(a + (double)lss_max["n_evnt_1"],
               1/(b + log2*(double)lss_max["tot_obst_1"]));
      m_pp_max(j, COL_RATIO) = m_pp_max(j, COL_LAMB1) - m_pp_max(j, COL_LAMB0);
    }
    // empirical posterior probability that ratio_lamb > 0
    ugt1 = arma::find(m_pp_max.col(COL_RATIO) > 0);
    ppos_max_ratio_gt1(i) =  (double)ugt1.n_elem / (double)post_draw;

    // reset to original state ready for the next posterior draw
    d.col(COL_EVTT) = arma::vec(d_orig.col(0));
    d.col(COL_CEN) = arma::vec(d_orig.col(1));
    d.col(COL_OBST) = arma::vec(d_orig.col(2));
    d.col(COL_REASON) = arma::vec(d_orig.col(3));
    d.col(COL_IMPUTE) = arma::vec(d_orig.col(4));
    d.col(COL_REFTIME) = arma::vec(d_orig.col(5));

  }

  double ppn = arma::mean(ppos_int_ratio_gt1);
  double ppmax = arma::mean(ppos_max_ratio_gt1);

  INFO(Rcpp::Rcout, idxsim, "clin: ppn " << ppn << " ppmax " << ppmax);

  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("ppn") = ppn,
                                      Rcpp::Named("ppmax") = ppmax,
                                      Rcpp::Named("lss_post") = lss_post,
                                      Rcpp::Named("lss_int") = lss_int,
                                      Rcpp::Named("lss_max") = lss_max,
                                      Rcpp::Named("uimpute") = uimpute,
                                      Rcpp::Named("m") = m,
                                      Rcpp::Named("m_pp_int") = m_pp_int,
                                      Rcpp::Named("m_pp_max") = m_pp_max,
                                      Rcpp::Named("d1") = d1,
                                      Rcpp::Named("d2") = d2,
                                      Rcpp::Named("d3") = d3,
                                      Rcpp::Named("ppos_int_ratio_gt1") = ppos_int_ratio_gt1,
                                      Rcpp::Named("reslist1") = reslist1,
                                      Rcpp::Named("reslist2") = reslist2);

  // Rcpp::List ret = Rcpp::List::create(Rcpp::Named("ppmax") = 0);

  return ret;





}



// [[Rcpp::export]]
Rcpp::List rcpp_clin_set_state(arma::mat& d, const int look,
                              const double ref_time,
                              const Rcpp::List& cfg){

  // this updates the state of d in place
  // and provides sufficient stats.

  int mylook = look - 1;
  double fudge = 0.0001;

  int n_evnt_0 = 0;
  int n_evnt_1 = 0;
  double tot_obst_0 = 0;
  double tot_obst_1 = 0;

  double a = (double)cfg["prior_gamma_a"];
  double b = (double)cfg["prior_gamma_b"];

  Rcpp::NumericVector looks = cfg["looks"];
  Rcpp::NumericVector months = cfg["interimmnths"];
  Rcpp::List lsuffstat;

  arma::uvec utmp;

  // set censoring and event times up to current enrolled
  // these kids were all enrolled prior to the current look
  for(int sub_idx = 0; sub_idx < looks[mylook]; sub_idx++){

    if(d(sub_idx, COL_ACCRT) + d(sub_idx, COL_EVTT) <= ref_time &&
      d(sub_idx, COL_AGE) + d(sub_idx, COL_EVTT) <= (double)cfg["max_age_fu_months"]){
      // observed event
      // dont impute
      d(sub_idx, COL_CEN) = 0;
      d(sub_idx, COL_OBST) = d(sub_idx, COL_EVTT);
      d(sub_idx, COL_REASON) = 1;
      d(sub_idx, COL_IMPUTE) = 0;
      d(sub_idx, COL_REFTIME) = ref_time;

      if(d(sub_idx, COL_TRT) == 0) {
        n_evnt_0 += 1;
      } else {
        n_evnt_1 += 1;
      }

    } else if (d(sub_idx, COL_ACCRT) + d(sub_idx, COL_EVTT) <= ref_time &&
      d(sub_idx, COL_AGE) + d(sub_idx, COL_EVTT) > (double)cfg["max_age_fu_months"]){
      // censor at max age
      // dont impute
      d(sub_idx, COL_CEN) = 1;
      d(sub_idx, COL_OBST) = (double)cfg["max_age_fu_months"] - d(sub_idx, COL_AGE);
      d(sub_idx, COL_REASON) = 2;
      d(sub_idx, COL_IMPUTE) = 0;
      d(sub_idx, COL_REFTIME) = ref_time;
    } else if (d(sub_idx, COL_ACCRT) + d(sub_idx, COL_EVTT) > ref_time &&
      ref_time - d(sub_idx, COL_ACCRT) <=
                  (double)cfg["max_age_fu_months"] - d(sub_idx, COL_AGE)){
      // censor at mnth - accrual (t2)
      // impute
      d(sub_idx, COL_CEN) = 1;
      d(sub_idx, COL_OBST) = ref_time - d(sub_idx, COL_ACCRT) ;
      d(sub_idx, COL_REASON) = 3;
      d(sub_idx, COL_IMPUTE) = 1;
      d(sub_idx, COL_REFTIME) = ref_time;
    } else { // mnth - accrual >= max age - age at accrual
      // censor at max age
      // dont impute
      d(sub_idx, COL_CEN) = 1;
      d(sub_idx, COL_OBST) = (double)cfg["max_age_fu_months"] - d(sub_idx, COL_AGE) ;
      d(sub_idx, COL_REASON) = 4;
      d(sub_idx, COL_IMPUTE) = 0;
      d(sub_idx, COL_REFTIME) = ref_time;
    }

    if(d(sub_idx, COL_TRT) == 0) {
      tot_obst_0 = tot_obst_0 + d(sub_idx, COL_OBST);
    } else {
      tot_obst_1 = tot_obst_1 + d(sub_idx, COL_OBST);
    }
  }

  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("n_evnt_0") = n_evnt_0,
                           Rcpp::Named("tot_obst_0") = tot_obst_0,
                           Rcpp::Named("n_evnt_1") = n_evnt_1,
                           Rcpp::Named("tot_obst_1") = tot_obst_1,
                           Rcpp::Named("ref_time") = ref_time);

  return ret;
}




// immunological endpoint



// [[Rcpp::export]]
Rcpp::List rcpp_immu(const arma::mat& d, const Rcpp::List& cfg, const int look){

  Rcpp::NumericVector looks_target = cfg["looks_target"];
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

    // how many records did we observe in total (assumes balance)
    int nobs = rcpp_n_obs(d, look, looks, months, (float)cfg["sero_info_delay"]);

    // how many successes in each arm?
    lnsero = rcpp_lnsero(d, nobs);

    // posterior at this interim
    arma::mat m = arma::zeros((int)cfg["post_draw"] , 3);
    rcpp_immu_interim_post(d, m, nobs, (int)cfg["post_draw"], lnsero);

    // therefore how many do we need to impute assuming that we
    // were enrolling at the 50 per interim rate?
    nimpute1 = looks_target[mylook] - nobs;

    // if nimpute > 0 then do the ppos calc
    double post1gt0 = 0;
    if(nimpute1 > 0){
      // predicted prob of success at interim
      pp1 = rcpp_immu_ppos_test(d, m, look, nobs, nimpute1,(int)cfg["post_draw"],lnsero, cfg);
    } else {
      // else compute the posterior prob that delta > 0
      arma::uvec tmp = arma::find(m.col(COL_DELTA) > 0);
      post1gt0 = (double)tmp.n_elem / (double)cfg["post_draw"];
    }

    // predicted prob of success at nmaxsero
    // if nimpute2 == 0 then we are at nmaxsero with no information delay so just report
    // the posterior prob that delta is gt 0 (post1gt0) which has already been computed above.
    nimpute2 = (int)cfg["nmaxsero"] - nobs;
    if(nimpute2 > 0){
      pp2 = rcpp_immu_ppos_test(d, m, look, nobs, nimpute2, (int)cfg["post_draw"],lnsero, cfg);
    }


    // assess posterior
    double mean_delta =  arma::mean(m.col(COL_DELTA));
    double sd_delta =  arma::stddev(m.col(COL_DELTA));
    double lwr = mean_delta - 1.96 * sd_delta;
    double upr = mean_delta + 1.96 * sd_delta;
    mean_delta = round(mean_delta * 1000) / 1000;
    lwr = round(lwr * 1000) / 1000;
    upr = round(upr * 1000) / 1000;

    ret = Rcpp::List::create(Rcpp::Named("ppos_n") = nimpute1 > 0 ? (double)pp1["ppos"] : post1gt0,
                             Rcpp::Named("ppos_max") = nimpute2 > 0 ? (double)pp2["ppos"] : post1gt0,
                             Rcpp::Named("nimpute1") = nimpute1,
                             Rcpp::Named("nimpute2") = nimpute2,
                             Rcpp::Named("delta") = mean_delta,
                             Rcpp::Named("lwr") = lwr,
                             Rcpp::Named("upr") = upr);

  }

  return ret;
}


// [[Rcpp::export]]
int rcpp_n_obs(const arma::mat& d,
               const int look,
               const Rcpp::NumericVector looks,
               const Rcpp::NumericVector months,
               const double info_delay){

  // set look to zero (first element of array)
  int mylook = look - 1;
  double obs_to_month = months[mylook] - info_delay;
  int nobs = 0;
  int flooraccrt = 0;
  float fudge = 0.0001;

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
      // DBG(Rcpp::Rcout, "(Greater than) ID " << d(i, COL_ID) << " ACCRT "
      //                                          << d(i, COL_ACCRT) << " at i = "
      //                                          << i << " nobs = " << nobs);
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
                                  const int look,
                                  const int nobs,
                                  const int nimpute,
                                  const int post_draw,
                                  const Rcpp::List& lnsero,
                                  const Rcpp::List& cfg){

  int mylook = look - 1;
  int n_sero_ctl = 0;
  int n_sero_trt = 0;
  int win = 0;
  arma::vec t0 = arma::zeros(post_draw);
  arma::vec t1 = arma::zeros(post_draw);
  arma::vec delta1 = arma::zeros(post_draw);
  arma::vec postprobdelta_gt0 = arma::zeros(post_draw);

  Rcpp::NumericVector post_sero_win_thresh = cfg["post_sero_win_thresh"];

  arma::vec n_gt0 = arma::zeros(post_draw);

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
    n_gt0(i) = (double)tmp.n_elem;
    postprobdelta_gt0(i) =  (double)tmp.n_elem / (double)post_draw;
    if(postprobdelta_gt0(i) > post_sero_win_thresh[mylook]){
      win++;
    }
  }

  double ppos = (double)win / (double)post_draw;

  DBG(Rcpp::Rcout, "immu pp impute " << nimpute << " num win " << win << " ppos " << ppos <<
    " post thresh for win " << post_sero_win_thresh[mylook] );


  Rcpp::List res = Rcpp::List::create(Rcpp::Named("ppos") = ppos,
                                      Rcpp::Named("postprobdelta_gt0") = postprobdelta_gt0);

  return res;

}



// [[Rcpp::export]]
Rcpp::List rcpp_immu_ppos_test(const arma::mat& d,
                                  const arma::mat& m,
                                  const int look,
                                  const int nobs,
                                  const int nimpute,
                                  const int post_draw,
                                  const Rcpp::List& lnsero,
                                  const Rcpp::List& cfg){

  int mylook = look - 1;
  int n_sero_ctl = 0;
  int n_sero_trt = 0;
  int win = 0;
  arma::vec t0 = arma::zeros(post_draw);
  arma::vec t1 = arma::zeros(post_draw);
  arma::vec delta1 = arma::zeros(post_draw);
  arma::vec postprobdelta_gt0 = arma::zeros(post_draw);

  Rcpp::NumericVector post_sero_win_thresh = cfg["post_sero_win_thresh"];

  arma::vec n_gt0 = arma::zeros(post_draw);

  int ntarget = nobs + nimpute;

  // create 1000 phony interims conditional on our current understanding
  // of theta0 and theta1.
  for(int i = 0; i < post_draw; i++){

    // This is a view of the total draws at a sample size of nobs + nimpute
    n_sero_ctl = lnsero["n_sero_ctl"] + R::rbinom((nimpute/2), m(i, COL_THETA0));
    n_sero_trt = lnsero["n_sero_trt"] + R::rbinom((nimpute/2), m(i, COL_THETA1));

    double a = n_sero_trt;
    double b = 1 + (ntarget/2) - n_sero_trt;

    double c = n_sero_ctl;
    double d = 1 + (ntarget/2) - n_sero_ctl;

    double m1 = a / (a + b);
    double v1 = a*b / (std::pow(a + b, 2.0) * (a + b + 1));

    double m2 = c / (c + d);
    double v2 = c*d / (std::pow(c + d, 2.0) * (c + d + 1));

    double z = (m1 - m2) / pow(v1 + v2, 0.5);

    postprobdelta_gt0(i) = R::pnorm(z, 0.0, 1.0, 1, 0);

    if(postprobdelta_gt0(i) > post_sero_win_thresh[mylook]){
      win++;
    }
  }

  double ppos = (double)win / (double)post_draw;

  DBG(Rcpp::Rcout, "immu pp impute " << nimpute << " num win " << win << " ppos " << ppos <<
    " post thresh for win " << post_sero_win_thresh[mylook] );

  Rcpp::List res = Rcpp::List::create(Rcpp::Named("ppos") = ppos,
                                      Rcpp::Named("postprobdelta_gt0") = postprobdelta_gt0);

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
arma::vec rcpp_gamma(const int n, const double a, const double b) {

  arma::vec v = arma::zeros(n);

  for(int i = 0; i < n; i++){
    v(i) = R::rgamma(a, b);
  }

  return v;
}







