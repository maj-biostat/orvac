

# To test, type test_file('util.R')


compute_sero_delta <- function(p_ctl, p_trt){
  
  delta <- (p_trt - p_ctl)/ (1 - p_ctl)
  delta
  
}


sim_cfg <- function(cfgfile = "cfg1.yaml", opt = NULL){
  
  tt <- tryCatch(configtmp <- read.config(file = cfgfile),
                 error=function(e) e, 
                 warning=function(w) w)
  ifelse(is(tt,"warning") | is(tt,"error"),"Configuration Warning/Error. 
         Please ensure configuration file has terminating empty line.",
         "Configuration File Loaded OK")
  
  flog.info("Using configuration: %s.", cfgfile)

  l <- list()
  
  thres <- unlist(tt$logging[[3]][1])
  tl <- list("WARN" = 4,
             "TRACE" = 9,
             "DEBUG" = 8,
             "INFO" = 6,
             "FATAL" = 1)
  whichappender <- unlist(tt$logging[[1]][1])
  flog.info("  logging to: %s appender.", whichappender)
  logfile <- unlist(tt$logging[[2]][1])
  flog.info("  logfile   : %s", logfile)
  appender <- NULL
  if(whichappender == "file"){
    flog.appender(appender.file(file.path(getwd(), "logs", logfile)), name='ROOT')
    
  }
  if(whichappender == "console"){
    flog.appender(appender.console(), name='ROOT')
  }
  flog.threshold(tt$logging[[3]][1])



  l$flog_appender <- whichappender
  l$flog_logfile <- logfile
  
  
  # not used
  l$desc <- tt$desc
  l$outfile <- tt$outfile
  
  # for return fields
  l$field_names <- c("idxsim",
                   "look",
                   "n_obs",
                   "ss_immu",
                   "ss_clin",
                   "n_max",
                   "n_max_sero",
                   
                   "i_ppos_n", "i_ppos_max", "i_post_prop_ctl", "i_post_prop_trt",  # immunological
                   "i_delta_mean", "i_delta_lwr_95", "i_delta_upr_95",          # immunological
                   
                   "c_ppos_n","c_ppos_max","c_post_lambda_ctl","c_post_lambda_trt", # clinical
                   "c_lambda_mean", "c_lambda_lwr_95", "c_lambda_upr_95",       # clinical
                   
                   "stop_v_samp",
                   "stop_i_fut",
                   "stop_c_fut",
                   "stop_c_sup")
  
  
  # "stop_i_sup",
  
  l$immu_rtn_names <- c("ppos_n", "ppos_max", "mean_post_prop_ctl", "mean_post_prop_trt",
                        "delta", "delta_lwr_95", "delta_upr_95")
  l$clin_rtn_names <- c("ppos_n", "ppos_max", "mean_post_lambda_ctl", "mean_post_lambda_trt", 
                        "delta", "delta_lwr_95", "delta_upr_95")


  # Sim control variables
  l$idsim <- tt$idsim
  l$nsims <- tt$nsims 
  l$seed <- tt$seed
  
  # interims
  l$nstart <- tt$nstart
  l$nstop <- tt$nstop 
  l$interim_period <- tt$interim_period 
  l$people_per_interim_period <- tt$people_per_interim_period 
  l$nmaxsero <- tt$nmaxsero  
  l$nstartclin <- tt$nstartclin
  
  l$looks <- seq(from = l$nstart, to=l$nstop, by = l$people_per_interim_period)
  # ensure max is exactly 1000
  if(max(l$looks) < l$nstop){
    l$looks <- c(l$looks, l$nstop)
  }
  l$nlooks <- length(l$looks)
  
  # accrual - it takes months_per_person months to recruite one person
  l$months_per_person <- l$interim_period / l$people_per_interim_period
  
  l$months_to_nstart <- l$months_per_person * l$nstart
  
  l$interimmnths <- seq(from = l$months_to_nstart , 
                        to= length(l$looks) * l$interim_period, 
                        by = l$interim_period)
  if(length(l$interimmnths) < length(l$looks)){
    l$interimmnths <- c(l$interimmnths, max(l$interimmnths)+l$interim_period)
  }

  stopifnot(length(l$interimmnths) == length(l$looks))
  stopifnot(max(l$looks) == l$nstop)
  
  # data generation

  
  
  # not currently used
  l$trtallocprob <- tt$trtallocprob   
  l$remoteprob <- tt$remoteprob      
  
  # used to simulate age of participants - used in the tte modelling 
  l$age_months_lwr <- tt$age_months_lwr 
  l$age_months_upr <- tt$age_months_upr 
  l$age_months_mean <- tt$age_months_mean 
  l$age_months_sd <- tt$age_months_sd
  
  l$max_age_fu_months <- tt$max_age_fu_months
  



  # seroconversion control variables
  # note - need to build utility functions to compute these
  l$baselineprobsero <- tt$baselineprobsero 
  l$trtprobsero <- tt$trtprobsero
  
  
  l$deltaserot3 <- compute_sero_delta(l$baselineprobsero, l$trtprobsero)
    
    
    
  # l$deltaremoteserot3 <- tt$deltaremoteserot3 

  # for significance testing of a win in the pp section
  l$post_sero_thresh <- tt$post_sero_thresh 
  



  # time to event control variables

  
  # time to event control variables
  # exponential
  # rates obtained from formula for med surv time (which are in months)
  # control log(2)/25 
  # treatment log(2)/30 - log(2)/25 
  # log(2)/25  = 0.027773
  # log(2)/30  = 0.023105
  
  # note - interp of survreg is fuckt 
  # see https://www.ms.uky.edu/~mai/Rsurv.pdf (top of page 4)
  # coefs aka. mu are related to rate via
  # mu = log(1/rate) 
  # i.e. rate = 1/exp(mu)
  
  # ignore the following - it relates to lognormal
  # for lognormal median is actually exp(mu) where mu represents
  # the mean parameter to a normal distribution, X for which
  # the lognormal dist is Y ~ exp(X) because the surv func
  # S(t) = 1 - CDF((ln(t) - μ) / σ) when equal to 0.5 gives
  # t_med = exp(μ).
  # So:
  # exp(3.4) ~= 30 months median survival for the control group
  # i.e. the time when surv prob equals 0.5 is about 30 months
  # exp(3.65) ~= 38 months
  #

  l$ctl_med_tte <- tt$ctl_med_tte
  l$trt_med_tte <- tt$trt_med_tte
  
  l$b0tte <- log(2)/l$ctl_med_tte 
  l$b1tte <- (log(2)/l$trt_med_tte) - (log(2)/l$ctl_med_tte)
  
  # l$b0tte <- tt$b0tte 
  # l$b1tte <- tt$b1tte  # trt effect
  # l$b2tte <- tt$b2tte # remote
  # l$b3tte <- tt$b3tte # baseline serot2 status
  
  # Ideally would like to incorporate serot2 status
  l$btte <- c(l$b0tte, l$b1tte)
  
  l$post_tte_thresh <- tt$post_tte_thresh
  
  l$ftte <- tt$ftte
  l$ttemodfile <- tt$ttemodfile
  
  # thresholds for interim decisions
  l$rule1_sero_dens_lwr <- tt$rule1_sero_dens_lwr
  l$rule1_sero_ppos_sup_thresh <- tt$rule1_sero_ppos_sup_thresh
  l$rule1_sero_pp_fut_thresh <- tt$rule1_sero_pp_fut_thresh
  l$rule1_sero_ppos_fut_thresh <- tt$rule1_sero_ppos_fut_thresh
  l$rule1_tte_postthresh  <- tt$rule1_tte_postthresh
  l$rule1_tte_ppthresh <- tt$rule1_tte_ppthresh
  l$rule1_tte_pposthresh  <- tt$rule1_tte_pposthresh
  l$rule2_tte_postthresh  <- tt$rule2_tte_postthresh
  l$rule3_sero_postthresh <- tt$rule3_sero_postthresh
  l$rule3_tte_postthresh  <- tt$rule3_tte_postthresh
  l$rule3_sero_ppos_thresh <- tt$rule3_sero_ppos_thresh

  
  # bayesian model control parameters
  
  # number of posterior draws to make 
  # relevant to both posterior and post predictive sections
  l$post_draw <- tt$post_draw
  
  # no longer using mcmc
  # mcmc
  l$mcmcchains <- tt$mcmcchains
  l$mcmciter <- tt$mcmciter
  l$mcmcburnin <- tt$mcmcburnin
  l$mcmcthin <- tt$mcmcthin
  l$mcmcadapt <- tt$mcmcadapt
  
  l$mcmc_gchains <- tt$mcmc_gchains
  l$mcmc_giter <- tt$mcmc_giter
  l$mcmc_gburnin <- tt$mcmc_gburnin
  l$mcmc_gthin <- tt$mcmc_gthin
  l$mcmc_gadapt <- tt$mcmc_gadapt
  l$mcmc_giterkeep <- seq(from = l$mcmc_gburnin+1,
                          to = l$mcmc_giter,
                          by = l$mcmc_gthin)
  
  l$mcmc_nim_chains <- tt$mcmc_nim_chains
  l$mcmc_nim_iter <- tt$mcmc_nim_iter
  l$mcmc_nim_burnin <- tt$mcmc_nim_burnin
  l$mcmc_nim_thin <- tt$mcmc_nim_thin
  l$mcmciterfin <- c(l$mcmciter - l$mcmcburnin) / l$mcmcthin
  
  
  if(opt$use){
    
    cat("Updating config.\n")
    flog.info("Updating configuration values based on command line arguments: %s", paste0(opt, collapse = " "))

    if(!is.null(opt$logfile)){
      l$flog_logfile <- opt$logfile

      flog.info("Updated logfile: %s.", l$flog_logfile)
      
      if(whichappender == "file"){
        cat(paste0("New logfile ", file.path(getwd(), "logs", l$flog_logfile), "\n"))
        flog.appender(appender.file(file.path(getwd(), "logs", l$flog_logfile)), name='ROOT')
      }
      flog.info("Updating configuration values based on command line arguments: %s", paste0(opt, collapse = " "))
    }
    
    if(!is.null(opt$idsim)){
      l$idsim <- opt$idsim
      flog.info("Updated idsim: %s.", l$idsim)
    }
    
    if(!is.null(opt$nsims)){
      l$nsims <- opt$nsims
      flog.info("Updated nsims: %s.", l$nsims)
    }
    
    if(!is.null(opt$seed)){
      l$seed <- opt$seed
      flog.info("Updated seed: %s.", l$seed)
    }

    if(!is.null(opt$accrual)){
      l$people_per_interim_period <- opt$accrual
      flog.info("Updated people_per_interim_period: %s.", l$people_per_interim_period)
      
      l$looks <- seq(from = l$nstart, to=l$nstop, by = l$people_per_interim_period)
      # ensure max is exactly 1000
      if(max(l$looks) < l$nstop){
        l$looks <- c(l$looks, l$nstop)
      }
      flog.info("Updated looks: %s.", paste0(l$looks, collapse = ", "))
      l$nlooks <- length(l$looks)
      flog.info("Updated nlooks: %s.", l$nlooks)
      
      # accrual - it takes months_per_person months to recruite one person
      l$months_per_person <- l$interim_period / l$people_per_interim_period
      flog.info("Updated months_per_person: %s.", l$months_per_person)
      l$months_to_nstart <- l$months_per_person * l$nstart
      flog.info("Updated months_to_nstart: %s.", l$months_to_nstart)
      
      l$interimmnths <- seq(from = l$months_to_nstart , 
                            to= length(l$looks) * l$interim_period, 
                            by = l$interim_period)
      
      if(length(l$interimmnths) < length(l$looks)){
        l$interimmnths <- c(l$interimmnths, max(l$interimmnths)+l$interim_period)
      }
      
      stopifnot(length(l$interimmnths) == length(l$looks))
      stopifnot(max(l$looks) == l$nstop)
      
      
      flog.info("Updated interimmnths: %s.", paste0(l$interimmnths, collapse = ", "))
    }
    
    
    if(!is.null(opt$basesero)){
      l$baselineprobsero <- opt$basesero
      flog.info("Updated baselineprobsero: %s.", l$baselineprobsero)

      l$deltaserot3 <- compute_sero_delta(l$baselineprobsero, l$trtprobsero)
      flog.info("Updated deltaserot3: %s.", l$deltaserot3)
    }
    
    if(!is.null(opt$trtprobsero)){
      l$trtprobsero <- opt$trtprobsero
      flog.info("Updated trtprobsero: %s.", l$trtprobsero)
      
      l$deltaserot3 <- compute_sero_delta(l$baselineprobsero, l$trtprobsero)
      flog.info("Updated deltaserot3: %s.", l$deltaserot3)
      
    }
    
    if(!is.null(opt$basemediantte)){
      l$ctl_med_tte <- opt$basemediantte
      flog.info("Updated ctl_med_tte: %s.", l$ctl_med_tte)
      
      l$b0tte <- log(2)/l$ctl_med_tte 
      l$b1tte <- (log(2)/l$trt_med_tte) - l$b0tte
      flog.info("Updated b0tte: %s.", l$b0tte)
      flog.info("Updated b1tte: %s.", l$b1tte)
    }
    
    if(!is.null(opt$trtmedtte)){
      l$trt_med_tte <- opt$trtmedtte
      flog.info("Updated trt_med_tte: %s.", l$trt_med_tte)

      l$b0tte <- log(2)/l$ctl_med_tte 
      l$b1tte <- (log(2)/l$trt_med_tte) - l$b0tte
      flog.info("Updated b0tte: %s.", l$b0tte)
      flog.info("Updated b1tte: %s.", l$b1tte)
    }
    
    # write.csv(opt, "test.csv")
 
  }
  
  
  
  
  

  return(l)
}




den_plot1 <- function(a, b, xlim = NULL){
  plot(density(a), xlim = xlim)
  lines(density(b), col = "red")
  abline(v = mean(a))
  abline(v = mean(b), col = "red")
  #legend("topleft", "(x,y)")
  legend("topleft", legend=c("a", "b"),
         col=c("black", "red"), lty=c(1,1), cex=0.8)
}

print_immu_res <- function(m){
  
  df <- as.data.frame(m)
  
  df <- df[, c("idxsim",
               "look", 
               "n_obs",
               "ss_immu",
               "n_max",
               "n_max_sero",
               "i_ppos_n", "i_ppos_max", "i_post_prop_ctl", "i_post_prop_trt",   # immunological
               "stop_v_samp",
               "stop_i_fut")]
  
  df
}

print_clin_res <- function(m){
  
  df <- as.data.frame(m)
  
  df <- df[, c("idxsim",
               "look", 
               "n_obs",
               "ss_clin",
               "n_max",
               "n_max_sero",
               "c_ppos_n","c_ppos_max","c_post_lambda_ctl","c_post_lambda_trt",  # clinical
               "stop_c_fut",
               "stop_c_sup")]
  
  df
}


stop_immu <- function(stop_ven_samp,
                      stop_immu_fut,
                      stop_clin_fut,
                      stop_clin_sup, n_look, n_maxsero){
  
  if(stop_ven_samp){
    return(T)
  }
  
  if(n_look > n_maxsero){
    return(T)
  }
  
  if(stop_clin(stop_immu_fut,
               stop_clin_fut,
               stop_clin_sup)){
    return(T)
  }

  return(F)
}


stop_clin <- function(stop_immu_fut,
                      stop_clin_fut,
                      stop_clin_sup){
  
  if(stop_immu_fut){
    return(T)
  }
  
  if(stop_clin_fut){
    return(T)
  }
  
  if(stop_clin_sup){
    return(T)
  }
  
  return(F)
}

surv_lnorm <-function(tte = 1:100, mu = 3, sig = 1){
  n <- length(tte)
  s <- numeric(n)
  surv <- function(x, mu, sig){
    1 - pnorm(log(tte[x]) - mu)/sig
  }
  s <- unlist(lapply(1:n, surv, mu, sig))
  s
}

mu_of_lognormal <- function(m, upsilon){
  log((m^2) / sqrt(upsilon + (m^2)))
}
sig_of_lognormal <- function(m, upsilon){
  sqrt(log((upsilon/m^2) + 1))
}

# The mean and variance of a lognormal distribution created
# with values mu and sig will be:
m_of_lognormal <- function(mu, sig){
  exp(mu + ((sig^2)/2))
}
upsilon_of_lognormal <- function(mu, sig){
  exp(2*mu + sig^2) * (exp(sig^2) - 1)
}

muformediansurvtime <- function(medtime, sig){
  log(medtime) - sig * qnorm(0.5) 
}

init <- function(){
  ggplot2::theme_set(theme_bw())
  ggplot2::theme_update(legend.position="bottom")
  ggplot2::theme_update(legend.title=element_blank())
  # See http://ggplot2.tidyverse.org/reference/theme.html
  ggplot2::theme_update(text=element_text(size=12,  family="sans"))
  ggplot2::theme_update(axis.text.x=element_text(size=10,  family="sans"))
  ggplot2::theme_update(axis.text.y=element_text(size=10,  family="sans"))
  f.sep <- .Platform$file.sep
}

# logit to p
inv_logit <- function(x){
  return(exp(x)/(1+exp(x)))
}
# p to logit
logit <- function(p){
  return(log(p/(1-p)))
}
prob_to_odd <- function(x){
  return(x/(1-x))
}
odd_to_prob <- function(x){
  return(x/(1+x))
}

get_null_sero <- function(cfgfile = "cfg1.yaml"){
  
  fs0 <- readRDS("fso.RDS")
  
  if(!is.null(fs0)){
    return(fs0)
  }
  
  cfg <- sim_cfg(cfgfile)
  dt1 = gen_dat(cfg)
  
  # obtain data and fit models to save in global scope (prevents recompile)
  dat <- make_standata(formula(cfg$fsero), 
                       data = dt1,
                       family = bernoulli())
  
  fs0 <- stan(file=cfg$seromodfile, 
              data=dat, iter=10, chains=1)
  
  saveRDS(fs0, "fso.RDS")
  
  return(fs0)
}


get_null_tte <- function(cfgfile = "cfg1.yaml"){
  
  ft0 <- readRDS("fto.RDS")
  
  if(!is.null(ft0)){
    return(ft0)
  }
  
  cfg <- sim_cfg(cfgfile)
  dt1 = gen_dat(cfg)
  
  # obtain data and fit models to save in global scope (prevents recompile)
  dat <- make_standata(formula(cfg$ftte), 
                       data = dt1,
                       family = lognormal())
  
  ft0 <- stan(file=cfg$ttemodfile, 
              data=dat, iter=10, chains=1)
  
  saveRDS(ft0, "fto.RDS")
  
  return(ft0)
}



med_surv <- function(eta){
  
  log(2)/eta
  
}

stan_model_workflow <- function(){
  
  # https://cran.r-project.org/web/packages/rstan/vignettes/rstan.html
  # A single call to stan performs all three steps, but they can also be
  # executed one by one (see the help pages for stanc, stan_model, and
  # sampling),
  
  stancode <- 'data {real y_mean;} parameters {real y;} model {y ~ normal(y_mean,1);}'
  mod <- stan_model(model_code = stancode)
  fit <- sampling(mod, data = list(y_mean = 0))
  fit2 <- sampling(mod, data = list(y_mean = 5))
  
}

save_demo_dat <- function(){
  
  cfg <- sim_cfg("cfg1.yaml")
  dt1 = gen_dat(cfg)
  
  saveRDS(dt1, "mydat.RDS")
}

write_mod1 <- function(){
  
  cfg <- sim_cfg()
  dt1 = gen_dat(cfg)
  
  dt1$serot3[1:5] <- NA
  
  mf <- bf(cfg$fsero)
  mp <- get_prior(mf,
                  data = dt1,
                  family = bernoulli())
  
  blm0 <- brms::brm(mf, 
                    data = dt1,
                    family = bernoulli(), 
                    prior = mp,
                    iter = 10,
                    chains = 1, 
                    save_model = "brm1.stan",
                    control = list(max_treedepth = 10))
  # summary(blm0, waic = TRUE)
}


write_mod2 <- function(){
  # , data = dt1, family = lognormal()
  # 
  # 
  
  cfg <- sim_cfg()
  dt1 = gen_dat(cfg)
  mf <- bf(evtt | cens(cen) ~ trt + remote)
  mp <- get_prior(mf,
                  data = dt1,
                  family = lognormal())
  
  blm0 <- brms::brm(mf, 
                    data = dt1,
                    family = lognormal(), 
                    prior = mp,
                    iter = 10,
                    chains = 1, 
                    save_model = "brm2.stan",
                    control = list(max_treedepth = 10))
 
}


check_res <- function(){
  
  
  dfdur <- readRDS("2018_12_19_1633/duration.RDS")
  dfres <- readRDS("2018_12_19_1633/dfresults.RDS")
  dfwarn <- readRDS("2018_12_19_1633/warn.RDS")

  # fek
 
}








plot_surv <- function(e3){
  dfnew = data.frame(trt = 0:1, remote = 0)
  predict(e3, newdata = dfnew, type=c("response"), p = 0.5, se = T)
  
  plot(predict(e3, newdata=list(trt=0, remote = 0),
               type="quantile",
               p=seq(.01,.99,by=.01)),
       seq(.99,.01,by=-.01), col="red", type = "l", xlim = c(0, 100))
  
  lines(predict(e3, newdata=list(trt=1, remote = 0),
                type="quantile",
                p=seq(.01,.99,by=.01)),
        seq(.99,.01,by=-.01), col="blue", type = "l", xlim = c(0, 100))
}



jags_init <- function(d, omit_const = T){
  
  # simulated event times for the rest of the trial
  d$Y <- copy(d$evtt)
  
  total_acc_time <- max(d$accrt) + 1
  
  # If the simulated event time plus the accrual time > total_acc_time then censor
  d$cen <- ifelse(d$Y + (d$accrt/cfg$dayspermonth) >  
                               total_acc_time/cfg$dayspermonth, 1, 0)
  d$Y <- ifelse(d$cen, NA, d$Y)
  d$Y_cen <- ifelse(d$cen, (total_acc_time - d$accrt)/cfg$dayspermonth, 0)
  d$Y_all <- ifelse(!is.na(d$Y), d$Y, d$Y_cen)
  
  dat <- list(Y = d$Y,
               Y_cen = d$Y_cen,
               cen = as.logical(d$cen),
               trt = d$trt,
               remote = d$remote)
  if(omit_const == F){
    dat$N = length(dat$Y)
  }
  
  
  init_Y <- rep(NA, length(d$cen))
  init_Y[dat$cen] <- dat$Y_cen[dat$cen]+1
  myinits <- list( 'Y' = init_Y,
          'tau' = runif(1),
          'b0' = rnorm(1),
          'b1' = rnorm(1))
  
  return(list(inits = myinits,
              dat = dat))
}

