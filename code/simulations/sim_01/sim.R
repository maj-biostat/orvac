





gen_dat <- function(cfg, n_override = NULL){

  n <- ifelse(is.null(n_override), cfg$nstop, n_override)

  dt <- data.table(id = numeric(n),
                   trt = numeric(n),
                   remote = numeric(n),
                   accrt = numeric(n),      # time of randomisation
                   age_months = numeric(n), # age at time of randomisation
                   trial_months_at_max_age = numeric(n),
                   serot2 = rep(NA, n),
                   serot3 = rep(NA, n),
                   probt3 = numeric(n),
                   evtt = numeric(n),
                   # cent = numeric(n), # time at which censored
                   cent_at_interim = rep(NA, n),
                   cent_at_max = rep(NA, n),
                   obst = rep(NA, n),
                   cen = rep(NA, n)) # censoring variable (0 = not censored, 1 = censored)

  dt[,id := 1:n]
  dt[,trt := rep(0:1, length.out = n)]
  dt[,remote := rbinom(n, 1, prob = cfg$remoteprob)]

  # age is at time of randomisation, i.e. at their accrual time
  dt[,age_months := rtruncnorm(n,
                               cfg$age_months_lwr,
                               cfg$age_months_upr,
                               cfg$age_months_mean,
                               cfg$age_months_sd)]


  dt[,accrt := (id - 1) * cfg$months_per_person]

  dt[,trial_months_at_max_age := cfg$max_age_fu_months - age_months + accrt]

  # seroconversion
  dt$serot2[1:cfg$nmaxsero] = rbinom(cfg$nmaxsero, 1, prob = cfg$baselineprobsero)
  dt$serot3[1:cfg$nmaxsero] = copy(dt$serot2[1:cfg$nmaxsero])

  setkey(dt, id, trt, serot3)
  dt[.(1:cfg$nmaxsero, 1, 0), probt3:= cfg$deltaserot3 * trt  ]
  n_t3_no_sero <- nrow(dt[serot2 == 0 & !is.na(serot2) & trt == 1])
  # n_t3_no_sero

  dt[serot2 == 0 & !is.na(serot2) & trt == 1, serot3 := rbinom(n_t3_no_sero, size = 1, prob = probt3)]

  # tte - measured in months - see cfg
  # tte is from zero origin, needs to be offset to account for accrual time.
  mu <- cfg$b0tte + dt$trt * cfg$b1tte
  dt[, evtt := stats::rexp(n, rate = mu)]

  # censoring time at interim are all different so not set here.
  # for convenience/consistency, censoring at max also set at time of interim

  # dt %>%
  #   dplyr::filter(!is.na(serot3)) %>%
  #   dplyr::group_by(trt) %>%
  #   dplyr::summarise(mean = mean(serot3))



  dt
}






# NEED TO CONSIDER DELAYED RESPONSE DUE TO TIME WINDOW THAT TESTS WILL NEED TO COME BACK.

# this model looks at a dichotomous outcome via conjugate beta-binmoial
model_immu <- function(d, cfg, look, idx){
  
  delta <- NA
  ppos_n <- NA
  delta_mean <- NA
  delta_wald_lwr_95 <- NA
  delta_wald_upr_95 <- NA
  ppos_max <- NA
  
  interim_month <- cfg$interimmnths[look]
  flog.info("model_immu_2 current simulation %s, look %s, interim_month = %s", idx, look, interim_month)


  # Both the interim analysis and the final are based on predictive prob.
  # While cfg$looks[look] can be > cfg$nmaxsero we only make the assessment on cfg$nmaxsero.
  if(cfg$looks[look] >= cfg$nmaxsero){
    dat_ctl <- immu_dat(d, cfg, look, n_target = cfg$nmaxsero, info_delay = 0, trt_status = 0)
    dat_trt <- immu_dat(d, cfg, look, n_target = cfg$nmaxsero, info_delay = 0, trt_status = 1)
    # final analysis for seroconversion looks at posterior
    delta <- dat_trt$theta - dat_ctl$theta
    ppos_n <- mean(delta > 0)
    
    delta_mean <- mean(delta)
    delta_wald_lwr_95 <- delta_mean + qnorm(0.025) * sd(delta)
    delta_wald_upr_95 <- delta_mean + qnorm(0.975) * sd(delta)
    
    p <- c(ppos_n, ppos_max, mean(dat_ctl$theta), mean(dat_trt$theta), 
           delta_mean, delta_wald_lwr_95, delta_wald_upr_95)
    
  } else {
    # get relevant data for each arm
    # n_max and n_look are divided by 2 in the called function
    dat_ctl <- immu_dat(d, cfg, look, n_target = cfg$looks[look], info_delay = cfg$sero_info_delay, trt_status = 0)
    dat_trt <- immu_dat(d, cfg, look, n_target = cfg$looks[look], info_delay = cfg$sero_info_delay, trt_status = 1)
    # this is the posterior based on the observed seroconversion events
    tctl <- mean(dat_ctl$theta)
    ttrt <- mean(dat_trt$theta)
    
    delta <- dat_trt$theta - dat_ctl$theta
    delta_mean <- mean(delta)
    delta_wald_lwr_95 <- delta_mean + qnorm(0.025) * sd(delta)
    delta_wald_upr_95 <- delta_mean + qnorm(0.975) * sd(delta)
    
    # interim analysis based on predictive prob due to info delay
    ppos_n <- immu_pred_prob(dat_trt, dat_ctl, cfg, n_target = cfg$looks[look])
    
    dat_ctl <- immu_dat(d, cfg, look, n_target = cfg$nmaxsero, info_delay = cfg$sero_info_delay, trt_status = 0)
    dat_trt <- immu_dat(d, cfg, look, n_target = cfg$nmaxsero, info_delay = cfg$sero_info_delay, trt_status = 1)
    ppos_max <- immu_pred_prob(dat_trt, dat_ctl, cfg, n_target = cfg$nmaxsero)
    
    p <- c(ppos_n, ppos_max, tctl, ttrt, 
           delta_mean, delta_wald_lwr_95, delta_wald_upr_95)
  } 

  names(p) <- cfg$immu_rtn_names
  
  flog.info("model_immu_2 current simulation %s, look %s, result %s", idx, look,
            paste0("(", paste(p, collapse = ", "), ")"))
  return(p)
}


# function to return relevant data by treatment group
# information delay is built into the calcs
# n_target is what n_obs (total obs) would be at this interim if we had no delayed data
immu_dat <- function(d, cfg, look, n_target, info_delay, trt_status = 0){
  
  accrualtimes <- d[d$accrt < (cfg$interimmnths[look] - info_delay) & 
                      trt == trt_status & 
                      id <= cfg$nmaxsero, accrt]
  
  # currently we can only see up to n_obs_grp observations for this group -
  # i.e. there are n_obs_grp that have been randomised to this arm
  n_obs_grp <- length(accrualtimes)
  n_seroconverted <- sum(d$serot3[d$trt == trt_status][1:n_obs_grp])
  n_impute <- floor(n_target/2) - n_obs_grp

  # draw from _posterior_ then the dgp distribution (rbinom) to impute the 'unobserved' values
  # these two lines are equivalent to making draws from the posterior predictive dist.
  # theta is posterior probability of seroconversion
  theta <- rbeta(n = cfg$post_draw,
                 shape1 = 1 + n_seroconverted,
                 shape2 = 1 + n_obs_grp - n_seroconverted)
  
 
  # immu_impute
  dat_interims <- NULL
  if(n_impute > 0){
    pp <- lapply(theta, rbinom, n = n_impute, size = 1)
    stopifnot(sum(unlist(lapply(pp, anyNA))) == 0)
    
    # combine the observed and pp
    combine_dat <- function(x){
      y_rep_x <- c(d$serot3[d$trt == trt_status][1:n_obs_grp], pp[[x]])
      
      n_sero_by_max_x <- sum(y_rep_x)
      return(list(y_rep_x = y_rep_x, n_sero_by_max_x = n_sero_by_max_x))
    }
    dat_interims <- lapply(1:length(pp), combine_dat)
  }

  
  return(list(n_obs_grp = n_obs_grp,
              n_impute = n_impute,
              n_max_grp = n_obs_grp + n_impute,
              n_seroconverted = n_seroconverted,
              theta = theta,
              dat_interims = dat_interims))
}




# compute difference over all simulated trials
immu_pred_prob <- function(dat_trt, dat_ctl, cfg, n_target){
  
  
  # compute difference distribution
  diff_prop <- function(x){
    
    # look at the complete trial for futility
    # x = x + 1
    n_sero_by_max_trt = dat_trt$dat_interims[[x]]$n_sero_by_max_x
    n_sero_by_max_ctl = dat_ctl$dat_interims[[x]]$n_sero_by_max_x
    
    # first look at the interim for superiority
    prop_trt <- rbeta(n = cfg$post_draw,
                      shape1 = 1 + n_sero_by_max_trt,
                      shape2 = 1 + floor(n_target/2) - n_sero_by_max_trt)

    prop_ctl <- rbeta(n = cfg$post_draw,
                      shape1 = 1 + n_sero_by_max_ctl,
                      shape2 = 1 + floor(n_target/2) - n_sero_by_max_ctl)
    
    # the distribution of the difference will be approximately normal
    diff <- prop_trt - prop_ctl
    
    # den_plot1(prop_trt, prop_ctl, xlim = c(0,1))
    # plot(density(diff))
    
    # do we see a treatment effect?
    win <- mean(diff > 0) > cfg$post_sero_thresh
    # win
    
    return(win)
  }
  
  res <- unlist(lapply(1:cfg$post_draw, diff_prop))
  
  # predictive probability for:
  # superiority (at current look) and
  # futility (assuming trial expends all resources)
  ppos_max <- mean(res, na.rm = T)
  
  ppos_max
  
}





model_clin <- function(d, cfg, look, idx){

  interim_month <- cfg$interimmnths[look]
  
  stopifnot(!is.na(interim_month) && !is.null(interim_month) && interim_month > 0)
  n_obs <- cfg$looks[look]
  flog.info("model_clin current simulation %s, look %s, interim_month %s", idx, look, cfg$interimmnths[look])

  # get relevant data for each arm
  # n_max and n_look are divided by 2 in the called function
  dat_ctl <- clin_dat(d, cfg, look, trt_status = 0)
  dat_trt <- clin_dat(d, cfg, look, trt_status = 1)

  # inference at the interim is based on the ratio of lambda_trt / lambda_ctl.
  # both distributions are gammas with integer shape parameter therefore we can transform
  # the variables to scaled chisqu and then the ratio is distributed as F (see page 316 bayesian ideas).
  # but we can also do this empircally by sampling directly to get a relative median because
  # median_1 = log(2)/lamda_1
  # median_2 = log(2)/lamda_2
  # median_1/median_2 = (log(2)/lambda_1)/(log(2)/lamda_2)
  #                   = lambda_2 / lambda_1
  # i.e. just remember to compute the right way around.
 
  # if there is an effect 'at the 97.5%' level then ppos_n
  # will be > 0.975
  delta <- dat_ctl$lambda / dat_trt$lambda
  ppos_n <- mean(delta > 1 , na.rm = T)

  # Wald isn't great but it is quick, the alternative is to compute a density 
  # kernel and then obtain the quantiles from that.
  delta_mean <- mean(delta)
  delta_wald_lwr_95 <- delta_mean + qnorm(0.025) * sd(delta)
  delta_wald_upr_95 <- delta_mean + qnorm(0.975) * sd(delta)
  
  # not applicable for the final analysis
  ppos_max <- NA
  if(look < length(cfg$looks)){
    ppos_max <- clin_pred_prob(dat_trt, dat_ctl, cfg)
  } 
  p <- c(ppos_n, ppos_max, mean(dat_ctl$lambda), mean(dat_trt$lambda), 
         delta_mean, delta_wald_lwr_95, delta_wald_upr_95)
  
  names(p) <- cfg$clin_rtn_names

  flog.info("model_clin current simulation %s, look %s, n_obs %s, result %s",
            idx, look, n_obs,
            paste0("(", paste(p, collapse = ", "), ")"))

  return(p)
}




# Extract the relevant datasets for this interim.
# Also initiates the call to create simulated data from this interim to max trial size.
clin_dat <- function(d, cfg, look, trt_status = 0){
  
  interim_month <- cfg$interimmnths[look]

  # < otherwise you get one extra obs than required
  accrualtimes <- d[d$accrt < interim_month & trt == trt_status, accrt]
  # currently we can only see up to n_obs_grp observations for this group -
  n_obs_grp <- length(accrualtimes)
  # we need to know the age at interim as we only fu for a max period of 36 months
  age_at_interim <- d[d$accrt < interim_month & trt == trt_status, age_months] + interim_month - accrualtimes 
  
  # if event time + age_at_rand > 36 then censor because we only follow to 36 months
  # note - evtt (event time) is the time from randomisation to the event
  cen_current <- ifelse(d$evtt[d$trt == trt_status][1:n_obs_grp] +
                          d$age_months[d$trt == trt_status][1:n_obs_grp] >
                          cfg$max_age_fu_months , 1, 0)

  if (look < length(cfg$looks)){
    # if the evtt + time_of_rand > month of interim then censor because the event
    # happens some after this interim
    cen_current <- ifelse(d$evtt[d$trt == trt_status][1:n_obs_grp] + accrualtimes >
                            interim_month , 1, cen_current)
    # what is their censoring time?
    cent_current <- pmin(cfg$max_age_fu,
                         cfg$interimmnths[look] - d$accrt[d$trt == trt_status][1:n_obs_grp])
  } else {
    # and now use the final_analysis_month as the reference point for censoring
    cen_current <- ifelse(d$evtt[d$trt == trt_status][1:n_obs_grp] + accrualtimes >
                            cfg$final_analysis_month , 1, cen_current)
    # assuming that individual i was censored, what is their censoring time?
    cent_current <- pmin(cfg$max_age_fu,
                         cfg$final_analysis_month - d$accrt[d$trt == trt_status][1:n_obs_grp])
  }
  
  # consolidate the evtt and the censoring time at interim into obst
  obst_current <- ifelse(cen_current == 1, cent_current, d$evtt[d$trt == trt_status][1:n_obs_grp])
  
  # this is the data we currently see in this trt_status group
  d_current <- cbind(idx = 1:n_obs_grp,
                     interim_month = interim_month,
                     age_at_rand = d$age_months[d$trt == trt_status][1:n_obs_grp],
                     age_at_interim = age_at_interim,
                     accrt = accrualtimes,
                     evtt = d$evtt[d$trt == trt_status][1:n_obs_grp],
                     cen_current = cen_current,
                     cent_current = cent_current,
                     obst_current = obst_current)
  # number of uncensored obs and total of time to events
  n_uncen <- sum(1-cen_current)
  sum_obst <- sum(obst_current)
  
  # project event times to the end of the trial
  n_impute <- floor(cfg$nstop/2) - n_obs_grp
  stopifnot(2*(n_impute + n_obs_grp) <= cfg$nstop)
  
  # posterior for lambda
  lambda <- rgamma(cfg$post_draw, shape = 1 + n_uncen, rate = 0.01 + sum_obst)
  
  dat_interims <- NULL
  if (look < length(cfg$looks)){
    dat_interims <- clin_impute(d, lambda, trt_status, 
                                n_impute, n_obs_grp, 
                                n_uncen, sum_obst, 
                                cfg)
  }
  
  return(list(n_obs_grp = n_obs_grp,
              n_impute = n_impute,
              n_max_grp = n_obs_grp + n_impute,
              n_uncen = n_uncen,
              sum_obst = sum_obst,
              lambda = lambda,
              dat_interims = dat_interims))
}






# Creates an imputed dataset (from interim to end of trial) for each posterior 
# draw in lambda
clin_impute <- function(d, lambda, trt_status, n_impute, n_obs_grp, 
                        n_uncen, sum_obst, cfg){
  
  # now to simulate many datasets that reflect possible scenarios
  # at the end of the trial assuming accrual remains at a constant rate
  
  # draw from posterior predictive to impute the 'unobserved' values
  # using gamma(1,1) prior
  pp <- lapply(lambda, rexp, n = n_impute)
  stopifnot(sum(unlist(lapply(pp, anyNA))) == 0)
  
  # combine the observed and pp
  combine_dat <- function(x){
    evtt_rep <- c(d$evtt[d$trt == trt_status][1:n_obs_grp], pp[[x]])
    return(evtt_rep)
  }
  evtt_rep <- lapply(1:length(pp), combine_dat)
  
  # imputed datasets - only do this if not on last trial
  hypoth_interim <- function(x){
    
    accrt_x <- d$accrt[d$trt == trt_status][1:(n_obs_grp + n_impute)]
    evtt_rep_x <- evtt_rep[[x]]
    
    cen_x <- ifelse(evtt_rep_x +
                      d$age_months[d$trt == trt_status][1:(n_obs_grp + n_impute)] >
                      cfg$max_age_fu_months , 1, 0)
    
    # our projection here is to the final analysis hence we compare to the 
    # final analysis month (which is the month when the youngest (at the last interim) reaches 36 months)
    cen_x <- ifelse(evtt_rep_x + accrt_x > cfg$final_analysis_month , 1, cen_x)
    cent_x <- pmin(cfg$max_age_fu, cfg$final_analysis_month - accrt_x)
    obst_x <- ifelse(cen_x == 1, cent_x, evtt_rep_x)
    
    m <- cbind(accrt_x = accrt_x,
               evtt_rep_x = evtt_rep_x,
               cen_x = cen_x,
               cent_x = cent_x,
               obst_x = obst_x)
    
    n_uncen_x = sum(1 - cen_x)
    sum_obst_x = sum(obst_x)
    
    stopifnot(length(accrt_x) == length(evtt_rep_x))
    stopifnot(length(accrt_x) == length(cen_x))
    stopifnot(length(accrt_x) == length(cent_x))
    stopifnot(length(accrt_x) == length(obst_x))
    
    # Optimisation - remove m
    return(list(m = m, n_uncen_x = n_uncen_x, sum_obst_x = sum_obst_x))
  }
  
  dat_interims <- lapply(1:length(evtt_rep), hypoth_interim)
  dat_interims
}









# This computes the predictive probabilities for assessing futility.
# Conditional imputed tte data from this interim to the end of the trial,
# did we see a difference in median surv time in more than 5% of the datasets?
clin_pred_prob <- function(dat_trt, dat_ctl, cfg){
  
  
  # compute difference distribution
  diff_med <- function(x){
    
    # look at the complete trial for futility
    n_uncen_trt = dat_trt$dat_interims[[x]]$n_uncen_x
    sum_obst_trt = dat_trt$dat_interims[[x]]$sum_obst_x
    
    n_uncen_ctl = dat_ctl$dat_interims[[x]]$n_uncen_x
    sum_obst_ctl = dat_ctl$dat_interims[[x]]$sum_obst_x
    
    # futility
    lamb_trt <- rgamma(cfg$post_draw, shape = 1 + n_uncen_trt, rate = 0.01 + sum_obst_trt)
    lamb_ctl <- rgamma(cfg$post_draw, shape = 1 + n_uncen_ctl, rate = 0.01 + sum_obst_ctl)
    
    # the distribution of the difference can be approximate as normal
    # see earlier narrative comment
    ratio_medians <- lamb_ctl / lamb_trt
    
    # do we see a treatment effect?
    win <- mean(ratio_medians > 1) > cfg$post_tte_thresh
    # win
    
    return(win)
  }
  
  res <- unlist(lapply(1:cfg$post_draw, diff_med))
  
  # predictive probability for:
  # superiority (at current look) and
  # futility (assuming trial expends all resources)
  ppos_max <- mean(res, na.rm = T)
  
  ppos_max
  
}


