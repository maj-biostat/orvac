





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
model_immu_2 <- function(d, cfg, look, idx){

  interim_month <- cfg$interimmnths[look]

  flog.info("model_immu_2 current simulation %s, look %s, interim_month = %s", idx, look, interim_month)

  # function to return relevant data by treatment group
  do_dat <- function(trt_status = 0){

    accrualtimes <- d[d$accrt < interim_month & trt == trt_status & id <= cfg$nmaxsero, accrt]

    # currently we can only see up to n_obs_grp observations for this group -
    # i.e. there are n_obs_grp that have been randomised to this arm
    n_obs_grp <- length(accrualtimes)

    n_seroconverted <- sum(d$serot3[d$trt == trt_status][1:n_obs_grp])
    n_impute <- floor(cfg$nmaxsero/2) - n_obs_grp

    # draw from posterior then the dgp distribution (rbinom) to impute the 'unobserved' values
    # these two lines are equivalent to making draws from the posterior predictive dist.

    # theta is probability of seroconversion
    theta <- rbeta(n = cfg$post_draw,
                   shape1 = 1 + n_seroconverted,
                   shape2 = 1 + n_obs_grp - n_seroconverted)
    pp <- lapply(theta, rbinom, n = n_impute, size = 1)

    # horrendous fix to address NA draws from binomial. there must be a prettier way...
    pp_incna <- which(unlist(lapply(pp, anyNA)) > 0)
    trytimes <- 0
    while(length(pp_incna) > 0 || trytimes > 100){
      trytimes <- trytimes + 1
      tmp <- lapply(lambda[pp_incna], rexp, n = n_impute)
      for(i in 1:length(pp_incna)){
        flog.debug("model_immu_2 attempting to fill %s NAs in pp", length(pp_incna))
        pp[[pp_incna[i]]] <- tmp[[i]]
      }
      pp_incna <- which(unlist(lapply(pp, anyNA)) > 0)
    }
    stopifnot(!is.null(which(unlist(lapply(pp, anyNA)) > 0)))



    # combine the observed and pp
    combine_dat <- function(x){
      y_rep_x <- c(d$serot3[d$trt == trt_status][1:n_obs_grp], pp[[x]])

      n_sero_by_max_x <- sum(y_rep_x)
      return(list(y_rep_x = y_rep_x, n_sero_by_max_x = n_sero_by_max_x))
    }
    dat_interims <- lapply(1:length(pp), combine_dat)

    return(list(n_obs_grp = n_obs_grp,
                n_impute = n_impute,
                n_max_grp = n_obs_grp + n_impute,
                n_seroconverted = n_seroconverted,
                theta = theta,
                dat_interims = dat_interims))
  }

  # get relevant data for each arm
  # n_max and n_look are divided by 2 in the called function
  dat_ctl <- do_dat(0)
  dat_trt <- do_dat(1)

  # den_plot1(dat_trt$theta, dat_ctl$theta)
  # plot(density(dat_trt$theta - dat_ctl$theta))

  # difference between proportions from the posteriors.
  # if there is an effect 'at the 97.5%' level then ppos_n
  # will be > 0.975
  ppos_n <- mean((dat_trt$theta - dat_ctl$theta) > 0)
  # plot(density(dat_trt$theta - dat_ctl$theta))



  # compute difference distribution
  diff_prop <- function(x){

    # look at the complete trial for futility
    # x = x + 1
    n_sero_by_max_trt = dat_trt$dat_interims[[x]]$n_sero_by_max_x
    n_sero_by_max_ctl = dat_ctl$dat_interims[[x]]$n_sero_by_max_x

    # first look at the interim for superiority
    prop_trt <- rbeta(n = cfg$post_draw,
                   shape1 = 1 + n_sero_by_max_trt,
                   shape2 = 1 + floor(cfg$nmaxsero/2) - n_sero_by_max_trt)

    prop_ctl <- rbeta(n = cfg$post_draw,
                   shape1 = 1 + n_sero_by_max_ctl,
                   shape2 = 1 + floor(cfg$nmaxsero/2) - n_sero_by_max_ctl)

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
  # hist(res)
  ppos_max <- mean(res)

  delta <- dat_trt$theta - dat_ctl$theta
  delta_mean <- mean(delta)
  delta_wald_lwr_95 <- delta_mean + qnorm(0.025) * sd(delta)
  delta_wald_upr_95 <- delta_mean + qnorm(0.975) * sd(delta)

  p <- c(ppos_n, ppos_max, mean(dat_ctl$theta), mean(dat_trt$theta), 
         delta_mean, delta_wald_lwr_95, delta_wald_upr_95)

  names(p) <- cfg$immu_rtn_names
  
  flog.info("model_immu_2 current simulation %s, look %s, result %s", idx, look,
            paste0("(", paste(p, collapse = ", "), ")"))
  return(p)
}













model_clin_1 <- function(d, cfg, look, idx){

  interim_month <- cfg$interimmnths[look]
  stopifnot(!is.na(interim_month) && !is.null(interim_month) && interim_month > 0)

  n_obs <- cfg$looks[look]

  # currently we can only see up to nobs observations - there are n_obs that have been randomised
  # this is total obs in both arms
  # n_obs <- which(d$accrt == max(d$accrt[d$accrt < interim_month]))
  #
  # # n_max (in total, i.e. n_trt + n_ctl) dependent on accrual rate -
  # n_max = min(max(cfg$looks), floor(  max(cfg$interimmnths) / cfg$months_per_person  ))

  flog.info("model_clin_1 current simulation %s, look %s, interim_month %s", idx, look, interim_month)

  # function to return relevant data
  do_dat <- function(trt_status = 0){

    # note - check the <. I think this is correct otherwise you get one extra obs than required
    accrualtimes <- d[d$accrt < interim_month & trt == trt_status, accrt]

    # currently we can only see up to n_obs_grp observations for this group -
    # i.e. there are n_obs_grp that have been randomised to this arm
    n_obs_grp <- length(accrualtimes)

    # we need to know the age at interim as we only fu for a max period of 36 months
    age_at_interim <- interim_month - accrualtimes + d[d$accrt < interim_month & trt == trt_status, age_months]

    # if event time + age_at_rand > 36 then censor because we only follow to 36 months
    # note - evtt (event time) is the time from randomisation to the event
    cen_current <- ifelse(d$evtt[d$trt[1:n_obs] == trt_status][1:n_obs_grp] +
                            d$age_months[d$trt[1:n_obs] == trt_status][1:n_obs_grp] >
                            cfg$max_age_fu_months , 1, 0)

    # if the evtt + time_of_rand > month of interim then censor because the event
    # happens some after this interim
    cen_current <- ifelse(d$evtt[d$trt[1:n_obs] == trt_status][1:n_obs_grp] + accrualtimes >
                            interim_month , 1, cen_current)

    # what is their censoring time?
    cent_current <- pmin(cfg$max_age_fu,
                         cfg$interimmnths[look] - d$accrt[d$trt[1:n_obs] == trt_status][1:n_obs_grp])



    # consolidate the evtt and the censoring time at interim into obst
    obst_current <- ifelse(cen_current == 1, cent_current,
                           d$evtt[d$trt[1:n_obs] == trt_status][1:n_obs_grp])

    # this is the data we currently see in this trt_status group
    # everything is in months
    d_current <- cbind(idx = 1:n_obs_grp,
                       interim_month = interim_month,
                       age_at_rand = d$age_months[d$trt[1:n_obs] == trt_status][1:n_obs_grp],
                       age_at_interim = age_at_interim,
                       accrt = accrualtimes,
                       evtt = d$evtt[d$trt[1:n_obs] == trt_status][1:n_obs_grp],
                       cen_current = cen_current,
                       cent_current = cent_current,
                       obst_current = obst_current)

    # number of uncensored obs and total of time to events
    n_uncen <- sum(1-cen_current)
    sum_obst <- sum(obst_current)

    # now to simulate many datasets that reflect possible scenarios
    # at the end of the trial assuming accrual remains at a constant rate

    # project event times to the end of the trial
    n_impute <- floor(cfg$nstop/2) - n_obs_grp

    stopifnot(2*(n_impute + n_obs_grp) <= cfg$nstop)

    # draw from posterior predictive to impute the 'unobserved' values
    # using gamma(1,1) prior
    lambda <- rgamma(cfg$post_draw, shape = 1 + n_uncen, rate = 1 + sum_obst)
    pp <- lapply(lambda, rexp, n = n_impute)

    # horrendous fix to address NA draws from exponential. there must be a prettier way...
    pp_incna <- which(unlist(lapply(pp, anyNA)) > 0)
    trytimes <- 0
    while(length(pp_incna) > 0 | trytimes > 100){
      trytimes <- trytimes + 1
      tmp <- lapply(lambda[pp_incna], rexp, n = n_impute)
      for(i in 1:length(pp_incna)){
        flog.debug("attempting to fill %s NAs in pp", length(pp_incna))
        pp[[pp_incna[i]]] <- tmp[[i]]
      }
      pp_incna <- which(unlist(lapply(pp, anyNA)) > 0)
    }
    stopifnot(!is.null(which(unlist(lapply(pp, anyNA)) > 0)))



    # combine the observed and pp
    combine_dat <- function(x){
      evtt_rep <- c(d$evtt[d$trt[1:n_obs] == trt_status][1:n_obs_grp], pp[[x]])
      return(evtt_rep)
    }
    evtt_rep <- lapply(1:length(pp), combine_dat)



    flog.info("simulation %s, look %s n_obs_grp = %s", idx, look, n_obs_grp)
    flog.info("simulation %s, look %s interim_month = %s", idx, look,
      interim_month)
    # flog.info("simulation %s, look %s accrualtimes = %s", idx, look,
    #   paste0(accrualtimes, collapse = ", "))


    stopifnot(n_obs_grp > 0)


    return(list(n_obs_grp = n_obs_grp,
                n_impute = n_impute,
                n_max_grp = n_obs_grp + n_impute,
                n_uncen = n_uncen,
                sum_obst = sum_obst,
                lambda = lambda,
                evtt_rep = evtt_rep))
  }



  hypoth_interim <- function(x, trt_status, dat){

    accrt_x <- d$accrt[d$trt == trt_status][1:(dat$n_obs_grp + dat$n_impute)]
    evtt_rep_x <- dat$evtt_rep[[x]]
    cen_x <- ifelse(evtt_rep_x +
                      d$age_months[d$trt == trt_status][1:(dat$n_obs_grp + dat$n_impute)] >
                      cfg$max_age_fu_months , 1, 0)

    # if the evtt + time of rand > month of interim then censor because the event
    # happens some after this interim
    cen_x <- ifelse(evtt_rep_x + accrt_x > max(cfg$interimmnths) , 1, cen_x)
    cent_x <- pmin(cfg$max_age_fu, max(cfg$interimmnths) - accrt_x)
    obst_x <- ifelse(cen_x == 1, cent_x, evtt_rep_x)


    if(length(accrt_x) != length(evtt_rep_x) ||
      length(accrt_x) != length(cen_x) ||
      length(accrt_x) != length(cent_x) ||
      length(accrt_x) != length(obst_x)){

      flog.info("hypoth_interim: dat$n_obs_grp = %s", dat$n_obs_grp)
      flog.info("hypoth_interim: dat$n_impute = %s", dat$n_impute)
      flog.info("hypoth_interim: len = %s, accrt_x = %s", length(accrt_x), paste(accrt_x, collapse = ", ") )
      flog.info("hypoth_interim: len = %s,  evtt_rep_x = %s", length(evtt_rep_x), paste(evtt_rep_x, collapse = ", ") )
      flog.info("hypoth_interim: len = %s,  cen_x = %s", length(cen_x), paste(cen_x, collapse = ", ") )
      flog.info("hypoth_interim: len = %s,  cent_x = %s", length(cent_x), paste(cent_x, collapse = ", ") )
      flog.info("hypoth_interim: len = %s,  obst_x = %s", length(obst_x), paste(obst_x, collapse = ", ") )

      saveRDS(list(type = "ERROR", where = "hypoth_interim" ,
              d=d, cfg = cfg,
              dat = dat,
              n_obs_grp = dat$n_obs_grp,
              n_impute = dat$n_impute,
          accrt_x = accrt_x,
          evtt_rep_x = evtt_rep_x,
          cen_x = cen_x,
          cent_x = cent_x,
          obst_x = obst_x), paste0("dbg/dbg-", format(Sys.time(), "%Y-%m-%d-%H-%M"), ".RDS"))
    }

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

    return(list(m = m, n_uncen_x = n_uncen_x, sum_obst_x = sum_obst_x))
  }





  # get relevant data for each arm
  # n_max and n_look are divided by 2 in the called function
  dat_ctl <- tryCatch({
    do_dat(0)
  }, error = function(err) {
    flog.info("CATCH ERROR model_clin_1 do_dat(0)  i = %s look = %s, err = %s", idx, look, err)
    saveRDS(list(type = "ERROR", d=d, cfg = cfg, err = err), paste0("dbg/dbg-", format(Sys.time(), "%Y-%m-%d-%H-%M"), ".RDS"))
    stop("Stopped on do_dat 0 error.")
  }, warning=function(cond) {
    flog.info("CATCH WARNING model_clin_1 do_dat(0) i = %s look = %s, err = %s", idx, look, cond)
    saveRDS(list(type = "WARNING", d=d, cfg = cfg, err = cond), paste0("dbg/dbg-", format(Sys.time(), "%Y-%m-%d-%H-%M"), ".RDS"))
    stop("Stopped on do_dat 0 warning")
  })



  dat_ctl <- tryCatch({
    dat_ctl$dat_interims <- lapply(1:length(dat_ctl$evtt_rep), hypoth_interim, trt_status = 0, dat = dat_ctl)
    dat_ctl
  }, error = function(err) {
    flog.info("CATCH ERROR model_clin_1 hypoth_interim(0)  i = %s look = %s, err = %s", idx, look, err)
    saveRDS(list(type = "ERROR", d=d, dat_ctl = dat_ctl, cfg = cfg, err = err), paste0("dbg/dbg-", format(Sys.time(), "%Y-%m-%d-%H-%M"), ".RDS"))
    stop("Stopped on hypoth_interim 0 error.")

  }, warning=function(cond) {
    flog.info("CATCH WARNING model_clin_1 hypoth_interim(0) i = %s look = %s, err = %s", idx, look, cond)
    saveRDS(list(type = "WARNING", d=d, dat_ctl = dat_ctl, cfg = cfg, err = cond), paste0("dbg/dbg-", format(Sys.time(), "%Y-%m-%d-%H-%M"), ".RDS"))
    stop("Stopped on hypoth_interim 0 error.")
  })





  dat_trt <- tryCatch({
    do_dat(1)
  }, error = function(err) {
    flog.info("CATCH ERROR model_clin_1 do_dat(1)  i = %s look = %s, err = %s", idx, look, err)
    saveRDS(list(type = "ERROR", d=d, cfg = cfg, err = err), paste0("dbg/dbg-", format(Sys.time(), "%Y-%m-%d-%H-%M"), ".RDS"))
    stop("Stopped on do_dat 1 error.")
  }, warning=function(cond) {
    flog.info("CATCH WARNING model_clin_1 do_dat(1) i = %s look = %s, err = %s", idx, look, cond)
    saveRDS(list(type = "WARNING", d=d, cfg = cfg, err = cond), paste0("dbg/dbg-", format(Sys.time(), "%Y-%m-%d-%H-%M"), ".RDS"))
    stop("Stopped on do_dat 1 warning")
  })



  dat_trt <- tryCatch({
    dat_trt$dat_interims <- lapply(1:length(dat_trt$evtt_rep), hypoth_interim, trt_status = 1, dat = dat_trt)
    dat_trt
  }, error = function(err) {
    flog.info("CATCH ERROR model_clin_1 hypoth_interim(1) i = %s look = %s, err = %s", idx, look, err)
    saveRDS(list(type = "ERROR", d=d, dat_trt = dat_trt, cfg = cfg, err = err), paste0("dbg/dbg-", format(Sys.time(), "%Y-%m-%d-%H-%M"), ".RDS"))
    stop("Stopped on hypoth_interim 1 error.")

  }, warning=function(cond) {
    flog.info("CATCH WARNING model_clin_1 hypoth_interim(1) i = %s look = %s, err = %s", idx, look, cond)
    saveRDS(list(type = "WARNING", d=d, dat_trt = dat_trt, cfg = cfg, err = cond), paste0("dbg/dbg-", format(Sys.time(), "%Y-%m-%d-%H-%M"), ".RDS"))
    stop("Stopped on hypoth_interim 1 warning")
  })


  # ratio of medians from the posterior estimates for lambda
  # note the ratio is intentionally ctl/trt due to the relationship between
  # lambda and median surv time
  # if there is an effect 'at the 97.5%' level then ppos_n
  # will be > 0.975
  ppos_n <- mean((dat_ctl$lambda / dat_trt$lambda)>1 , na.rm = T)

  # inference at the interim is based on the ratio of lambda_trt / lambda_ctl.
  # both distributed are gammas with integer shape parameter therefore we can transform
  # the variables to scaled chisqu and then the ratio is distributed as F (see page 316 bayesian ideas).
  # but we can also do this empircally by sampling directly to get a relative median because
  # median_1 = log(2)/lamda_1
  # median_2 = log(2)/lamda_2
  # median_1/median_2 = (log(2)/lambda_1)/(log(2)/lamda_2)
  #                   = lambda_2 / lambda_1
  # i.e. just remember to compute the right way around.

  # compute difference distribution
  diff_med <- function(x){

    # look at the complete trial for futility
    n_uncen_trt = dat_trt$dat_interims[[x]]$n_uncen_x
    sum_obst_trt = dat_trt$dat_interims[[x]]$sum_obst_x

    n_uncen_ctl = dat_ctl$dat_interims[[x]]$n_uncen_x
    sum_obst_ctl = dat_ctl$dat_interims[[x]]$sum_obst_x

    # futility
    lamb_trt <- rgamma(cfg$post_draw, shape = 1 + n_uncen_trt, rate = 1 + sum_obst_trt)
    lamb_ctl <- rgamma(cfg$post_draw, shape = 1 + n_uncen_ctl, rate = 1 + sum_obst_ctl)

    # the distribution of the difference can be approximate as normal
    # see earlier narrative comment
    ratio_medians <- lamb_ctl / lamb_trt

    # do we see a treatment effect?
    win <- mean(ratio_medians > 1) > cfg$post_tte_thresh
    # win

    return(win)
  }


  res <- tryCatch({
    unlist(lapply(1:cfg$post_draw, diff_med))
  }, error = function(err) {
    flog.info("CATCH ERROR model_clin_1 diff_med i = %s look = %s, err = %s", idx, look, err)
    saveRDS(list(type = "ERROR", d=d, cfg = cfg, err = err), paste0("dbg/dbg-", format(Sys.time(), "%Y-%m-%d-%H-%M"), ".RDS"))
    stop("Stopped on diff_med error.")
  }, warning=function(cond) {
    flog.info("CATCH WARNING model_clin_1 diff_med i = %s look = %s, err = %s", idx, look, cond)
    saveRDS(list(type = "WARNING", d=d, cfg = cfg, err = cond), paste0("dbg/dbg-", format(Sys.time(), "%Y-%m-%d-%H-%M"), ".RDS"))
    stop("Stopped on diff_med warning")
  })

  # predictive probability for:
  # superiority (at current look) and
  # futility (assuming trial expends all resources)
  ppos_max <- mean(res, na.rm = T)

  p <- tryCatch({
    
    # Wald isn't great but it is quick, the alternative is to compute a density 
    # kernel and then obtain the quantiles from that.
    delta <- dat_ctl$lambda / dat_trt$lambda
    delta_mean <- mean(delta)
    delta_wald_lwr_95 <- delta_mean + qnorm(0.025) * sd(delta)
    delta_wald_upr_95 <- delta_mean + qnorm(0.975) * sd(delta)
    
    p <- c(ppos_n, ppos_max, mean(dat_ctl$lambda), mean(dat_trt$lambda), 
           delta_mean, delta_wald_lwr_95, delta_wald_upr_95)
    names(p) <- cfg$clin_rtn_names
    p
  }, error = function(err) {
    flog.info("CATCH ERROR model_clin_1 creating return value i = %s look = %s, err = %s", idx, look, err)
    saveRDS(list(type = "ERROR", d=d, cfg = cfg, err = err), paste0("dbg/dbg-", format(Sys.time(), "%Y-%m-%d-%H-%M"), ".RDS"))
    stop("Stopped on creating return value error.")
  }, warning=function(cond) {
    flog.info("CATCH WARNING model_clin_1 creating return value i = %s look = %s, err = %s", idx, look, cond)
    saveRDS(list(type = "WARNING", d=d, cfg = cfg, err = cond), paste0("dbg/dbg-", format(Sys.time(), "%Y-%m-%d-%H-%M"), ".RDS"))
    stop("Stopped on creating return value warning")
  })


  flog.info("model_clin_1 current simulation %s, look %s, n_obs %s, result %s",
            idx, look, n_obs,
            paste0("(", paste(p, collapse = ", "), ")"))

  return(p)
}















