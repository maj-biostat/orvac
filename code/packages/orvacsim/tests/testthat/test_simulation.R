  library(testthat)
library(orvacsim)
library(data.table)
source("../../../../simulations/sim_01/util.R")
library(truncnorm)
library(survival)
library(boot)
library(rstan)
library(eha)


context("data generation")




test_that("dgp for sero correct", {



  cfg <- readRDS("cfg-example.RDS")

  # cfg$baselineprobsero
  # cfg$trtprobsero
  # compute_sero_delta(cfg$baselineprobsero, cfg$trtprobsero)

  d <- rcpp_dat(cfg)
  colnames(d) <- dnames

  m <- matrix(0, nrow = 1000, ncol = 2)

  for(i in 1:1000){
    d <- rcpp_dat(cfg)
    d <- d[1:cfg$nmaxsero,]

    m[i,1] <- mean(d[d[,COL_TRT] == 0, COL_SEROT3])
    m[i,2] <- mean(d[d[,COL_TRT] == 1, COL_SEROT3])

  }

  expect_equal(cfg$baselineprobsero, mean(m[,1]),  tolerance = cfg$baselineprobsero * 0.01)
  expect_equal(cfg$trtprobsero, mean(m[,2]),  tolerance = cfg$trtprobsero * 0.01)

})



test_that("retrieves correct number obs - information delay", {

  cfg <- readRDS("cfg-example.RDS")


  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(d)
  names(d2) <- dnames
  # View(d2)
  look <- 1
  info_delay <- 0
  expect_equal(rcpp_n_obs(d, look, cfg$looks, cfg$interimmnths, info_delay), cfg$looks[1])
  look <- 32
  expect_equal(rcpp_n_obs(d, look, cfg$looks, cfg$interimmnths, info_delay), max(cfg$looks))

  info_delay <- 1
  look <- 1
  expect_equal(rcpp_n_obs(d, look, cfg$looks, cfg$interimmnths, info_delay),
    (cfg$interimmnths[1]-info_delay)/cfg$months_per_person)

  # inexact match
  info_delay <- 0.5
  look <- 1
  expect_equal(rcpp_n_obs(d, look, cfg$looks, cfg$interimmnths, info_delay), 64)

  look <- 32
  expect_equal(rcpp_n_obs(d, look, cfg$looks, cfg$interimmnths, info_delay), 994)

})


test_that("correct sero counts", {

  cfg <- readRDS("cfg-example.RDS")

  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(copy(d))
  names(d2) <- dnames
  # View(d2)
  look <- 1
  info_delay <- 0
  nobs <- rcpp_n_obs(d, look, cfg$looks, cfg$interimmnths, info_delay)

  expect_equal(nobs, cfg$looks[1])

  d3 <- d2[1:nobs,]

  nseroctl <- sum(d3$serot3[d3$trt == 0])
  nserotrt <- sum(d3$serot3[d3$trt == 1])

  expect_equal(as.numeric(unlist(rcpp_lnsero(d, nobs))), c(nseroctl, nserotrt))

  # inexact match
  info_delay <- 0.5
  look <- 1
  nobs <- rcpp_n_obs(d, look, cfg$looks, cfg$interimmnths, info_delay)

  expect_equal(nobs, 64)
  d3 <- d2[1:nobs,]

  nseroctl <- sum(d3$serot3[d3$trt == 0])
  nserotrt <- sum(d3$serot3[d3$trt == 1])

  expect_equal(as.numeric(unlist(rcpp_lnsero(d, nobs))), c(nseroctl, nserotrt))


})


test_that("dgp for tte correct", {

  # note that cpp version of rexp uses the scale parameterisation for
  # the exponential distribution

  # ie
  # f = (1/b) exp (- x/b) rather than f = (lamb) exp (-lamb x)

  cfg <- readRDS("cfg-example.RDS")

  look <- 10

  lamb0 <- log(2)/cfg$ctl_med_tte
  lamb1 <- log(2)/cfg$trt_med_tte

  nsim <- 1000
  m <- matrix(0, nrow = nsim, ncol = 2)
  v <- matrix(0, nrow = nsim, ncol = 2)

  for(i in 1:nsim){
    d <- rcpp_dat(cfg)

    m[i,1] <- median(d[d[,COL_TRT] == 0, COL_EVTT])
    m[i,2] <- median(d[d[,COL_TRT] == 1, COL_EVTT])

    v[i,1] <- var(d[d[,COL_TRT] == 0, COL_EVTT])
    v[i,2] <- var(d[d[,COL_TRT] == 1, COL_EVTT])
  }

  expect_equal(cfg$ctl_med_tte,  median(m[,1]), tolerance = 0.3)
  expect_equal(cfg$trt_med_tte,  median(m[,2]), tolerance = 0.3)


  expect_equal(1/(lamb0^2), mean(v[,1]),  tolerance = 10)
  expect_equal(1/(lamb1^2), mean(v[,2]),  tolerance = 10)

})


test_that("all data big and small", {

  look <- 1

  cfg <- readRDS("cfg-example.RDS")
  d <- rcpp_dat(cfg)

  # assignment is otherwise by reference.
  d2 <- as.data.frame(copy(d))
  colnames(d2) <- dnames

  # this updates the age, evtt, fu1, fu2, cen and obst from observations
  # that come after the current look. rcpp_dat_small is used in generating
  # simulated datasets based on the posterior estimates for lamb0 and lamb1
  rcpp_dat_small(d, cfg, look, l0 = cfg$b0tte, l1 = cfg$b0tte + cfg$b1tte)

  d3 <- as.data.frame(copy(d))
  colnames(d3) <- dnames

  cfg$looks[look]

  # the records should stay the same up to cfg$looks[look]
  expect_equal(d2$evtt[1:cfg$looks[look]], d3$evtt[1:cfg$looks[look]], tolerance = 0.1)

  startidx <- cfg$looks[look]+1
  endidx <- max(cfg$looks)
  v2 <- d2$evtt[startidx:endidx]
  v3 <- d3$evtt[startidx:endidx]

  # and the event times after  cfg$looks[look] should all be different
  expect_false(isTRUE(all.equal(v2, v3)))

  # summary statistics of the new evtt should be similar to old
  nsim <- 1000
  mymeds <- matrix(0, nrow = nsim, ncol = 2)
  myvars <- matrix(0, nrow = nsim, ncol = 2)
  for(i in 1:nsim){
    d <- rcpp_dat(cfg)
    d2 <- as.data.frame(copy(d))
    colnames(d2) <- dnames
    rcpp_dat_small(d, cfg, look, l0 = cfg$b0tte, l1 = cfg$b0tte + cfg$b1tte)
    d3 <- as.data.frame(copy(d))
    colnames(d3) <- dnames
    startidx <- cfg$looks[look]+1
    endidx <- max(cfg$looks)
    v2 <- d2$evtt[startidx:endidx]
    v3 <- d3$evtt[startidx:endidx]
    mymeds[i, ] <- c(median(v2), median(v3))
    myvars[i, ] <- c(var(v2), var(v3))
  }

  # hist(mymeds[,1]-mymeds[,2])
  # hist(myvars[,1]-myvars[,2])

  expect_equal(mean(mymeds[,1] - mymeds[,2]), 0, tolerance = 1)
  expect_equal(mean(myvars[,1]-myvars[,2]), 0, tolerance = 20)

})

test_that("visit times", {

  cfg <- readRDS("cfg-example.RDS")

  set.seed(4343)
  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(d[1:5,])
  names(d2) <- dnames

  # test assuming at first interim analysis (look = 1)

  d2[1, "accrt"] = 0.2
  d2[1, "age"] = 6
  d2[1, "fu1"] = 0.5
  d2[1, "fu2"] = 2
  look <- 1
  idxcpp <- 0

  cfg$interimmnths[look]
  visits <- rcpp_visits(as.matrix(d2), idxcpp, look, cfg)

  # correct length
  expect_equal(length(visits), 3)
  # correct values
  expect_equal(visits[1], d2[1, "accrt"] + d2[1, "fu1"])
  expect_equal(visits[2], d2[1, "accrt"] + d2[1, "fu2"])
  expect_lt(max(visits), cfg$interimmnths[1])

  # test assuming at last interim analysis (look = 32)

  d2[1, "accrt"] = 0.2
  d2[1, "age"] = 6
  d2[1, "fu1"] = 0.5
  d2[1, "fu2"] = 2
  look <- 32

  cfg$interimmnths[look]
  visits2 <- rcpp_visits(as.matrix(d2), idxcpp, look, cfg)

  expect_equal(visits2[1], d2[1, "accrt"] + d2[1, "fu1"])
  expect_equal(visits2[2], d2[1, "accrt"] + d2[1, "fu2"])
  expect_lt(max(visits), cfg$max_age_fu_months)

  # test assuming at mid interim analysis (look = 7)

  d2[1, "accrt"] = 0.2
  d2[1, "age"] = 6
  d2[1, "fu1"] = 0.5
  d2[1, "fu2"] = 2
  look <- 7

  cfg$interimmnths[look]
  visits <- rcpp_visits(as.matrix(d2), idxcpp, look, cfg)

  # correct length - maximum visits is the number of visits
  # such that (cfg$visit_lwr * 5 ) + 6 < 36 and then add two for fu1 and fu2
  # add one to do a lt check
  expect_lt(length(visits), length(visits2))
  # correct values
  expect_equal(visits[1], d2[1, "accrt"] + d2[1, "fu1"])
  expect_equal(visits[2], d2[1, "accrt"] + d2[1, "fu2"])
  expect_lt(max(visits) + d2[1, "age"], cfg$max_age_fu_months + 0.0001)
  expect_lt(max(visits), cfg$interimmnths[look] + 0.0001)

  # test 4 - test zero visits returned

  d2[1, "accrt"] = 6.9
  d2[1, "age"] = 6
  d2[1, "fu1"] = 0.5
  d2[1, "fu2"] = 2
  look <- 1
  idxcpp <- 0

  cfg$interimmnths[look]
  visits <- rcpp_visits(as.matrix(d2), idxcpp, look, cfg)

  expect_equal(length(visits), 0)


  # retest at max visit - late enrollers

  set.seed(4343)
  d <- rcpp_dat(cfg)
  startidx <- nrow(d) - 5
  endidx <- nrow(d)
  idx <- endidx - startidx + 1
  d2 <- as.data.frame(d[startidx:endidx,])
  names(d2) <- dnames
  d2$age[idx] = 6
  d2[idx, "fu1"] = 0.5
  d2[idx, "fu2"] = 2
  d2[idx, "evtt"] = 3
  look <- 32


  cfg$interimmnths[look]
  visits2 <- rcpp_visits(as.matrix(d2), idx - 1, look, cfg)
  visits2

  # correct length - maximum visits is the number of visits
  # such that (cfg$visit_lwr * 5 ) + 6 < 36 and then add two for fu1 and fu2
  # add one to do a lt check
  expect_equal(length(visits2), 3)
  # correct values
  expect_equal(visits2[1], d2[idx, "accrt"] + d2[idx, "fu1"])
  expect_equal(visits2[2], d2[idx, "accrt"] + d2[idx, "fu2"])
  expect_lte(max(visits2) + d2[idx, "age"] - d2[idx, "accrt"], cfg$max_age_fu_months + 1)


})

test_that("censoring at interim", {


  # set.seed(4343)
  cfg <- readRDS("cfg-example.RDS")
  #cfg$b0tte <- log(2)/35
  cfg$b1tte <- 0
  d <- rcpp_dat(cfg)
  look <- 31
  cfg$interimmnths[look]

  for(i in 1:nrow(d)){

    # i = i + 1
    (visits <- rcpp_visits(d, i-1, look, cfg))
    (cens <- rcpp_cens_interim(d, visits, i-1, look, cfg))

    d[i, COL_CEN] <- cens$cen
    d[i, COL_OBST] <- cens$obst

  }

  expect_equal(max(d[1:cfg$looks[look], COL_OBST] + d[1:cfg$looks[look], COL_AGE]),
    cfg$max_age_fu_months, tolerance = 0.5)

  # suffstats
  # e1 <- d[1:cfg$looks[look], COL_EVTT]
  # c1 <- rep(0, length(e1))
  # e2 <- d[1:cfg$looks[look], COL_OBST]
  # c2 <- d[1:cfg$looks[look], COL_CEN]
  # nobs0 <- sum(1-c1)
  # nobs1 <- sum(1-c2)
  # tott0 <- sum(e1)
  # tott1 <- sum(e2)
  # rg0 <- rgamma(1000, shape = 0.01 + nobs0, rate = 0.01 + tott0)
  # rg1 <- rgamma(1000, shape = 0.01 + nobs1, rate = 0.01 + tott1)
  # rat <- rg0/rg1
  # hist(rat)

  # tbd could test ...
  # km1 <- survfit(Surv(e, 1-c)~g, data=df)
  # plot(km1)
  # quantile(km1)


})

test_that("censoring at final", {


  cfg <- readRDS("cfg-example.RDS")
  #cfg$b0tte <- log(2)/35
  cfg$b1tte <- 0
  d <- rcpp_dat(cfg)
  look <- length(cfg$looks)
  cfg$interimmnths[look]
  cfg$max_age_fu_months <- 36
  # cfg$visit_lwr <- 1.8
  # cfg$visit_lwr <- 2.2

  for(i in 1:nrow(d)){
    visits <- rcpp_visits(d, i-1, look, cfg)
    cens <- rcpp_cens_final(d, visits, i-1, look, cfg)
    d[i, COL_CEN] <- cens$cen
    d[i, COL_OBST] <- cens$obst
  }

  median(d[, COL_EVTT])
  median(d[, COL_OBST])

})

test_that("censoring", {

  cfg <- readRDS("cfg-example.RDS")

  set.seed(4343)
  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(d[1:5,])
  names(d2) <- dnames

  # test -  not censored gives correct indicator and tte

  d2[1, "accrt"] = 0.2
  d2[1, "age"] = 6
  d2[1, "evtt"] = 4
  d2[1, "fu1"] = 0.5
  d2[1, "fu2"] = 2
  look <- 1
  idxcpp <- 0

  cfg$interimmnths[look]
  visits <- rcpp_visits(as.matrix(d2), idxcpp, look, cfg)
  cens <- rcpp_cens(as.matrix(d2), visits, idxcpp, look, 0, cfg)

  expect_equal(cens$cen, 0)
  expect_equal(cens$obst, d2[1, "evtt"] , tolerance = 0.001)
  expect_lt(cens$obst, cfg$max_age_fu_months)


  # test - censored because too old event though event occurred before max visit

  d2[1, "accrt"] = 0.2
  d2[1, "age"] = 10
  d2[1, "evtt"] = 28
  d2[1, "fu1"] = 0.5
  d2[1, "fu2"] = 2
  look <- 10
  idxcpp <- 0

  cfg$interimmnths[look]
  visits <- rcpp_visits(as.matrix(d2), idxcpp, look, cfg)
  visits

  visits[length(visits)] <- d2[1, "evtt"] + d2[1, "accrt"] + 0.5

  cens <- rcpp_cens(as.matrix(d2), visits, idxcpp, look, 0, cfg)

  expect_equal(cens$cen, 1)
  expect_lte(cens$obst, cfg$max_age_fu_months)
  # i contrived the age so do not test

  # test - no visits yet

  d2[1, "accrt"] = 6.8
  d2[1, "age"] = 6
  d2[1, "evtt"] = 8
  d2[1, "fu1"] = 0.5
  d2[1, "fu2"] = 2
  look <- 1
  idxcpp <- 0

  cfg$interimmnths[look]
  (visits <- rcpp_visits(as.matrix(d2), idxcpp, look, cfg))

  (cens <- rcpp_cens(as.matrix(d2), visits, idxcpp, look, 0, cfg))

  expect_equal(cens$cen, 1)
  expect_equal(cens$obst, cfg$interimmnths[look] - d2[1, "accrt"], tolerance = 0.001)
  expect_lt(cens$obst, cfg$max_age_fu_months)

  # test - no visits yet and hasn't even enrolled

  d2[1, "accrt"] = 7.2
  d2[1, "age"] = 6
  d2[1, "evtt"] = 8
  d2[1, "fu1"] = 0.5
  d2[1, "fu2"] = 2
  look <- 1
  idxcpp <- 0

  cfg$interimmnths[look]
  (visits <- rcpp_visits(as.matrix(d2), idxcpp, look, cfg))

  (cens <- rcpp_cens(as.matrix(d2), visits, idxcpp, look, 0, cfg))

  expect_false(!is.na(cens$cen))
  expect_false(!is.na(cens$obst))

  # test - censored event is after last visit

  d2[1, "accrt"] = 98
  d2[1, "age"] = 6
  d2[1, "evtt"] = 200
  d2[1, "fu1"] = 0.5
  d2[1, "fu2"] = 2
  look <- length(cfg$looks)
  idxcpp <- 0

  cfg$interimmnths[look]
  (visits <- rcpp_visits(as.matrix(d2), idxcpp, look, cfg))

  cens <- rcpp_cens(as.matrix(d2), visits, idxcpp, look, 0, cfg)
  cens

  expect_equal(cens$cen, 1)
  expect_lt(cens$obst, cfg$max_age_fu_months)


  # test -  censored because whlie event is b4 max visit, age is already > 36 at max visit

  d2[1, "accrt"] = 33
  d2[1, "age"] = 36
  d2[1, "evtt"] = 200
  d2[1, "fu1"] = 0.5
  d2[1, "fu2"] = 2
  look <- 10
  idxcpp <- 0

  cfg$interimmnths[look]
  (visits <- rcpp_visits(as.matrix(d2), idxcpp, look, cfg))

  cens <- rcpp_cens(as.matrix(d2), visits, idxcpp, look, 0, cfg)
  cens

  expect_equal(cens$cen, 1)
  expect_equal(cfg$max_age_fu_months - d2[1, "age"], cens$obst, tolerance = 0.001)


})

test_that("setting obst and censor status", {

  cfg <- readRDS("cfg-example.RDS")

  # test - are the correct number of obs times and censor status set at first interim

  set.seed(4343)
  look <- 1
  cfg$interimmnths[look]
  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(copy(d))
  colnames(d2) <- dnames

  rcpp_clin_set_obst(d, cfg, look, 0)

  d3 <- as.data.frame(copy(d))
  colnames(d3) <- dnames

  sum(!is.na(d3$obst))
  cfg$looks[look]

  expect_equal(sum(!is.na(d3$obst)), cfg$looks[look], tolerance = 0.01)
  expect_equal(sum(!is.na(d3$cen)), cfg$looks[look], tolerance = 0.01)

  # test - are the correct number of obs times and censor status set at last interim

  look <- length(cfg$looks)
  cfg$interimmnths[look]
  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(copy(d))
  colnames(d2) <- dnames

  rcpp_clin_set_obst(d, cfg, look, 0)

  d3 <- as.data.frame(copy(d))
  colnames(d3) <- dnames

  sum(!is.na(d3$obst)) ;   cfg$looks[look]

  expect_equal(sum(!is.na(d3$obst)), cfg$looks[look], tolerance = 0.01)
  expect_equal(sum(!is.na(d3$cen)), cfg$looks[look], tolerance = 0.01)

  # test - are the correct number of obs times and censor status set at mid interim

  look <- 12
  cfg$interimmnths[look]
  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(copy(d))
  colnames(d2) <- dnames

  rcpp_clin_set_obst(d, cfg, look, 0)

  d3 <- as.data.frame(copy(d))
  colnames(d3) <- dnames

  sum(!is.na(d3$obst)) ;   cfg$looks[look]

  expect_equal(sum(!is.na(d3$obst)), cfg$looks[look], tolerance = 0.01)
  expect_equal(sum(!is.na(d3$cen)), cfg$looks[look], tolerance = 0.01)




})

context("clinical endpoint framework")

test_that("clinical endpoint posterior", {

  cfg <- readRDS("cfg-example.RDS")


  # test - estimated median tte is biased in both arms due to censoring mechanism
  # but ratio is preserved with lower than 2% bias
  look <- 32
  nsim <- 1000
  v <- numeric(nsim)

  for(i in 1:nsim){
    # i = i + 1
    d <- rcpp_dat(cfg)
    d2 <- as.data.frame(copy(d))
    colnames(d2) <- dnames

    # get sufficieitn stats
    (lsuffstat1 <- rcpp_clin_set_obst(d, cfg, look, 1))

    d2 <- as.data.frame(copy(d))
    colnames(d2) <- dnames

    #lsuffstat2 <- rcpp_clin_set_obst(d, cfg, look)

    # obtain posterior based on current look
    m <- matrix(0, nrow = cfg$post_draw, ncol = 3)
    rcpp_clin_interim_post(m,
                           lsuffstat1$n_uncen_0, lsuffstat1$tot_obst_0,
                           lsuffstat1$n_uncen_1, lsuffstat1$tot_obst_1,
                           cfg$post_draw, cfg);
   # plot_tte_hist(m)
    v[i] <- mean(m[, 1]/m[, 2])

    # should be (log(2)/(cfg$b0tte + cfg$b1tte))/ (log(2)/cfg$b0tte)
  }
  # hist(v)
  # abline(v = (log(2)/(cfg$b0tte + cfg$b1tte))/ (log(2)/cfg$b0tte), col = cbp[2], lwd = 2)
  # abline(v = mean(v), col = cbp[3], lwd = 2, lty = 3)

  expect_equal(mean(v), (log(2)/(cfg$b0tte + cfg$b1tte))/ (log(2)/cfg$b0tte), tolerance = 0.02)



})


test_that("clin tte posterior check against mcmc", {


  cfg <- readRDS("cfg-example.RDS")

  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(copy(d))
  colnames(d2) <- dnames
  look <- 32
  lsuffstat <- rcpp_clin_set_obst(d, cfg, look, 1)
  d2 <- as.data.frame(copy(d)); colnames(d2) <- dnames
  # plot_tte_hist_dat(d2$obst, d2$trt, cfg$looks[look])

  m <- matrix(0, nrow = cfg$post_draw, ncol = 3)
  rcpp_clin_interim_post(m,
                         lsuffstat$n_uncen_0, lsuffstat$tot_obst_0,
                         lsuffstat$n_uncen_1, lsuffstat$tot_obst_1,
                         cfg$post_draw, cfg);
  # plot_tte_hist(m)

  sdat <- list(
    ## Number of event individuals
    Nobs = sum(d2$cen == 0),
    ## Number of censored individuals
    Ncen = sum(d2$cen == 1),
    ## Times for event individuals
    yobs = d2$obst[d2$cen == 0],
    ## Times for censored individuals
    ycen = d2$obst[d2$cen == 1],
    ## Covariates for event individuals as a matrix
    Xobs_bg = d2$trt[d2$cen == 0],
    ## Covariates for censored individuals as a matrix
    Xcen_bg = d2$trt[d2$cen == 1],
    alpha0 = 1,
    beta0 = 30
  )

  s2 <- rstan::stan(file = "expo_2.stan",
                    data = sdat, chains = 1)

  beta <- as.matrix(s2)
  beta[,2] <- beta[,1] + beta[,2]
  beta[,3] <- beta[,1]/beta[,2]
  # plot_tte_hist(beta)


  expect_equal(mean(m[,1]), mean(beta[,1]), tolerance = 0.02 * 0.05)
  expect_equal(mean(m[,2]), mean(beta[,2]), tolerance = 0.02 * 0.05)
  expect_equal(mean(m[,3]), mean(beta[,3]), tolerance = 1.2 * 0.05)

})

test_that("clin tte data ppos", {


  cfg <- readRDS("cfg-example.RDS")

  # no effect
  cfg$b1tte <- 0
  d <- rcpp_dat(cfg)
  look <- 10

  m_clin_res <- rcpp_clin_opt(d, cfg, look)

  # ppn is just the post prob that ratio is > 1
  expect_lt(m_clin_res$ppn, cfg$post_tte_sup_thresh[look])
  expect_lt(m_clin_res$lwr, 1)
  expect_gt(m_clin_res$upr, 1)

  # expect_lt(m_clin_res$ppmax, cfg$pp_tte_fut_thresh)

  # effect present

  cfg <- readRDS("cfg-example.RDS")

  cfg$b0tte
  cfg$b1tte <- -0.006
  d <- rcpp_dat(cfg)
  look <- 20

  m_clin_res <- rcpp_clin_opt(d, cfg, look)

  # ppn is just the post prob that ratio is > 1
  expect_gt(m_clin_res$ppn, cfg$post_tte_sup_thresh[look])
  expect_gt(m_clin_res$lwr, 1)
  expect_gt(m_clin_res$upr, 1)






})


context("immu endpoint framework")


test_that("immu model - estimated post and diff", {


  cfg <- readRDS("cfg-example.RDS")
  look = 7

  # look at average results.
  res <- matrix(0, nrow = 1000, ncol = 3)

  for(i in 1:1000){
    d <- rcpp_dat(cfg)

    info_delay <- 0
    look <- 32
    nobs <- rcpp_n_obs(d, look, cfg$looks, cfg$interimmnths, info_delay)
    lnsero <- rcpp_lnsero(d, nobs)

    # call this instead.
    m <- matrix(0, nrow = cfg$post_draw, ncol = 3)
    rcpp_immu_interim_post(d, m, nobs, cfg$post_draw, lnsero)

    res[i,1] <- mean(m[, COL_THETA0 ])
    res[i,2] <- mean(m[, COL_THETA1 ])
    res[i,3] <- mean(m[, COL_DELTA ])
  }


  # got to be within 1%
  expect_equal(cfg$baselineprobsero, mean(res[,1]),  cfg$baselineprobsero * 1/100)
  expect_equal(cfg$trtprobsero, mean(res[,2]),  cfg$trtprobsero * 1/100)

  # got to be within 1%
  diff_true <- cfg$trtprobsero - cfg$baselineprobsero
  diff_mean_est <- mean(res[,2] - res[,1])
  expect_equal(diff_true, diff_mean_est,  diff_true * 1/100)

})


test_that("dotrial", {

  doglm <- function(df){
    pci <- function(data, indices) {
      d <- data[indices,] # allows boot to select sample
      fit <- glm(serot3 ~ trt, data = d, family = binomial)
      p <- predict(fit, type = "response", newdata = data.frame(trt = 0:1))
      return(p[2]-p[1])
    }
    res <- boot(data=df, statistic=pci, R=1000)
    return(res)
  }

  doexp <- function(df){
    rci <- function(data, indices) {
      d <- data[indices,] # allows boot to select sample
      # aft - in aft (expon), +ve param est imply increasing surv time and longer exp durn
      fit1 <- survreg(Surv(obst, cen == 0) ~ trt, data = d, dist = "exponential")
      summary(fit1)
      r <- exp(coef(fit1)[2])
      # # ph
      # fit1 <- phreg(Surv(obst, cen == 0) ~ trt, data = df, dist = "weibull", shape = 1)
      # summary(fit1)
      return(r)
    }
    res <- boot(data=df, statistic=rci, R=1000)
    return(res)
  }



  cfg <- readRDS("cfg-example.RDS")

  cfg$baselineprobsero
  cfg$trtprobsero <- 0.55
  cfg$deltaserot3 <- compute_sero_delta(cfg$baselineprobsero, cfg$trtprobsero)

  l <- rcpp_dotrial(1, cfg, TRUE)

  df <- as.data.frame(l$d)
  names(df) <- dnames
  b1 <- boot.ci(doglm(df[1:cfg$nmaxsero, ]), type="norm")

  expect_equal(b1$normal[2], l$i_lwr, tolerance = 0.2)
  expect_equal(b1$normal[3], l$i_upr, tolerance = 0.2)


  b1 <- boot.ci(doexp(df), type="norm")

  expect_equal(b1$normal[2], l$c_lwr, tolerance = 0.2)
  expect_equal(b1$normal[3], l$c_upr, tolerance = 0.2)



})




