library(testthat)
library(orvacsim)
library(data.table)




context("data generation")


test_that("use pass by reference if possible", {
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  
  mydim <- 1000
  nsim <- 1000
  start <- proc.time()
  for(i in 1:nsim){
    m1 <- matrix(0, ncol = mydim, nrow = mydim)
    rcpp_test_1(m1)
  }
  end <- proc.time()
  duration1 <- end - start
  
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  mydim <- 1000
  nsim <- 1000
  start <- proc.time()
  for(i in 1:nsim){
    m1 <- matrix(0, ncol = mydim, nrow = mydim)
    m2 <- rcpp_test_2(m1)
  }
  end <- proc.time()
  duration2 <- end - start
  
  
  summary(duration1)
  summary(duration2)
  
  expect_lt(duration1[3],  duration2[3])
  # By reference is MUCH quicker!!!
  #
  # >   summary(duration1)
  # user  system elapsed 
  # 2.556   0.681   3.239 
  # >   summary(duration2)
  # user  system elapsed 
  # 13.974   1.792  15.791 
  
})


test_that("cpp quicker than data.table", {
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  source("sim.R")
  library(truncnorm)
  
  d <- readRDS("tests/dat-example.RDS")
  cfg <- readRDS("tests/cfg-example.RDS")
  
  n_obs_grp <- 170
  look <- 10
  
  expect_equal(cfg$looks, c(70L, 100L, 130L, 160L, 190L, 220L, 250L, 280L, 310L, 340L, 
                            370L, 400L, 430L, 460L, 490L, 520L, 550L, 580L, 610L, 640L, 670L, 
                            700L, 730L, 760L, 790L, 820L, 850L, 880L, 910L, 940L, 970L, 1000L))
  
  expect_equal(cfg$interimmnths, c(7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49, 
                                   52, 55, 58, 61, 64, 67, 70, 73, 76, 79, 82, 85, 88, 91, 94, 97, 
                                   100))
  
  expect_equal(cfg$interimmnths[look], 34)
  expect_equal(cfg$months_per_person, 0.1)
  
  # cfg$baselineprobsero
  # cfg$trtprobsero
  # compute_sero_delta(cfg$baselineprobsero, cfg$trtprobsero)
  start <- proc.time()
  for(i in 1:1000){
    
    d <- gen_dat(cfg)
    
  }
  end <- proc.time()
  duration1 <- end - start
  summary(duration1)
  
  start <- proc.time()
  for(i in 1:1000){
    
    d <- rcpp_dat(cfg)
    
  }
  end <- proc.time()
  duration2 <- end - start
  summary(duration2)
  
  expect_lt(duration2[3],  duration1[3])
  
  
})


test_that("all datas big and small", {
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  source("sim.R")

  look <- 1
  
  cfg <- readRDS("tests/cfg-example.RDS")
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
  nsim <- 100
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
    v2 <- d2$evtt[startidx:endidx] + d2$age[startidx:endidx]
    v3 <- d3$evtt[startidx:endidx] + d3$age[startidx:endidx]
    mymeds[i, ] <- c(median(v2), median(v3))
    myvars[i, ] <- c(var(v2), var(v3))
  }
  
  hist(mymeds[,1]-mymeds[,2])
  hist(myvars[,1]-myvars[,2])
  
  expect_equal(mean(mymeds[,1] - mymeds[,2]), 0, tolerance = 3)
  expect_equal(mean(myvars[,1]-myvars[,2]), 0, tolerance = 20)

})


test_that("dgp for sero correct", {
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  
  d <- readRDS("tests/dat-example.RDS")
  cfg <- readRDS("tests/cfg-example.RDS")
  
  n_obs_grp <- 170
  look <- 10
  

  expect_equal(cfg$looks, c(70L, 100L, 130L, 160L, 190L, 220L, 250L, 280L, 310L, 340L, 
                            370L, 400L, 430L, 460L, 490L, 520L, 550L, 580L, 610L, 640L, 670L, 
                            700L, 730L, 760L, 790L, 820L, 850L, 880L, 910L, 940L, 970L, 1000L))
  
  expect_equal(cfg$interimmnths, c(7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49, 
                                   52, 55, 58, 61, 64, 67, 70, 73, 76, 79, 82, 85, 88, 91, 94, 97, 
                                   100))
  
  expect_equal(cfg$interimmnths[look], 34)
  expect_equal(cfg$months_per_person, 0.1)
  
  # cfg$baselineprobsero
  # cfg$trtprobsero
  # compute_sero_delta(cfg$baselineprobsero, cfg$trtprobsero)
  
  m <- matrix(0, nrow = 1000, ncol = 2)
  
  for(i in 1:1000){
    d <- rcpp_dat(cfg)
    d <- d[1:cfg$nmaxsero,]
    
    m[i,1] <- mean(d[d[,COL_TRT] == 0, COL_SEROT3])
    m[i,2] <- mean(d[d[,COL_TRT] == 1, COL_SEROT3])
    
  }
  
  expect_equal(abs(cfg$baselineprobsero - mean(m[,1])), 0,  tolerance = cfg$baselineprobsero * 0.01)
  expect_equal(abs(cfg$trtprobsero - mean(m[,2])), 0,  tolerance = cfg$trtprobsero * 0.01)
  
  

  # tests 
  # - repeat sampling of stochastic elements with summary stats
  # balance for tretament aloc
  # linear accrual
  
  
  
})


test_that("dgp for tte correct", {
  
  # note that cpp version of rexp uses the scale parameterisation for 
  # the exponential distribution
  
  # ie
  # f = (1/b) exp (- x/b) rather than f = (lamb) exp (-lamb x)
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  
  d <- readRDS("tests/dat-example.RDS")
  cfg <- readRDS("tests/cfg-example.RDS")
  
  n_obs_grp <- 170
  look <- 10
  
  
  stopifnot(cfg$looks == c(70L, 100L, 130L, 160L, 190L, 220L, 250L, 280L, 310L, 340L, 
                            370L, 400L, 430L, 460L, 490L, 520L, 550L, 580L, 610L, 640L, 670L, 
                            700L, 730L, 760L, 790L, 820L, 850L, 880L, 910L, 940L, 970L, 1000L))
  
  stopifnot(cfg$interimmnths == c(7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49, 
                                   52, 55, 58, 61, 64, 67, 70, 73, 76, 79, 82, 85, 88, 91, 94, 97, 
                                   100))
  
  stopifnot(cfg$interimmnths[look] == 34)
  stopifnot(cfg$months_per_person ==  0.1)
  
  stopifnot(cfg$ctl_med_tte == 30)
  stopifnot(cfg$trt_med_tte == 35)
  
  
  lamb0 <- log(2)/cfg$ctl_med_tte
  lamb1 <- log(2)/cfg$trt_med_tte
  
  nsim <- 1000
  m <- matrix(0, nrow = nsim, ncol = 2)
  v <- matrix(0, nrow = nsim, ncol = 2)
  
  for(i in 1:nsim){
    d <- rcpp_dat(cfg)
    
    m[i,1] <- median(d[d[,COL_TRT] == 0, COL_EVTT] + d[d[,COL_TRT] == 0, COL_AGE])
    m[i,2] <- median(d[d[,COL_TRT] == 1, COL_EVTT] + d[d[,COL_TRT] == 1, COL_AGE])
    
    v[i,1] <- var(d[d[,COL_TRT] == 0, COL_EVTT] + d[d[,COL_TRT] == 0, COL_AGE])
    v[i,2] <- var(d[d[,COL_TRT] == 1, COL_EVTT] + d[d[,COL_TRT] == 1, COL_AGE])
  }
  
  expect_equal(cfg$ctl_med_tte,  median(m[,1]), tolerance = 0.3)
  expect_equal(cfg$trt_med_tte,  median(m[,2]), tolerance = 0.3)
  
  
  expect_equal(1/(lamb0^2), mean(v[,1]),  tolerance = 10)
  expect_equal(1/(lamb1^2), mean(v[,2]),  tolerance = 10)
  
  
  # tests 
  # - repeat sampling of stochastic elements with summary stats
  # balance for tretament aloc
  # linear accrual
  
  
  
})


test_that("retrieves correct number obs", {
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  source("sim.R")
  library(truncnorm)
  
  d <- readRDS("tests/dat-example.RDS")
  cfg <- readRDS("tests/cfg-example.RDS")
  
  n_obs_grp <- 170
  look <- 10
  
  expect_equal(cfg$looks, c(70L, 100L, 130L, 160L, 190L, 220L, 250L, 280L, 310L, 340L, 
                            370L, 400L, 430L, 460L, 490L, 520L, 550L, 580L, 610L, 640L, 670L, 
                            700L, 730L, 760L, 790L, 820L, 850L, 880L, 910L, 940L, 970L, 1000L))
  
  expect_equal(cfg$interimmnths, c(7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49, 
                                   52, 55, 58, 61, 64, 67, 70, 73, 76, 79, 82, 85, 88, 91, 94, 97, 
                                   100))
  
  expect_equal(cfg$interimmnths[look], 34)
  expect_equal(cfg$months_per_person, 0.1)
  

  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(d)
  names(d2) <- dnames
  # View(d2)
  look <- 1
  info_delay <- 0
  expect_equal(rcpp_n_obs(d, look, cfg$looks, cfg$interimmnths, info_delay), 70)
  look <- 32
  expect_equal(rcpp_n_obs(d, look, cfg$looks, cfg$interimmnths, info_delay), 1000)
  
  info_delay <- 1
  look <- 1
  expect_equal(rcpp_n_obs(d, look, cfg$looks, cfg$interimmnths, info_delay), 60)
  
  # inexact match
  info_delay <- 0.5
  look <- 1
  expect_equal(rcpp_n_obs(d, look, cfg$looks, cfg$interimmnths, info_delay), 64)
  
  look <- 32
  expect_equal(rcpp_n_obs(d, look, cfg$looks, cfg$interimmnths, info_delay), 994)
  
})


test_that("correct sero counts", {
  
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  source("sim.R")
  library(truncnorm)
  
  d <- readRDS("tests/dat-example.RDS")
  cfg <- readRDS("tests/cfg-example.RDS")
  
  n_obs_grp <- 170
  look <- 10
  
  expect_equal(cfg$looks, c(70L, 100L, 130L, 160L, 190L, 220L, 250L, 280L, 310L, 340L, 
                            370L, 400L, 430L, 460L, 490L, 520L, 550L, 580L, 610L, 640L, 670L, 
                            700L, 730L, 760L, 790L, 820L, 850L, 880L, 910L, 940L, 970L, 1000L))
  
  expect_equal(cfg$interimmnths, c(7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49, 
                                   52, 55, 58, 61, 64, 67, 70, 73, 76, 79, 82, 85, 88, 91, 94, 97, 
                                   100))
  
  expect_equal(cfg$interimmnths[look], 34)
  expect_equal(cfg$months_per_person, 0.1)
  
  
  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(copy(d))
  names(d2) <- dnames
  # View(d2)
  look <- 1
  info_delay <- 0
  nobs <- rcpp_n_obs(d, look, cfg$looks, cfg$interimmnths, info_delay)
  
  expect_equal(nobs, 70)
  
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




context("clinical endpoint framework")


test_that("visit times", {
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  cfg <- readRDS("tests/cfg-example.RDS")
  
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
  
  
  
  # odd visit

  # cfg <- readRDS("tests/cfg-example.RDS")
  # d <- readRDS("tests/oddvisit.RDS")
  # d2 <- as.data.frame(copy(d))
  # names(d2) <- dnames
  # look = length(cfg$looks)
  # 
  # 
  # i <- 2
  # visits <- rcpp_visits(d, i-1, look, cfg)
  
})


test_that("censoring at final", {
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  

  # set.seed(4343)
  cfg <- readRDS("tests/cfg-example.RDS")
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
    cat(paste0(visits, "\n"))
    cens <- rcpp_cens_final(d, visits, i-1, look, cfg)
    
    d[i, COL_CEN] <- cens$cen
    d[i, COL_OBST] <- cens$obst
    
  }
  median(d[, COL_EVTT])
  median(d[, COL_OBST])
  m <- cbind(d[, COL_EVTT], d[, COL_OBST], 0)
  plot_tte_meds_hist(m, 30, 30)
  prop.table(table(d[, COL_CEN]))
  
  # i remain unconvince about the feasibility of a 36 month fu when your 
  # median tte is 30
  
})


test_that("censoring at interim", {
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  
  
  # set.seed(4343)
  cfg <- readRDS("tests/cfg-example.RDS")
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
    
    d2 <- as.data.frame(d)
    names(d2) <- dnames
    d2[i,]
  }
  median(d[1:cfg$looks[look], COL_EVTT])
  median(d[1:cfg$looks[look], COL_OBST])
  m <- cbind(d[, COL_EVTT], d[, COL_OBST], 0)
  plot_tte_meds_hist(m, 30, 30)
  prop.table(table(d[, COL_CEN]))
  
  # i remain unconvince about the feasibility of a 36 month fu when your 
  # median tte is 30
  
})


test_that("censoring", {
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  cfg <- readRDS("tests/cfg-example.RDS")
  

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
  
  cens <- rcpp_cens(as.matrix(d2), visits, idxcpp, look, cfg)
  
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
  
  cens <- rcpp_cens(as.matrix(d2), visits, idxcpp, look, cfg)
  
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
  
  (cens <- rcpp_cens(as.matrix(d2), visits, idxcpp, look, cfg))
  
  expect_equal(cens$cen, 1)
  expect_equal(cens$obst, d2[1, "age"] + cfg$interimmnths[look] - d2[1, "accrt"], tolerance = 0.001)
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
  
  (cens <- rcpp_cens(as.matrix(d2), visits, idxcpp, look, cfg))
  
  expect_false(!is.na(cens$cen))
  expect_false(!is.na(cens$obst)) 
  
  
  
  # test - censored event is after last visit
  
  d2[1, "accrt"] = 98
  d2[1, "age"] = 6
  d2[1, "evtt"] = 1000
  d2[1, "fu1"] = 0.5
  d2[1, "fu2"] = 2
  look <- length(cfg$looks)
  idxcpp <- 0
  
  cfg$interimmnths[look]
  (visits <- rcpp_visits(as.matrix(d2), idxcpp, look, cfg))
  
  cens <- rcpp_cens(as.matrix(d2), visits, idxcpp, look, cfg)
  cens
  
  expect_equal(cens$cen, 1)
  expect_lt(cens$obst, cfg$max_age_fu_months)
  
  
  # test -  censored because event is after max visit but age is already > 36 at max visit
  
  d2[1, "accrt"] = 33
  d2[1, "age"] = 36
  d2[1, "evtt"] = 1000
  d2[1, "fu1"] = 0.5
  d2[1, "fu2"] = 2
  look <- 10
  idxcpp <- 0
  
  cfg$interimmnths[look]
  (visits <- rcpp_visits(as.matrix(d2), idxcpp, look, cfg))
  
  cens <- rcpp_cens(as.matrix(d2), visits, idxcpp, look, cfg)
  cens
  
  expect_equal(cens$cen, 1)
  expect_equal(cfg$max_age_fu_months, cens$obst, tolerance = 0.001)
  

  # test - sample statistics
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  cfg <- readRDS("tests/cfg-example.RDS")
  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(copy(d))
  names(d2) <- dnames
  look = length(cfg$looks)
  
  nsim <- 1000
  mest <- matrix(NA, nrow = nsim, ncol = 2)
  
  for(j in 1:nsim){
    
    d <- rcpp_dat(cfg)
    m <- matrix(NA, nrow = nrow(d), ncol = 2)
    
    for(i in 1:nrow(d)){
      
      visits <- rcpp_visits(d, i-1, look, cfg)
      cens1 <- rcpp_cens(d, visits, i-1, look, cfg)
      
      
      m[i, ] <- c(d[i, COL_EVTT] + d[i, COL_AGE], cens1$obst)
      
    }
    mest[j, ] <- apply(m, 2, median, na.rm = T)
  }
  
  
  hist(mest[,1])
  hist(mest[,2])
  
  
  
  
  
})


test_that("setting obst and censor status", {
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  cfg <- readRDS("tests/cfg-example.RDS")
  
  

  # test - are the correct number of obs times and censor status set at first interim
  
  set.seed(4343)
  look <- 1
  cfg$interimmnths[look]
  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(copy(d))
  colnames(d2) <- dnames
  
  rcpp_clin_set_obst(d, cfg, look)
  
  d3 <- as.data.frame(copy(d))
  colnames(d3) <- dnames
  
  sum(!is.na(d3$obst))
  cfg$looks[look]
  
  expect_equal(sum(!is.na(d3$obst)), cfg$looks[look], tolerance = 0.01)
  expect_equal(sum(!is.na(d3$cen)), cfg$looks[look], tolerance = 0.01)
  
  
  
  # test - 
  
  look <- length(cfg$looks)
  cfg$interimmnths[look]
  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(copy(d))
  colnames(d2) <- dnames
  
  d2$accrt <- 0
  d2$age <- truncnorm::rtruncnorm(nrow(d2), 
                                  cfg$age_months_lwr, cfg$age_months_upr,
                                  cfg$age_months_mean, cfg$age_months_sd)
  d2$evtt <- rexp(nrow(d2), cfg$b0tte) - d2$age
  cfg$max_age_fu_months <- 36
  d2$obst <- ifelse(d2$evtt + d2$age > cfg$max_age_fu_months, cfg$max_age_fu_months, d2$evtt + d2$age)
  d2$cen <- ifelse(d2$evtt + d2$age > cfg$max_age_fu_months, 1, 0)
  m <- cbind(d2$evtt + d2$age, d2$obst, 0)
  plot_tte_meds_hist(m, log(2)/cfg$b0tte, log(2)/cfg$b0tte)
  prop.table(table(d2$cen))
  
  #
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  cfg <- readRDS("tests/cfg-example.RDS")
  
  look <- length(cfg$looks)
  cfg$interimmnths[look]
  cfg$b1tte <- 0
  
  #cfg$max_age_fu_months <- 100
  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(copy(d))
  colnames(d2) <- dnames
  #head(d, 5)
  
  #set.seed(4343)
  # (visit <- rcpp_visits(d, 2, 32, cfg))
  # rcpp_cens(d, visit, 2, 32, cfg)
  
  # set.seed(4343)
  # (vis <- rcpp_visits(as.matrix(d), 3, look, cfg))
  # rcpp_cens(d, vis, 3, look, cfg)
  
  lsuf <- rcpp_clin_set_obst(d, cfg, look)
  head(d, 5)
  d2 <- as.data.frame(copy(d))
  colnames(d2) <- dnames
  head(d2, 5)
  median(d2$obst)
  prop.table(table(d2$cen))
  
  m <- cbind(d2$evtt + d2$age, d2$obst, 0)
  plot_tte_meds_hist(m, log(2)/cfg$b0tte, log(2)/cfg$b0tte)
  
  
  
  
  
  
  
  
  
  rcpp_clin_set_obst(d, cfg, look)
  
  d3 <- as.data.frame(copy(d))
  colnames(d3) <- dnames
  
  sum(!is.na(d3$obst))
  cfg$looks[look]
  
  expect_equal(sum(!is.na(d3$obst)), cfg$looks[look], tolerance = 0.01)
  expect_equal(sum(!is.na(d3$cen)), cfg$looks[look], tolerance = 0.01)
  
  nsim <- 10
  for(i in 1:nsim){
    
    
    
    
  }
  
  
  
  
  # test - are the correct number of obs times and censor status set at last interim
  
  look <- length(cfg$looks)
  cfg$interimmnths[look]
  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(copy(d))
  colnames(d2) <- dnames
  
  rcpp_clin_set_obst(d, cfg, look)
  
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
  
  rcpp_clin_set_obst(d, cfg, look)
  
  d3 <- as.data.frame(copy(d))
  colnames(d3) <- dnames
  
  sum(!is.na(d3$obst)) ;   cfg$looks[look]
  
  expect_equal(sum(!is.na(d3$obst)), cfg$looks[look], tolerance = 0.01)
  expect_equal(sum(!is.na(d3$cen)), cfg$looks[look], tolerance = 0.01)
 
  
  # sample stats for obst times (that include censor times give unbiased view)
  cfg <- readRDS("tests/cfg-example.RDS")
  look <- length(cfg$looks)
  nsim <- 10
  m <- matrix(0, ncol = 6, nrow = nsim)
  for(i in 1:nsim){
    
    d <- rcpp_dat(cfg)
    d2 <- as.data.frame(copy(d))
    colnames(d2) <- dnames
    rcpp_clin_set_obst(d, cfg, look)
    d3 <- as.data.frame(copy(d))
    colnames(d3) <- dnames
    
    d2_ctl <- d2[d2$trt == 0,]
    d2_trt <- d2[d2$trt == 1,]
    
    d3_ctl <- d3[d3$trt == 0,]
    d3_trt <- d3[d3$trt == 1,]
    
    
    v2_ctl <- d2_ctl$evtt[1:(0.5*nrow(d2))] + d2_ctl$age[1:(0.5*nrow(d2))]
    v2_trt <- d2_trt$evtt[1:(0.5*nrow(d2))] + d2_trt$age[1:(0.5*nrow(d2))]
    
    v3_ctl <- d3_ctl$obst[1:(0.5*nrow(d2))]
    v3_trt <- d3_trt$obst[1:(0.5*nrow(d2))]
    
    propcensctl = sum(d3_ctl$cen, na.rm = T) / (0.5*nrow(d2))
    propcenstrt = sum(d3_ctl$cen, na.rm = T) / (0.5*nrow(d2))
    
    m[i, ] <- c(median(v2_ctl), median(v2_trt), 
                median(v3_ctl), median(v3_trt),
                propcensctl, propcenstrt)
  }
  m
  
})


test_that("clinical endpoint posterior", {
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  cfg <- readRDS("tests/cfg-example.RDS")
  
  
  
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
    (lsuffstat1 <- rcpp_clin_set_obst(d, cfg, look))
    
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
  hist(v)
  abline(v = (log(2)/(cfg$b0tte + cfg$b1tte))/ (log(2)/cfg$b0tte), col = cbp[2], lwd = 2)
  abline(v = mean(v), col = cbp[3], lwd = 2, lty = 3)
  
  expect_equal(mean(v), (log(2)/(cfg$b0tte + cfg$b1tte))/ (log(2)/cfg$b0tte), tolerance = 0.02)
  

  
  
  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(copy(d))
  colnames(d2) <- dnames
  
  # get sufficieitn stats
  (lsuffstat1 <- rcpp_clin_set_obst(d, cfg, look))
  
  d2 <- as.data.frame(copy(d))
  colnames(d2) <- dnames
  
  #lsuffstat2 <- rcpp_clin_set_obst(d, cfg, look)
  
  # obtain posterior based on current look 
  m <- matrix(0, nrow = cfg$post_draw, ncol = 3)
  rcpp_clin_interim_post(m, 
                         lsuffstat1$n_uncen_0, lsuffstat1$tot_obst_0,
                         lsuffstat1$n_uncen_1, lsuffstat1$tot_obst_1,
                         cfg$post_draw, cfg);
  
  # colMeans(m)
  
  par(mfrow = c(2, 2))
  hist(log(2)/m[, 1], probability = T, main = "")
  abline(v = log(2)/cfg$b0tte, col = "red", lwd = 2)
  hist(log(2)/m[, 2], probability = T, main = "")
  abline(v = log(2)/(cfg$b0tte+cfg$b1tte), col = "red", lwd = 2)
  hist(m[, 1]/m[, 2], probability = T, main = "")
  abline(v = (log(2)/(cfg$b0tte + cfg$b1tte))/ (log(2)/cfg$b0tte), col = "red", lwd = 2)
  plot(c(0, 10), c(0, 10))
  legend(0, 5, legend=c("true med", "sample med"),
         col=c("red", "blue"), lty=1:2, cex=0.8)
  par(mfrow = c(1, 1))
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # test - are the correct number of obs times and censor status set at first interim
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  set.seed(4343)
  look <- 32
  cfg$interimmnths[look]
  cfg$b1tte <- 0
  
  cfg$max_age_fu_months <- 36
  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(copy(d))
  colnames(d2) <- dnames
  # about 30
  median(d2$evtt)
  
  (lsuffstat <- rcpp_clin_set_obst(d, cfg, look))
  d2 <- as.data.frame(copy(d))
  colnames(d2) <- dnames
  
  sum(1-d2$cen)
  sum(d2$obst)
  median(d2$obst)
  d2$testcen <- ifelse(d2$evtt + d2$age > 48, 1, 0)
  d2$testobst <- ifelse(d2$testcen == 1, 48, d2$evtt)
  sum(1-d2$testcen)
  median(d2$testobst)
  
  m <- matrix(0, nrow = cfg$post_draw, ncol = 3)
  rcpp_clin_interim_post(m, 
                         lsuffstat$n_uncen_0, lsuffstat$tot_obst_0,
                         lsuffstat$n_uncen_1, lsuffstat$tot_obst_1,
                         cfg$post_draw, cfg);
  apply(m , 2, function(x)log(2)/median(x))
  plot_tte_hist(m)
  
  # x0 <- rgamma(1000, 1000, lsuffstat$tot_obst_0)
  # x1 <- rgamma(1000, 1000, lsuffstat$tot_obst_1)
  # m2 <- cbind(x0, x1, x0/x1)
  # 
  # plot_tte_hist(m2)
  

  
  library(survival)
  library(spBayesSurv)
  library(coda)
  
  km_fit <- survfit(Surv(obst, 1-cen) ~ trt, data = d2)
  plot(km_fit)
  abline(v = 30)
  abline(v = 35)
  
  
  nburn=500; nsave=300; nskip=0;
  # Note larger nburn, nsave and nskip should be used in practice.
  mcmc=list(nburn=nburn, nsave=nsave, nskip=nskip, ndisplay=1000);
  prior = list(M=10, r0=1);
  # Fit the Cox PH model
  res1 = indeptCoxph(formula = Surv(obst, 1-cen)~trt, data=d2, 
                     prior=prior, mcmc=mcmc);
  sfit1=summary(res1); sfit1;
  
  par(mfrow = c(1,1))
  tgrid = seq(1e-10,60,0.1);
  xpred = data.frame(trt=c(0,1)); 
  plot(res1, xnewdata=xpred, tgrid=tgrid);
  
  
  
  
  
  
  
  
  
  
  # test stan
  
  
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  library(rstan)
  source("util.R")
  cfg <- readRDS("tests/cfg-example.RDS")
  
  look = 32
  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(copy(d))
  colnames(d2) <- dnames
  (lsuffstat1 <- rcpp_clin_set_obst(d, cfg, look))
  d2 <- as.data.frame(copy(d)); colnames(d2) <- dnames
  d3 <- d2[d2$trt == 0,]
  
  sdat <- list(
    ## Number of event individuals
    Nobs = sum(d3$cen == 0),
    ## Number of censored individuals
    Ncen = sum(d3$cen == 1),
    ## Times for event individuals
    yobs = d3$obst[d3$cen == 0],
    ## Times for censored individuals
    ycen = d3$obst[d3$cen == 1]
  )
  
  s1 <- rstan::stan(file = "expo.stan",
              data = sdat)
  
  
  
  
  
  
  look <- 32
  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(copy(d))
  colnames(d2) <- dnames
  (lsuffstat1 <- rcpp_clin_set_obst(d, cfg, look))
  d2 <- as.data.frame(copy(d))
  colnames(d2) <- dnames
  stan_weibull_survival_model_data <-
    list(
      ## Number of event individuals
      Nobs = sum(d2$cen == 0),
      ## Number of censored individuals
      Ncen = sum(d2$cen == 1),
      ## Number of covariates
      M_bg = 1,
      ## Times for event individuals
      yobs = d2$obst[d2$cen == 0],
      ## Times for censored individuals
      ycen = d2$obst[d2$cen == 1],
      ## Covariates for event individuals as a matrix
      Xobs_bg = matrix(d2$trt[d2$cen == 0]),
      ## Covariates for censored individuals as a matrix
      Xcen_bg = matrix(d2$trt[d2$cen == 1])
    )
  stan_weibull_survival_model_data
  
  
 
  rstan::stan(file = "weibull.stan",
              data = stan_weibull_survival_model_data)
  
})


# both produce unreliable parameter estimates. this is because we 
# look at the data at discrete times and therefore observe lower
# numbers of events than we would if we observed all the data at the time
# of each interim rather than at some arbitrary earlier point.
test_that("clinical endpoint post - conjugate versus stan", {
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  library(rstan)
  source("util.R")
  cfg <- readRDS("tests/cfg-example.RDS"); look = 32
  cfg$max_age_fu_months <- 45
  cfg$use_alt_censoring <- 1
  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(copy(d))
  colnames(d2) <- dnames
  plot_tte_hist_dat(d2$evtt, d2$trt, cfg$looks[look])
  (lsuffstat <- rcpp_clin_set_obst(d, cfg, look))
  d2 <- as.data.frame(copy(d)); colnames(d2) <- dnames
  plot_tte_hist_dat(d2$obst, d2$trt, cfg$looks[look])
  
  m <- matrix(0, nrow = cfg$post_draw, ncol = 3)
  rcpp_clin_interim_post(m, 
                         lsuffstat$n_uncen_0, lsuffstat$tot_obst_0,
                         lsuffstat$n_uncen_1, lsuffstat$tot_obst_1,
                         cfg$post_draw, cfg);
  plot_tte_hist(m)
  
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
  plot_tte_hist(beta)
  
  
})


test_that("clinical endpoint post - does the posterior make any sense?", {
  
  
  # test - use alternative censoring to see if it improves stability in being able to recover params
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  library(rstan)
  source("util.R")
  cfg <- readRDS("tests/cfg-example.RDS"); look = 31
  cfg$max_age_fu_months <- 36
  cfg$use_alt_censoring <- 1
  
  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(copy(d))
  colnames(d2) <- dnames
  #plot_tte_hist_dat(d2$evtt, d2$trt, cfg$looks[look])
  (lsuffstat <- rcpp_clin_set_obst(d, cfg, look))
  d2 <- as.data.frame(copy(d)); colnames(d2) <- dnames
 
  # 
  d2$cen2 <- NA
  d2$cen2 <- ifelse(d2$evtt > cfg$max_age_fu_months &
                      as.numeric(rownames(d2)) <= cfg$looks[look], 1, d2$cen2)
  d2$cen2 <- ifelse(d2$evtt <= cfg$max_age_fu_months &
                      as.numeric(rownames(d2)) <= cfg$looks[look], 0, d2$cen2)
  d2$obst2 <- NA
  d2$obst2 <- ifelse(!is.na(d2$cen2) & d2$cen2 == 0, d2$evtt, d2$obst2)
  d2$obst2 <- ifelse(!is.na(d2$cen2) & d2$cen2 == 1, 
                     pmin(cfg$max_age_fu_months, d2$age + d2$evtt ), d2$obst2)
  hist(d2$obst2)
  
  sum(1-d2$cen2, na.rm = T)
  sum(1-d2$cen, na.rm = T)
  
  sum(d2$obst2, na.rm = T)
  sum(d2$obst, na.rm = T)
  
  lsuff <- tte_suff_stats(obst = d2$obst2, 
                 trt = d2$trt, 
                 cen = d2$cen2, 
                 n = cfg$looks[look])
  
  m <- matrix(0, nrow = cfg$post_draw, ncol = 3)
  
  m[,1] <- rgamma(1000, cfg$prior_gamma_a + lsuff$n_uncen_0,   cfg$prior_gamma_b + lsuff$tot_obst_0)
  m[,2] <- rgamma(1000, cfg$prior_gamma_a + lsuff$n_uncen_1,  cfg$prior_gamma_b + lsuff$tot_obst_1)
  m[,3] <- m[,1]/m[,2]
  plot_tte_hist(m)
  
  m <- matrix(0, nrow = cfg$post_draw, ncol = 3)
  rcpp_clin_interim_post(m, 
                         lsuff$n_uncen_0, lsuff$tot_obst_0,
                         lsuff$n_uncen_1, lsuff$tot_obst_1,
                         cfg$post_draw, cfg);
  plot_tte_hist(m)
  
  m <- matrix(0, nrow = cfg$post_draw, ncol = 3)
  rcpp_clin_interim_post(m, 
                         lsuffstat$n_uncen_0, lsuffstat$tot_obst_0,
                         lsuffstat$n_uncen_1, lsuffstat$tot_obst_1,
                         cfg$post_draw, cfg);
  
  plot_tte_hist(m)
  
  # run the above as many times as you like. anecdotably it looks ok.
  

  
  
  
  
  # sim to examine whether the alternative censoring aligns with naive censoring
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  library(rstan)
  source("util.R")
  cfg <- readRDS("tests/cfg-example.RDS"); look = 31
  cfg$max_age_fu_months <- 36
  cfg$use_alt_censoring <- 1
  nsim <- 1000
  mr <- matrix(0, ncol = 6, nrow = nsim)
  
  for(i in 1:nsim){
    
    d <- rcpp_dat(cfg)
    (lsuffstat <- rcpp_clin_set_obst(d, cfg, look))
    d2 <- as.data.frame(copy(d)); colnames(d2) <- dnames
    # 
    d2$cen2 <- NA
    d2$cen2 <- ifelse(d2$evtt > cfg$max_age_fu_months &
                        as.numeric(rownames(d2)) <= cfg$looks[look], 1, d2$cen2)
    d2$cen2 <- ifelse(d2$evtt <= cfg$max_age_fu_months &
                        as.numeric(rownames(d2)) <= cfg$looks[look], 0, d2$cen2)
    d2$obst2 <- NA
    d2$obst2 <- ifelse(!is.na(d2$cen2) & d2$cen2 == 0, d2$evtt, d2$obst2)
    d2$obst2 <- ifelse(!is.na(d2$cen2) & d2$cen2 == 1, 
                       pmin(cfg$max_age_fu_months, d2$age + d2$evtt ), d2$obst2)

    lsuff <- tte_suff_stats(obst = d2$obst2, 
                            trt = d2$trt, 
                            cen = d2$cen2, 
                            n = cfg$looks[look])
    
    post_ctl <- rgamma(1000, cfg$prior_gamma_a + lsuff$n_uncen_0,   cfg$prior_gamma_b + lsuff$tot_obst_0)
    post_trt <- rgamma(1000, cfg$prior_gamma_a + lsuff$n_uncen_1,  cfg$prior_gamma_b + lsuff$tot_obst_1)
    post_ratio <- post_ctl / post_trt
    
    post_ctl <- mean(post_ctl)
    post_trt <- mean(post_trt)
    post_ratio <- mean(post_ratio)

    m <- matrix(0, nrow = cfg$post_draw, ncol = 3)
    rcpp_clin_interim_post(m, 
                           lsuffstat$n_uncen_0, lsuffstat$tot_obst_0,
                           lsuffstat$n_uncen_1, lsuffstat$tot_obst_1,
                           cfg$post_draw, cfg);

    post_ctl2 <- mean(m[,1])
    post_trt2 <- mean(m[,2])
    post_ratio2 <- mean(m[,3])
    
    mr[i,] <- c(post_ctl, post_trt, post_ratio, post_ctl2, post_trt2, post_ratio2)
  }
  
  pdf(file = "test-alt-censoring.pdf")
  par(mfrow = c(1, 3))
  hist(mr[,1], prob = T)
  lines(density(mr[,4]))
  hist(mr[,2], prob = T)
  lines(density(mr[,5]))
  hist(mr[,3], prob = T)
  lines(density(mr[,6]))
  par(mfrow = c(1, 1))
  dev.off()
  
  res <- log(2)/colMeans(mr)
  
  expect_equal(res[1], 30, tolerance = 0.5)
  expect_equal(res[2], 35, tolerance = 0.5)
  expect_equal(res[4], 30, tolerance = 0.5)
  expect_equal(res[5], 35, tolerance = 0.5)

  
  # sim to look at orvac censoring alignment
  
  
  
  # sim to examine whether the alternative censoring aligns with naive censoring
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  library(rstan)
  source("util.R")
  cfg <- readRDS("tests/cfg-example.RDS"); look = 31
  cfg$max_age_fu_months <- 36
  cfg$use_alt_censoring <- 0
  nsim <- 1000
  mr <- matrix(0, ncol = 6, nrow = nsim)
  
  for(i in 1:nsim){
    
    d <- rcpp_dat(cfg)
    (lsuffstat <- rcpp_clin_set_obst(d, cfg, look))
    d2 <- as.data.frame(copy(d)); colnames(d2) <- dnames
    # 
    d2$cen2 <- NA
    d2$cen2 <- ifelse(d2$evtt > cfg$max_age_fu_months &
                        as.numeric(rownames(d2)) <= cfg$looks[look], 1, d2$cen2)
    d2$cen2 <- ifelse(d2$evtt <= cfg$max_age_fu_months &
                        as.numeric(rownames(d2)) <= cfg$looks[look], 0, d2$cen2)
    d2$obst2 <- NA
    d2$obst2 <- ifelse(!is.na(d2$cen2) & d2$cen2 == 0, d2$evtt, d2$obst2)
    d2$obst2 <- ifelse(!is.na(d2$cen2) & d2$cen2 == 1, 
                       pmin(cfg$max_age_fu_months, d2$age + d2$evtt ), d2$obst2)
    
    lsuff <- tte_suff_stats(obst = d2$obst2, 
                            trt = d2$trt, 
                            cen = d2$cen2, 
                            n = cfg$looks[look])
    
    post_ctl <- rgamma(1000, cfg$prior_gamma_a + lsuff$n_uncen_0,   cfg$prior_gamma_b + lsuff$tot_obst_0)
    post_trt <- rgamma(1000, cfg$prior_gamma_a + lsuff$n_uncen_1,  cfg$prior_gamma_b + lsuff$tot_obst_1)
    post_ratio <- post_ctl / post_trt
    
    post_ctl <- mean(post_ctl)
    post_trt <- mean(post_trt)
    post_ratio <- mean(post_ratio)
    
    m <- matrix(0, nrow = cfg$post_draw, ncol = 3)
    rcpp_clin_interim_post(m, 
                           lsuffstat$n_uncen_0, lsuffstat$tot_obst_0,
                           lsuffstat$n_uncen_1, lsuffstat$tot_obst_1,
                           cfg$post_draw, cfg);
    
    post_ctl2 <- mean(m[,1])
    post_trt2 <- mean(m[,2])
    post_ratio2 <- mean(m[,3])
    
    mr[i,] <- c(post_ctl, post_trt, post_ratio, post_ctl2, post_trt2, post_ratio2)
  }
  
  pdf(file = "test-orvac-censoring.pdf")
  par(mfrow = c(1, 3))
  hist(mr[,1], prob = T)
  lines(density(mr[,4]))
  hist(mr[,2], prob = T)
  lines(density(mr[,5]))
  hist(mr[,3], prob = T)
  lines(density(mr[,6]))
  par(mfrow = c(1, 1))
  dev.off()
  
  res <- log(2)/colMeans(mr)
  
  expect_equal(res[1], 30, tolerance = 0.5)
  expect_equal(res[2], 35, tolerance = 0.5)
  expect_equal(res[4], 30, tolerance = 0.5)
  expect_equal(res[5], 35, tolerance = 0.5)
  
})

test_that("clinical endpoint ppos - does the posterior make any sense?", {
  

  library(testthat)
  library(orvacsim)
  library(data.table)
  library(survival)
  source("util.R")
  cfg <- readRDS("tests/cfg-example.RDS"); look = 31
  cfg$max_age_fu_months <- 36
  cfg$use_alt_censoring <- 0

  d <- rcpp_dat(cfg)
  lsuffstat <- rcpp_clin_set_obst(d, cfg, look)
  d2 <- as.data.frame(copy(d)); colnames(d2) <- dnames
  m <- matrix(0, nrow = cfg$post_draw, ncol = 3)
  rcpp_clin_interim_post(m, 
                         lsuffstat$n_uncen_0, lsuffstat$tot_obst_0,
                         lsuffstat$n_uncen_1, lsuffstat$tot_obst_1,
                         cfg$post_draw, cfg);
  plot_tte_hist(m)
  
  d_new <- copy(d)
  (nimpute = max(cfg$looks) - cfg$looks[look])
  lres <- rcpp_clin_interim_ppos(d_new, m, nimpute, look, cfg)
  str(lres)
  mean(lres$pvalue < 0.05)
  
  # get data
  
  # hazard is equal to lambda, => hazard ratio = 
  # compute_exp_rate(39)/compute_exp_rate(30) because:
  # compute_exp_rate(35)/compute_exp_rate(30) because:
  # h(x) = f(x)/S(x) = f(x)/(1-F(x))
  # f(x) = lambda exp(-lambda * x)
  # F(x) = 1 - exp(-lambda x)
  # S(x) = 1 - (1 - exp(-lambda x)) = exp(-lambda x)
  # h(x) = lambda exp(-lambda * x) / exp(-lambda x) = lambda
  
  nsim <- 1000
  hr <- numeric(nsim)
  
  for(i in 1:nsim){
    # i = i + 1
    d <- rcpp_dat(cfg)
    d2 <- as.data.frame(copy(d))
    colnames(d2) <- dnames
    # get sufficieitn stats
    (lsuffstat1 <- rcpp_clin_set_obst(d, cfg, look))
    rcpp_logrank(d, look, cfg)
    d2 <- as.data.frame(copy(d))
    colnames(d2) <- dnames
    s <- summary(cph1 <- coxph(Surv(obst, 1-cen) ~ trt, data = d2))
    hr[i] <- s$coefficients[2]
  }
  
  
  
  
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  library(survival)
  source("util.R")
  i = 0
  cfg <- readRDS("tests/cfg-example.RDS")
  look <- 8
  nsim <- 500
  pp1 <- numeric(nsim)
  pp2 <- numeric(nsim)
  
  for(i in 1:nsim){
    # i = i + 1
    d <- rcpp_dat(cfg)
    d2 <- as.data.frame(copy(d))
    colnames(d2) <- dnames
    # get sufficieitn stats
    (lsuffstat1 <- rcpp_clin_set_obst(d, cfg, look))
    
    # how many obs in total?
    cfg$looks[look]
    # obtain posterior based on current look 
    m <- matrix(0, nrow = cfg$post_draw, ncol = 3)
    rcpp_clin_interim_post(m, 
                           lsuffstat1$n_uncen_0, lsuffstat1$tot_obst_0,
                           lsuffstat1$n_uncen_1, lsuffstat1$tot_obst_1,
                           cfg$post_draw, cfg);
    plot_tte_hist(m)
    
    d_new <- copy(d)
    (nimpute = max(cfg$looks) - cfg$looks[look])
    lres <- rcpp_clin_interim_ppos(d_new, m, nimpute, look, cfg)
    str(lres)
    mean(lres$pvalue < 0.05)
    
    pp1[i] <- lres$ppos
    pp2[i] <- mean(lres$pvalue < 0.05)
  }
  
  
  
  
  
  
  
  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(copy(d))
  colnames(d2) <- dnames
  
  nsim <- nrow(m)
  win <- numeric(nsim)
  for(i in 1:nsim){
    
    # i = i + 1
    d_new <- copy(d)
    rcpp_dat_small(d_new, cfg, look, m[i, 1], m[i, 2])
    d2 <- as.data.frame(copy(d_new))
    colnames(d2) <- dnames
    tail(d2, 5)
    (lsuffstat <- rcpp_clin_set_obst(d_new, cfg, length(cfg$looks)))
    d2 <- as.data.frame(copy(d_new))
    colnames(d2) <- dnames
    tail(d2, 5)
    
    lrfit1 <- survdiff(Surv(obst, 1-cen) ~ trt, data = d2)
    
    p1 <- pchisq(lrfit1$chisq, 1, lower.tail = F)
    

    lrfit2 <- LogRank2(z1 = d2$obst[d2$trt == 0],
             delta1 = 1-d2$cen[d2$trt == 0],
             z2 = d2$obst[d2$trt == 1],
             delta2 = 1-d2$cen[d2$trt == 1])
    
    p2 <- lrfit2$ApproxPvalue2side
    win[i] <- p2<0.05
    
    expect_equal(p1, p2, tolerance = 0.0001)
    #c(p1, p2)
    #plot(survfit(Surv(obst, (1-cen)) ~ trt, data = d2))
  }
  d_new <- copy(d)
  lres <- rcpp_clin_interim_ppos(d_new, m, nimpute, look, cfg)
  str(lres)
  mean(win)
  
  
  lsuffstat <- rcpp_clin_set_obst(d_new, cfg, length(cfg$looks))

  m_new = matrix(0, ncol = 3, nrow = cfg$post_draw)
  rcpp_clin_interim_post(m_new,
                                 lsuffstat$n_uncen_0, lsuffstat$tot_obst_0,
                                 lsuffstat$n_uncen_1, lsuffstat$tot_obst_1,
                                 cfg$post_draw, cfg);
                         
  plot_tte_hist(m_new)
  
  
  
  
  # test
  
  
  nsim <- 1000
  ppmax <- numeric(nsim)
  for(i in 1:nsim){
    cfg <- readRDS("tests/cfg-example.RDS")
    look <- 8
    cfg$looks[look]
    cfg$interimmnths
    cfg$interimmnths[look]
    cfg$b1tte <- compute_exp_rate(42) - compute_exp_rate(30)
    # get data
    d <- rcpp_dat(cfg)
    (lsuffstat1 <- rcpp_clin_set_obst(d, cfg, look))
    d2 <- as.data.frame(copy(d))
    colnames(d2) <- dnames
    # obtain posterior based on current look 
    m <- matrix(0, nrow = cfg$post_draw, ncol = 3)
    rcpp_clin_interim_post(m, 
                           lsuffstat1$n_uncen_0, lsuffstat1$tot_obst_0,
                           lsuffstat1$n_uncen_1, lsuffstat1$tot_obst_1,
                           cfg$post_draw, cfg);
    d_new <- copy(d)
    (nimpute = max(cfg$looks) - cfg$looks[look])
    lres <- rcpp_clin_interim_ppos(d_new, m, nimpute, look, cfg)
    ppmax[i] <- lres$ppos
  }

  str(lres)
  # plot_tte_hist(m)
  # hist(lres$mean_ratio)
  #  
  
  
  library(Rcpp)
  cppFunction('int fx(){
    Rcpp::NumericVector v(10);
      return v.size();
  }')
  
})


# finish rcpp_clin loop

# test for rcpp_clin loop

# calling clin from main_3

# decision rules for main_3

# refactor sim.cpp to optimise speed - update rather than return.


tmptest <- function(){
  library(Rcpp)
  cppFunction("NumericVector callrexp(int n, double mean) { 
              NumericVector x(n);
              for(int i = 0; i < n; i ++){
              x[i] = R::rexp(mean);
              }
              return(x); }")
  
  
  cppFunction("NumericVector callrgamma(int n, double shape, double scale) { 
              return(rgamma(n, shape, scale)); }")
  
  
  # generate tte data using expo with mean beta (reciprocal of rate lambda)
  x1 <- callrexp(1000, 1/cfg$b0tte)
  x2 <- rexp(1000, cfg$b0tte)
  
  hist(x1, prob = T); lines(density(x2), col = "red")
  
  # censor based on max fu of 36 months
  c1 <- ifelse(x1 > 36, 1, 0)
  c2 <- ifelse(x2 > 36, 1, 0)
  
  x1 <- ifelse(x1 > 36, 36, x1)
  x2 <- ifelse(x2 > 36, 36, x2)
  hist(x1, prob = T); lines(density(x2), col = "red")
  
  # sufficient stats
  (nuncen1 <- sum(1-c1))
  (nuncen2 <- sum(1-c2))
  (obst1 <- sum(x1))
  (obst2 <- sum(x2))
  
  # generate posterior
  hist(rgamma(1000, 1, 50))
  abline(v = cfg$b0tte, col = "red", lwd = 2)
  b0post1 <- rgamma(1000, 1 + nuncen1, 50 + obst1)
  b0post2 <- rgamma(1000, 1 + nuncen2, 50 + obst2)
  
  # hist(callrgamma(1000, 1, 1/50))
  # abline(v = cfg$b0tte, col = "red", lwd = 2)
  
  b0post1b <- callrgamma(1000, 1 + nuncen1, 1/(50 + obst1))
  b0post2b <- callrgamma(1000, 1 + nuncen2, 1/(50 + obst2))
  
  par(mfrow = c(1, 2))
  hist(log(2)/b0post1, prob = T)
  lines(density(log(2)/b0post1b), col = "red")
  abline(v = log(2)/cfg$b0tte, col = "blue", lwd = 2)
  hist(log(2)/b0post2, prob = T)
  lines(density(log(2)/b0post2b), col = "red")
  abline(v = log(2)/cfg$b0tte, col = "blue", lwd = 2)
  par(mfrow = c(1, 1))
}



test_that("clin tte data - all", {
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  cfg <- readRDS("tests/cfg-example.RDS")
  
  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(d)
  names(d2) <- dnames
  look <- 32
  cfg$interimmnths[look]
  
  
  lres <- rcpp_clin(d, cfg, look)
  str(lres)
  lres$ppos_n
  lres$ppos_max
  
  
  
})


test_that("clin tte data - ppos", {
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  cfg <- readRDS("tests/cfg-example.RDS")
  
  d <- rcpp_dat(cfg)
  
  look <- 10
  cfg$interimmnths[look]
  cfg$looks[look]
  
  # add obst and censoring status based on current look and compute posterior
  l <- rcpp_clin(d, cfg, look)
  d2 <- as.data.frame(copy(d))
  names(d2) <- dnames
  
  expect_equal(sum(!is.na(d2$obst)), cfg$looks[look])
  
  l2 <- rcpp_clin_interim_ppos(d, l$posterior, cfg$nstop - cfg$looks[look], look, cfg)
  
  
  d3 <- as.data.frame(copy(d))
  names(d3) <- dnames
  
  
  nsim <- 10000
  m <- matrix(0, ncol = 3, nrow = nsim)
  
  for(i in 1:nsim){
    
    l <- rcpp_clin(d, cfg, look)
    
    m[i, 1] <- mean(l$posterior[, 1])
    m[i, 2] <- mean(l$posterior[, 2])
    m[i, 3] <- mean(l$posterior[, 3])
    
    # Think there is bias due to the discrete observation pattern
    # x <- rgamma(1000, 1 + l$n_uncen_0, 0.03 + l$tot_obst_0)
    # hist(log(2)/x)
    # abline(v = 30, col = "red")
    # abline(v = median(log(2)/x), col = "green",  lwd = 3)
    # median(log(2)/x)
    # 
    # x <- rgamma(1000, 1 + l$n_uncen_1, 0.03 + l$tot_obst_1)
    # hist(log(2)/x)
    # abline(v = 35, col = "red")
    # abline(v = median(log(2)/x), col = "green",  lwd = 3)
    # median(log(2)/x)
  }
  
  # this highlights the bias in each group but fortunately the
  # ratio looks ok (just)
  # l <- rcpp_clin(d, cfg, look)
  # plot_tte_hist(l$posterior)
  
  rat <- cfg$trt_med_tte / cfg$ctl_med_tte
  expect_lt(abs(mean(m[i, 3]) - rat), rat * 0.15)
  
  
  
  
  
})

test_that("clin tte data - all", {
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  cfg <- readRDS("tests/cfg-example.RDS")
  
  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(d)
  names(d2) <- dnames
  look <- 32
  cfg$interimmnths[look]
  
  cfg$post_draw <- 2000
  
  nsim <- 10000
  m <- matrix(0, ncol = 3, nrow = nsim)
  
  for(i in 1:nsim){
    
    l <- rcpp_clin(d, cfg, look)
    
    m[i, 1] <- mean(l$posterior[, 1])
    m[i, 2] <- mean(l$posterior[, 2])
    m[i, 3] <- mean(l$posterior[, 3])
    
    # Think there is bias due to the discrete observation pattern
    # x <- rgamma(1000, 1 + l$n_uncen_0, 0.03 + l$tot_obst_0)
    # hist(log(2)/x)
    # abline(v = 30, col = "red")
    # abline(v = median(log(2)/x), col = "green",  lwd = 3)
    # median(log(2)/x)
    # 
    # x <- rgamma(1000, 1 + l$n_uncen_1, 0.03 + l$tot_obst_1)
    # hist(log(2)/x)
    # abline(v = 35, col = "red")
    # abline(v = median(log(2)/x), col = "green",  lwd = 3)
    # median(log(2)/x)
  }
  
  # this highlights the bias in each group but fortunately the
  # ratio looks ok (just)
  # l <- rcpp_clin(d, cfg, look)
  # plot_tte_hist(l$posterior)
  
  rat <- cfg$trt_med_tte / cfg$ctl_med_tte
  expect_lt(abs(mean(m[i, 3]) - rat), rat * 0.15)
  
  
  
  
  
})




context("immu endpoint framework")


test_that("immu model - estimated post and diff", {
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  
  d <- readRDS("tests/dat-example.RDS")
  cfg <- readRDS("tests/cfg-example.RDS")
  
  # n_obs_grp <- 170
  # look <- 10
  
  look = 7
  
  # ensure setup is as expected
  expect_equal(cfg$looks, c(70L, 100L, 130L, 160L, 190L, 220L, 250L, 280L, 310L, 340L, 
                            370L, 400L, 430L, 460L, 490L, 520L, 550L, 580L, 610L, 640L, 670L, 
                            700L, 730L, 760L, 790L, 820L, 850L, 880L, 910L, 940L, 970L, 1000L))
  
  expect_equal(cfg$interimmnths, c(7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49, 
                                   52, 55, 58, 61, 64, 67, 70, 73, 76, 79, 82, 85, 88, 91, 94, 97, 
                                   100))
  
  expect_equal(cfg$interimmnths[look], 25)
  expect_equal(cfg$months_per_person, 0.1)
  
  # look at average results.
  res <- matrix(0, nrow = 1000, ncol = 3)
  
  for(i in 1:1000){
    d <- rcpp_dat(cfg)
    
    info_delay <- 0
    look <- 32
    nobs <- rcpp_n_obs(d, look, cfg$looks, cfg$interimmnths, info_delay)
    lnsero <- rcpp_lnsero(d, nobs)

    # call this instead.
    m <- rcpp_immu_interim_post(d, nobs, cfg$post_draw, lnsero)

    res[i,1] <- mean(m[, COL_THETA0 ])
    res[i,2] <- mean(m[, COL_THETA1 ])
    res[i,3] <- mean(m[, COL_DELTA ])
  }
  
  # hist(res[, 1 ])
  # hist(res[, 2 ])
  # hist(res[, 3 ])
  
  # got to be within 1%
  expect_lt(abs(cfg$baselineprobsero - mean(res[,1])),  cfg$baselineprobsero * 1/100)
  expect_lt(abs(cfg$trtprobsero - mean(res[,2])),  cfg$trtprobsero * 1/100)
  
  # got to be within 3%
  diff_true <- cfg$trtprobsero - cfg$baselineprobsero
  diff_mean_est <- mean(res[,2] - res[,1])
  expect_lt(abs( diff_true -   diff_mean_est  ),  diff_true * 3/100)
  
  # got to be within 2%
  expect_lt(abs(cfg$baselineprobsero - median(res[,1])), cfg$baselineprobsero * 2/100)
  expect_lt(abs(cfg$trtprobsero - median(res[,2])), cfg$trtprobsero * 2/100)
  
  # got to be within 3%
  diff_true <- cfg$trtprobsero - cfg$baselineprobsero
  diff_mean_est <- median(res[,2] - res[,1])
  expect_lt(abs( diff_true -   diff_mean_est  ),  diff_true * 3/100)

})


test_that("immu model - posterior predictive", {
  

  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  
  d <- readRDS("tests/dat-example.RDS")
  cfg <- readRDS("tests/cfg-example.RDS")

  run_loop <- function(cfg){
    v1 <- numeric(1000)
    v2 <- numeric(1000)
    for(i in 1:1000){
      
      d <- rcpp_dat(cfg)
      info_delay <- 2
      look <- 6
      nobs <- rcpp_n_obs(d, look, cfg$looks, cfg$interimmnths, info_delay)
      lnsero <- rcpp_lnsero(d, nobs)
      m <- rcpp_immu_interim_post(d, nobs, cfg$post_draw, lnsero)
      nimpute1 <- cfg$looks[look] - nobs
      ppos_n = rcpp_immu_interim_ppos(d, m, nobs, nimpute1, cfg$post_draw, lnsero, cfg)
      v1[i] <- ppos_n[[1]]
      v2[i] <- ppos_n[[2]]
    }
    return(list(v1, v2))
  }
  
    
  # ensure setup is as expected
  expect_equal(cfg$looks, c(70L, 100L, 130L, 160L, 190L, 220L, 250L, 280L, 310L, 340L, 
                            370L, 400L, 430L, 460L, 490L, 520L, 550L, 580L, 610L, 640L, 670L, 
                            700L, 730L, 760L, 790L, 820L, 850L, 880L, 910L, 940L, 970L, 1000L))
  
  expect_equal(cfg$interimmnths, c(7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49, 
                                   52, 55, 58, 61, 64, 67, 70, 73, 76, 79, 82, 85, 88, 91, 94, 97, 
                                   100))
  
  expect_equal(cfg$months_per_person, 0.1)
  
  cfg$baselineprobsero <- 0.4
  cfg$trtprobsero <- 0.5
  cfg$deltaserot3 <- compute_sero_delta(cfg$baselineprobsero, cfg$trtprobsero)
  cfg$deltaserot3
  v <- run_loop(cfg)
  hist(v[[1]])
  hist(v[[2]])
  # got to be within 1%
  true_diff <- cfg$trtprobsero - cfg$baselineprobsero
  expect_lt(abs(mean(v[[2]]) - true_diff),  1/100)
  
 
  cfg$baselineprobsero <- 0.4
  cfg$trtprobsero <- 0.4
  cfg$deltaserot3 <- compute_sero_delta(cfg$baselineprobsero, cfg$trtprobsero)
  cfg$deltaserot3
  v <- run_loop(cfg)
  hist(v[[1]])
  hist(v[[2]])
  true_diff <- cfg$trtprobsero - cfg$baselineprobsero
  expect_lt(abs(mean(v[[2]]) - true_diff),  1/100)
  mean(v[[1]])
  
  cfg$baselineprobsero <- 0.4
  cfg$trtprobsero <- 0.6
  cfg$deltaserot3 <- compute_sero_delta(cfg$baselineprobsero, cfg$trtprobsero)
  cfg$deltaserot3
  vnew <- run_loop(cfg)
  hist(vnew[[1]])
  hist(vnew[[2]])
  true_diff <- cfg$trtprobsero - cfg$baselineprobsero
  expect_lt(abs(mean(vnew[[2]]) - true_diff),  1/100)
  mean(vnew[[1]])
  
  
  
})


test_that("immu model - rcpp_immu call", {
  
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  
  d <- readRDS("tests/dat-example.RDS")
  cfg <- readRDS("tests/cfg-example.RDS")
  
  # ensure setup is as expected
  #                            1   2    3     4      5     6     7
  expect_equal(cfg$looks, c(70L, 100L, 130L, 160L, 190L, 220L, 250L, 280L, 310L, 340L, 
                            370L, 400L, 430L, 460L, 490L, 520L, 550L, 580L, 610L, 640L, 670L, 
                            700L, 730L, 760L, 790L, 820L, 850L, 880L, 910L, 940L, 970L, 1000L))
  
  expect_equal(cfg$interimmnths, c(7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49, 
                                   52, 55, 58, 61, 64, 67, 70, 73, 76, 79, 82, 85, 88, 91, 94, 97, 
                                   100))
  
  expect_equal(cfg$months_per_person, 0.1)
  
  cfg$baselineprobsero <- 0.4
  cfg$trtprobsero <- 0.5
  cfg$deltaserot3 <- compute_sero_delta(cfg$baselineprobsero, cfg$trtprobsero)
  cfg$deltaserot3
  cfg$sero_info_delay

  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(d)
  names(d2) <- dnames
  look <- 6
  l <- rcpp_immu(d, cfg, look)
  expect_equal(length(l), 0)
  
  hist(l$posterior[, COL_DELTA])
  abline(v = 0.1, col = "red")
 
  str(l)  
  
  
  
  
  
  
  # test for zero length when beyond nmaxsero
  
  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(d)
  names(d2) <- dnames
  look <- 7
  l <- rcpp_immu(d, cfg, look)
  expect_equal(length(l), 0)
  
  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(d)
  names(d2) <- dnames
  look <- 8
  l <- rcpp_immu(d, cfg, look)
  expect_equal(length(l), 0)
  
})








test_that("clin tte data - tmp", {
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  cfg <- readRDS("tests/cfg-example.RDS")
  
  set.seed(4343)
  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(d)
  names(d2) <- dnames
  look <- 32
  cfg$interimmnths[look]
  # saveRDS(d, "tests\tmp2.RDS")
  
  l <- rcpp_clin(d, cfg, look)
  
  d3 <- as.data.frame(l$d)
  names(d3) <- dnames
  
  sum(d3$obst[d3$trt == 0], na.rm = T)
  sum(d3$obst[d3$trt == 1], na.rm = T)
  
  # rcpp_visits(d, 1, look = 32, cfg)
  
  # tte for second record doesn't look right.
  
  addmargins(table(d3$trt, d3$cen, useNA = "always"))
  l$n_uncen_0
  l$n_uncen_1
  hist(log(2)/l$posterior[,1])
  hist(log(2)/l$posterior[,2])
  
  
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  cfg <- readRDS("tests/cfg-example.RDS")
  
  set.seed(4343)
  n <- 1000
  a <- 9
  b <- 1/0.5
  
  
  x <- seq(from =0.1, to = 16, length.out = 100)
  y <- dgamma(x, a, b)
  
  
  y2 <- rcpp_gamma(n, a, 1/b)
  
  hist(y2, probability = T)
  lines(x, y)
  
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
  
  library(eha)
  library(survival)
  library(boot)
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  source("sim.R")
  library(truncnorm)
  
  cfg <- readRDS("tests/cfg-example.RDS")
  
  cfg$baselineprobsero
  cfg$trtprobsero <- 0.55
  cfg$deltaserot3 <- compute_sero_delta(cfg$baselineprobsero, cfg$trtprobsero)
  
  l <- rcpp_dotrial(cfg)
  
  df <- as.data.frame(l$d)
  names(df) <- dnames
  boot.ci(doglm(df[1:cfg$nmaxsero, ]), type="norm")
  boot.ci(doexp(df), type="norm")
 
  
  
  # instead of final analysis - immu posterior: mean theta0 0.387026 mean theta1 0.489018 mean delta 0.101992 n delta gt0 950 post_prob_gt0 0.95
  # report lwr and upr 95% ci
  
  
  
  
  
})




