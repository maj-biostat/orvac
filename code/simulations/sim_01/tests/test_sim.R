library(testthat)
library(orvacsim)
library(data.table)




context("data generation")


test_that("data creation - by reference, by value", {
  
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
  
  
  # By reference is MUCH quicker!!!
  #
  # >   summary(duration1)
  # user  system elapsed 
  # 2.556   0.681   3.239 
  # >   summary(duration2)
  # user  system elapsed 
  # 13.974   1.792  15.791 
  
})




test_that("data creation - cpp quicker than data.table", {
  
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


test_that("data creation - all datas big and small", {
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  source("sim.R")

  cfg <- readRDS("tests/cfg-example.RDS")
  d <- rcpp_dat(cfg)
  
  # assignment is otherwise by reference.
  d2 <- copy(d)

  # this updates the evtt, fu1, fu2, cen and obst from observations 
  # on or after look
  rcpp_dat_small(d, cfg, look = 1, l0 = 0.02, l1 = 0.015)
  expect_false(isTRUE(all.equal(d2, d)))
  
  
  
  
})


test_that("data creation - dgp for sero looks correct", {
  
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
  
  expect_lt(abs(cfg$baselineprobsero - mean(m[,1])),  0.001)
  expect_lt(abs(cfg$trtprobsero - mean(m[,2])),  0.001)
  
  expect_lt(abs(cfg$baselineprobsero - median(m[,1])),  0.001)
  expect_lt(abs(cfg$trtprobsero - median(m[,2])),  0.001)


  # tests 
  # - repeat sampling of stochastic elements with summary stats
  # balance for tretament aloc
  # linear accrual
  
  
  
})


test_that("data creation - dgp for tte looks correct", {
  
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
  
  m <- matrix(0, nrow = 1000, ncol = 2)
  v <- matrix(0, nrow = 1000, ncol = 2)
  
  for(i in 1:1000){
    d <- rcpp_dat(cfg)
    
    m[i,1] <- mean(d[d[,COL_TRT] == 0, COL_EVTT])
    m[i,2] <- mean(d[d[,COL_TRT] == 1, COL_EVTT])
    
    v[i,1] <- var(d[d[,COL_TRT] == 0, COL_EVTT])
    v[i,2] <- var(d[d[,COL_TRT] == 1, COL_EVTT])
  }
  
  expect_lt(abs(1/lamb0 - mean(m[,1])),    mean(m[,1]) * 0.01)
  expect_lt(abs(1/lamb1 - mean(m[,2])),    mean(m[,2]) * 0.01)
  
  expect_lt(abs(1/(lamb0^2) - mean(v[,1])),    mean(v[,1]) * 0.01)
  expect_lt(abs(1/(lamb1^2) - mean(v[,2])),    mean(v[,2]) * 0.01)
  
  
  
  
  # tests 
  # - repeat sampling of stochastic elements with summary stats
  # balance for tretament aloc
  # linear accrual
  
  
  
})


test_that("data creation - retrieves correct number obs", {
  
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


test_that("data creation - correct sero counts", {
  
  
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
  nobs <- rcpp_n_obs(d, look, cfg$looks, cfg$interimmnths, info_delay)
  
  stopifnot(nobs == 70)
  d3 <- d2[1:nobs,]
  
  nseroctl <- sum(d3$serot3[d3$trt == 0])
  nserotrt <- sum(d3$serot3[d3$trt == 1])

  expect_equal(rcpp_lnsero(d, nobs), list(nseroctl, nserotrt))
  
  # inexact match
  info_delay <- 0.5
  look <- 1
  nobs <- rcpp_n_obs(d, look, cfg$looks, cfg$interimmnths, info_delay)
  
  stopifnot(nobs == 64)
  d3 <- d2[1:nobs,]
  
  nseroctl <- sum(d3$serot3[d3$trt == 0])
  nserotrt <- sum(d3$serot3[d3$trt == 1])
  
  expect_equal(rcpp_lnsero(d, nobs), list(nseroctl, nserotrt))
  
  
})


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


test_that("clin tte data - visit times", {
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  cfg <- readRDS("tests/cfg-example.RDS")

  set.seed(4343)
  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(d[1:5,])
  names(d2) <- dnames
  
  d2[1, "accrt"] = 0.2
  d2[1, "age"] = 6
  d2[1, "evtt"] = 4
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
  
  
  d2[1, "accrt"] = 0.2
  d2[1, "age"] = 6
  d2[1, "evtt"] = 4
  d2[1, "fu1"] = 0.5
  d2[1, "fu2"] = 2
  look <- 32
  
  cfg$interimmnths[look]
  visits2 <- rcpp_visits(as.matrix(d2), idxcpp, look, cfg)
  
  # correct length - maximum visits is the number of visits
  # such that (cfg$visit_lwr * 5 ) + 6 < 36 and then add two for fu1 and fu2
  # add one to do a lt check
  expect_lt(length(visits2), 8)
  # correct values
  expect_equal(visits2[1], d2[1, "accrt"] + d2[1, "fu1"])
  expect_equal(visits2[2], d2[1, "accrt"] + d2[1, "fu2"])
  expect_lt(max(visits2) + d2[1, "age"], cfg$max_age_fu_months + 0.000001)
  
  
  d2[1, "accrt"] = 0.2
  d2[1, "age"] = 6
  d2[1, "evtt"] = 4
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
  expect_lt(max(visits) + d2[1, "age"], cfg$max_age_fu_months + 0.000001)
  expect_lt(max(visits), cfg$interimmnths[look] + 0.000001)
  
  
  
  # test 4 - zero visits
  
  d2[1, "accrt"] = 6.9
  d2[1, "age"] = 6
  d2[1, "evtt"] = 1000
  d2[1, "fu1"] = 0.5
  d2[1, "fu2"] = 2
  look <- 1
  idxcpp <- 0
  
  cfg$interimmnths[look]
  visits <- rcpp_visits(as.matrix(d2), idxcpp, look, cfg)
  
  expect_equal(length(visits), 0)
})



test_that("clin tte data - censoring", {
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  cfg <- readRDS("tests/cfg-example.RDS")
  
  set.seed(4343)
  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(d[1:5,])
  names(d2) <- dnames
  
  # test 1 - not censored
  
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
  expect_equal(cens$obst, d2[1, "evtt"])
  
  
  
  # test 2 - censored due to evtt being after interim look
  
  d2[1, "accrt"] = 0.2
  d2[1, "age"] = 6
  d2[1, "evtt"] = 8
  d2[1, "fu1"] = 0.5
  d2[1, "fu2"] = 2
  look <- 1
  idxcpp <- 0
  
  cfg$interimmnths[look]
  visits <- rcpp_visits(as.matrix(d2), idxcpp, look, cfg)
  
  cens <- rcpp_cens(as.matrix(d2), visits, idxcpp, look, cfg)
  cens

  expect_equal(cens$cen, 1)
  expect_equal(cens$obst, max(visits) - d2[1, "accrt"])
  
  
  
  # test 3 - censored due to evtt being after age
  
  d2[1, "accrt"] = 0.2
  d2[1, "age"] = 6
  d2[1, "evtt"] = 1000
  d2[1, "fu1"] = 0.5
  d2[1, "fu2"] = 2
  look <- length(cfg$looks)
  idxcpp <- 0
  
  cfg$interimmnths[look]
  visits <- rcpp_visits(as.matrix(d2), idxcpp, look, cfg)
  
  cens <- rcpp_cens(as.matrix(d2), visits, idxcpp, look, cfg)
  cens

  expect_equal(cens$cen, 1)
  expect_lt(max(visits) - d2[1, "accrt"] + d2[1, "age"], cfg$max_age_fu_months)
  
  
  
  
  
  
  
  # test 4 - zero visits
  
  d2[1, "accrt"] = 6.9
  d2[1, "age"] = 6
  d2[1, "evtt"] = 1000
  d2[1, "fu1"] = 0.5
  d2[1, "fu2"] = 2
  look <- 1
  idxcpp <- 0
  
  cfg$interimmnths[look]
  visits <- rcpp_visits(as.matrix(d2), idxcpp, look, cfg)
  
  expect_equal(length(visits), 0)
  
  cens <- rcpp_cens(as.matrix(d2), visits, idxcpp, look, cfg)
  cens
  
  
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



plot_tte_hist <- function(m){
  
  par(mfrow = c(2, 2))
  
  med0 <- log(2)/m[, 1]
  med1 <- log(2)/m[, 2]
  
  hist(med0, probability = T, main = "")
  abline(v = 30, col = "red", lwd = 2)
  abline(v = median(med0), col = "blue", lwd = 2)
  hist(med1, probability = T, main = "")
  abline(v = 35, col = "red", lwd = 2)
  abline(v = median(med1), col = "blue", lwd = 2)
  hist(m[, 3], probability = T, main = "")
  abline(v = 35/30, col = "red", lwd = 2)
  abline(v = median(med1/med0) , col = "blue", lwd = 2)
  
  par(mfrow = c(1, 1))
}

test_gammy <- function(){
  set.seed(4343)
  n <- 1000
  a <- 1
  b <- 10
  
  hist(rgamma(n, a, b))
  
  
  x <- seq(from = 0.0, to = 1.5, length.out = 1000)
  y <- dgamma(x, shape = a, rate = b)
  
  
  plot(x, y, type = "l")
  
  # 0.009902103
  
  test <- rgamma(1000, a, b)
  mean(test)
  hist(test)
  
  
  y2 <- rcpp_gamma(n, a, 1/b)

  hist(y2, probability = T)
  lines(x, y, col = "red", lwd = 3)
  
  # 48 36 
  # 2070 2400
  
  c <- 1/b
  y2 <- rcpp_gamma(n, a, c)
  
  hist(y2, probability = T)  

  
  y3 <- rcpp_gamma(n, a + 48, c / (1 + c * 2000))
  hist(y3, probability = T)  
  
  y3 <- rcpp_gamma(n, a + 36, c / (1 + c * 2400))
  hist(y3, probability = T)   
  
  
  hist(rcpp_gamma(1000, 1 + 48, (1/20) / (1 + (1/20) * 2000)))
    
}




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







