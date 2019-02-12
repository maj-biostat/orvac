library(testthat)
library(orvacsim)
library(data.table)




context("data generation")



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




test_that("immu model - ", {
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  
  d <- readRDS("tests/dat-example.RDS")
  
  
  cfg <- readRDS("tests/cfg-example.RDS")
  
  # n_obs_grp <- 170
  # look <- 10
  
  
  
  # look = 2
  # # ensure setup is as expected
  # cfg$looks <- c(70L, 100L, 130L)
  # cfg$interimmnths <- c(7, 10, 13)
  
  
  # look at average results.
  res <- matrix(0, nrow = 1000, ncol = 3)
  
  #for(i in 1:1000){
    # d <- rcpp_dat(cfg)
    # d <- d[1:130,]
    
    # saveRDS(d, "tests/tmp.RDS")
    d <- readRDS("tests/tmp.RDS")
    
    d2 <- rcpp_clin(d, cfg, cfg$nlooks)
  
    d2 <- as.data.frame(d2)
    names(d2) <- dnames
  #}
  
  
  
})


