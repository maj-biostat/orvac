library(testthat)
library(orvacsim)
library(data.table)




#context("censoring")

test_that("data creation - dgp for sero looks correct", {
  
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("util.R")
  
  d <- readRDS("tests/dat-example.RDS")
  cfg <- readRDS("tests/cfg-example.RDS")
  
  n_obs_grp <- 170
  look <- 10
  
  trt_status <- 0
  istart <- 1
  iend <- 10
  
  cfg$baselineprobsero
  cfg$trtprobsero
  
  compute_sero_delta(cfg$baselineprobsero, cfg$trtprobsero)
  
  m <- matrix(0, nrow = 1000, ncol = 2)
  
  for(i in 1:1000){
    d <- rcpp_gen_dat(cfg)
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
    d <- rcpp_gen_dat(cfg)
    m <- rcpp_mod_immu(d, cfg, look)
    
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
  
  # got to be within 2%
  diff_true <- cfg$trtprobsero - cfg$baselineprobsero
  diff_mean_est <- mean(res[,2] - res[,1])
  expect_lt(abs( diff_true -   diff_mean_est  ),  diff_true * 2/100)
  
  # got to be within 1%
  expect_lt(abs(cfg$baselineprobsero - median(res[,1])), cfg$baselineprobsero * 1/100)
  expect_lt(abs(cfg$trtprobsero - median(res[,2])), cfg$trtprobsero * 1/100)
  
  # got to be within 2%
  diff_true <- cfg$trtprobsero - cfg$baselineprobsero
  diff_mean_est <- median(res[,2] - res[,1])
  expect_lt(abs( diff_true -   diff_mean_est  ),  diff_true * 2/100)

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
    d <- rcpp_gen_dat(cfg)
    m <- rcpp_mod_immu(d, cfg, look)
    
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
  
  # got to be within 2%
  diff_true <- cfg$trtprobsero - cfg$baselineprobsero
  diff_mean_est <- mean(res[,2] - res[,1])
  expect_lt(abs( diff_true -   diff_mean_est  ),  diff_true * 2/100)
  
  # got to be within 1%
  expect_lt(abs(cfg$baselineprobsero - median(res[,1])), cfg$baselineprobsero * 1/100)
  expect_lt(abs(cfg$trtprobsero - median(res[,2])), cfg$trtprobsero * 1/100)
  
  # got to be within 2%
  diff_true <- cfg$trtprobsero - cfg$baselineprobsero
  diff_mean_est <- median(res[,2] - res[,1])
  expect_lt(abs( diff_true -   diff_mean_est  ),  diff_true * 2/100)
  
})














test_that("running mean with constant x or position", {
  
  d <- readRDS("tests/dat-example.RDS")
  cfg <- readRDS("tests/cfg-example.RDS")
  
  n_obs_grp <- 170
  look <- 10
  
  trt_status <- 0
  istart <- 1
  iend <- 10
  
  mynames <- names(d)
  data.table(rcpp_censoring(as.matrix(d), look, trt_status, iend, 10, 6))
  test <- data.table(rcpp_censoring(as.matrix(d), look, trt_status, iend, 10, 6))
  names(test) <- c(mynames, "currentage")
  
  n <- 100
  x <- rnorm(n)
  pos <- rep(0, n)
  expect_equal( runningmean(pos, x, window=1), rep(mean(x), n) )
  
  mu <- mean(x)
  x <- rep(mu, n)
  pos <- runif(n, 0, 5)
  expect_equal( runningmean(pos, x, window=1), x)
  
})