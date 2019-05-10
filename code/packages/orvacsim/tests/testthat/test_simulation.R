# setwd("~/Documents/orvac/code/packages/orvacsim/tests/testthat")
# setwd("/home/mjones/Documents/orvac/code/packages/orvacsim/tests/testthat")
# setwd("/home/mark/Documents/orvac/code/packages/orvacsim/tests/testthat")

library(testthat)
library(ggplot2)
library(orvacsim)
library(data.table)
source("../../../../simulations/sim_01/util.R")
library(truncnorm)
library(survival)
library(boot)
library(rstan)
# library(eha)
library(rbenchmark)



test_that("dotrial", {



  setwd("~/Documents/orvac/code/packages/orvacsim/tests/testthat")
  #  setwd("/home/mark/Documents/orvac/code/packages/orvacsim/tests/testthat")

  library(testthat)
  library(orvacsim)
  library(data.table)
  source("../../../../simulations/sim_01/util.R")
  library(truncnorm)
  library(survival)
  library(boot)
  library(rbenchmark)

  cfg <- readRDS("cfg-example.RDS")

  cfg$baselineprobsero
  cfg$trtprobsero <- 0.55
  cfg$deltaserot3 <- compute_sero_delta(cfg$baselineprobsero, cfg$trtprobsero)

  set.seed(1)
  l <- rcpp_dotrial(1, cfg, TRUE)
  set.seed(2)
  l <- rcpp_dotrial(1, cfg, TRUE)
  set.seed(3)
  l <- rcpp_dotrial(1, cfg, TRUE)

  cfg$baselineprobsero
  cfg$trtprobsero <- 0.5
  cfg$deltaserot3 <- compute_sero_delta(cfg$baselineprobsero, cfg$trtprobsero)
  set.seed(4)
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


test_that("experimentation with clinical main loop", {


  setwd("~/Documents/orvac/code/packages/orvacsim/tests/testthat")
  #  setwd("/home/mark/Documents/orvac/code/packages/orvacsim/tests/testthat")

  library(testthat)
  library(orvacsim)
  library(data.table)
  source("../../../../simulations/sim_01/util.R")
  library(truncnorm)
  library(survival)
  library(boot)
  library(rbenchmark)

  cfg <- readRDS("cfg-example.RDS")

  # test - are the correct number of obs times and censor status set at first interim

  # question
  # is ppn less variable with larger post_draw?

  look <- 10
  simid <- 2
  nsim <- 1000
  c(cfg$interimmnths[look], cfg$looks[look])

  post_draw <- c(5000)
  pgt1 <- matrix(NA, nrow = nsim, ncol = length(post_draw))
  i = j = 1

  for(i in 1:length(post_draw)){
    cfg$post_draw <- post_draw[i]
    for(j in 1:nsim){
      d <- rcpp_dat(cfg)
      test <- rcpp_clin(d, cfg, look, j)
      pgt1[j, i] <- test$ppn
    }
  }

  plot(density(pgt1[,1]))
  lines(density(pgt1[,2]))
  lines(density(pgt1[,3]))

  colMeans(pgt1)
})


test_that("visualisation with clinical main loop", {



  mysummary <- function(clin_res){
    message("Observed ", paste(round(unlist(clin_res$lss_post),1), sup = " "))
    message("Interim  ", paste(round(unlist(clin_res$lss_int),1), sup = " "))
    message("Maximum  ", paste(round(unlist(clin_res$lss_max),1), sup = " "))
    message("pp intrm ", clin_res$ppn)
    message("pp max   ", clin_res$ppmax)
  }


  cfg <- readRDS("cfg-example.RDS")

  look <- 10
  simid <- 1
  c(cfg$interimmnths[look], cfg$looks[look])
  cfg$post_draw <- 10000

  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(copy(d))
  colnames(d2) <- dnames

  test <- rcpp_clin(d, cfg, look, simid)

  poster <- data.frame(l0_int = test$m[,1],
                       l1_int = test$m[,2],
                       l0_max = test$m[,1],
                       l1_max = test$m[,2])
  poster <- poster %>%
    tidyr::gather(key = "param", value = "val")

  pp_post1 <- data.frame(id = rep(1:cfg$post_draw, each = cfg$post_draw),
                               l0_int = test$int_post[,1],
                                l1_int = test$int_post[,2],
                         l0_max = test$max_post[,1],
                         l1_max = test$max_post[,2])

  pp_post1 <- pp_post1 %>%
    tidyr::gather(key = "param", value = "val", -id)

  idx <- sample(1:cfg$post_draw, size = 100, replace = F)
  pp_post1 <- pp_post1[pp_post1$id == idx,]
  ggplot(pp_post1, aes(x = val, group = factor(id)))+
    geom_line(stat = "density", alpha = 0.2)+
    geom_line(data = poster, aes(x = val),
              inherit.aes = F, stat = "density", col = "red")+
    ggtitle(paste0(cfg$post_draw)) +
    facet_wrap(~param)


#
  poster <- data.frame(ratio_int = test$m[,1]/test$m[,2],
                       ratio_max = test$m[,1]/test$m[,2])

  poster <- poster %>%
    tidyr::gather(key = "param", value = "val")

  pp_post1 <- data.frame(id = rep(1:cfg$post_draw, each = cfg$post_draw),
                         ratio_int = test$int_post[,1]/test$int_post[,2],
                         ratio_max = test$max_post[,1]/test$max_post[,2])

  pp_post1 <- pp_post1 %>%
    tidyr::gather(key = "param", value = "val", -id)

  idx <- sample(1:cfg$post_draw, size = 100, replace = F)
  pp_post1 <- pp_post1[pp_post1$id == idx,]
  ggplot(pp_post1, aes(x = val))+
    geom_line(stat = "density")+
    geom_line(stat = "density", aes(group = factor(id)), alpha = 0.2)+
    geom_line(data = poster, aes(x = val),
              inherit.aes = F, stat = "density", col = "red")+
    ggtitle(paste0(cfg$post_draw)) +
    facet_wrap(~param)

})


test_that("examine rcpp_clin_set_state updates at interim plus fu", {

  setwd("~/Documents/orvac/code/packages/orvacsim/tests/testthat")
  #  setwd("/home/mark/Documents/orvac/code/packages/orvacsim/tests/testthat")

  library(testthat)
  library(orvacsim)
  library(data.table)
  source("../../../../simulations/sim_01/util.R")
  library(truncnorm)
  library(survival)
  library(boot)
  library(rbenchmark)


  cfg <- readRDS("cfg-example.RDS")

  set.seed(1)
  simid <- 1
  look <- 2
  cfg$interimmnths[look]
  cfg$post_draw <- 1000
  d <- rcpp_dat(cfg)
  dat <- as.data.frame(copy(d))
  colnames(dat) <- dnames

  # observed event
  dat$accrt[1] <- 5
  dat$age[1] <- 6
  dat$evtt[1] <- 12
  obs_summary(dat, 1, look, cfg)
  dat <- as.matrix(dat)
  rcpp_clin_set_state(as.matrix(dat), look, 36, cfg)
  dat <- as.data.frame(copy(dat))
  colnames(dat) <- dnames
  obs_summary(dat, 1, look, cfg)
  expect_equal(dat$obst[1], 12)
  expect_equal(dat$cen[1], 0)
  expect_equal(dat$impute[1], 0)
  expect_equal(dat$reason[1], 1)
  expect_equal(dat$reftime[1], 35)

  # event after reftime, subj > max fu age
  dat$accrt[2] <- 5
  dat$age[2] <- 8
  dat$evtt[2] <- 32
  obs_summary(dat, 2, look, cfg)
  dat <- as.matrix(dat)
  rcpp_clin_set_state(as.matrix(dat), look, 36, cfg)
  dat <- as.data.frame(copy(dat))
  colnames(dat) <- dnames
  obs_summary(dat, 2, look, cfg)
  expect_equal(dat$obst[2], 28)
  expect_equal(dat$cen[2], 1)
  expect_equal(dat$impute[2], 1)
  expect_equal(dat$reason[2], 3)
  expect_equal(dat$reftime[2], 33)



  # max n - observed
  look <- length(cfg$looks)
  cfg$interimmnths[look]

  dat$accrt[1000] <- cfg$interimmnths[look]
  dat$age[1000] <- 10
  dat$evtt[1000] <- 12
  obs_summary(dat, 1000, look, cfg)
  dat <- as.matrix(dat)
  rcpp_clin_set_state(as.matrix(dat), look, 36, cfg)
  dat <- as.data.frame(copy(dat))
  colnames(dat) <- dnames
  obs_summary(dat, 1000, look, cfg)
  expect_equal(dat$obst[1000], 12)
  expect_equal(dat$cen[1000], 0)
  expect_equal(dat$impute[1000], 0)
  expect_equal(dat$reason[1000], 1)
  expect_equal(dat$reftime[1000], 126)

  # max n - censored
  look <- length(cfg$looks)
  cfg$interimmnths[look]

  dat$accrt[1000] <- cfg$interimmnths[look]
  dat$age[1000] <- 10
  dat$evtt[1000] <- 32
  obs_summary(dat, 1000, look, cfg)
  dat <- as.matrix(dat)
  rcpp_clin_set_state(as.matrix(dat), look, 36, cfg)
  dat <- as.data.frame(copy(dat))
  colnames(dat) <- dnames
  obs_summary(dat, 1000, look, cfg)
  expect_equal(dat$obst[1000], 26)
  expect_equal(dat$cen[1000], 1)
  expect_equal(dat$impute[1000], 1)
  expect_equal(dat$reason[1000], 3)
  expect_equal(dat$reftime[1000], 126)




})

test_that("examine rcpp_clin_set_state updates at interim", {

  cfg <- readRDS("cfg-example.RDS")

  set.seed(1)
  simid <- 1
  look <- 10
  cfg$post_draw <- 1000
  d <- rcpp_dat(cfg)
  dat <- as.data.frame(copy(d))
  colnames(dat) <- dnames
  # observed event
  dat$age[1] <- 8
  dat$evtt[1] <- 4
  obs_summary(dat, 1, look, cfg)
  # event prior to reftime but sub too old
  dat$accrt[2] <- 10
  dat$age[2] <- 35
  dat$evtt[2] <- 2
  obs_summary(dat, 2, look, cfg)
  # event after reftime and age not above max
  dat$accrt[3] <- 33
  dat$age[3] <- 6
  dat$evtt[3] <- 5
  obs_summary(dat, 3, look, cfg)
  # event after reftime and age above max prior to reftime
  dat$accrt[4] <- 10
  dat$age[4] <- 34
  dat$evtt[4] <- 30
  obs_summary(dat, 3, look, cfg)

  dat <- as.matrix(dat)
  rcpp_clin_set_state(as.matrix(dat), look, 0, cfg)
  dat <- as.data.frame(copy(dat))
  colnames(dat) <- dnames
  obs_summary(dat, 1, look, cfg)
  obs_summary(dat, 2, look, cfg)
  obs_summary(dat, 3, look, cfg)
  obs_summary(dat, 4, look, cfg)

  expect_equal(dat$obst[1], 4)
  expect_equal(dat$obst[2], 1)
  expect_equal(dat$obst[3], 1)
  expect_equal(dat$obst[4], 2)

  expect_equal(dat$cen[1], 0)
  expect_equal(dat$cen[2], 1)
  expect_equal(dat$cen[3], 1)
  expect_equal(dat$cen[4], 1)

  expect_equal(dat$impute[1], 0)
  expect_equal(dat$impute[2], 0)
  expect_equal(dat$impute[3], 1)
  expect_equal(dat$impute[4], 0)

  expect_equal(dat$reason[1], 1)
  expect_equal(dat$reason[2], 2)
  expect_equal(dat$reason[3], 3)
  expect_equal(dat$reason[4], 4)

  expect_equal(dat$reftime[1], 34)
  expect_equal(dat$reftime[2], 34)
  expect_equal(dat$reftime[3], 34)
  expect_equal(dat$reftime[4], 34)


  expect_equal(dat$reftime[cfg$looks[look]], 34)

})


test_that("visit times", {


  cfg <- readRDS("cfg-example.RDS")

  set.seed(4343)
  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(d[1:5,])
  names(d2) <- dnames

  # enrolled at time 0.2 and age 6 months with fu1 = 0.5 after enrolling
  # and fu2=2 time units after enrolling implies fu1 happes at 0.7,
  # and fu2 happens at 2.2
  # the surveillance then happens every six months after fu2 implying the
  # next visit will be at 6.2 then 12.2 etc.
  d2[1, "accrt"] = 60
  d2[1, "evtt"] = 10
  d2[1, "age"] = 6
  d2[1, "fu1"] = 0.5
  d2[1, "fu2"] = 2
  look <- 20
  idxcpp <- 0

  c(cfg$interimmnths[look], cfg$looks[look], cfg$looks_target[look])
  cfg$interimmnths
  cfg$looks
  cfg$looks_target
  d2
  # without followup past the time of the interim, which is month 64
  visits <- rcpp_clin_visits(as.matrix(d2), idxcpp, look, cfg, 0)
  visits
  # __with__ followup i.e. up to 36 months past the current interim
  visits <- rcpp_clin_visits(as.matrix(d2), idxcpp, look, cfg, 1)
  visits

  # set.seed(4343)
  # d <- rcpp_dat(cfg)
  # d2 <- as.data.frame(d[1:5,])
  # names(d2) <- dnames
  #
  # # test assuming at first interim analysis (look = 1)
  #
  # d2[1, "accrt"] = 0.2
  # d2[1, "age"] = 6
  # d2[1, "fu1"] = 0.5
  # d2[1, "fu2"] = 2
  # look <- 1
  # idxcpp <- 0
  #
  # cfg$interimmnths[look]
  # cfg$looks
  # cfg$looks_target
  # visits <- rcpp_visits(as.matrix(d2), idxcpp, look, cfg, 0)
  #
  # # correct length
  # expect_equal(length(visits), 3)
  # # correct values
  # expect_equal(visits[1], d2[1, "accrt"] + d2[1, "fu1"])
  # expect_equal(visits[2], d2[1, "accrt"] + d2[1, "fu2"])
  # expect_lt(max(visits), cfg$interimmnths[1])
  #
  # # test assuming at last interim analysis (look = 32)
  #
  # d2[1, "accrt"] = 0.2
  # d2[1, "age"] = 6
  # d2[1, "fu1"] = 0.5
  # d2[1, "fu2"] = 2
  # look <- 32
  #
  # cfg$interimmnths[look]
  # visits2 <- rcpp_visits(as.matrix(d2), idxcpp, look, cfg, 0)
  #
  # expect_equal(visits2[1], d2[1, "accrt"] + d2[1, "fu1"])
  # expect_equal(visits2[2], d2[1, "accrt"] + d2[1, "fu2"])
  # expect_lt(max(visits), cfg$max_age_fu_months)
  #
  # # test assuming at mid interim analysis (look = 7)
  #
  # d2[1, "accrt"] = 0.2
  # d2[1, "age"] = 6
  # d2[1, "fu1"] = 0.5
  # d2[1, "fu2"] = 2
  # look <- 7
  #
  # cfg$interimmnths[look]
  # visits <- rcpp_visits(as.matrix(d2), idxcpp, look, cfg, 0)
  #
  # # correct length - maximum visits is the number of visits
  # # such that (cfg$visit_lwr * 5 ) + 6 < 36 and then add two for fu1 and fu2
  # # add one to do a lt check
  # expect_lt(length(visits), length(visits2))
  # # correct values
  # expect_equal(visits[1], d2[1, "accrt"] + d2[1, "fu1"])
  # expect_equal(visits[2], d2[1, "accrt"] + d2[1, "fu2"])
  # expect_lt(max(visits) + d2[1, "age"], cfg$max_age_fu_months + 0.0001)
  # expect_lt(max(visits), cfg$interimmnths[look] + 0.0001)
  #
  # # test 4 - test zero visits returned
  #
  # d2[1, "accrt"] = 6.9
  # d2[1, "age"] = 6
  # d2[1, "fu1"] = 0.5
  # d2[1, "fu2"] = 2
  # look <- 1
  # idxcpp <- 0
  #
  # cfg$interimmnths[look]
  # visits <- rcpp_visits(as.matrix(d2), idxcpp, look, cfg, 0)
  #
  # expect_equal(length(visits), 0)
  #
  #
  # # retest at max visit - late enrollers
  #
  # # if this is just the last interim then we expect zero visits
  # set.seed(4343)
  # d <- rcpp_dat(cfg)
  # startidx <- nrow(d) - 5
  # endidx <- nrow(d)
  # idx <- endidx - startidx + 1
  # d2 <- as.data.frame(d[startidx:endidx,])
  # names(d2) <- dnames
  # d2$age[idx] = 6
  # d2[idx, "fu1"] = 0.5
  # d2[idx, "fu2"] = 2
  # d2[idx, "evtt"] = 3
  # look <- 32
  #
  # cfg$interimmnths[look]
  # visits2 <- rcpp_visits(as.matrix(d2), idx - 1, look, cfg, 0)
  # visits2
  #
  # # correct length - maximum visits is the number of visits
  # # such that (cfg$visit_lwr * 5 ) + 6 < 36 and then add two for fu1 and fu2
  # # add one to do a lt check
  # expect_equal(length(visits2), 0)
  #
  # # however if this is the final analysis we want to wait until we have all
  # # that have enrolled
  # cfg$interimmnths[look]
  # dofinal = 1
  # visits2 <- rcpp_visits(as.matrix(d2), idx - 1, look, cfg, dofinal)
  # visits2
  #
  # # correct length - maximum visits is the number of visits
  # # such that (cfg$visit_lwr * 5 ) + 6 < 36 and then add two for fu1 and fu2
  # # add one to do a lt check
  # expect_equal(length(visits2), 3)
  # # correct values
  # expect_equal(visits2[1], d2[idx, "accrt"] + d2[idx, "fu1"])
  # expect_equal(visits2[2], d2[idx, "accrt"] + d2[idx, "fu2"])
  # expect_lte(max(visits2) + d2[idx, "age"] - d2[idx, "accrt"], cfg$max_age_fu_months + 1)




})

test_that("censoring at interim", {


  cfg <- readRDS("cfg-example.RDS")

  # censored - because no visit
  set.seed(4343)
  d <- rcpp_dat(cfg)
  # just keep row 1 to 5
  d2 <- as.data.frame(d[1:5,])
  names(d2) <- dnames
  d2[1, "accrt"] = 63.8
  d2[1, "evtt"] = 10
  d2[1, "age"] = 6
  d2[1, "fu1"] = 0.5
  d2[1, "fu2"] = 2
  look <- 20
  idxcpp <- 0

  c(cfg$interimmnths[look], cfg$looks[look], cfg$looks_target[look])
  cfg$interimmnths
  cfg$looks
  cfg$looks_target
  d2
  withfollowup <- 0
  visits <- rcpp_clin_visits(as.matrix(d2), idxcpp, look, cfg, withfollowup)
  visits


  cens_status <- rcpp_clin_cens(as.matrix(d2), idxcpp, look, cfg, visits, withfollowup)
  expect_equal(cens_status$cen, 1)
  expect_equal(cens_status$obst, 0.2)




  # observed event
  set.seed(4343)
  d <- rcpp_dat(cfg)
  # just keep row 1 to 5
  d2 <- as.data.frame(d[1:5,])
  names(d2) <- dnames
  # enrolled at time 60 with an event 1 month later
  # age at event < 36 months
  d2[1, "accrt"] = 60
  d2[1, "evtt"] = 1
  d2[1, "age"] = 6
  d2[1, "fu1"] = 0.5
  d2[1, "fu2"] = 2
  look <- 20
  idxcpp <- 0

  c(cfg$interimmnths[look], cfg$looks[look], cfg$looks_target[look])
  cfg$interimmnths
  cfg$looks
  cfg$looks_target
  d2
  withfollowup <- 0
  # has a visit at 60.5 and 62
  visits <- rcpp_clin_visits(as.matrix(d2), idxcpp, look, cfg, withfollowup)
  visits

  cens_status <- rcpp_clin_cens(as.matrix(d2), idxcpp, look, cfg, visits, withfollowup)
  cens_status
  expect_equal(cens_status$cen, 0)
  expect_equal(cens_status$obst, 1)




  # censored - event before max visits but too old
  # event observed at time 5 + 30 when kid is 40 months
  # set.seed(4343)
  # d <- rcpp_dat(cfg)
  # # just keep row 1 to 5
  # d2 <- as.data.frame(d[1:5,])
  # names(d2) <- dnames
  # d2[1, "accrt"] = 5
  # d2[1, "evtt"] = 30
  # d2[1, "age"] = 10
  # d2[1, "fu1"] = 0.5
  # d2[1, "fu2"] = 2
  # look <- 28
  # idxcpp <- 0
  # c(cfg$interimmnths[look], cfg$looks[look], cfg$looks_target[look])
  # # cfg$interimmnths
  # # cfg$looks
  # # cfg$looks_target
  # d2
  # visits <- rcpp_visits(as.matrix(d2), idxcpp, look, cfg, 0)
  # # contrive last visit
  # visits[7] <- 36
  # visits
  #
  # cens_status <- rcpp_cens_interim(as.matrix(d2), visits, idxcpp, look, F, cfg)
  # cens_status
  # expect_equal(cens_status$cen, 1)
  # expect_equal(cens_status$obst, 26)
  #
  #
  #
  #
  #
  #
  # # censored - event after visit and too old
  # # event observed at time 5 + 30 when kid is 40 months
  # set.seed(4343)
  # d <- rcpp_dat(cfg)
  # # just keep row 1 to 5
  # d2 <- as.data.frame(d[1:5,])
  # names(d2) <- dnames
  # d2[1, "accrt"] = 5
  # d2[1, "evtt"] = 32
  # d2[1, "age"] = 10
  # d2[1, "fu1"] = 0.5
  # d2[1, "fu2"] = 2
  # look <- 28
  # idxcpp <- 0
  # c(cfg$interimmnths[look], cfg$looks[look], cfg$looks_target[look])
  # # cfg$interimmnths
  # # cfg$looks
  # # cfg$looks_target
  # d2
  # visits <- rcpp_visits(as.matrix(d2), idxcpp, look, cfg, 0)
  # visits
  #
  # cens_status <- rcpp_cens_interim(as.matrix(d2), visits, idxcpp, look, F, cfg)
  # cens_status
  # expect_equal(cens_status$cen, 1)
  # expect_equal(cens_status$obst, 26)
  #
  #
  #
  #
  #
  # # censored - event after visit and age < max age
  # # event observed at time 5 + 30 when kid is 40 months
  # set.seed(4343)
  # d <- rcpp_dat(cfg)
  # # just keep row 1 to 5
  # d2 <- as.data.frame(d[1:5,])
  # names(d2) <- dnames
  # d2[1, "accrt"] = 5
  # d2[1, "evtt"] = 32
  # d2[1, "age"] = 10
  # d2[1, "fu1"] = 0.5
  # d2[1, "fu2"] = 2
  # look <- 4
  # idxcpp <- 0
  # c(cfg$interimmnths[look], cfg$looks[look], cfg$looks_target[look])
  # # cfg$interimmnths
  # # cfg$looks
  # # cfg$looks_target
  # d2
  # visits <- rcpp_visits(as.matrix(d2), idxcpp, look, cfg, 0)
  # visits
  #
  # cens_status <- rcpp_cens_interim(as.matrix(d2), visits, idxcpp, look, F, cfg)
  # cens_status
  # expect_equal(cens_status$cen, 1)
  # expect_equal(cens_status$obst, 11)
  #
  #
  #
  #
  # # observed event
  # set.seed(4343)
  # d <- rcpp_dat(cfg)
  # # just keep row 1 to 5
  # d2 <- as.data.frame(d[1:5,])
  # names(d2) <- dnames
  # d2[1, "accrt"] = 60
  # d2[1, "evtt"] = 1
  # d2[1, "age"] = 6
  # d2[1, "fu1"] = 0.5
  # d2[1, "fu2"] = 2
  # look <- 20
  # idxcpp <- 0
  #
  # c(cfg$interimmnths[look], cfg$looks[look], cfg$looks_target[look])
  # cfg$interimmnths
  # cfg$looks
  # cfg$looks_target
  # d2
  # visits <- rcpp_visits(as.matrix(d2), idxcpp, look, cfg, 0)
  # visits
  #
  # cens_status <- rcpp_cens_interim(as.matrix(d2), visits, idxcpp, look, T, cfg)
  # cens_status
  # expect_equal(cens_status$cen, 0)
  # expect_equal(cens_status$obst, 1)

})

test_that("censoring at final", {



  cfg <- readRDS("cfg-example.RDS")

  # observed event
  set.seed(4343)
  d <- rcpp_dat(cfg)
  # just keep row 1 to 5
  d2 <- as.data.frame(d[1:5,])
  names(d2) <- dnames
  d2[1, "accrt"] = 63.8
  d2[1, "evtt"] = 10
  d2[1, "age"] = 6
  d2[1, "fu1"] = 0.5
  d2[1, "fu2"] = 2
  look <- 32
  idxcpp <- 0

  c(cfg$interimmnths[look], cfg$looks[look], cfg$looks_target[look])
  cfg$interimmnths
  cfg$looks
  cfg$looks_target
  d2

  withfollowup <- 1
  visits <- rcpp_clin_visits(as.matrix(d2), idxcpp, look, cfg, withfollowup)
  visits

  cens_status <- rcpp_clin_cens(as.matrix(d2), idxcpp, look, cfg, visits, withfollowup)
  expect_equal(cens_status$cen, 0)
  expect_equal(cens_status$obst, 10)





  # # event after last visit but kid too old
  # set.seed(4343)
  # d <- rcpp_dat(cfg)
  # # just keep row 1 to 5
  # d2 <- as.data.frame(d[999:1000,])
  # names(d2) <- dnames
  # d2[1, "accrt"] = 100
  # d2[1, "evtt"] = 32
  # d2[1, "age"] = 10
  # d2[1, "fu1"] = 0.5
  # d2[1, "fu2"] = 2
  # look <- 32
  # idxcpp <- 0
  #
  # c(cfg$interimmnths[look], cfg$looks[look], cfg$looks_target[look])
  # cfg$interimmnths
  # cfg$looks
  # cfg$looks_target
  # d2
  # visits <- rcpp_visits(as.matrix(d2), idxcpp, look, cfg, 1)
  # visits
  #
  # cens_status <- rcpp_cens_final(as.matrix(d2),visits,idxcpp,look,cfg)
  # expect_equal(cens_status$cen, 1)
  # expect_equal(cens_status$obst, 36)



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

  setwd("/home/mark/Documents/orvac/code/packages/orvacsim/tests/testthat")

library(testthat)
library(orvacsim)
library(data.table)
source("../../../../simulations/sim_01/util.R")
library(truncnorm)
library(survival)
library(boot)
library(rstan)
# library(eha)
library(rbenchmark)

  cfg <- readRDS("cfg-example.RDS")

  # test - are the correct number of obs times and censor status set at first interim

  set.seed(4343)
  look <- 7
  cfg$interimmnths[look]
  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(copy(d))
  colnames(d2) <- dnames

  c(cfg$interimmnths[look], cfg$looks[look], cfg$looks_target[look])

  with_followup <- 0

  idxcpp <- 0
  cfg$interimmnths[look]

  suffstat <- rcpp_clin_set_state(d, look, cfg, with_followup)
  suffstat

  #rcpp_clin_set_obst(d, cfg, look, F, F)

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

  rcpp_clin_set_obst(d, cfg, look, F, F)

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

  rcpp_clin_set_obst(d, cfg, look, F, F)

  d3 <- as.data.frame(copy(d))
  colnames(d3) <- dnames

  sum(!is.na(d3$obst)) ;   cfg$looks[look]

  expect_equal(sum(!is.na(d3$obst)), cfg$looks[look], tolerance = 0.01)
  expect_equal(sum(!is.na(d3$cen)), cfg$looks[look], tolerance = 0.01)





  # what happens at the final?

  look <- 32
  cfg$interimmnths[look]
  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(copy(d))
  colnames(d2) <- dnames

  rcpp_clin_set_obst(d, cfg, look, T, F)


})

context("clinical endpoint framework")



test_that("clinical endpoint posterior", {

  cfg <- readRDS("cfg-example.RDS")


  # test - estimated median tte is biased in both arms due to censoring mechanism
  # but ratio is preserved with lower than 2% bias
  look <- 28
  nsim <- 1000
  mres <- matrix(0, ncol = 3, nrow = nsim)
  cfg$post_draw <- 5000

  log(2)/(cfg$b0tte)
  log(2)/(cfg$b0tte + cfg$b1tte)

  for(i in 1:nsim){
    # i = i + 1
    d <- rcpp_dat(cfg)
    d2 <- as.data.frame(copy(d))
    colnames(d2) <- dnames

    # get sufficieitn stats
    (lsuffstat1 <- rcpp_clin_set_obst(d, cfg, look, F, F))

    d2 <- as.data.frame(copy(d))
    colnames(d2) <- dnames

    #lsuffstat2 <- rcpp_clin_set_obst(d, cfg, look)

    # obtain posterior based on current look
    m <- matrix(0, nrow = cfg$post_draw, ncol = 3)
    rcpp_clin_interim_post(m,
                           lsuffstat1$n_uncen_0, lsuffstat1$tot_obst_0,
                           lsuffstat1$n_uncen_1, lsuffstat1$tot_obst_1,
                           cfg$post_draw, cfg);

    mres[i,] <- colMeans(m)
  }

  # par(mfrow = c(2, 2))
  # plot(density(mres[,1]))
  # plot(density(mres[,2]))
  # plot(density(mres[,3]))
  # par(mfrow = c(1, 1))

  expect_equal(mres[1], cfg$b0tte, tolerance = 0.005)
  expect_equal(mres[2], cfg$b0tte + cfg$b1tte , tolerance = 0.005)

})


test_that("clin tte posterior check against mcmc", {


  cfg <- readRDS("cfg-example.RDS")

  d <- rcpp_dat(cfg)
  d2 <- as.data.frame(copy(d))
  colnames(d2) <- dnames
  look <- 32
  lsuffstat <- rcpp_clin_set_obst(d, cfg, look, T, F)
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


test_that("clin tte data ppos basic", {

  setwd("/home/mark/Documents/orvac/code/packages/orvacsim/tests/testthat")
  library(testthat)
  library(orvacsim)
  library(data.table)
  source("../../../../simulations/sim_01/util.R")
  library(truncnorm)
  library(survival)
  library(rbenchmark)
  library(configr)

  cfg <- readRDS("cfg-example.RDS")
  #cfg <- readRDS("cfg-50.RDS")

  look <- 12
  c(cfg$interimmnths[look], cfg$looks[look], cfg$looks_target[look])
  cfg$interimmnths
  cfg$looks
  cfg$looks_target

  cfg$ctl_med_tte <- 30
  cfg$trt_med_tte <- 35

  cfg$b0tte <- log(2)/cfg$ctl_med_tte
  cfg$b1tte <- (log(2)/cfg$trt_med_tte) - (log(2)/cfg$ctl_med_tte)

  #cfg$post_draw <- 1000


  #set.seed(5)
  nsim <- 500
  lres <- list()
  ppos_int <- numeric(nsim)
  ppos_fin <- numeric(nsim)
  for(i in 1:nsim){
    d <- rcpp_dat(cfg)
    lres[[i]] <- rcpp_clin_opt(d, cfg, look, i)
    ppos_int[i] <- lres[[i]]$ppn
    ppos_fin[i] <- lres[[i]]$ppmax
  }


  # save the results
  # l30 <- list(ppos_int = ppos_int, ppos_fin = ppos_fin, lres = lres)
  # l50 <- list(ppos_int = ppos_int, ppos_fin = ppos_fin, lres = lres)
  # run with different accrual
  # compare

  par(mfrow = c(2, 1))
  hist(l30$ppos_int, prob = T)
  abline(v = mean(l30$ppos_int))
  hist(l50$ppos_int, prob = T)
  abline(v = mean(l50$ppos_int))



  dp1 <- test$pp_int
  dp2 <- test$pp_fin
  dmm1 <- as.data.frame(test$m_post)
  dmm2 <- as.data.frame(test$m_newint)
  dmm3 <- as.data.frame(test$m_newfin)

  hist(dp1)
  hist(dp2)
  hist(dmm1[,3], xlim = c(0.3, 3))
  hist(dmm2[,3], xlim = c(0.3, 3))
  hist(dmm3[,3], xlim = c(0.3, 3))
  #

})

test_that("clin tte data ppos", {

  setwd("/home/mjones/Documents/orvac/code/packages/orvacsim/tests/testthat")

  library(testthat)
  library(orvacsim)
  library(data.table)
  source("../../../../simulations/sim_01/util.R")
  library(truncnorm)
  library(survival)
  library(boot)
  library(rstan)
  # library(eha)
  library(rbenchmark)

  cfg <- readRDS("cfg-example.RDS")

  # no effect
  cfg$b1tte <- 0
  cfg$post_draw <- 2000
  d <- rcpp_dat(cfg)
  look <- 6

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




context("data generation")




test_that("dgp for sero correct", {

  cfg <- readRDS("cfg-example.RDS")

  d <- rcpp_dat(cfg)

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
  look <- 2
  expect_equal(rcpp_n_obs(d, look, cfg$looks, cfg$interimmnths, info_delay), 100)
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


test_that("do immu trial", {

  cfg <- readRDS("cfg-example.RDS")
  cfg$sero_info_delay = 1
  d <- rcpp_dat(cfg)

  look <- 2
  nobs <- rcpp_n_obs(d, look, cfg$looks, cfg$interimmnths, cfg$sero_info_delay)

  lnsero <- rcpp_lnsero(d, nobs)

  m <- matrix(0, ncol = 3, nrow = cfg$post_draw)
  rcpp_immu_interim_post(d, m, nobs, cfg$post_draw, lnsero);

  cfg$looks[look]
  cfg$interimmnths
  nimpute1 <- cfg$looks_target[look] - nobs;

  # onst arma::mat& d,
  # const arma::mat& m,
  # const int look,
  # const int nobs,
  # const int nimpute,
  # const int post_draw,
  # const Rcpp::List& lnsero,
  # const Rcpp::List& cfg

  pp1 <- rcpp_immu_ppos_test(d, m, look, nobs, nimpute1, cfg$post_draw,lnsero, cfg);

  rcpp_immu(d, cfg, 1)

  expect_equal(rcpp_immu(d, cfg, 1)$nimpute1, 10)
  expect_equal(rcpp_immu(d, cfg, 2)$nimpute1,
               cfg$looks_target[2] - cfg$looks[2] + 10 )



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







test_that("immu model - compare post to approx", {

  cfg <- readRDS("cfg-example.RDS")
  cfg$post_draw <- 1000
  cfg$baselineprobsero
  cfg$trtprobsero <- 0.45
  cfg$deltaserot3 <- compute_sero_delta(cfg$baselineprobsero, cfg$trtprobsero)

  d <- rcpp_dat(cfg)
  df <- as.data.frame(d)
  colnames(df) <- dnames

  m <- matrix(0, ncol  = 3, nrow = cfg$post_draw)
  look <- 5
  nobs <- rcpp_n_obs(d, look, cfg$looks, cfg$interimmnths, cfg$sero_info_delay)

  # seroconversion events in each arm
  lnsero <- rcpp_lnsero(d, nobs)
  m <- matrix(0, ncol  = 3, nrow = cfg$post_draw)
  rcpp_immu_interim_post(d, m, nobs, cfg$post_draw, lnsero)
  nimpute1 <- cfg$looks[look] - nobs

  # how many trials are successful? this is what we rely on to do the
  # stop venous sampling check.
  pp1a <- rcpp_immu_interim_ppos(d, m, look, nobs, nimpute1, cfg$post_draw, lnsero, cfg)
  dens1a <- density(pp1a$postprobdelta_gt0)
  pp1b <- rcpp_test_immu_ppos(d, m, look, nobs, nimpute1, cfg$post_draw, lnsero, cfg)
  dens1b <- density(pp1b$postprobdelta_gt0)

  c(pp1a$ppos, pp1b$ppos)

  plot(dens1a)
  lines(dens1b, col = "red")
  #



})


test_that("immu model - main call", {


  cfg <- readRDS("cfg-example.RDS")
  cfg$post_draw <- 1000
  cfg$baselineprobsero
  cfg$trtprobsero <- 0.4
  cfg$deltaserot3 <- compute_sero_delta(cfg$baselineprobsero, cfg$trtprobsero)

  expect_equal(cfg$nmaxsero, 250)


  #rcpp_testing(cfg, 10)

  # generate data
  d <- rcpp_dat(cfg)
  df <- as.data.frame(d)
  colnames(df) <- dnames
  nobs <- 272
  m <- matrix(0, ncol  = 3, nrow = cfg$post_draw)



  # profile
  d <- rcpp_dat(cfg)
  look <- 7
  nobs <- 242
  lnsero <- rcpp_lnsero(d, nobs)
  library(rbenchmark)
  benchmark(rcpp_dat(cfg),
            rcpp_n_obs(d, 7, cfg$looks, cfg$interimmnths, cfg$sero_info_delay),
            rcpp_immu_interim_post(d, m, nobs, cfg$post_draw, lnsero),
            rcpp_test_immu_ppos(d, m, look, nobs, nimpute1, cfg$post_draw, lnsero, cfg),
            rcpp_immu_interim_ppos(d, m, look, nobs, nimpute1, cfg$post_draw, lnsero, cfg),
            columns=c("test", "elapsed", "relative"),
            order="relative", replications=100)

  #


  # simulate a trial - just doing immu endpoint.
  runtrial <- function(nsims){
    wins <- numeric(nsims)
    for(i in 1:nsims){
      # generate data
      d <- rcpp_dat(cfg)
      # doesn't matter that this goes all the way up to looks,
      for(look in 1:length(cfg$looks)){
        if(cfg$looks[look] > cfg$nmaxsero){
          break
        }
        # look minus 1 for c++ world is done inside method
        nobs <- rcpp_n_obs(d, look, cfg$looks, cfg$interimmnths, cfg$sero_info_delay)
        expect_lt(nobs, cfg$looks[look])
        # seroconversion events in each arm
        lnsero <- rcpp_lnsero(d, nobs);
        m <- matrix(0, ncol  = 3, nrow = cfg$post_draw)
        rcpp_immu_interim_post(d, m, nobs, cfg$post_draw, lnsero)
        nimpute1 <- cfg$looks[look] - nobs
        expect_gt(nimpute1, 0)
        # how many trials are successful? this is what we rely on to do the
        # stop venous sampling check.
        pp1 <- rcpp_test_immu_ppos(d, m, look, nobs, nimpute1, cfg$post_draw, lnsero, cfg)
        if(pp1$ppos > cfg$pp_sero_sup_thresh){
          wins[i] <- 1
          break
        }
      }
    }
    return(wins)
  }

  mywin <- runtrial(1000)
  mean(mywin)

  #











  # this is what is used to do the futility check
  nimpute2 <- cfg$nmaxsero - nobs
  pp2 <- rcpp_immu_interim_ppos(d, m, look, nobs, nimpute2, cfg$post_draw, lnsero, cfg)



  m_immu_res <- rcpp_immu(d, cfg, 1);




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





