


library(configr)
library(data.table)
library(futile.logger)
library(survival)
library(doParallel)
# library(doSNOW)
# library(doMC)
library(foreach)
library(truncnorm)
#library(gibby)
source("util.R")
source("sim.R")


# Load config

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  cfgfile = "cfg1.yaml"
} else if (length(args)==1) {
  # default output file
  cfgfile = args[1]
}

cfg <- sim_cfg(cfgfile)

# dummy dataset
dt1 = gen_dat(cfg)


interval_display <- 5




# Initiate cluster
cl <- makeCluster(parallel::detectCores() - 2, outfile="")
registerDoParallel(cl)
# registerDoSNOW(cl)
# registerDoMC(parallel::detectCores() - 2)

# ntasks <- 100
# pb <- tkProgressBar(max=ntasks)
# progress <- function(n) setTkProgressBar(pb, n)
# opts <- list(progress=progress)

start <- proc.time()

# parallelised  .combine = rbind,

packs <- c("data.table", "futile.logger", "configr", "survival", "foreach", "truncnorm")



set.seed(666)

results <- foreach(i = 1:cfg$nsims,
                   .errorhandling = 'pass',
                   .packages=packs,
                   #.options.snow=opts,
                   .combine = 'rbind'
                   ) %dopar%{

  # sample size and venous sampling status
  ss_clin = 0
  ss_imm = 0
  stop_ven_samp = 0
  dt1 = gen_dat(cfg)
  
  # dt1 = gen_dat(cfg); d = copy(dt1)
  
  flog.appender(appender.file(cfg$flog_logfile), name='ROOT')
  
   
  foreach(look = 1:cfg$nlooks,
          .combine = 'rbind'
          ) %dopar% {
    
    
    # easiest way to accomodate the 'whichever earlier' rule: 
    # loop over every day up to iterim[loo]. 
    # check n_obs and month of trial at each iteration. 
    # if n_obs hits a interim then do modelling.
    

    # idx = i = look = 1
    # idx = i = look = 5
    # idx = i = look = 6
    # idx = i = 1; look = 15

    lr1 <- NULL
    m_immu_res <- m_clin_res <- NULL
    
    # idisplay <- ifelse(i %% interval_display == 0 && look == 1, 1, 0)

    # Take a copy of the trial data so that we can fiddle with it.
    d = copy(dt1)
    
    # how many observations can we see? this is conditional on the accural rate
    interim_month <- cfg$interimmnths[look]
    n_obs <- which(d$accrt == max(d$accrt[d$accrt < interim_month]))
    
    # what is the maximum number of obs we will see if accrual remains at the current rate
    n_max <- min(max(cfg$looks), floor(  max(cfg$interimmnths) / cfg$months_per_person  )) 
    # paranoid android:
    n_max <- ifelse(n_max > max(cfg$looks), max(cfg$looks), n_max) 

    # immunological model
    # Are we still venous sampling?
    if(stop_ven_samp != 1 & cfg$looks[look] <= cfg$nmaxsero){
      
      m_immu_res <- tryCatch({
        model_immu_2(d, cfg, look, n_obs, n_max, i)
      }, error = function(err) {
        print(paste("model_immu_2 err:  ",err))
        print(sys.calls())
        flog.info("CATCH ERROR model_immu_2 i = %s look = %s", i, look)
        flog.info(sys.calls())
      })
    }
    
    # clinical model 
    
    # start clinical only when we have 250 kids and then every 3 months afterwards
    
    if(cfg$looks[look] >= cfg$nstartclin){
      
      m_clin_res <- tryCatch({
        model_clin_1(d, cfg, look, n_obs, n_max, i)
      }, error = function(err) {
        print(paste("model_clin_1 err:  ",err))
        print(sys.calls())
        flog.info("CATCH ERROR model_clin_1 i = %s look = %s", i, look)
        flog.info(sys.calls())
      })
    }
      
    # do various tests on contents of m_immu_res  m_clin_res
    
    
    
    
    # to ensure that we can populate the return list with something
    if(exists("m_immu_res") & !is.null(m_immu_res)){
      immu_res = m_immu_res
    } else {
      immu_res = c(NA, NA)
    }
    
    if(exists("m_clin_res") & !is.null(m_clin_res)){
      clin_res = m_clin_res
    } else {
      clin_res = c(NA, NA)
    }

    # update control variables
    lr1 <- c(idxsim = i,
             look = look, 
             n_obs = n_obs,
             n_max = n_max,
             m_immu_res = immu_res,
             m_clin_res = clin_res)
    
    names(lr1) <- c("simidx", "look", 
                    "n_obs", "n_max", 
                    "immu_ppos_n", "immu_ppos_max", 
                    "clin_ppos_n", "clin_ppos_max")
    
    return(lr1)

  }

}




end <- proc.time() 

duration <- end - start
flog.info("proc.time:" )
flog.info("user    %s", round(duration[1], 2) )
flog.info("system  %s", round(duration[2], 2) )
flog.info("elapsed %s", round(duration[3], 2) )


# 
# ppos_sup <- results[,1]
# 
# 
# mycombine <- function(x){
#   
# }
# 
# mtest <- matrix(unlist(results), ncol = 2, byrow = T)
# class(mtest)
# test <-colMeans(mtest)
# names(test) <- c("ppos_at_look_immu", "ppos_at_max_immu")


stopCluster(cl)


# results <- mclapply( 1:cfg$nsims, FUN=mcmc2_pp, d = d, cfg = cfg )
# 
# 
# 
# flog.info("main.R finishing NOW!!! %s")
# flog.info("time taken %s", difftime(Sys.time() - start))
# 
# 
# saveRDS(warnings(), "warn.RDS")
# saveRDS(results, "results.RDS")
# 
# dfresults <- do.call(rbind, lapply(results, function(x){return(x$result)}))
# dfresults
# 
# 
# saveRDS(list(duration = duration, start = start, end = end), "duration.RDS")
# saveRDS(dfresults, "dfresults.RDS")


# cmcmc$run(5000)
# n_samp <- as.matrix(cmcmc$mvSamples)
# t(apply(as.matrix(n_samp), 2, quantile, probs = c(0.025, 0.5, 0.975)))

# This is the pattern
# cnmod$resetData()
# dt1 = gen_dat(cfg)
# jdat <- jags_init(d = copy(dt1), omit_const = T)
# cnmod$setData(jdat$dat)
# cnmod$setInits(jdat$inits)
# cmcmc$run(5000)
# n_samp <- as.matrix(cmcmc$mvSamples)
# t(apply(as.matrix(n_samp), 2, quantile, probs = c(0.025, 0.5, 0.975)))

# sim <- function(a, b, c) return(c(sim = a,
#                                   look = b,
#                                   value = c))

# results <- foreach(i = 1:cfg$nsims,
#                    .errorhandling = 'pass',
#                    .packages=packs,
#                    .combine = 'rbind') %dopar%{
#                      
#   thisisatest = 1
#   
#   foreach(look = 1:cfg$nlooks,
#          .combine = 'rbind') %dopar% {
#            
#            sim(i, look, rnorm(1))
#          }
#                                   
#                      
# }
#   
# results
