

library(ggplot2)
library(dplyr)
library(configr)
library(data.table)
library(futile.logger)
library(survival)
library(doParallel)
library(foreach)
library(truncnorm)
library(beepr)
library(optparse)
# c++ bits
library(orvacsim)

source("util.R")
source("sim.R")

ggplot2::theme_set(theme_bw())

# parse command line args
option_list <- list(
  make_option(c("-f", "--cfgfile"), type = "character", default = "cfg1.yaml",
              help = "config file name", metavar = "character"),
  make_option(c("-o", "--use"), type = "logical", default = FALSE,
              help = "override config file with command line settings",
              metavar = "logical"),
  make_option(c("-i", "--idsim"), type = "character", default = FALSE,
              help = "label for current simulation", metavar = "character"),
  make_option(c("-l", "--logfile"), type = "character", default = NULL,
              help = "log file name", metavar = "character"),
  make_option(c("-n", "--nsims"), type = "integer", default = NULL,
              help = "number of simulations", metavar = "integer"),
  make_option(c("-s", "--seed"), type = "integer", default = NULL,
              help = "random seed", metavar = "integer"),
  make_option(c("-a", "--accrual"), type = "integer", default = NULL,
              help = "accrual rate i.e. people_per_interim_period",
              metavar = "integer"),
  make_option(c("-d", "--delay"), type = "double", default = NULL,
              help = "seroconversion information delay",
              metavar = "double"),
  make_option(c("-b", "--basesero"), type = "double", default = NULL,
              help = "baseline seroconversion prob", metavar = "double"),
  make_option(c("-p", "--trtprobsero"), type = "double", default = NULL,
              help = "trt seroconversion prob", metavar = "double"),
  make_option(c("-m", "--basemediantte"), type = "double", default = NULL,
              help = "baseline median time to med attendance (months)",
              metavar = "double"),
  make_option(c("-t", "--trtmedtte"), type = "double", default = NULL,
              help = "treatment arm median time to med attendance (months)",
              metavar = "double")
);

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

if (length(opt) < 3){
  print_help(opt_parser)
  stop("At a minimum you need to specify the config file.\n", call. = FALSE)
}

debug = F
if(!debug){
  
  if (opt$use == T && length(opt) <= 3){
    print_help(opt_parser)
    stop("If you are going to override you need to specify further command line args.\n", call.=FALSE)
  }
  cfg <- sim_cfg(opt$cfgfile, opt)
  
} else {
  cfgfile = "cfg1.yaml"
  cfg <- sim_cfg("cfg1.yaml", opt)
}


# Initiate cluster
cl <- NA
if(!debug){
  cl <- makeCluster(parallel::detectCores() - 2, outfile="")
  # cl <- makeCluster(3, outfile="")
  registerDoParallel(cl)
} else {
  registerDoSEQ()
}

nworkers <- getDoParWorkers()
flog.info("number of workers %s", nworkers)

# start the clock
start <- proc.time()
packs <- c("data.table", "futile.logger", "configr", "survival", 
           "foreach", "truncnorm", "beepr", "orvacsim")


# demon seed
flog.info("Setting random seed: %s.", cfg$seed)
set.seed(cfg$seed)

results <- foreach(i = 1:cfg$nsims,
                   .errorhandling = 'pass',
                   .packages=packs,
                   #.options.snow=opts,
                   .combine = 'rbind'
                   ) %dopar%{

                     
                     
                     
  # i = 1
  #flog.info("Starting trial: i = %s", i)
  set.seed(cfg$seed + i)
  # sample size and venous sampling status
  ss_clin <- 0
  ss_immu <- 0

  # stop indicators
  stop_ven_samp <- 0
  stop_immu_fut <- 0
  stop_clin_fut <- 0
  stop_clin_sup <- 0
  inconclusive <- 0

  dt1 <- rcpp_dat(cfg)
  
  # Precedence:
  # TRACE, DEBUG, INFO, WARN, ERROR, FATAL
  flog.appender(appender.file(file.path(getwd(), "logs", cfg$flog_logfile)), name='ROOT')
  flog.threshold(DEBUG)
  

  # individual trial processing is 'embarrassingly parallel' hence the foreach dopar loop
  dotrial <- function(look){
    
    #flog.info("Started look, trial sim = %s, look = %s", i, look)

    # idx = i = look = 1
    # idx = i = 2; look = look + 1; look
    # idx = i = look = 5
    # idx = i = look = 10
    # idx = i = 1; look = 28
    # idx = i = 1; look = 33
    # look = look + 1

    m_immu_res <- NULL
    m_clin_res <- NULL

    # note! do not reset stop_i to zero!!!
    # stop_i_this is a (bit of a) cludge so that we do not get confusing output
    # and can reduce unnecessary sim processing
    trial_state <- list(stop_ven_samp = 0,
                        stop_immu_fut = 0,
                        stop_clin_fut = 0,
                        stop_clin_sup = 0,
                        inconclusive = 0)

    # immunological model
    if (!stop_immu(stop_ven_samp,
                  stop_immu_fut,
                  stop_clin_fut,
                  stop_clin_sup,
                  cfg$looks[look],
                  cfg$nmaxsero)){

      m_immu_res <- tryCatch({
        rcpp_immu(dt1, cfg, look)
        # m_immu_res <- rcpp_immu(dt1, cfg, look)
        # model_immu(d, cfg, look, i)
      }, error = function(err) {
        flog.fatal("CATCH ERROR model_immu_2 err = %s \n sim = %s look = %s", err, i, look)
        flog.fatal(sys.calls())
        stop("Stopped in main loop model_immu_2 error")
      }, warning=function(cond) {
        flog.fatal("CATCH WARNING model_immu_2 err = %s \n sim = %s look = %s", cond, i, look)
        flog.fatal(sys.calls())
        stop("Stopped in main loop model_immu_2 warning")
      })

      # rule 1 - "futility test"
      # we do not make superiority assessments based on the immu endpoint. 
      # analysis for making decision is the same in the interim and at max sero
      # namely we look at the ppos based on simulated trials. 
      if (m_immu_res$ppos_max < cfg$rule1_sero_pp_fut_thresh){
       flog.info("Immunological ep futile: ppos_max = %s threshold %s, sim = %s look = %s ", 
                 m_immu_res$ppos_max, cfg$rule1_sero_pp_fut_thresh, i, look)
        stop_immu_fut <<- 1
        trial_state$stop_immu_fut <- 1
      }

      # rule 3 - "stop venous sampling"
      if (m_immu_res$ppos_n > cfg$rule3_sero_pp_sup_thresh && !trial_state$stop_immu_fut){
        flog.info("Immunological ep stopped sampling: ppos_n = %s threshold %s, sim = %s look = %s", 
                  m_immu_res$ppos_n, cfg$rule3_sero_pp_sup_thresh, i, look)
        stop_ven_samp <<- 1
        trial_state$stop_ven_samp <- 1
      }
      
      ss_immu <- cfg$looks[look]
    }

    
    
    if (!stop_clin(stop_immu_fut,
                   stop_clin_fut,
                   stop_clin_sup) & cfg$looks[look] >= cfg$nstartclin){
      
      m_clin_res <- tryCatch({
        rcpp_clin(dt1, cfg, look)
      }, error = function(err) {
        flog.fatal("CATCH ERROR model_clin err = %s \n sim = %s look = %s", err, i, look)
        flog.fatal(sys.calls())
        stop("Stopped in main loop model_clin error")
      }, warning=function(cond) {
        flog.fatal("CATCH WARNING model_clin err = %s \n sim = %s look = %s", cond, i, look)
        flog.fatal(sys.calls())
        stop("Stopped in main loop model_clin warning")
      })
      
      # ppmax will be NA at max looks since that is the final analysis and we just
      # look at the posterior rather than the predictive probability
      if (!is.na(m_clin_res$ppmax) && m_clin_res$ppmax < cfg$rule1_tte_pp_fut_thresh){
        flog.info("Clin ep futile: ppmax = %s threshold %s (pgt1 = %.3f l0 = %.3f l1 = %.3f ratio = %.3f n_uncen_0 = %s n_uncen_1 = %s), sim = %s look = %s", 
                  m_clin_res$ppmax, cfg$rule1_tte_pp_fut_thresh, 
                  m_clin_res$pgt1, m_clin_res$l0, m_clin_res$l1,  m_clin_res$ratio,
                  m_clin_res$n_uncen_0, m_clin_res$n_uncen_1,
                  i, look)
        stop_clin_fut <<- 1
        trial_state$stop_clin_fut <- 1
      }
      
      if (m_clin_res$ppn > cfg$rule2_tte_pp_sup_thresh[look] &&
          !trial_state$stop_clin_fut){
        flog.info("Clin ep superior: ppn = %s threshold %s (pgt1 = %.3f l0 = %.3f l1 = %.3f ratio = %.3f n_uncen_0 = %s n_uncen_1 = %s), sim = %s look = %s", 
                  m_clin_res$ppn, cfg$rule2_tte_pp_sup_thresh[look], 
                  m_clin_res$pgt1, m_clin_res$l0, m_clin_res$l1,  m_clin_res$ratio, 
                  m_clin_res$n_uncen_0, m_clin_res$n_uncen_1,
                  i, look)
        stop_clin_sup <<- 1
        trial_state$stop_clin_sup <- 1
      }
      
      if (look == length(cfg$looks) && stop_clin_fut == 0 && stop_clin_sup == 0){
        flog.info("Clin ep inconclusive: ppmax = %s threshold %s (pgt1 = %.3f l0 = %.3f l1 = %.3f ratio = %.3f n_uncen_0 = %s n_uncen_1 = %s), sim = %s look = %s", 
                  m_clin_res$ppmax, cfg$rule1_tte_pp_fut_thresh, 
                  m_clin_res$pgt1, m_clin_res$l0, m_clin_res$l1,  m_clin_res$ratio,
                  m_clin_res$n_uncen_0, m_clin_res$n_uncen_1,
                  i, look)
        inconclusive <<- 1
        trial_state$inconclusive <- 1
      }
      
      
      # we are not concerned with the zeros.
      ss_clin <- cfg$looks[look]
    }
    
    

    if (!is.null(warnings())){
      flog.warn("main loop warnings - current simulation %s, look %s \n warnings %s", i, look, warnings())
      # clears any warnings.
      assign("last.warning", NULL, envir = baseenv())
    }

    # update control variables

    lr <- tryCatch({
      
      # ensures all columns always present
      lr <- list(idxsim = i,
                 look = look,
                 n_obs = cfg$looks[look],
                 ss_immu = ss_immu,
                 n_max = cfg$nstop,
                 n_max_sero = cfg$nmaxsero,
                 stop_v_samp = trial_state$stop_ven_samp,
                 stop_i_fut = trial_state$stop_immu_fut,
                 stop_c_fut = trial_state$stop_clin_fut,
                 stop_c_sup = trial_state$stop_clin_sup,
                 inconclusive = trial_state$inconclusive,
                 i_pposn = NA,
                 i_pposmax = NA,
                 c_ppn = NA,
                 c_ppmax = NA,
                 c_ppmax_mean_ratio = NA,
                 c_ppmax_sd_ratio = NA
                 )
      
      if (exists("m_immu_res") && !is.null(m_immu_res)){
        lr$i_pposn <- m_immu_res$ppos_n
        lr$i_pposmax <- m_immu_res$ppos_max
      }
      
      if (exists("m_clin_res") && !is.null(m_clin_res)){
        lr$c_ppn <- m_clin_res$ppn
        lr$c_ppmax <- m_clin_res$ppmax
        lr$c_ppmax_mean_ratio <- m_clin_res$ppmax_mean_ratio
        lr$c_ppmax_sd_ratio <- m_clin_res$ppmax_sd_ratio
      }
      
      lr

      }, error = function(err) {
        flog.fatal("CATCH ERROR err = %s", err)
        flog.fatal(sys.calls())
        lrerr <- rep(NA, length(cfg$field_names))
        return(lrerr)
      }, warning=function(cond) {
        flog.fatal("CATCH WARNING err = %s", cond)
        flog.fatal(sys.calls())
        lrerr <- rep(NA, length(cfg$field_names))
        return(lrerr)
      })
    
    if(i%%100 == 0){
      flog.info("Completed look, trial: sim = %s, look = %s", i, look)
    }
    

    return(unlist(lr))
  }

  res <- do.call(rbind, lapply(1:cfg$nlooks, dotrial))
  
  # flog.info("Finished trial: sim = %s", i)
  return(res)
}

results <- as.data.frame(results)
end <- proc.time()

duration <- end - start
flog.info("proc.time:" )
flog.info("user    %s", round(duration[1], 2) )
flog.info("system  %s", round(duration[2], 2) )
flog.info("elapsed %s", round(duration[3], 2) )

beepr::beep()

w <- warnings()

rdsfilename <- paste0("out/res-",format(Sys.time(), "%Y-%m-%d-%H-%M-%S"), ".RDS")
flog.info("saving rdsfilename : %s", rdsfilename )
saveRDS(list(results=results, cfg = cfg, warnings = w), rdsfilename)
assign("last.warning", NULL, envir = baseenv())

if(!debug){
  stopCluster(cl)
}

flog.info("Done. Cluster stopped. Exiting now." )
