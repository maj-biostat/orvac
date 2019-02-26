

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

# parse command line args/config
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


# initiate clusters
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







# main parallel loop
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
  
  # Precedence:
  # TRACE, DEBUG, INFO, WARN, ERROR, FATAL
  flog.appender(appender.file(file.path(getwd(), "logs", cfg$flog_logfile)), name='ROOT')
  flog.threshold(DEBUG)
  
  res <- rcpp_dotrial(i, cfg)

  # lr <- list(idxsim = NA,
  #             look = NA,
  #             ss_immu = NA,
  #             ss_clin = NA,
  #             stop_v_samp = NA,
  #             stop_i_fut = NA,
  #             stop_c_fut = NA,
  #             stop_c_sup = NA,
  #             i_final = NA,
  #             c_final = NA,
  #             i_ppn = NA,
  #             i_ppmax = NA,
  #             c_ppn = NA,
  #             c_ppmax = NA,
  #             i_mean = NA,
  #             i_lwr = NA,
  #             i_upr = NA,
  #             c_mean = NA,
  #             c_lwr = NA,
  #             c_upr = NA)
  # 
  # lr$idxsim <- res$idxsim
  # lr$look <- res$look
  # lr$ss_immu <- res$ss_immu
  # lr$ss_clin <- res$ss_clin
  # lr$stop_v_samp <- res$stop_v_samp
  # lr$stop_i_fut <- res$stop_i_fut
  # lr$stop_c_fut <- res$stop_c_fut
  # lr$stop_c_sup <- res$stop_c_sup
  # lr$i_final <- res$i_final
  # lr$c_final <- res$c_final
  # lr$i_ppn <- res$i_ppn
  # lr$i_ppmax <- res$i_ppmax
  # lr$c_ppn <- res$c_ppn
  # lr$c_ppmax <- res$c_ppmax
  # lr$i_mean <- res$i_mean
  # lr$i_lwr <- res$i_lwr
  # lr$i_upr <- res$i_upr
  # lr$c_mean <- res$c_mean
  # lr$c_lwr <- res$c_lwr
  # lr$c_upr <- res$c_upr
  

  # flog.info("Finished trial: sim = %s", i)
  return(unlist(res))
}









# save results to file
results <- as.data.frame(results)
rownames(results) <- NULL
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
