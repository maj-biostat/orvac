

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


starttime <- Sys.time()




# main parallel loop
flog.info("Setting random seed: %s.", cfg$seed)
set.seed(cfg$seed)

results <- foreach(i = 1:cfg$nsims,
                   .errorhandling = 'pass',
                   .packages=packs
                   #.options.snow=opts,
                   ) %dopar%{
    
  # i = 1
  #flog.info("Starting trial: i = %s", i)
  set.seed(cfg$seed + i)
  
  # Precedence:
  # TRACE, DEBUG, INFO, WARN, ERROR, FATAL
  flog.appender(appender.file(file.path(getwd(), "logs", cfg$flog_logfile)), name='ROOT')
  flog.threshold(DEBUG)
  
  res <- rcpp_dotrial(i, cfg, FALSE)

  # flog.info("Finished trial: sim = %s", i)
  return(res)
}





# save results to file

dfres1 <- data.frame()
dfres2 <- data.frame()

for(i in 1:length(results)){

  myv <- unlist(results[[i]][1:25])
  nm <- names(myv)
  dfres1 <- rbind(dfres1, myv)
  colnames(dfres1) <- nm

}

endtime <- Sys.time()

end <- proc.time()

duration <- end - start
flog.info("proc.time:" )
flog.info("user    %s", round(duration[1], 2) )
flog.info("system  %s", round(duration[2], 2) )
flog.info("elapsed %s", round(duration[3], 2) )
beepr::beep()
w <- warnings()
# rdsfilename <- paste0("out/list-",format(Sys.time(), "%Y-%m-%d-%H-%M-%S"), ".RDS")
# flog.info("saving rdsfilename : %s", rdsfilename )
# saveRDS(results, rdsfilename)
rdsfilename <- paste0("out/res-",format(Sys.time(), "%Y-%m-%d-%H-%M-%S"), ".RDS")
flog.info("saving rdsfilename : %s", rdsfilename )
saveRDS(list(results=dfres1, cfg = cfg, warnings = w,
             starttime = starttime, endtime = endtime,
             duration = difftime(endtime, starttime, units = "hours")),
        rdsfilename)
assign("last.warning", NULL, envir = baseenv())

if(!debug){
  stopCluster(cl)
}

flog.info("Done. Cluster stopped. Exiting now." )


