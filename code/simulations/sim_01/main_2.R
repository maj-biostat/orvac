


library(configr)
library(data.table)
library(futile.logger)
library(survival)
library(doParallel)
library(foreach)
library(truncnorm)
library(beepr)
library(optparse)
source("util.R")
source("sim.R")

# parse command line args
option_list <- list(
  make_option(c("-f", "--cfgfile"), type = "character", default = NULL,
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
    stop("If you are going to override you need to specify further
  command line args.\n", call.=FALSE)
  }
  cfg <- sim_cfg(opt$cfgfile, opt)
  
} else {
  
  cfgfile = "cfg1.yaml"
  cfg <- sim_cfg("cfg1.yaml", opt)
}



# dummy dataset
dt1 <- gen_dat(cfg)


# Initiate cluster
cl <- NA
if(!debug){
  cl <- makeCluster(parallel::detectCores() - 2, outfile="")
  registerDoParallel(cl)
} else {
  registerDoSEQ()
}

nworkers <- getDoParWorkers()
flog.info("number of workers %s", nworkers)

# start the clock
start <- proc.time()
packs <- c("data.table", "futile.logger", "configr", "survival", "foreach", "truncnorm", "beepr")


# demon seed
flog.info("Setting random seed: %s.", cfg$seed)
set.seed(cfg$seed)

results <- foreach(i = 1:cfg$nsims,
                   .errorhandling = 'pass',
                   .packages=packs,
                   #.options.snow=opts,
                   .combine = 'rbind'
                   ) %dopar%{

  #
  set.seed(cfg$seed + i)
  # sample size and venous sampling status
  ss_clin <- ss_immu <- 0

  # stop indicators
  stop_ven_samp <- 0
  stop_immu_fut <- 0
  stop_clin_fut <- 0
  stop_clin_sup <- 0
  inconclusive <- 0

  dt1 <- gen_dat(cfg)

  flog.appender(appender.file(file.path(getwd(), "logs", cfg$flog_logfile)), name='ROOT')
  # flog.threshold(DEBUG)

  # individual trial processing is 'embarrassingly parallel' hence the foreach dopar loop
  dotrial <- function(look){

    # idx = i = look = 1
    # idx = i = 2; look = look + 1; look
    # idx = i = look = 5
    # idx = i = look = 8
    # idx = i = 1; look = 28
    # idx = i = 1; look = 33

    m_immu_res <- m_clin_res <- NULL

    # note! do not reset stop_i to zero!!!
    # stop_i_this is a (bit of a) cludge so that we do not get confusing output
    # and can reduce unnecessary sim processing
    trial_state <- list(stop_ven_samp = 0,
                                  stop_immu_fut = 0,
                                  stop_clin_fut = 0,
                                  stop_clin_sup = 0,
                                 inconclusive = 0)


    # Take a copy of the trial data so that we can fiddle with it.
    d <- copy(dt1)

    # how many observations can we see? this is conditional on the accural rate
    interim_month <- cfg$interimmnths[look]
    n_obs <- cfg$looks[look]
    n_max <- cfg$nstop


    # immunological model
    if (!stop_immu(stop_ven_samp,
                  stop_immu_fut,
                  stop_clin_fut,
                  stop_clin_sup,
                  cfg$looks[look],
                  cfg$nmaxsero)){
      m_immu_res <- tryCatch({
        model_immu(d, cfg, look, i)
      }, error = function(err) {
        flog.info("CATCH ERROR model_immu_2 err = %s \n i = %s look = %s", err, i, look)
        flog.info(sys.calls())
        stop("Stopped in main loop model_immu_2 error")
      }, warning=function(cond) {
        flog.info("CATCH WARNING model_immu_2 err = %s \n i = %s look = %s", cond, i, look)
        flog.info(sys.calls())
        stop("Stopped in main loop model_immu_2 warning")
      })

      # rule 1 - "futility test"
      # at the nmaxsero (which I am assuming is the final analysis for this endpoint) ppos_max is NA
      if (!is.na(m_immu_res["ppos_max"]) && m_immu_res["ppos_max"] < cfg$rule1_sero_pp_fut_thresh){
        
        if(look >= cfg$nstartclin && stop_clin_fut){
          flog.info("Both immu and clin futile: ppos_max = %s threshold %s, i = %s look = %s", 
                    m_immu_res["ppos_max"], cfg$rule1_sero_pp_fut_thresh, i, look)
          stop_immu_fut <<- 1
          trial_state$stop_immu_fut <- 1
        } else {
          flog.info("Immu futile: ppos_max = %s threshold %s, i = %s look = %s", 
                    m_immu_res["ppos_max"], cfg$rule1_sero_pp_fut_thresh, i, look)
          stop_immu_fut <<- 1
          trial_state$stop_immu_fut <- 1
        }
        
        
      }

      # rule 3 - "stop venous sampling"
      if (m_immu_res["ppos_n"] > cfg$rule3_sero_pp_sup_thresh &&
          !trial_state$stop_immu_fut){
        stop_ven_samp <<- 1
        trial_state$stop_ven_samp <- 1
      }
      
      ss_immu <- n_obs
    }

    # clinical model
    # start clinical only when we have 250 kids and then every next look

    if (!stop_clin(stop_immu_fut,
                  stop_clin_fut,
                  stop_clin_sup) & cfg$looks[look] >= cfg$nstartclin){

      m_clin_res <- tryCatch({
        model_clin(d, cfg, look, i)
      }, error = function(err) {
        flog.info("CATCH ERROR model_clin err = %s \n i = %s look = %s", err, i, look)
        flog.info(sys.calls())
        stop("Stopped in main loop model_clin error")
      }, warning=function(cond) {
        flog.info("CATCH WARNING model_clin err = %s \n i = %s look = %s", cond, i, look)
        flog.info(sys.calls())
        stop("Stopped in main loop model_clin warning")
      })

      # ppos_max will be NA at max looks since that is the final analysis and we just
      # look at the posterior rather than the predictive probability
      if (!is.na(m_clin_res["ppos_max"]) && m_clin_res["ppos_max"] < cfg$rule1_tte_pp_fut_thresh){
        flog.info("Treatment futile: ppos_max = %s threshold %s, i = %s look = %s", 
                  m_clin_res["ppos_max"], cfg$rule1_tte_pp_fut_thresh, i, look)
        stop_clin_fut <<- 1
        trial_state$stop_clin_fut <- 1
      }

      if (m_clin_res["ppos_n"] > cfg$rule2_tte_pp_sup_thresh[look] &&
          !trial_state$stop_clin_fut){
        flog.info("Treatment superior: ppos_n = %s threshold %s, i = %s look = %s", 
                  m_clin_res["ppos_n"], cfg$rule2_tte_pp_sup_thresh[look], i, look)
        stop_clin_sup <<- 1
        trial_state$stop_clin_sup <- 1
      }
      
      if (look == length(cfg$looks) && stop_clin_fut == 0 && stop_clin_sup == 0){
        inconclusive <<- 1
        trial_state$inconclusive <- 1
      }
      

      # we are not concerned with the zeros.
      ss_clin <- n_obs
    }


    # do various tests on contents of m_immu_res  m_clin_res

    # to ensure that we can populate the return list with something
    if (exists("m_immu_res") & !is.null(m_immu_res)){
      immu_res = m_immu_res
    } else {
      immu_res = rep(NA, length(cfg$immu_rtn_names))
      names(immu_res) <- cfg$immu_rtn_names
    }

    if (exists("m_clin_res") & !is.null(m_clin_res)){
      clin_res = m_clin_res
    } else {
      clin_res = rep(NA, length(cfg$clin_rtn_names))
      names(clin_res) <- cfg$clin_rtn_names
    }

    if (!is.null(warnings())){
      flog.info("main loop warnings - current simulation %s, look %s \n warnings %s", i, look, warnings())
      # clears any warnings.
      assign("last.warning", NULL, envir = baseenv())
    }

    # update control variables

    lr <- tryCatch({
      lr <- c(idxsim = i,
            look = look,
            n_obs = n_obs,
            ss_immu = ss_immu,
            ss_clin = ss_clin,
            n_max = n_max,
            n_max_sero = cfg$nmaxsero,
            immu_res = immu_res,
            clin_res = clin_res,
            stop_v_samp = trial_state$stop_ven_samp,
            stop_i_fut = trial_state$stop_immu_fut,
            stop_c_fut = trial_state$stop_clin_fut,
            stop_c_sup = trial_state$stop_clin_sup,
            inconclusive = trial_state$inconclusive)

        names(lr) <- cfg$field_names
        lr

      }, error = function(err) {
        flog.info("CATCH ERROR err = %s \n
                                 i = %s look = %s n_obs = %s \n
                                 ss_immu = %s \n
                                 ss_clin = %s \n
                                 n_max = %s \n
                                 n_max_sero = %s \n
                                 immu_res = %s \n
                                 clin_res = %s \n
                                 stop_v_samp = %s \n
                                 stop_i_fut = %s \n
                                 stop_c_fut = %s \n
                                 stop_c_sup = %s \n
                                  inconclusive = %s", 
                  err, i, look, n_obs, 
                  ss_immu, ss_clin, 
                  n_max, cfg$nmaxsero, 
                  paste0(immu_res, collapse = ", "), 
                  paste0(clin_res, collapse = ", "), 
                  trial_state$stop_ven_samp, 
                  trial_state$stop_immu_fut, 
                  trial_state$stop_clin_fut, 
                  trial_state$stop_clin_sup,
                  trial_state$inconclusive)
        flog.info(sys.calls())
        lrerr <- rep(NA, length(cfg$field_names))
        return(lrerr)
      }, warning=function(cond) {
        flog.info("CATCH WARNING err = %s \n
                                 i = %s look = %s n_obs = %s \n
                                 ss_immu = %s \n
                                 ss_clin = %s \n
                                 n_max = %s \n
                                 n_max_sero = %s \n
                                 immu_res = %s \n
                                 clin_res = %s \n
                                 stop_v_samp = %s \n
                                 stop_i_fut = %s \n
                                 stop_c_fut = %s \n
                                 stop_c_sup = %s \n
                                  inconclusive = %s", 
                  cond, i, look, n_obs, 
                  ss_immu, ss_clin, 
                  n_max, cfg$nmaxsero, 
                  paste0(immu_res, collapse = ", "), 
                  paste0(clin_res, collapse = ", "), 
                  trial_state$stop_ven_samp, 
                  trial_state$stop_immu_fut, 
                  trial_state$stop_clin_fut, 
                  trial_state$stop_clin_sup,
                  trial_state$inconclusive)
        flog.info(sys.calls())
        lrerr <- rep(NA, length(cfg$field_names))
        return(lrerr)
      })

    return(lr)
  }

  res <- do.call(rbind, lapply(1:cfg$nlooks, dotrial))

  return(res)
}


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
