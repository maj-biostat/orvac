


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
              help = "accrual rate - people_per_interim_period",
              metavar = "integer"),
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

if (opt$use == T && length(opt) <= 3){
  print_help(opt_parser)
  stop("If you are going to override you need to specify further
  command line args.\n", call.=FALSE)
}

cfg <- sim_cfg(opt$cfgfile, opt)

# cfgfile = "cfg1.yaml"
# cfg <- sim_cfg("cfg1.yaml", opt)
# cfg$seed <- 4565


# dummy dataset
dt1 <- gen_dat(cfg)


# Initiate cluster
cl <- makeCluster(parallel::detectCores() - 2, outfile="")
registerDoParallel(cl)
# to move back to sequential
# registerDoSEQ()
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

  # sample size and venous sampling status
  ss_clin <- ss_immu <- 0

  # stop indicators
  stop_ven_samp <- 0
  stop_immu_fut <- 0
  stop_clin_fut <- 0
  stop_clin_sup <- 0


  dt1 <- gen_dat(cfg)

  flog.appender(appender.file(file.path(getwd(), "logs", cfg$flog_logfile)), name='ROOT')
  # flog.threshold(DEBUG)

  # individual trial processing is 'embarrassingly parallel' hence the foreach dopar loop
  dotrial <- function(look){

    # idx = i = look = 1
    # idx = i = 2; look = look + 1; look
    # idx = i = look = 5
    # idx = i = look = 6
    # idx = i = 1; look = 33

    m_immu_res <- m_clin_res <- NULL

    # note! do not reset stop_i to zero!!!
    # stop_i_this is a (bit of a) cludge so that we do not get confusing output
    # and can reduce unnecessary sim processing
    stop_i_this_iteration <- list(stop_ven_samp = 0,
                                  stop_immu_fut = 0,
                                  stop_clin_fut = 0,
                                  stop_clin_sup = 0)


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
        model_immu_2(d, cfg, look, i)
      }, error = function(err) {
        flog.info("CATCH ERROR model_immu_2 err = %s", err)
        flog.info("CATCH ERROR model_immu_2 i = %s look = %s", i, look)
        flog.info("CATCH ERROR model_immu_2 sys.calls follow")
        flog.info(sys.calls())
        stop("Stopped in main loop model_immu_2 error")
      }, warning=function(cond) {
        flog.info("CATCH WARNING model_immu_2 err = %s", cond)
        flog.info("CATCH WARNING model_immu_2 i = %s look = %s", i, look)
        flog.info("CATCH WARNING model_immu_2 sys.calls follow")
        flog.info(sys.calls())
        stop("Stopped in main loop model_immu_2 warning")
      })

      # rule 1 - "futility test"
      if (m_immu_res["ppos_max"] < cfg$rule1_sero_pp_fut_thresh){
        # ok, so I will (almost certainly) burn in hell for this, but it should be safe because of
        # the way that things are parallelised...
        stop_immu_fut <<- 1
        stop_i_this_iteration$stop_immu_fut <- 1
      }

      # rule 3 - "stop venous sampling"
      if (m_immu_res["ppos_n"] > cfg$rule3_sero_ppos_thresh){
        # ok, so I will (almost certainly) burn in hell for this, but it should be safe because of
        # the way that things are parallelised...
        stop_ven_samp <<- 1
        stop_i_this_iteration$stop_ven_samp <- 1
      }

      ss_immu <- n_obs
    }

    # clinical model
    # start clinical only when we have 250 kids and then every next look

    if (!stop_clin(stop_immu_fut,
                  stop_clin_fut,
                  stop_clin_sup) & cfg$looks[look] >= cfg$nstartclin){

      m_clin_res <- tryCatch({
        model_clin_1(d, cfg, look, i)
      }, error = function(err) {
        flog.info("CATCH ERROR model_clin_1 err = %s", err)
        flog.info("CATCH ERROR model_clin_1 i = %s look = %s", i, look)
        flog.info("CATCH ERROR model_clin_1 sys.calls follow")
        flog.info(sys.calls())
        stop("Stopped in main loop model_clin_1 error")
      }, warning=function(cond) {
        flog.info("CATCH WARNING model_clin_1 err = %s", cond)
        flog.info("CATCH WARNING model_clin_1 i = %s look = %s", i, look)
        flog.info("CATCH WARNING model_clin_1 sys.calls follow")
        flog.info(sys.calls())
        stop("Stopped in main loop model_clin_1 warning")
      })

      if (m_clin_res["ppos_max"] < cfg$rule1_sero_pp_fut_thresh){
        # ok, so I will (almost certainly) burn in hell for this, but it should be safe because of
        # the way that things are parallelised...
        stop_clin_fut <<- 1
        stop_i_this_iteration$stop_clin_fut <- 1
      }

      if (m_clin_res["ppos_n"] > cfg$rule2_tte_postthresh){
        # ok, so I will (almost certainly) burn in hell for this, but it should be safe because of
        # the way that things are parallelised...
        stop_clin_sup <<- 1
        stop_i_this_iteration$stop_clin_sup <- 1
      }
      # we are not concerned with the zeros.
      ss_clin <- n_obs
    }


    # do various tests on contents of m_immu_res  m_clin_res

    # to ensure that we can populate the return list with something
    if (exists("m_immu_res") & !is.null(m_immu_res)){
      immu_res = m_immu_res
    } else {
      immu_res = rep(NA, 4)
    }

    if (exists("m_clin_res") & !is.null(m_clin_res)){
      clin_res = m_clin_res
    } else {
      clin_res = rep(NA, 4)
    }

    if (!is.null(warnings())){
      flog.info("main loop warnings - current simulation %s, look %s", i, look)
      flog.info("Warnings: %s.", warnings())
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
            stop_v_samp = stop_i_this_iteration$stop_ven_samp,
            stop_i_fut = stop_i_this_iteration$stop_immu_fut,
            stop_c_fut = stop_i_this_iteration$stop_clin_fut,
            stop_c_sup = stop_i_this_iteration$stop_clin_sup)

        names(lr) <- c("idxsim",
                    "look",
                   "n_obs",
                    "ss_immu",
                    "ss_clin",
                    "n_max",
                    "n_max_sero",
                    "i_ppos_n", "i_ppos_max", "i_post_prop_ctl", "i_post_prop_trt",   # immunological
                    "c_ppos_n","c_ppos_max","c_post_lambda_ctl","c_post_lambda_trt",  # clinical
                    "stop_v_samp",
                    "stop_i_fut",
                    "stop_c_fut",
                    "stop_c_sup")
        lr

      }, error = function(err) {
        flog.info("CATCH ERROR err = %s", err)
        flog.info("CATCH ERROR i = %s look = %s", i, look)
        flog.info("CATCH ERROR n_obs = %s", n_obs)
        flog.info("CATCH ERROR ss_immu = %s", ss_immu)
        flog.info("CATCH ERROR ss_clin = %s", ss_clin)
        flog.info("CATCH ERROR n_max = %s", n_max)
        flog.info("CATCH ERROR n_max_sero = %s", cfg$nmaxsero)
        flog.info("CATCH ERROR immu_res = %s", paste0(immu_res, collapse = ", "))
        flog.info("CATCH ERROR clin_res = %s", paste0(clin_res, collapse = ", "))
        flog.info("CATCH ERROR stop_v_samp = %s", stop_i_this_iteration$stop_ven_samp)
        flog.info("CATCH ERROR stop_i_fut = %s", stop_i_this_iteration$stop_immu_fut)
        flog.info("CATCH ERROR stop_c_fut = %s", stop_i_this_iteration$stop_clin_fut)
        flog.info("CATCH ERROR stop_c_sup = %s", stop_i_this_iteration$stop_clin_sup)
        flog.info("CATCH ERROR system calls follow now:")
        flog.info(sys.calls())
        lrerr <- rep(NA, 19)
        return(lrerr)
      }, warning=function(cond) {
        flog.info("CATCH WARNING cond = %s", cond)
        flog.info("CATCH WARNING i = %s look = %s", i, look)
        flog.info("CATCH WARNING n_obs = %s", n_obs)
        flog.info("CATCH WARNING ss_immu = %s", ss_immu)
        flog.info("CATCH WARNING ss_clin = %s", ss_clin)
        flog.info("CATCH WARNING n_max = %s", n_max)
        flog.info("CATCH WARNING n_max_sero = %s", cfg$nmaxsero)
        flog.info("CATCH WARNING immu_res = %s", paste0(immu_res, collapse = ", "))
        flog.info("CATCH WARNING clin_res = %s", paste0(clin_res, collapse = ", "))
        flog.info("CATCH WARNING stop_v_samp = %s", stop_i_this_iteration$stop_ven_samp)
        flog.info("CATCH WARNING stop_i_fut = %s", stop_i_this_iteration$stop_immu_fut)
        flog.info("CATCH WARNING stop_c_fut = %s", stop_i_this_iteration$stop_clin_fut)
        flog.info("CATCH WARNING stop_c_sup = %s", stop_i_this_iteration$stop_clin_sup)
        flog.info("CATCH WARNING system calls follow now:")
        flog.info(sys.calls())
        lrerr <- rep(NA, 19)
        return(lrerr)
      })


    # round(lr, 4)

    # at the moment creating the whole trial then post processing but could
    # terminate early if using a while or for structure, i.e.
    # if(stop_immu_fut | stop_clin_fut | stop_clin_sup){
    #   return(lr1)
    # }

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



dothis <- FALSE
if (dothis){
  head(print_immu_res(results), 25)
  tail(print_immu_res(results), 10)
  head(print_clin_res(results), 25)
  tail(print_clin_res(results), 10)
}

w <- warnings()

rdsfilename <- paste0("out/res-",format(Sys.time(), "%Y-%m-%d-%H-%M-%S"), ".RDS")
flog.info("saving rdsfilename : %s", rdsfilename )
saveRDS(list(results=results, cfg = cfg, warnings = w), rdsfilename)
assign("last.warning", NULL, envir = baseenv())

stopCluster(cl)
flog.info("Done. Cluster stopped. Exiting now." )
