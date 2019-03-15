
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


cat("Done")