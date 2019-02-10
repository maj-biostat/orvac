

library(doParallel)
library(foreach)
library(data.table)

test_func <- function(single_clust = T, nsim = 100){
  
  cl <- NA
  if(single_clust){
    registerDoSEQ()
  } else {
    cl <- makeCluster(parallel::detectCores() - 2, outfile="")
    registerDoParallel(cl)
  }
 
  n <- 10000
  looks <- 10
  beta0 <- 10
  beta1 <- 5
  e_sd <- 2
  
  # start clock
  start <- proc.time()
  
  
  # Parallel loop
  packs <- c("data.table", "foreach")
  results <- foreach(i = 1:nsim,
                     .errorhandling = 'pass',
                     .packages=packs,
                     #.options.snow=opts,
                     .combine = 'rbind'
  ) %dopar%{
    
    dotrial <- function(look){
      
      ninterim <- (n/looks) + (look-1)*(n/looks)
      
      dt <- data.table(trt = rep(0:1, each = ninterim/2))
      dt[, y:= beta0 + beta1 * trt + rnorm(ninterim, 0, e_sd)]
      
      lm1 <- lm(y ~ trt, data = dt)
      
      return(confint(lm1)["trt", ])
    }
    
    res <- do.call(rbind, lapply(1:looks, dotrial))
    
    if(i%%100 == 0){
      cat(paste0("i = ", i, "\n"))
    }
    
    return(res)
  }
  
  saveRDS(results, "results.RDS")
  
  if(!single_clust) stopCluster(cl)
  
  # end clock
  end <- proc.time()
  
  duration <- end - start
  return(duration)
}

nsim <- 1000
duration <- test_func(single_clust = T, nsim)

cat(paste0("\nuser is the CPU time spent by the current process\nsystem CPU time” gives the CPU time spent by the kernel\nelapsed is the wall clock time taken to execute the function\n
============================================================\n\n"))

cat(sprintf(" user    %s\n system  %s\n elapsed %s\n", 
            round(duration[1], 2), 
            round(duration[2], 2), 
            round(duration[3], 2) ))

cat(sprintf(" user time per sim   %s\n\n", 
            round(duration[1], 2) / nsim))


nsim <- 10000
duration <- test_func(single_clust = F, nsim)

cat(paste0("\nuser is the CPU time spent by the current process\nsystem CPU time” gives the CPU time spent by the kernel\nelapsed is the wall clock time taken to execute the function\n
============================================================\n\n"))

cat(sprintf(" user    %s\n system  %s\n elapsed %s\n", 
            round(duration[1], 2), 
            round(duration[2], 2), 
            round(duration[3], 2) ))

cat(sprintf(" elapsed time per sim   %s\n\n", 
            round(duration[3], 2) / nsim))


