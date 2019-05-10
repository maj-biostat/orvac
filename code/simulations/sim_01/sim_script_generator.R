library(dplyr)

nsim <- 1000  # 10000
accrual <- c(30, 50) 
info_delay <- c(0.5)
p0 <- c(0.1, 0.4, 0.7)
pdelta <- c(0.05, 0.1, 0.15) # 0
t0 <- c(20, 35, 50)
tdelta <-  c(5, 10, 15) # 0

d <- expand.grid(accrual, info_delay, 
                 p0, pdelta, 
                 t0, tdelta)
names(d) <- c("accural", "delay", "p0", "pdelta", "t0", "tdelta")
d$p1 <- d$p0 + d$pdelta
d$t1 <- d$t0 + d$tdelta
d <- d[, c(-4, -6)]
d <- d[, c("accural", "delay", "p0", "p1", "t0", "t1")]



idx_start <- 1000000
idx_end <-   9999999

fileConn<-file("run_sim.sh")

out <- c("#!/bin/bash")

seeds <- sample(idx_start:idx_end, nrow(d))
for(i in 1:nrow(d)){
  
  out <- c(out, paste0("/usr/bin/Rscript main_4.R -f cfg1.yaml -o T -i SIM", seeds[i],
                 " -l log_", seeds[i], ".log",
                 " -n ", nsim,
                 " -s ", seeds[i],
                 " -a ", d$accural[i],
                 " -d ", d$delay[i],
                 " -b ", d$p0[i],
                 " -p ", d$p1[i],
                 " -m ", d$t0[i],
                 " -t ", d$t1[i]))
  
  
}

writeLines(out, fileConn)
close(fileConn)







