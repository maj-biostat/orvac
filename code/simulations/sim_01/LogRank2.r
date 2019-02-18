LogRank2 <- function(z1,delta1,z2,delta2,k=1)
{
# This function compute the 2 sample log-rank test statistic
# for the censored data. 
# Input: z1 is observed times (may be censored), a vector of length n. 
#        delta1 is censoring indicator, 0/1's means censored/not-censored.
#        z1 must have same length as delta1. 
#        z2 and delta2 are similar, they can have length m. 
#        k is a constant (alternative specific) try to achieve better power.
#############################################################################
   t1 <- z1[delta1 == 1]
   t2 <- z2[delta2 == 1]
   risk1t1 <- psum( t( outer(z1, t1, FUN=">=") ) )
   risk1t2 <- psum( t( outer(z1, t2, FUN=">=") ) )
   risk2t1 <- psum( t( outer(z2, t1, FUN=">=") ) )
   risk2t2 <- psum( t( outer(z2, t2, FUN=">=") ) )
   
    # z1 = d2$obst[d2$trt == 0];    z1 <- z1[!is.na(z1)]
    # delta1 = 1-d2$cen[d2$trt == 0]; delta1 <- delta1[!is.na(delta1)]
    # 
    # z2 = d2$obst[d2$trt == 1];    z2 <- z2[!is.na(z2)]
    # delta2 = 1-d2$cen[d2$trt == 1]; delta2 <- delta2[!is.na(delta2)]
    


   num <- sum(risk2t1/(risk1t1+k*risk2t1))-sum(risk1t2/(risk1t2+k*risk2t2))
   var2 <- sum(risk1t1*risk2t1/(risk1t1+k*risk2t1)^2) + 
                       sum(risk1t2*risk2t2/(risk1t2+k*risk2t2)^2)

# As far as SAS Wilcoxson test inside the lifetest,  when there is no tie 
#    numGH <- sum(risk2t1)-sum(risk1t2)                    This version 
#    varGH <- sum(risk1t1*risk2t1) + sum(risk1t2*risk2t2)  agrees with SAS.
#    varGH1 <- sum(risk2t1^2) + sum(risk1t2^2)             Not this version
# Ref. for use this version of var est.?  Gill 1981 (V2) (4.1.21) 
# This agrees with Splus when there is no tie.

   temp2 <- num/sqrt(var2)
   pval <- 1-pchisq((temp2)^2, df=1)

   list(VarA2=var2, Logrank=temp2, ApproxPvalue2side=pval)
} 
####another possibility of replacing psum() below is to use rowsum()
####to make a psum2().  Which one is faster??
psum <- function (..., na.rm = FALSE)  {
  args <- do.call("cbind", list(...))
  nargs <- ncol(args)
#  if (na.rm) {
#    args[is.na(args)] <- 0
#  }
  as.vector(args %*% rep(1, nargs))
}

# In R, 
# psum2 and psum3 only work for matrix, not a single number?? The
# problem is in the rowsum(), try to give a dimname to single name?
# In Splus2000 it is OK. (but slow)
psum2 <- function(mat) {
         mat <- as.matrix(mat)
         mat2 <- cbind(0, mat)
         rowsum( t(mat2), c(0, rep(1, ncol(mat))), reorder=F)[2,]
}

####For Splus2000, the function rowsum() changed: not allow a reorder=
####so the following version, psum3, is for Splus2000
psum3 <- function(mat) {
         mat <- as.matrix(mat) 
         mat2 <- cbind(0, mat)
         as.vector(rowsum( t(mat2), c(0, rep(1, ncol(mat))))[2,])
}
