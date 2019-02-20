/*  Variable naming:
 obs       = observed
 cen       = (right) censored
 N         = number of samples
 tau       = scale parameter
*/
data {
  int<lower=0> Nobs;
  int<lower=0> Ncen;
  vector[Nobs] yobs;
  vector[Ncen] ycen;
}

parameters {
  real rate;
}

model {
  //yobs ~ weibull(alpha, exp(-(mu)/alpha));
  
  target += exponential_lpdf( yobs | rate);
  target += exponential_lccdf(ycen | rate);

  rate ~ normal(0.0, 1.0);
}

