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
  vector[Nobs] Xobs_bg;
  vector[Ncen] Xcen_bg;
  real<lower=0> alpha0; // shape hyperpar for gamma
  real<lower=0> beta0; // scale hyperpar for gamma
}

parameters {
  real<lower=0> b0;   // baseline rate
  real b1;          // coefficients
}

model {

  target += exponential_lpdf( yobs | b0 + Xobs_bg * b1  );
  target += exponential_lccdf(ycen | b0 + Xcen_bg * b1  );

  b0 ~ gamma(alpha0, beta0);
  b1 ~ normal(0, 1);
}

