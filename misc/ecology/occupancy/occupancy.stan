/**
 * Occupancy model; Single season, no covariates
 *
 * Translated from Marc Kery's WinBUGS model:
 * http://www.fisheriesstockassessment.com/TikiWiki/tiki-index.php?page=BUGS+Occupancy
 *
 * Based on: 
 * J. Andrew Royle and Marc Kery. 2007. A Bayesian state-space
 * formulation of dynamic occupancy models.  Ecology 88(7),
 * 1813--1823.
 * URL: http://www.fisheriesstockassessment.com/files/RoyleKery2007.pdf
 */

data {
  int<lower=0> R;
  int<lower=0> T;
  int<lower=0,upper=1> y[R,T];
}
parameters {
  real<lower=0,upper=1> psi1;
  real<lower=0,upper=1> p;
}
model {
  // local variables to avoid recomputing log(psi1) and log(1 - psi1)
  real log_psi1;
  real log1m_psi1;
  log_psi1 = log(psi1);
  log1m_psi1 = log1m(psi1);

  // priors
  psi1 ~ uniform(0,1);
  p ~ uniform(0,1);
  
  // likelihood
  for (r in 1:R) {
    if (sum(y[r]) > 0)
      target += log_psi1 + bernoulli_lpmf(y[r] | p);
    else
      target += log_sum_exp(log_psi1 + bernoulli_lpmf(y[r] | p),
			    log1m_psi1);
  }
}
