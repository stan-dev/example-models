data {
  int<lower=0> nobs; // Number of observed data
  int<lower=0> nmis; // Number of missing data
  int<lower=0> nyear; // Number of years
  int<lower=0> nsite; // Number of sites
  array[nobs] int<lower=0> obs; // Observed counts
  array[nobs] int<lower=0> obsyear; // Years in observed data
  array[nobs] int<lower=0> obssite; // Sites in observed data
  array[nmis] int<lower=0> misyear; // Years in missing data
  array[nmis] int<lower=0> missite; // Sites in missing data
  array[nyear, nsite] int<lower=0, upper=1> first; // First-year observer?
  array[nyear] real year; // Year
  array[nyear, nsite] int<lower=0> newobs; // Observers
  int<lower=0> nnewobs;
}
parameters {
  real mu; // Overall intercept
  real beta1; // Overall trend
  real beta2; // First-year observer effect
  vector[nsite] alpha; // Random site effects
  real<lower=0, upper=3> sd_alpha;
  vector[nyear] eps; // Random year effects
  real<lower=0, upper=1> sd_eps;
  vector[nnewobs] gamma; // Random observer effects
  real<lower=0, upper=1> sd_gamma;
}
transformed parameters {
  matrix[nyear, nsite] log_lambda;
  
  for (j in 1 : nsite) {
    for (i in 1 : nyear) {
      log_lambda[i, j] = mu + beta1 * year[i] + beta2 * first[i, j]
                         + alpha[j] + gamma[newobs[i, j]] + eps[i];
    }
  }
}
model {
  // Priors
  mu ~ normal(0, 10);
  beta1 ~ normal(0, 10);
  beta2 ~ normal(0, 10);
  
  alpha ~ normal(0, sd_alpha);
  //  sd_alpha ~ uniform(0, 3); // Implicitly defined
  
  eps ~ normal(0, sd_eps);
  //  sd_eps ~ uniform(0, 1);   // Implicitly defined
  
  gamma ~ normal(0, sd_gamma);
  //  sd_gamma ~ uniform(0, 1); // Implicitly defined
  
  // Likelihood
  for (i in 1 : nobs) {
    obs[i] ~ poisson_log(log_lambda[obsyear[i], obssite[i]]);
  }
}
generated quantities {
  array[nmis] int<lower=0> mis;
  
  for (i in 1 : nmis) {
    mis[i] = poisson_log_rng(log_lambda[misyear[i], missite[i]]);
  }
}
