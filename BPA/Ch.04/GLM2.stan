data {
  int<lower=0> nobs;            // Number of observed data
  int<lower=0> nmis;            // Number of missing data
  int<lower=0> nyear;           // Number of years
  int<lower=0> nsite;           // Number of sites
  int<lower=0> obs[nobs];       // Observed counts
  int<lower=0> obsyear[nobs];   // Years in observed data
  int<lower=0> obssite[nobs];   // Sites in observed data
  int<lower=0> misyear[nmis];   // Years in missing data
  int<lower=0> missite[nmis];   // Sites in missing data
}

parameters {
  vector[nsite] alpha;          // Site effects
  vector[nyear - 1] eps2;       // Year effects (year > 1)
}

transformed parameters {
  vector[nyear] eps;            // Year effects
  matrix[nyear, nsite] log_lambda;

  eps[1] = 0;
  eps[2:nyear] = eps2[1:nyear - 1];
  log_lambda = rep_matrix(alpha', nyear)
             + rep_matrix(eps, nsite);
}

model {
  // Priors
  alpha ~ normal(0, 10);
  eps2 ~ normal(0, 10);

  // Likelihood
  for (i in 1:nobs)
    obs[i] ~ poisson_log(log_lambda[obsyear[i], obssite[i]]);
}

generated quantities {
  int<lower=0> mis[nmis];

  for (i in 1:nmis)
    mis[i] = poisson_log_rng(log_lambda[misyear[i], missite[i]]);
}
