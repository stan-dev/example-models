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
  vector[nsite] alpha;          // Random site effects
  real mu_alpha;
  real<lower=0,upper=5> sd_alpha;
}

transformed parameters {
  matrix[nyear, nsite] log_lambda;

  log_lambda = rep_matrix(alpha', nyear);
}

model {
  // Priors
  alpha ~ normal(mu_alpha, sd_alpha);
  mu_alpha ~ normal(0, 10);
  //  sd_alpha ~ uniform(0, 5);  // Implicitly defined

  for (i in 1:nobs)
    obs[i] ~ poisson_log(log_lambda[obsyear[i], obssite[i]]);
}

generated quantities {
  int<lower=0> mis[nmis];

  for (i in 1:nmis)
    mis[i] = poisson_log_rng(log_lambda[misyear[i], missite[i]]);
}
