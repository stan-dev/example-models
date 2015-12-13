data {
  int<lower=0> nobs;
  int<lower=0> nmis;
  int<lower=0> nyear;
  int<lower=0> nsite;
  int<lower=0> obs[nobs];
  int<lower=0> obsyear[nobs];
  int<lower=0> obssite[nobs];
  int<lower=0> misyear[nmis];
  int<lower=0> missite[nmis];
}

parameters {
  real alpha[nsite];            // site effects
  real eps2[nyear - 1];	        // year effects (year > 1)
}

transformed parameters {
  real eps[nyear];	        // year effects
  real log_lambda[nyear, nsite];

  eps[1] <- 0;
  for (i in 2:nyear)
    eps[i] <- eps2[i - 1];
  for (i in 1:nyear)
    for (j in 1:nsite)
       log_lambda[i, j] <- alpha[j] + eps[i];
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
    mis[i] <- poisson_log_rng(log_lambda[misyear[i], missite[i]]);
}