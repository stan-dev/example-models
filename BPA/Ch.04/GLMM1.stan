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
  real alpha[nsite];            // Random site effects
  real mu_alpha;
  real<lower=0> sd_alpha;
}

transformed parameters {
  real log_lambda[nyear, nsite];

  for (i in 1:nyear)
    for (j in 1:nsite)
      log_lambda[i, j] <- alpha[j];
}

model {
  // Priors
  for (j in 1:nsite) {
    alpha[j] ~ normal(mu_alpha, sd_alpha);
  }
  mu_alpha ~ normal(0, 10);
  sd_alpha ~ uniform(0, 5);

  for (i in 1:nobs)
    obs[i] ~ poisson_log(log_lambda[obsyear[i], obssite[i]]);
}

generated quantities {
  int<lower=0> mis[nmis];

  for (i in 1:nmis)
    mis[i] <- poisson_log_rng(log_lambda[misyear[i], missite[i]]);
}
