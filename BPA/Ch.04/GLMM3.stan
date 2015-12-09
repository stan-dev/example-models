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
  int<lower=0,upper=1> first[nyear, nsite];
}

parameters {
  real mu;                      // Overall mean
  real alpha[nsite];            // Random site effects
  real<lower=0> sd_alpha;
  real beta2;                   // First-year observer effect
  real eps[nyear];              // Random year effects
  real<lower=0> sd_eps;
}

transformed parameters {
  real log_lambda[nyear, nsite];

  for (i in 1:nyear)
    for (j in 1:nsite)
      log_lambda[i, j] <- mu + beta2 * first[i, j] +
                          alpha[j] + eps[i];
}

model {
  // Priors
  mu ~ normal(0, 10);
  beta2 ~ normal(0, 10);

  alpha ~ normal(0, sd_alpha);
  sd_alpha ~ uniform(0, 5);

  eps ~ normal(0, sd_eps);
  sd_eps ~ uniform(0, 5);

  // Likelihood
  for (i in 1:nobs)
    obs[i] ~ poisson_log(log_lambda[obsyear[i], obssite[i]]);
}

generated quantities {
  int<lower=0> mis[nmis];

  for (i in 1:nmis)
    mis[i] <- poisson_log_rng(log_lambda[misyear[i], missite[i]]);
}
