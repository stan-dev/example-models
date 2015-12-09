data {
  int<lower=0> nobs;
  int<lower=0> nmis;
  int<lower=0> nyear;
  int<lower=0> nsite;
  int<lower=0> nnewobs;
  int<lower=0> obs[nobs];
  int<lower=0> obsyear[nobs];
  int<lower=0> obssite[nobs];
  int<lower=0> misyear[nmis];
  int<lower=0> missite[nmis];
  int<lower=0,upper=1> first[nyear, nsite];
  real year[nyear];
  int<lower=0> newobs[nyear, nsite];
}

parameters {
  real mu;                      // Overall intercept
  real beta1;                   // Overall trend 
  real beta2;                   // First-year observer effect
  real alpha[nsite];            // Random site effects
  real<lower=0> sd_alpha;
  real eps[nyear];              // Random year effects
  real<lower=0> sd_eps;
  real gamma[nnewobs];          // Random observer effects
  real<lower=0> sd_gamma;
}

transformed parameters {
  real log_lambda[nyear, nsite];

  for (i in 1:nyear)
    for (j in 1:nsite)
      log_lambda[i, j] <- mu + beta1 * year[i] + beta2 * first[i, j] +
                          alpha[j] + gamma[newobs[i, j]] + eps[i];
}

model {
  // Priors
  mu ~ normal(0, 10);
  beta1 ~ normal(0, 10);
  beta2 ~ normal(0, 10);

  alpha ~ normal(0, sd_alpha);
  sd_alpha ~ uniform(0, 3);

  eps ~ normal(0, sd_eps);
  sd_eps ~ uniform(0, 1);

  gamma ~ normal(0, sd_gamma); 
  sd_gamma ~ uniform(0, 1);

  // Likelihood
  for (i in 1:nobs)
    obs[i] ~ poisson_log(log_lambda[obsyear[i], obssite[i]]);
}

generated quantities {
  int<lower=0> mis[nmis];

  for (i in 1:nmis)
    mis[i] <- poisson_log_rng(log_lambda[misyear[i], missite[i]]);
}
