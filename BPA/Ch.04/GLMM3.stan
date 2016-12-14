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
  int<lower=0,upper=1> first[nyear, nsite];     // First-year observer?
}

parameters {
  real mu;                      // Overall mean
  vector[nsite] alpha;          // Random site effects
  real<lower=0,upper=5> sd_alpha;
  real beta2;                   // First-year observer effect
  vector[nyear] eps;            // Random year effects
  real<lower=0,upper=5> sd_eps;
}

transformed parameters {
  matrix[nyear, nsite] log_lambda;

  for (j in 1:nsite)
    for (i in 1:nyear)
      log_lambda[i, j] = mu + beta2 * first[i, j]
                       + alpha[j] + eps[i];
}

model {
  // Priors
  mu ~ normal(0, 10);
  beta2 ~ normal(0, 10);

  alpha ~ normal(0, sd_alpha);
  //  sd_alpha ~ uniform(0, 5);  // Implicitly defined

  eps ~ normal(0, sd_eps);
  //  sd_eps ~ uniform(0, 5);    // Implicitly defined

  // Likelihood
  for (i in 1:nobs)
    obs[i] ~ poisson_log(log_lambda[obsyear[i], obssite[i]]);
}

generated quantities {
  int<lower=0> mis[nmis];

  for (i in 1:nmis)
    mis[i] = poisson_log_rng(log_lambda[misyear[i], missite[i]]);
}
