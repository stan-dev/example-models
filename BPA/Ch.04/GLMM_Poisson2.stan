data {
  int<lower=0> nsite;           // Number of populations
  int<lower=0> nyear;           // Number of years
  int<lower=0> C[nyear, nsite]; // Counts
  real year[nyear];             // Year covariate
}

parameters {
  real mu;
  real alpha[nsite];
  real eps[nyear];
  real beta[3];
  real sd_alpha;
  real sd_year;
}

transformed parameters {
  real log_lambda[nyear, nsite];

  for (i in 1:nyear)
    for (j in 1:nsite)
    // Linear predictor including random site and random year effects
      log_lambda[i, j] <- alpha[j] + beta[1] * year[i] +
                          beta[2] * year[i]^2 + beta[3] * year[i]^3 +
                          eps[i];
}

model {
  // Priors

  // Random site effects
  alpha ~ normal(mu, sd_alpha);

  // Hyperparameter 1
  mu ~ normal(0, 10);

  // Hyperparameter 2
  sd_alpha ~ uniform(0, 2);

  beta ~ normal(0, 10);

  // Hyperparameter 3
  sd_year ~ uniform(0, 1);

  // Likelihood
  for (i in 1:nyear) {
    // Random year effects
    eps[i] ~ normal(0, sd_year);
    for (j in 1:nsite) {
      // Distribution for random part
      // Link function
      C[i, j] ~ poisson_log(log_lambda[i, j]);
    } #j
  } #i
}
