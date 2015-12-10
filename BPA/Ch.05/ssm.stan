data {
  int<lower=0> T;
  real y[T];
}

parameters {
  real<lower=0> mean_lambda;    // Mean growth rate
  real<lower=0> sigma_proc;     // SD of state process
  real<lower=0> sigma_obs;      // SD of observation process
  real<lower=0> lambda[T-1];
  real<lower=0> N_est1;         // Initial population size
}

transformed parameters {
  real<lower=0> N_est[T];

  N_est[1] <- N_est1;
  // State process
  for (t in 1:(T - 1))
    N_est[t + 1] <- N_est[t] * lambda[t];
}

model {
  // Priors
  N_est1 ~ uniform(0, 500);
  mean_lambda ~ uniform(0, 10);
  sigma_proc ~ uniform(0, 10);
  sigma_obs ~ uniform(0, 100);

  // Likelihood
  lambda ~ normal(mean_lambda, sigma_proc);

  // Observation process
  y ~ normal(N_est, sigma_obs);
}

generated quantities {
  real<lower=0> sigma2_obs;
  real<lower=0> sigma2_proc;

  sigma2_obs <- square(sigma_obs);
  sigma2_proc <- square(sigma_proc);
}
