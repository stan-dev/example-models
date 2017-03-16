data {
  int<lower=0> T;
  vector[T] y;
}

parameters {
  real<lower=0,upper=10> mean_lambda;    // Mean growth rate
  real<lower=0,upper=10> sigma_proc;     // SD of state process
  real<lower=0,upper=100> sigma_obs;     // SD of observation process
  vector<lower=0>[T - 1] lambda;
  real<lower=0,upper=500> N_est1;        // Initial population size
}

transformed parameters {
  vector<lower=0>[T] N_est;

  N_est[1] = N_est1;
  // State process
  for (t in 1:(T - 1))
    N_est[t + 1] = N_est[t] * lambda[t];
}

model {
  // Priors are implicitly defined
  //  N_est1 ~ uniform(0, 500);
  //  mean_lambda ~ uniform(0, 10);
  //  sigma_proc ~ uniform(0, 10);
  //  sigma_obs ~ uniform(0, 100);

  // Likelihood
  lambda ~ normal(mean_lambda, sigma_proc);

  // Observation process
  y ~ normal(N_est, sigma_obs);
}

generated quantities {
  real<lower=0> sigma2_obs;
  real<lower=0> sigma2_proc;

  sigma2_obs = square(sigma_obs);
  sigma2_proc = square(sigma_proc);
}
