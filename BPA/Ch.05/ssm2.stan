data {
  int<lower=0> T;
  vector[T] y;
  int<lower=0> pyears;          // Number of years for prediction
}

parameters {
  real logN_est1;               // Initial log population size
  real mean_r;                  // Mean growth rate
  real<lower=0,upper=1> sigma_proc; // SD of state procesdatas
  real<lower=0,upper=1> sigma_obs;  // SD of observation process
  vector[T - 1] r;
}

transformed parameters {
  vector[T] logN_est;           // Log population size

  // State process
  logN_est[1] = logN_est1;
  for (t in 1:(T - 1))
    logN_est[t + 1] = logN_est[t] + r[t];
}

model {
  // Priors and constraints
  logN_est1 ~ normal(5.6, 10);
  mean_r ~ normal(1, sqrt(1000));
  //  sigma_proc ~ uniform(0, 1);
  //  sigma_obs ~ uniform(0, 1);

  // Likelihood
  r ~ normal(mean_r, sigma_proc);

  // Observation process
  y ~ normal(logN_est, sigma_obs);
}

generated quantities {
  // Population sizes on real scale
  real sigma2_proc;
  real sigma2_obs;
  vector[pyears] pr;
  vector[pyears] plogN_est;             // Predicted log population size
  vector<lower=0>[T + pyears] N_est;    // Population size

  sigma2_proc = square(sigma_proc);
  sigma2_obs = square(sigma_obs);

  pr[1] = normal_rng(mean_r, sigma_proc);
  plogN_est[1] = logN_est[T] + pr[1];
  for (t in 2:pyears) {
    pr[t] = normal_rng(mean_r, sigma_proc);
    plogN_est[t] = plogN_est[t - 1] + pr[t];
  }
  for (t in 1:T)
    N_est[t] = exp(logN_est[t]);
  for (t in 1:pyears)
    N_est[T + t] = exp(plogN_est[t]);
}
