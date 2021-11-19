data {
  int<lower=0> N; // sample size
  vector[N] y; // observed outcome
  vector[N] w; // treatment assigned
  real<lower=-1, upper=1> rho; // assumed correlation between the potential outcomes
}
parameters {
  real alpha; // intercept
  real tau; // super-population average treatment effect
  real<lower=0> sigma_c; // residual SD for the control
  real<lower=0> sigma_t; // residual SD for the treated
}
model {
  // PRIORS
  alpha ~ normal(0, 5);
  tau ~ normal(0, 5);
  sigma_c ~ normal(0, 5);
  sigma_t ~ normal(0, 5);
  
  // LIKELIHOOD
  y ~ normal(alpha + tau * w, sigma_t * w + sigma_c * (1 - w));
}
generated quantities {
  real tau_fs; // finite-sample average treatment effect  
  array[N] real y0; // potential outcome if W = 0
  array[N] real y1; // potential outcome if W = 1
  array[N] real tau_unit; // unit-level treatment effect
  for (n in 1 : N) {
    real mu_c = alpha;
    real mu_t = alpha + tau;
    if (w[n] == 1) {
      y0[n] = normal_rng(mu_c + rho * (sigma_c / sigma_t) * (y[n] - mu_t),
                         sigma_c * sqrt(1 - rho ^ 2));
      y1[n] = y[n];
    } else {
      y0[n] = y[n];
      y1[n] = normal_rng(mu_t + rho * (sigma_t / sigma_c) * (y[n] - mu_c),
                         sigma_t * sqrt(1 - rho ^ 2));
    }
    tau_unit[n] = y1[n] - y0[n];
  }
  tau_fs = mean(tau_unit);
}
