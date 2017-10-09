data {
  int<lower=0> N;
  int<lower=0> n_groups;
  int<lower=0> n_scenarios;
  int<lower=1,upper=n_groups> group_id[N];
  int<lower=1,upper=n_scenarios> scenario_id[N];
  vector[N] y;
}
parameters {
  vector[n_groups] gamma;
  vector[n_scenarios] delta;
  real mu;
  real<lower=0,upper=100> sigma_gamma;
  real<lower=0,upper=100> sigma_delta;
  real<lower=0,upper=100> sigma_y;
}
transformed parameters {
  vector[N] y_hat;

  for (i in 1:N)
    y_hat[i] = mu + gamma[group_id[i]] + delta[scenario_id[i]];
}
model {
  # mu ~ unif(-inf, +inf); # implicit
  gamma ~ normal(0, sigma_gamma);
  delta ~ normal(0, sigma_delta);
  y ~ normal(y_hat, sigma_y);
}
