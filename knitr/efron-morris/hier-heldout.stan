data {
  int<lower=0> N;                     // players
  int<lower=0> K;                     // first at-bats
  int<lower=0, upper=K> y[N];         // first hits
  int<lower=0> K_new[N];              // remaining at-bats
  int<lower=0> y_new[N];              // remaining hits
}
parameters {
  real<lower=0, upper=1> phi;         // population ability
  real<lower=1> kappa;                // population concentration
  vector<lower=0, upper=1>[N] theta;  // player ability
}
model {
  kappa ~ pareto(1, 1.5);                        // hyperprior
  theta ~ beta(phi * kappa, (1 - phi) * kappa);  // prior
  y ~ binomial(K, theta);                        // likelihood
}
generated quantities {
  real log_p_new;  // posterior predictive density of rest of season
  log_p_new <- 0;
  for (n in 1:N)
    log_p_new <- log_p_new + binomial_log(y_new[n], K_new[n], theta[n]);
}
