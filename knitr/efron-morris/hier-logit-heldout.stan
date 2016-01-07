data {
  int<lower=0> N;                // items
  int<lower=0> K;                // trials
  int<lower=0, upper=K> y[N];    // successes
  int<lower=0> K_new[N];         // new trials
  int<lower=0> y_new[N];         // new successes
}
parameters {
  real mu;                       // population mean of success log-odds
  real<lower=0> sigma;           // population sd of success log-odds
  vector[N] alpha;               // success log-odds
}
model {
  mu ~ normal(-1, 1);            // hyperprior
  sigma ~ normal(0, 1);          // hyperprior
  alpha ~ normal(0, 1);          // hierarchical prior
  y ~ binomial_logit(K, mu + sigma * alpha);  // likelihood
}
generated quantities {
  real log_p_new;  // posterior predictive density of new data
  log_p_new <- 0;
  for (n in 1:N)
    log_p_new <- log_p_new + binomial_log(y_new[n], K_new[n],
                                          inv_logit(mu + sigma * alpha[n]));
}
