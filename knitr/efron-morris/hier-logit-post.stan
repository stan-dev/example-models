data {
  int<lower=0> N;                // items
  int<lower=0> K;                // trials
  int<lower=0, upper=K> y[N];    // successes
}
parameters {
  real mu;                       // population mean of success log-odds
  real<lower=0> sigma;           // population sd of success log-odds
  vector[N] alpha;               // success log-odds
}
model {
  mu ~ normal(-1, 1);                         // hyperprior
  sigma ~ normal(0, 1);                       // hyperprior
  alpha ~ normal(0, 1);                       // hierarchical prior
  y ~ binomial_logit(K, mu + sigma * alpha);  // likelihood
}
generated quantities {
  int<lower=0, upper=K> y_rep[N];      // replications for existing items
  int<lower=0, upper=K> y_pop_rep[N];  // replications for simulated items
  for (n in 1:N)
    y_pop_rep[n] <- binomial_rng(K, inv_logit(mu + sigma * normal_rng(0, 1)));
  for (n in 1:N)
    y_rep[n] <- binomial_rng(K, inv_logit(mu + sigma * alpha[n]));
}
