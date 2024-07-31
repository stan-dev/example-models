#include data-blocks.stan
parameters {
  real mu; // population mean of success log-odds
  real<lower=0> sigma; // population sd of success log-odds
  vector<offset=mu, multiplier=sigma>[N] alpha_std; // success log-odds (standardized)
}
model {
  mu ~ normal(-1, 1); // hyperprior
  sigma ~ normal(0, 1); // hyperprior
  alpha_std ~ normal(mu, sigma); // prior (hierarchical)
  y ~ binomial_logit(K, alpha_std); // likelihood
}
generated quantities {
  vector[N] theta = inv_logit(alpha_std);
  #include gq-postpred.stan
  #include gq-ranking.stan

  array[N] int<lower=0> y_pop_rep; // replications for simulated items
  for (n in 1 : N) {
    y_pop_rep[n] = binomial_rng(K[n], inv_logit(normal_rng(mu, sigma)));
  }
}
