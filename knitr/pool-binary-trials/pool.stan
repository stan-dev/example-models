#include data-blocks.stan
parameters {
  real<lower=0, upper=1> phi; // chance of success (pooled)
}
model {
  y ~ binomial(K, phi); // likelihood
}
generated quantities {
  // define theta, per-player chance-of-success
  vector<lower=0, upper=1>[N] theta = rep_vector(phi, N);

  #include gq-postpred.stan
}
