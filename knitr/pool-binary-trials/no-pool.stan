#include data-blocks.stan
parameters {
  vector<lower=0, upper=1>[N] theta; // chance of success
}
model {
  y ~ binomial(K, theta); // likelihood
}
generated quantities {
  #include gq-postpred.stan
  #include gq-ranking.stan
}
