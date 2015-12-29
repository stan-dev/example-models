data {
  int<lower=0> N;                     // players
  int<lower=0> K;                     // first at-bats
  int<lower=0, upper=K> y[N];         // first hits
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
