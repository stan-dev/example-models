#include data-blocks.stan
parameters {
  real<lower=0, upper=1> phi; // population chance of success
  real<lower=1> kappa; // population concentration
  vector<lower=0, upper=1>[N] theta; // chance of success 
}
model {
  kappa ~ pareto(1, 1.5); // hyperprior
  theta ~ beta(phi * kappa, (1 - phi) * kappa); // prior
  y ~ binomial(K, theta); // likelihood
}
generated quantities {
  #include gq-postpred.stan
  #include gq-ranking.stan

  // replications for simulated items  
  array[N] int<lower=0> y_pop_rep; 
  for (n in 1 : N) {
    y_pop_rep[n] = binomial_rng(K[n],
                                beta_rng(phi * kappa, (1 - phi) * kappa));
  }
}
