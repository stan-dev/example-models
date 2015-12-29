/**
 * Model for Rat survival from Chapter 5 of Gelman et al., Bayesian
 * Data Analysis and Chapter 16 of Stan Manual.
 *
 * Implicit uniform prior on mu, and p(kappa) propto kappa^(-2.5)
 *
 * Prediction must use the RNG for proper predictive uncertainty.
 */
data {
  int<lower=0> N;                     // players
  int<lower=0> K1;                    // first at-bats
  int<lower=0,upper=K1> y1[N];        // first hits
  int<lower=0> K2[N];                 // remaining at-bats
}
parameters {
  real mu;                            // population ability mean
  real<lower=0> sigma;                // population ability sd
  vector[N] theta_std;                // player ability
}
transformed parameters {
  vector[N] theta;
  theta <- mu + sigma * theta_std;
}
model {
  sigma ~ normal(0, 1);               // weakly informative prior
  mu ~ normal(-1, 1);                 // weakly informative prior
  theta_std ~ normal(0, 1);           // hierarchical prior (vect)
  y1 ~ binomial_logit(K1, theta);     // likelihood (vect)
}
