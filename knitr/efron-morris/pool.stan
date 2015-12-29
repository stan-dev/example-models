data {
  int<lower=0> N;                // items
  int<lower=0> K;                // initial trials
  int<lower=0, upper=K> y[N];    // initial successes
}
parameters {
  real<lower=0, upper=1> theta;  // player ability
}
model {
  y ~ binomial(K, theta);        // likelihood
}
