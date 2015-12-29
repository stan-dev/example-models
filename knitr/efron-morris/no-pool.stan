data {
  int<lower=0> N;                     // items
  int<lower=0> K;                     // initial trials
  int<lower=0, upper=K> y[N];         // initial successes
}
parameters {
  vector<lower=0, upper=1>[N] theta;  // player ability
}
model {
  y ~ binomial(K, theta);             // likelihood
}
