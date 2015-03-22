data {
  int<lower=2> K;             // number of categories
  int<lower=2> D;             // number of predictors
  int<lower=0> N;             // number of observations
  matrix[N, D] x;             // predictors
  int<lower=1,upper=K> y[N];  // observations
}
parameters {
  matrix[K, D] beta;      // slopes
}
model {
  matrix[N, K] gamma;
  gamma <- x * beta';

  // prior
  to_vector(beta) ~ cauchy(0, 2.5);

  // likelihood
  for (n in 1:N)
    y[n] ~ categorical_logit(gamma[n]');
}
