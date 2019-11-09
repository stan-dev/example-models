/**
 * One layer RBM.
 */
data {
  int<lower=0> N;               // num train instances
  int<lower=0> M;               // num train predictors
  matrix[N, M] x;               // train predictors
  int<lower=2> K;               // num categories
  int<lower=1, upper=K> y[N];   // train category
  int<lower=1> J;               // num hidden units
}
parameters {
  matrix[M, J] alpha;
  matrix[J, K] beta;
}
model {
  matrix[K, N] v = (tanh(x * alpha) * beta)';

  // priors
  to_vector(alpha) ~ normal(0, 2);
  to_vector(beta) ~ normal(0, 2);

  // likelihood
  for (n in 1:N)
    y[n] ~ categorical_logit(v[ , n]);
}
