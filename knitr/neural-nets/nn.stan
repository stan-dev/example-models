functions {

  /**
   * Returns linear predictor for restricted Boltzman machine (RBM).
   * Assumes one-hidden layer with logistic sigmoid activation.
   *
   * @param x Predictors (N x M)
   * @param alpha First-layer weights (M x J)
   * @param beta Second-layer weights (J x (K - 1))
   * @return Linear predictor for output layer of RBM.
   */
  matrix rbm(matrix x, matrix alpha, matrix beta) {
    return tanh(x * alpha) * beta;
  }

}
data {
  int<lower=0> N;               // num train instances
  int<lower=0> M;               // num train predictors
  matrix[N, M] x;               // train predictors
  int<lower=2> K;               // num categories
  int<lower=1, upper=K> y[N];   // train category

  int<lower=1> J;               // num hidden units

  int<lower=0> Nt;              // num test instances
  matrix[Nt, M] xt;             // test predictors
  int<lower=1, upper=K> yt[N];  // test category
}
parameters {
  matrix[M, J] alpha;
  matrix[J, K] beta;
}
model {
  matrix[K, N] v = rbm(x, alpha, beta)';

  // priors
  to_vector(alpha) ~ normal(0, 2);
  to_vector(beta) ~ normal(0, 2);

  // likelihood
  for (n in 1:N)
    y[n] ~ categorical_logit(v[ , n]);
}


generated quantities {
  int<lower=1, upper=K> yt_sim[Nt];  // test observations
  real log_p_yt;
  real accuracy;
  matrix[K, Nt] v;

  v = rbm(xt, alpha, beta)';
  log_p_yt = 0;
  accuracy = 0;
  for (n in 1:Nt) {
    vector[K] theta;
    theta = softmax(col(v, n));
    log_p_yt = log_p_yt + categorical_lpmf(yt_sim | theta);
    yt_sim[n] = categorical_rng(theta);
    if (yt_sim[n] == yt[n])
      accuracy = accuracy + 1;
  }
  accuracy = accuracy / Nt;
}
