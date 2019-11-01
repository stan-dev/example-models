data {
  int<lower = 1> K;
  int<lower = 0> N[K];
  int<lower = 0> y[K];
}
parameters {
  vector<lower = 0, upper = 1>[K] theta;
}
model {
  y ~ binomial(N, theta);
}
generated quantities {
  int<lower = 0, upper = 1> is_best[K];
  for (k in 1:K)
    is_best[k] = (theta[k] >= best_prob);
}
