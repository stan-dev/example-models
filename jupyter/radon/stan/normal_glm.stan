/* normal linear regression model */
data {
  int<lower=0> N;  // number of observations
  int<lower=1> M;  // number of predictors
  matrix[N, M] xs;
  vector[N] y;
}
parameters {
  real alpha;
  vector beta[M];
  real<lower=0> sigma;
}
model {
  y ~ normal_glm_id(xs, alpha, beta, sigma);
  sigma ~ normal(0, 10);   // override default uniform(-inf, +inf)
}
