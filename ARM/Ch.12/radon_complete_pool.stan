data {
  int<lower=1> N;
  vector[N] x;
  vector[N] y;
}
transformed data {
  matrix[N, 1] cov = [x']';
}
parameters {
  real alpha;
  vector[1] beta;
  real<lower=0> sigma;
}
model {
  sigma ~ cauchy(0, 2.5);
  y ~ normal_id_glm(cov, alpha, beta, sigma);
}
