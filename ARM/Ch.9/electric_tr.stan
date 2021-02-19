data {
  int<lower=0> N;
  vector[N] post_test;
  vector[N] treatment;
}
transformed data {
  matrix[N,1] x = [treatment']';
}
parameters {
  real alpha;
  vector[2] beta;
  real<lower=0> sigma;
}
model {
  post_test ~ normal_id_glm(x, alpha, beta, sigma);
}
