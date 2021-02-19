data {
  int<lower=0> N;
  vector[N] post_test;
  vector[N] treatment;
  vector[N] pre_test;
}
transformed data {
  // interaction
  vector[N] inter = treatment .* pre_test;
  matrix[N,3] x = [treatment', pre_test', inter']';
}
parameters {
  real alpha;
  vector[3] beta;
  real<lower=0> sigma;
}
model {
  post_test ~ normal_id_glm(x, alpha, beta, sigma);
}
