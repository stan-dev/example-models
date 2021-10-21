data {
  int<lower=0> N;
  vector[N] incumbency_88;
  vector[N] vote_86;
  vector[N] vote_88;
}
transformed data {
  matrix[N, 2] x = [vote_86', incumbency_88']';
}
parameters {
  real alpha;
  vector[2] beta;
  real<lower=0> sigma;
}
model {
  vote_88 ~ normal_id_glm(x, alpha, beta, sigma);
}
