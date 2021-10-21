data {
  int<lower=0> N;
  vector[N] party;
  vector[N] score1;
  vector[N] x;
}
transformed data {
  vector[N] inter = party .* x;
  matrix[N, 3] cov = [party', x', inter']';
}
parameters {
  real alpha;
  vector[3] beta;
  real<lower=0> sigma;
}
model {
  score1 ~ normal_id_glm(cov, alpha, beta, sigma);
}
