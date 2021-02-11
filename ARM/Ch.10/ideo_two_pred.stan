data {
  int<lower=0> N; 
  vector[N] party;
  vector[N] score1;
  vector[N] x;
}
transformed data {
  matrix[N,2] cov = [party', x']';
}
parameters {
  real alpha;
  vector[2] beta;
  real<lower=0> sigma;
} 
model {
  score1 ~ normal_id_glm(cov, alpha, beta, sigma);
}
