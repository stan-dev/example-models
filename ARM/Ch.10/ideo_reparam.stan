data {
  int<lower=0> N; 
  vector[N] party;
  vector[N] score1;
  vector[N] z1; //z value for party 0, 0 otherwise 
  vector[N] z2; //z value for party 1, 0 otherwise 
}
transformed data {
  matrix[N,3] x = [party', z1', z2']';
}
parameters {
  real alpha;
  vector[3] beta;
  real<lower=0> sigma;
} 
model {
  score1 ~ normal_id_glm(x, alpha, beta, sigma);
}
