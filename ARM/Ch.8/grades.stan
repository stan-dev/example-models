data {
  int<lower=0> N;
  vector[N] final;
  vector[N] midterm;
}
transformed data {
  matrix[N, 1] x = [midterm']';
}
parameters {
  real alpha;
  vector[1] beta;
  real<lower=0> sigma;
}
model {
  final ~ normal_id_glm(x, alpha, beta, sigma);
}
