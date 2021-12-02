data {
  int<lower=0> N;
  vector[N] earnings;
  vector[N] height;
  vector[N] sex;
}
transformed data {
  vector[N] log_earnings = log(earnings);
  vector[N] male = 2 - sex;
  matrix[N, 2] x = [height', male']';
}
parameters {
  real alpha;
  vector[2] beta;
  real<lower=0> sigma;
}
model {
  log_earnings ~ normal_id_glm(x, alpha, beta, sigma);
}
