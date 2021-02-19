data {
  int<lower=0> N;
  vector[N] earn;
  vector[N] height;
  vector[N] male;
}
transformed data {           // log transformation
  vector[N] log_earn = log(earn);
  matrix[N,2] x = [height', male']';
}
parameters {
  real alpha;
  vector[2] beta;
  real<lower=0> sigma;
}
model {
  log_earn ~ normal_id_glm(x, alpha, beta, sigma);
}
