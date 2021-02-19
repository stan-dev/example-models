data {
  int<lower=0> N;
  vector[N] earn;
  vector[N] height;
}
transformed data {           // log 10 transformation
  vector[N] log10_earn = log10(earn);
  matrix[N,1] x = [height']';
}
parameters {
  real alpha;
  vector[1] beta;
  real<lower=0> sigma;
}
model {
  log10_earn ~ normal_id_glm(x, alpha, beta, sigma);
}
