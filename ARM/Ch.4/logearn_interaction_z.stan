data {
  int<lower=0> N;
  vector[N] earn;
  vector[N] height;
  vector[N] male;
}
transformed data {
  // log transformation
  vector[N] log_earn = log(earn);
  // standardization
  vector[N] z_height = (height - mean(height)) / sd(height);
  // interaction
  vector[N] inter = z_height .* male;
  matrix[N,3] x = [z_height', male', inter']';
}
parameters {
  real alpha;
  vector[3] beta;
  real<lower=0> sigma;
}
model {
  log_earn ~ normal_id_glm(x, alpha, beta, sigma);
}
