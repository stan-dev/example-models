data {
  int<lower=0> N;
  int<lower=0,upper=1> switched[N];
  vector[N] dist;
  vector[N] arsenic;
}
transformed data {
  // rescaling
  vector[N] dist100 = dist / 100.0;
  matrix[N,2] x = [dist100', arsenic']';
}
parameters {
  real alpha;
  vector[2] beta;
}
model {
  switched ~ bernoulli_logit_glm(x, alpha, beta);
}
