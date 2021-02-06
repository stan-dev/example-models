data {
  int<lower=0> N;
  int<lower=0,upper=1> switched[N];
  vector[N] dist;
}
transformed data {
  matrix[N,1] x = [dist']';
}
parameters {
  real alpha;
  vector[1] beta;
}
model {
  switched ~ bernoulli_logit_glm(x, alpha, beta);
}
