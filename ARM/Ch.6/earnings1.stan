data {
  int<lower=0> N;
  array[N] int<lower=0, upper=1> earn_pos;
  vector[N] height;
  vector[N] male;
}
transformed data {
  matrix[N, 2] x = [height', male']';
}
parameters {
  real alpha;
  vector[2] beta;
}
model {
  earn_pos ~ bernoulli_logit_glm(x, alpha, beta);
}
