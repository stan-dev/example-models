data {
  int<lower=0> N;
  array[N] int<lower=0, upper=1> switched;
  vector[N] dist;
  vector[N] arsenic;
  vector[N] educ;
}
transformed data {
  // rescaling
  vector[N] dist100 = dist / 100.0;
  vector[N] educ4 = educ / 4.0;
  matrix[N, 3] x = [dist100', arsenic', educ4']';
}
parameters {
  real alpha;
  vector[3] beta;
}
model {
  switched ~ bernoulli_logit_glm(x, alpha, beta);
}
