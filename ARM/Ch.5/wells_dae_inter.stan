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
  // interaction
  vector[N] inter = dist100 .* arsenic;
  matrix[N, 4] x = [dist100', arsenic', educ4', inter']';
}
parameters {
  real alpha;
  vector[4] beta;
}
model {
  switched ~ bernoulli_logit_glm(x, alpha, beta);
}
