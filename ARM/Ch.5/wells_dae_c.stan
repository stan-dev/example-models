data {
  int<lower=0> N;
  array[N] int<lower=0, upper=1> switched;
  vector[N] dist;
  vector[N] arsenic;
  vector[N] educ;
}
transformed data {
  vector[N] c_dist100 = (dist - mean(dist)) / 100.0;
  vector[N] c_arsenic = arsenic - mean(arsenic);
  vector[N] da_inter = c_dist100 .* c_arsenic;
  vector[N] educ4 = educ / 4.0;
  matrix[N, 4] x = [c_dist100', c_arsenic', da_inter', educ4']';
}
parameters {
  real alpha;
  vector[4] beta;
}
model {
  switched ~ bernoulli_logit_glm(x, alpha, beta);
}
