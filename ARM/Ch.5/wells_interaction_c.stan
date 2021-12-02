data {
  int<lower=0> N;
  array[N] int<lower=0, upper=1> switched;
  vector[N] dist;
  vector[N] arsenic;
}
transformed data {
  // centering
  vector[N] c_dist100 = (dist - mean(dist)) / 100.0;
  vector[N] c_arsenic = arsenic - mean(arsenic);
  // interaction
  vector[N] inter = c_dist100 .* c_arsenic;
  matrix[N, 3] x = [c_dist100', c_arsenic', inter']';
}
parameters {
  real alpha;
  vector[3] beta;
}
model {
  switched ~ bernoulli_logit_glm(x, alpha, beta);
}
