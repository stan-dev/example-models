data {
  int<lower=0> N;
  vector[N] dist;
  array[N] int<lower=0, upper=1> switc;
}
transformed data {
  vector[N] dist100 = dist / 100.0;
  matrix[N, 1] x = [dist100']';
}
parameters {
  real alpha;
  vector[1] beta;
}
model {
  switc ~ bernoulli(Phi(alpha + x * beta));
}
