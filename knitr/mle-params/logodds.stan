data {
  int<lower=0> N;
  array[N] int<lower=0, upper=1> y;
}
parameters {
  real alpha;
}
transformed parameters {
  real<lower=0, upper=1> theta = inv_logit(alpha);
}
model {
  theta ~ uniform(0, 1);
  y ~ bernoulli(alpha);
}
