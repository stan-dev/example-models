data {
  int<lower=0> N;
  int<lower=0, upper=1> y[N];
}
parameters {
  real alpha;
}
model {
  for (n in 1:N)
    y[n] ~ bernoulli(inv_logit(alpha));
}
generated quantities {
  real<lower=0, upper=1> theta;
  theta <- inv_logit(alpha);
}
