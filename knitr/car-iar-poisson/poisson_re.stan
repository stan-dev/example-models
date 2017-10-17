data {
  int<lower=0> N;
  int<lower=0> y[N];              // count outcomes
  vector[N] x;                    // predictor
}
parameters {
  real beta0;                // intercept
  real beta1;                // slope
  vector[N] theta;           // heterogeneous random effects
  real<lower=0> sigma;       // non-centered re variance 
}
model {
  y ~ poisson_log(beta0 + beta1 * x + theta * sigma);
  beta0 ~ normal(0, 5);
  beta1 ~ normal(0, 5);
  theta ~ normal(0, 1);
  sigma ~ normal(0, 5);
}
generated quantities {
  vector[N] mu = exp(beta0 + beta1 * x + theta * sigma);
}
