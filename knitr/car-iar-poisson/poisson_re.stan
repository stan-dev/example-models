data {
  int<lower=0> N;
  int<lower=0> y[N];              // count outcomes
  vector[N] x;                    // predictor
}
parameters {
  real beta0;                // intercept
  real beta1;                // slope
  real<lower=0> sigma;        // overall standard deviation
  vector[N] theta_std;       // standardized heterogeneous effects
}
transformed parameters {
  vector[N] theta = sigma * theta_std;
}
model {
  y ~ poisson_log(beta0 + beta1 * x + theta);
  beta0 ~ normal(0, 5);
  beta1 ~ normal(0, 5);
  theta_std ~ normal(0, 1);
  sigma ~ normal(0, 5);
}
generated quantities {
  vector[N] mu = exp(beta0 + beta1 * x + theta);
}
