data {
  int<lower=0> N;
  array[N] real x;
  array[N] real Y;
}
parameters {
  real<lower=0> sigma;
  real<lower=0> alpha;
  array[2] real beta;
  real<lower=min(x), upper=max(x)> x_change;
}
model {
  array[N] real mu;
  
  alpha ~ normal(0, 5);
  beta ~ normal(0, 5);
  sigma ~ cauchy(0, 5);
  
  for (n in 1 : N) {
    mu[n] = alpha
            + ((x[n] < x_change) ? beta[1] : beta[2]) * (x[n] - x_change);
  }
  
  Y ~ normal(mu, sigma);
}
