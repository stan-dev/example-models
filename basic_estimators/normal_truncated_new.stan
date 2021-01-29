data {
  real U;
  int<lower=1> N;
  vector<upper=U>[N] y;
}
parameters {
  real mu;
  real<lower=0, upper=2> sigma;
}
model {
  // Truncation not yet available with vectorised distributions
  for (n in 1:N)
    y[n] ~ normal(mu, sigma) T[,U];
}
