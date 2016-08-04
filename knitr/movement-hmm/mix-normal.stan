data {
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real<lower=0, upper=1> theta;
  vector<lower=0>[2] sigma;
  ordered[2] mu;
}
model {
  for (n in 1:N)
    target += log_mix(theta,
                      lognormal_lpdf(y[n] | mu[1], sigma[1]),
                      lognormal_lpdf(y[n] | mu[2], sigma[2]));
}                     
