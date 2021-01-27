transformed data {
  real<lower=0,upper=1> theta = 0.25;
  vector[2] mu = [0, 4]';
  vector<lower=0>[2] sigma = [0.5, 3]';
}
parameters {
  real y;
}
model {
  target += log_mix(theta,
                    normal_lpdf(y | mu[1], sigma[1]),
                    normal_lpdf(y | mu[2], sigma[2]))
}
