data {
  int<lower=1> K;
  int<lower=1> N;
  vector[N] y;
}
parameters {
  simplex[K] theta;
  simplex[K] mu_prop;
  real mu_loc;
  real<lower=0> mu_scale;
  vector<lower=0>[K] sigma;
}
transformed parameters {
  ordered[K] mu = mu_loc + mu_scale * cumulative_sum(mu_prop);
}
model {
  array[N] vector[K] ps;
  
  // prior
  mu_loc ~ cauchy(0, 5);
  mu_scale ~ cauchy(0, 5);
  sigma ~ cauchy(0, 5);
  
  // likelihood
  for (n in 1 : N) {
    for (k in 1 : K) {
      ps[n][k] = normal_lupdf(y[n] | mu[k], sigma[k]);
    }
  }
  
  target += log_mix(theta, ps);
}
