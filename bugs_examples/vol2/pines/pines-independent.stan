data {
  int<lower=0> N;
  vector[N] y;
  vector[N] x;
  vector[N] z;
}
transformed data {
  vector[N] y_std;
  vector[N] x_std;
  vector[N] z_std;
  
  y_std = (y - mean(y)) / sd(y);
  x_std = (x - mean(x)) / sd(x);
  z_std = (z - mean(z)) / sd(z);
  print("pines_independant");
}
parameters {
  real alpha;
  real beta;
  real gamma;
  real delta;
  vector<lower=0>[2] sigma;
}
model {
  array[2] vector[N] mu;
  
  alpha ~ normal(0, 10);
  beta ~ normal(0, 5);
  mu[1] = alpha + beta * x_std;
  
  gamma ~ normal(0, 10);
  delta ~ normal(0, 5);
  mu[2] = gamma + delta * z_std;
  
  sigma ~ cauchy(0, 5);
  
  // estimate separately
  y_std ~ normal(mu[1], sigma[1]);
  y_std ~ normal(mu[2], sigma[2]);
}
