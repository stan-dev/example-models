data {
  int<lower=0> N;
  vector[N] y;
  vector[N] x;
  vector[N] z;
}
transformed data {
  vector[N] y_std = (y - mean(y)) / sd(y);
  vector[N] x_std = (x - mean(x)) / sd(x);
  vector[N] z_std = (z - mean(z)) / sd(z);
  print("pines-3");
}
parameters {
  real alpha;
  real beta;
  real gamma;
  real delta;
  vector<lower=0>[2] tau;
}
transformed parameters {
  vector<lower=0>[2] sigma;
  for (i in 1 : 2) {
    sigma[i] = 1 / sqrt(tau[i]);
  }
}
model {
  alpha ~ normal(0, sqrt(1e6));
  beta ~ normal(0, sqrt(1e4));
  tau[1] ~ gamma(0.0001, 0.0001);
  
  y_std ~ normal(alpha + beta * x_std, sigma[1]);
  
  gamma ~ normal(0, 10);
  delta ~ normal(0, 5);
  tau[2] ~ gamma(0.0001, 0.0001);
  y_std ~ normal(gamma + delta * z_std, sigma[2]);
}
generated quantities {
  vector[2] log_py;
  real lambda;
  
  log_py[1] = log(0.9995)
              + normal_lpdf(y_std | alpha + beta * x_std, sigma[1])
              + normal_lpdf(alpha | 0, sqrt(1e6))
              + normal_lpdf(beta | 0, sqrt(1e4))
              + gamma_lpdf(tau[1] | 0.0001, 0.0001)
              + normal_lpdf(gamma | 0.0001, 0.0001);
  
  log_py[2] = log(0.0005)
              + normal_lpdf(y_std | gamma + delta * z_std, sigma[2])
              + normal_lpdf(gamma | 0, sqrt(1e6))
              + normal_lpdf(delta | 0, sqrt(1e4))
              + gamma_lpdf(tau[2] | 0.0001, 0.0001);
  
  lambda = softmax(log_py)[1];
}
