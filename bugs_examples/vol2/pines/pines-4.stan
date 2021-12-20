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
  print("pines-4");
}
parameters {
  real alpha;
  real beta;
  real gamma;
  real delta;
  vector<lower=0>[2] tau;
  real<lower=0, upper=1> lambda;
}
transformed parameters {
  vector<lower=0>[2] sigma;
  for (i in 1 : 2) {
    sigma[i] = 1 / sqrt(tau[i]);
  }
}
model {
  vector[2] log_p;
  
  log_p[1] = log(lambda)
             + normal_lpdf(y_std | alpha + beta * x_std, sigma[1])
             + normal_lpdf(alpha | 0, sqrt(1e6))
             + normal_lpdf(beta | 0, sqrt(1e4))
             + gamma_lpdf(tau[1] | 0.0001, 0.0001);
  // + normal_log(gamma,0,1/sqrt(400))
  // + normal_log(delta,1,1/sqrt(400))
  // + gamma_log(tau[2],46,4.5);
  
  log_p[2] = log1m(lambda)
             + normal_lpdf(y_std | gamma + delta * z_std, sigma[2])
             + normal_lpdf(gamma | 0, sqrt(1e6))
             + normal_lpdf(delta | 0, sqrt(1e4))
             + gamma_lpdf(tau[2] | 0.0001, 0.0001);
  // + normal_log(alpha,0,1/sqrt(400))
  // + normal_log(beta,1,1/sqrt(400))
  // + gamma_log(tau[1],46,4.5);
  
  target += log_sum_exp(log_p);
}
