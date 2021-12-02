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
  vector[2] log_py;
  
  for (i in 1 : 2) {
    sigma[i] = 1 / sqrt(tau[i]);
  }
  
  log_py[1] = log(lambda) + log(0.9995)
              + normal_lpdf(y_std | alpha + beta * x_std, sigma[1])
              + normal_lpdf(alpha | 0, sqrt(1e6))
              + normal_lpdf(beta | 0, sqrt(1e4))
              + gamma_lpdf(tau[1] | 0.0001, 0.0001)
              // pseudo-priors
              + normal_lpdf(gamma | 0, sqrt(1 / 400.0))
              + normal_lpdf(delta | 1, sqrt(1 / 400.0))
              + gamma_lpdf(tau[2] | 46, 4.5);
  
  log_py[2] = log(lambda) + log1m(0.0005)
              + normal_lpdf(y_std | gamma + delta * z_std, sigma[2])
              + normal_lpdf(gamma | 0, sqrt(1e6))
              + normal_lpdf(delta | 0, sqrt(1e4))
              + gamma_lpdf(tau[2] | 0.0001, 0.0001)
              // pseudo-priors
              + normal_lpdf(alpha | 0, sqrt(1 / 256.0))
              + normal_lpdf(beta | 1, sqrt(1 / 256.0))
              + gamma_lpdf(tau[1] | 30, 4.5);
}
model {
  target += log_sum_exp(log_py);
}
generated quantities {
  real pM2;
  pM2 = bernoulli_rng(softmax(log_py)[2]);
}
/*
 * p(lambda) = 1
 *
 * p(y|lambda,zeta1,zeta2) 
 *   = lambda * p1(y|zeta1) + (1 - lambda) * p2(y|zeta2)
 *
 * log p(y|lambda,zeta1,zeta2) 
 *   = log_sum_exp(log(lambda) + log(p1(y|zeta1)),
 *                 log1m(lambda) + log(p2(y|zeta2));
 *
 * p(lambda|y,zeta1,zeta2) propto p(y|lambda,zeta1,zet2) p(lambda)
 *
 */