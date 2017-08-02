// Sample from a Gaussian process using Stan's built in exponentiated quadratic
// covariance function.
// Fixed kernel hyperparameters: rho=1, alpha=1, sigma=sqrt(0.1)

data {
  int<lower=1> N;
  real x[N];
}
transformed data {
  matrix[N, N] K = cov_exp_quad(x, 1.0, 1.0);
  vector[N] mu = rep_vector(0, N);
  for (n in 1:N) 
    K[n, n] = K[n, n] + 0.1;
}
parameters {
  vector[N] y;
}
model {
  y ~ multi_normal(mu, K);
}
