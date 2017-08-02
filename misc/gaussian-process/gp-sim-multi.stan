// Sample from Gaussian process
// Fixed covar function: eta_sq=1, rho_sq=1, sigma_sq=0.1

data {
  int<lower=1> N;
  int<lower=1> D;
  vector[D] x[N];
}
transformed data {
  matrix[N, N] K = cov_exp_quad(x, 1.0, 1.0);
  vector[N] mu = rep_vector(0, N);
  for (i in 1:N) 
    K[i, i] = K[i, i] + 0.1;
}
parameters {
  vector[N] y;
}
model {
  y ~ multi_normal(mu, K);
}
