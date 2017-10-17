// Sample from a Gaussian process using Stan's built in exponentiated quadratic
// covariance function and the Cholesky parameterization of a latent-variable
// Gaussian process.  
// Fixed kernel hyperparameters: rho=1, alpha=1, sigma=sqrt(0.1)

data {
  int<lower=1> N;
  real x[N];
}
transformed data {
  vector[N] mu = rep_vector(0, N);
  matrix[N, N] L;
  {
    matrix[N, N] K = cov_exp_quad(x, 1.0, 1.0);
    for (n in 1:N) 
      K[n, n] = K[n, n] + 0.1;
    L = cholesky_decompose(K);
  }
}
parameters {
  vector[N] eta;
}
model {
  eta ~ normal(0, 1);
}
generated quantities {
  vector[N] y;
  y = mu + L * eta;
}
