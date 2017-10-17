// Fit the hyperparameters of a latent-variable Gaussian process with an
// exponentiated quadratic kernel and a Poisson likelihood

data {
  int<lower=1> N;
  real x[N];
  int<lower=0> y[N];
}
transformed data {
  real delta = 1e-9;
}
parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  real a;
  vector[N] eta;
}
model {
  vector[N] f;
  {
    matrix[N, N] L_K;
    matrix[N, N] K = cov_exp_quad(x, alpha, rho);
  
    // diagonal elements
    for (n in 1:N)
      K[n, n] = K[n, n] + delta;
    
    L_K = cholesky_decompose(K);
    f = L_K * eta;
  }
  
  rho ~ inv_gamma(5, 5);
  alpha ~ normal(0, 1);
  a ~ normal(0, 1);
  eta ~ normal(0, 1);

  y ~ poisson_log(a + f);
}
