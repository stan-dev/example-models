// Sample from a multi-output Gaussian process using Stan's built in
// exponentiated quadratic covariance function.
// Fixed kernel hyperparameters: rho=1, alpha=1, sigma=sqrt(0.1)
data {
  int<lower=1> N;
  int<lower=1> D;
  real x[N];
}
transformed data {
  real delta = 1e-9;
  real<lower=0> rho = 1.0;
  vector<lower=0>[D] alpha = rep_vector(1.0, D);
  real<lower=0> sigma = sqrt(0.1);
}
parameters {
  cholesky_factor_corr[D] L_Omega;
  matrix[N, D] eta;
  matrix[N, D] y;
}
model {
  matrix[N, D] f;
  {
    matrix[N, N] K = cov_exp_quad(x, 1.0, rho);
    matrix[N, N] L_K;

    // diagonal elements
    for (n in 1:N)
      K[n, n] = K[n, n] + delta; 

    L_K = cholesky_decompose(K);
    f = L_K * eta
        * diag_pre_multiply(alpha, L_Omega)';
  }

  L_Omega ~ lkj_corr_cholesky(3);
  to_vector(eta) ~ normal(0, 1);

  to_vector(y) ~ normal(to_vector(f), sigma);
}
generated quantities {
  matrix[D, D] Omega;
  Omega = L_Omega * L_Omega';
}
