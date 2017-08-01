// Fit a Gaussian process's hyperparameters
// for squared exponential prior
data {
  int<lower=1> N;
  int<lower=1> M;
  real x[N];
  matrix[N, M] y;
}
transformed data {
  real delta = 1e-9;
}
parameters {
  real<lower=0> rho;
  vector<lower=0>[M] alpha;
  real<lower=0> sigma;
  cholesky_factor_corr[M] L_Omega;
  matrix[N, M] eta;
}
model {
  matrix[N, M] f;
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

  rho ~ gamma(4, 4);
  alpha ~ normal(0, 1);
  sigma ~ normal(0, 1);
  L_Omega ~ lkj_corr_cholesky(3);
  to_vector(eta) ~ normal(0, 1);

  to_vector(y) ~ normal(to_vector(f), sigma);
}
