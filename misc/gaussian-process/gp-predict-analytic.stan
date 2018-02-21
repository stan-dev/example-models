// Fit the hyperparameters of a Gaussian process with an exponentiated
// quadratic kernel and analytically predict out-of-sample observations
functions{
    vector gp_pred_rng(real[] x2,
                       vector y1,
                       real[] x1,
                       real alpha,
                       real rho,
                       real sigma,
                       real delta) {
      int N1 = rows(y1);
      int N2 = size(x2);
      vector[N2] f2;
  
      {
        matrix[N1, N1] L_K;
        vector[N1] L_div_y1;
        matrix[N1, N2] k_x1_x2;
        matrix[N1, N2] L_div_k_x1_x2;
        vector[N2] f2_mu;
        matrix[N2, N2] cov_f2;
        matrix[N2, N2] diag_delta;
        matrix[N1, N1] K;
        K = cov_exp_quad(x1, alpha, rho);
        for (n in 1:N1)
          K[n, n] = K[n,n] + square(sigma);
        L_K = cholesky_decompose(K);
        L_div_y1 = mdivide_left_tri_low(L_K, y1);
        k_x1_x2 = cov_exp_quad(x1, x2, alpha, rho);
        L_div_k_x1_x2 = mdivide_left_tri_low(L_K, k_x1_x2);
        f2_mu = L_div_k_x1_x2' * L_div_y1;
        cov_f2 = cov_exp_quad(x2, alpha, rho) - L_div_k_x1_x2' * L_div_k_x1_x2;
        diag_delta = diag_matrix(rep_vector(delta,N2));
  
        f2 = multi_normal_rng(f2_mu, cov_f2 + diag_delta);
      }
    return f2;
  }
}
data {
  int<lower=1> N1;
  real x1[N1];
  vector[N1] y1;
  int<lower=1> N2;
  real x2[N2];
}
transformed data {
  vector[N1] mu = rep_vector(0, N1);
  real delta = 1e-9;
}
parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
}
model {
  matrix[N1, N1] L_K;
  {
    matrix[N1, N1] K = cov_exp_quad(x1, alpha, rho);
    real sq_sigma = square(sigma);
  
    // diagonal elements
    for (n1 in 1:N1)
      K[n1, n1] = K[n1, n1] + sq_sigma;
    
    L_K = cholesky_decompose(K);
  }
  
  rho ~ inv_gamma(5, 5);
  alpha ~ normal(0, 1);
  sigma ~ normal(0, 1);

  y1 ~ multi_normal_cholesky(mu, L_K);
}
generated quantities {
  vector[N2] f2;
  vector[N2] y2;

  f2 = gp_pred_rng(x2, y1, x1, alpha, rho, sigma, delta);
  for (n2 in 1:N2)
    y2[n2] = normal_rng(f2[n2], sigma);
}
