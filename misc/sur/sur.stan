data {
  int<lower=1> K;
  int<lower=1> J;
  int<lower=0> N;
  array[N] vector[J] x;
  array[N] vector[K] y;
}
parameters {
  matrix[K, J] beta;
  cholesky_factor_corr[K] L_Omega;
  vector<lower=0>[K] L_sigma;
}
model {
  array[N] vector[K] mu;
  matrix[K, K] L_Sigma;
  for (n in 1 : N) {
    mu[n] = beta * x[n];
  }
  L_Sigma = diag_post_multiply(L_Omega, L_sigma);
  
  to_vector(beta) ~ normal(0, 5);
  L_Omega ~ lkj_corr_cholesky(4);
  L_sigma ~ cauchy(0, 2.5);
  
  y ~ multi_normal_cholesky(mu, L_Sigma);
}
