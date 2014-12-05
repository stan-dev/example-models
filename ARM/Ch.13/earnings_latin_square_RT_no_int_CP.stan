data {
  int<lower=0> N; 
  int<lower=0> n_eth; 
  int<lower=0> n_age; 
  int<lower=0> n_eth_age; 
  int<lower=1,upper=n_eth> eth[N];
  int<lower=1,upper=n_age> age[N];
  int<lower=1,upper=n_eth_age> eth_age[N];
  row_vector[2] x[N];
  vector[N] y;
} 
parameters {
  vector[2] mu_01;
  matrix[n_eth,2] gamma_01_eth;
  vector<lower=0>[2] sigma_gamma_01_eth;
  real<lower=0> sigma_y;
  cholesky_factor_corr[2] Omega_eth;
}
transformed parameters {
  vector[N] y_hat;

  for (i in 1:N)
    y_hat[i] <- x[i] * gamma_01_eth[eth[i]]'; 
} 
model {
  mu_01 ~ normal(0,10);
  Omega_eth ~ lkj_corr_cholesky(2);
  sigma_gamma_01_eth ~ cauchy(1,5);

  for(j in 1:n_eth)
    gamma_01_eth[j] ~ multi_normal_cholesky(mu_01,diag_pre_multiply(sigma_gamma_01_eth,Omega_eth));

  sigma_y ~ cauchy(1,5);

  y ~ normal(y_hat, sigma_y);
}