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
  matrix[2,n_eth] gamma_01_eth_std;
  matrix[2,n_age] gamma_01_age_std;
  matrix[2,n_eth_age] gamma_01_eth_age_std;
  vector<lower=0>[2] sigma_gamma_01_eth;
  vector<lower=0>[2] sigma_gamma_01_age;
  vector<lower=0>[2] sigma_gamma_01_eth_age;
  real<lower=0> sigma_y;
  cholesky_factor_corr[2] Omega_eth;
  cholesky_factor_corr[2] Omega_age;
  cholesky_factor_corr[2] Omega_eth_age;
}
transformed parameters {
  matrix[n_eth,2] gamma_01_eth;
  matrix[n_age,2] gamma_01_age;
  matrix[n_eth_age,2] gamma_01_eth_age;
  vector[N] y_hat;

  gamma_01_eth <- (diag_pre_multiply(sigma_gamma_01_eth, Omega_eth) 
                   * gamma_01_eth_std)';

  gamma_01_age <- (diag_pre_multiply(sigma_gamma_01_age, Omega_age) 
                   * gamma_01_age_std)';
  gamma_01_eth_age <- (diag_pre_multiply(sigma_gamma_01_eth_age, Omega_eth_age) 
                   * gamma_01_eth_age_std)';

  for (i in 1:N)
    y_hat[i] <- x[i] * (mu_01 + gamma_01_eth[eth[i]]' 
                + gamma_01_age[age[i]]' + gamma_01_eth_age[eth_age[i]]');
} 
model {
  mu_01 ~ normal(0,10);

  Omega_eth ~ lkj_corr_cholesky(0.5);
  Omega_age ~ lkj_corr_cholesky(0.5);
  Omega_eth_age ~ lkj_corr_cholesky(0.5);

  sigma_gamma_01_eth ~ cauchy(1,5);
  sigma_gamma_01_age ~ cauchy(1,5);
  sigma_gamma_01_eth_age ~ cauchy(1,5);

  to_vector(gamma_01_eth_std) ~ normal(0,1);
  to_vector(gamma_01_age_std) ~ normal(0,1);
  to_vector(gamma_01_eth_age_std) ~ normal(0,1);

  sigma_y ~ cauchy(1,5);

  y ~ normal(y_hat, sigma_y);
}