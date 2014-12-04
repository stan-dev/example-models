data {
  int<lower=0> N; 
  int<lower=0> n_age; 
  int<lower=0> n_eth; 
  int<lower=0> n_eth_age; 
  int<lower=1,upper=n_age> age[N];
  int<lower=1,upper=n_eth> eth[N];
  int<lower=1,upper=n_eth_age> eth_age[N];
  row_vector[2] x[N];
  vector[N] y;
} 
parameters {
  vector[2] mu_01;
  matrix[2,n_eth] gamma_01_eth_std;
  matrix[2,n_age] gamma_01_age_std;
  matrix[2,n_eth_age] gamma_01_eth_age_std;
  real<lower=0> sigma_gamma_0_eth;
  real<lower=0> sigma_gamma_1_eth;
  real<lower=0> sigma_gamma_0_age;
  real<lower=0> sigma_gamma_1_age;
  real<lower=0> sigma_gamma_0_eth_age;
  real<lower=0> sigma_gamma_1_eth_age;
  real<lower=0> sigma_y;
  real<lower=-1,upper=1> rho_eth;
  real<lower=-1,upper=1> rho_age;
  real<lower=-1,upper=1> rho_eth_age;
//  real<lower=0,upper=1> rho_eth_raw;
//  real<lower=0,upper=1> rho_age_raw;
//  real<lower=0,upper=1> rho_eth_age_raw;
}
transformed parameters {
  matrix[n_eth,2] gamma_01_eth;
  matrix[n_age,2] gamma_01_age;
  matrix[n_eth_age,2] gamma_01_eth_age;
  vector[N] y_hat;
  matrix[2,2] cov_mat_eth;
  matrix[2,2] cov_mat_age;
  matrix[2,2] cov_mat_eth_age;
//  real<lower=-1,upper=1> rho_eth;
//  real<lower=-1,upper=1> rho_age;
//  real<lower=-1,upper=1> rho_eth_age;

  cov_mat_eth[1,1] <- pow(sigma_gamma_0_eth, 2);
  cov_mat_eth[2,2] <- pow(sigma_gamma_1_eth, 2);
  cov_mat_eth[2,1] <- rho_eth * sigma_gamma_0_eth * sigma_gamma_1_eth;
  cov_mat_eth[1,2] <- cov_mat_eth[2,1];

  cov_mat_age[1,1] <- pow(sigma_gamma_0_age, 2);
  cov_mat_age[2,2] <- pow(sigma_gamma_1_age, 2);
  cov_mat_age[2,1] <- rho_age * sigma_gamma_0_age * sigma_gamma_1_age;
  cov_mat_age[1,2] <- cov_mat_age[2,1];

  cov_mat_eth_age[1,1] <- pow(sigma_gamma_0_eth_age, 2);
  cov_mat_eth_age[2,2] <- pow(sigma_gamma_1_eth_age, 2);
  cov_mat_eth_age[2,1] <- rho_eth_age * sigma_gamma_0_eth_age 
                          * sigma_gamma_1_eth_age;
  cov_mat_eth_age[1,2] <- cov_mat_eth_age[2,1];

  gamma_01_eth <- (cholesky_decompose(cov_mat_eth) * gamma_01_eth_std)';
  gamma_01_age <- (cholesky_decompose(cov_mat_age) * gamma_01_age_std)';
  gamma_01_eth_age <- (cholesky_decompose(cov_mat_eth_age) * gamma_01_eth_age_std)';

//  rho_eth <- 2 * rho_eth_raw - 1;
//  rho_age <- 2 * rho_age_raw - 1;
//  rho_eth_age <- 2 * rho_eth_age_raw - 1;

  for (i in 1:N)
    y_hat[i] <- x[i] * (mu_01 + gamma_01_eth[eth[i]]' 
                + gamma_01_age[age[i]]' + gamma_01_eth_age[eth_age[i]]');
} 
model {
//  rho_eth_raw ~ beta(3.0/4.0,3.0/4.0);
//  rho_age_raw ~ beta(3.0/4.0,3.0/4.0);
//  rho_eth_age_raw ~ beta(3.0/4.0,3.0/4.0);

  mu_01 ~ normal(0,10);

  sigma_gamma_0_age ~ cauchy(1,5);
  sigma_gamma_1_age ~ cauchy(1,5);

  sigma_gamma_0_eth ~ cauchy(1,5);
  sigma_gamma_1_eth ~ cauchy(1,5);

  sigma_gamma_0_eth_age ~ cauchy(1,5);
  sigma_gamma_1_eth_age ~ cauchy(1,5);

  to_vector(gamma_01_eth) ~ normal(0,1);
  to_vector(gamma_01_age) ~ normal(0,1);
  to_vector(gamma_01_eth_age) ~ normal(0,1);

  sigma_y ~ cauchy(1,5);

  y ~ normal(y_hat, sigma_y);
}