// multi-level model for binomial data with 4 categorical predictors.
data {
  int<lower=1> N; // number of strata
  int<lower=1> N_age;
  int<lower=1> N_eth;
  int<lower=1> N_edu;

  // hyperparameters
  real<lower=0, upper=1> sens;
  real<lower=0, upper=1> spec;
}
parameters {
  real beta_0;
  real beta_sex_raw;
  real<lower=0> sigma_age, sigma_eth, sigma_edu;
  vector[N_age - 1] beta_age_raw;
  vector[N_eth - 1] beta_eth_raw;
  vector[N_edu - 1] beta_edu_raw;
}
transformed parameters {
  vector[2] beta_sex = [beta_sex_raw, -beta_sex_raw]';

  vector[N_age] beta_age = append_row(beta_age_raw, -sum(beta_age_raw));
  vector[N_eth] beta_eth = append_row(beta_eth_raw, -sum(beta_eth_raw));
  vector[N_edu] beta_edu = append_row(beta_edu_raw, -sum(beta_edu_raw));
}
model {
  // priors
  beta_0 ~ normal(0, 2.5);
  beta_sex ~ std_normal();
  // centered parameterization
  beta_age_raw ~ normal(0, sigma_age);
  beta_eth_raw ~ normal(0, sigma_eth);
  beta_edu_raw ~ normal(0, sigma_edu);
  sigma_eth ~ std_normal();
  sigma_age ~ std_normal();
  sigma_edu ~ std_normal();
}
