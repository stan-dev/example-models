// generate sample from model priors, (before seeing any data)
data {
  int<lower=1> N; // number of strata
  int<lower=1> N_age;
  int<lower=1> N_eth;
  int<lower=1> N_edu;
  // omit observational data
}
transformed data {
  // scaling factors for marginal variances of sum_to_zero_vectors
  // https://discourse.mc-stan.org/t/zero-sum-vector-and-normal-distribution/38296
  real s_age = sqrt(N_age * inv(N_age - 1));
  real s_eth = sqrt(N_eth * inv(N_eth - 1));
  real s_edu = sqrt(N_edu * inv(N_edu - 1));
}
parameters {
  real beta_0;
  real beta_sex;
  real<lower=0> sigma_age, sigma_eth, sigma_edu;
  sum_to_zero_vector[N_age] beta_age;
  sum_to_zero_vector[N_eth] beta_eth;
  sum_to_zero_vector[N_edu] beta_edu;
}
model {
  // omit likelihood
  // priors
  beta_0 ~ normal(0, 2.5);
  beta_sex ~ std_normal();
  sigma_eth ~ std_normal();
  sigma_age ~ std_normal();
  sigma_edu ~ std_normal();

  // centered parameterization
  // scale normal priors on sum_to_zero_vectors
  beta_age ~ normal(0, s_age * sigma_age);
  beta_eth ~ normal(0, s_eth * sigma_eth);
  beta_edu ~ normal(0, s_edu * sigma_edu);
}
