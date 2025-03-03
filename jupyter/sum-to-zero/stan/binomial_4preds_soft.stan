// multi-level model for binomial data with 4 categorical predictors.
data {
  int<lower=1> N; // number of strata
  int<lower=1> N_age;
  int<lower=1> N_eth;
  int<lower=1> N_edu;

  array[N] int<lower=0> pos_tests;
  array[N] int<lower=0> tests;
  array[N] int<lower=1, upper=2> sex;
  array[N] int<lower=1, upper=N_age> age;
  array[N] int<lower=1, upper=N_eth> eth;
  array[N] int<lower=1, upper=N_edu> edu;

  // hyperparameters
  real<lower=0, upper=1> sens;
  real<lower=0, upper=1> spec;
}
transformed data {
  real mean_sex = mean(sex);
  vector[N] sex_c = to_vector(sex) - mean_sex;
}
parameters {
  real beta_0;
  real beta_sex;
  real<lower=0> sigma_age, sigma_eth, sigma_edu;
  vector[N_age] beta_age;
  vector[N_eth] beta_eth;
  vector[N_edu] beta_edu;
}
transformed parameters {
  // non-standard link function
  vector[N] p =  inv_logit(beta_0 + beta_sex * sex_c + beta_age[age]
			   + beta_eth[eth] +  beta_edu[edu]);
  vector[N] p_sample = p * sens + (1 - p) * (1 - spec);
}
model {
  pos_tests ~ binomial(tests, p_sample);  // likelihood

  // priors
  beta_0 ~ normal(0, 2.5);
  beta_sex ~ std_normal();
  // centered parameterization
  beta_age ~ normal(0, sigma_age);
  beta_eth ~ normal(0, sigma_eth);
  beta_edu ~ normal(0, sigma_edu);
  sigma_eth ~ std_normal();
  sigma_age ~ std_normal();
  sigma_edu ~ std_normal();
  // soft sum-to-zero constraint
  sum(beta_age) ~ normal(0, 0.001 * N_age); 
  sum(beta_eth) ~ normal(0, 0.001 * N_eth);
  sum(beta_edu) ~ normal(0, 0.001 * N_edu);
}
generated quantities {
  real beta_intercept = beta_0 - mean_sex * beta_sex;
  array[N] int<lower=0>y_rep = binomial_rng(tests, p_sample);
}
