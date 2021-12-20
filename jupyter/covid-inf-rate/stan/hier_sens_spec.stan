/**
 * meta-analysis of disease infection rate studies
 * hierarchical model of sensitivity and specificity
 * non-centered parameterization
 */
data {
  int y_sample;
  int n_sample;
  int J_spec;
  array[J_spec] int y_spec;
  array[J_spec] int n_spec;
  int J_sens;
  array[J_sens] int y_sens;
  array[J_sens] int n_sens;
}
parameters {
  real<lower=0, upper=1> p;
  real mu_spec;
  real<lower=0> sigma_spec;
  real mu_sens;
  real<lower=0> sigma_sens;
  vector[J_spec] e_spec;
  vector[J_sens] e_sens;
}
transformed parameters {
  vector[J_spec] spec;
  vector[J_sens] sens;
  spec = inv_logit(mu_spec + sigma_spec * e_spec);
  sens = inv_logit(mu_sens + sigma_sens * e_sens);
}
model {
  real p_sample;
  p_sample = p * sens[1] + (1 - p) * (1 - spec[1]);
  y_sample ~ binomial(n_sample, p_sample);
  y_spec ~ binomial(n_spec, spec);
  y_sens ~ binomial(n_sens, sens);
  e_spec ~ normal(0, 1);
  e_sens ~ normal(0, 1);
  sigma_spec ~ normal(0, 1);
  sigma_sens ~ normal(0, 0.2);
}
