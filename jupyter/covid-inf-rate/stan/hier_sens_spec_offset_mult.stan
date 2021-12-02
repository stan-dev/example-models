/**
 * meta-analysis of disease infection rate studies
 * hierarchical model of sensitivity and specificity
 */
data {
  int y_sample; // for j = 1
  int n_sample; // for j = 1
  int J_spec;
  array[J_spec] int y_spec; // no samples for j > 1
  array[J_spec] int n_spec;
  int J_sens;
  array[J_sens] int y_sens; // no samples for j > 1
  array[J_sens] int n_sens;
}
parameters {
  real<lower=0, upper=1> p;
  real mu_logit_spec;
  real<lower=0> sigma_logit_spec;
  real mu_logit_sens;
  real<lower=0> sigma_logit_sens;
  vector<offset=mu_logit_spec, multiplier=sigma_logit_spec>[J_spec] logit_spec;
  vector<offset=mu_logit_sens, multiplier=sigma_logit_sens>[J_sens] logit_sens;
}
transformed parameters {
  vector<lower=0, upper=1>[J_spec] spec = inv_logit(logit_spec);
  vector<lower=0, upper=1>[J_sens] sens = inv_logit(logit_sens);
}
model {
  real p_sample = p * sens[1] + (1 - p) * (1 - spec[1]); // j = 1
  y_sample ~ binomial(n_sample, p_sample); // j = 1
  y_spec ~ binomial_logit(n_spec, logit_spec);
  y_sens ~ binomial_logit(n_sens, logit_sens);
  logit_spec ~ normal(mu_logit_spec, sigma_logit_spec);
  logit_sens ~ normal(mu_logit_sens, sigma_logit_sens);
  sigma_logit_spec ~ normal(0, 1);
  sigma_logit_sens ~ normal(0, 0.2);
}
