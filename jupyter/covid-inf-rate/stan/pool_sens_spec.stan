/**
 * meta-analysis over studies of disease infection rate 
 * for variation in sensitivity and specificity
 * pool data from all studies
 */
data {
  int y_sample;
  int n_sample;
  int y_spec;
  int n_spec;
  int y_sens;
  int n_sens;
}
parameters {
  real<lower = 0, upper = 1> p;
  real<lower = 0, upper = 1> spec;
  real<lower = 0, upper = 1> sens;
}
model {
  real p_sample;
  p_sample = p * sens + (1 - p) * (1 - spec);
  y_sample ~ binomial(n_sample, p_sample);
  y_spec ~ binomial(n_spec, spec);
  y_sens ~ binomial(n_sens, sens);
}
