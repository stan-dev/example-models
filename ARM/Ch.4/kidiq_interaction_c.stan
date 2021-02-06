data {
  int<lower=0> N;
  vector[N] kid_score;
  vector[N] mom_hs;
  vector[N] mom_iq;
}
transformed data {           // centered predictors
  vector[N] c_mom_hs = mom_hs - mean(mom_hs);
  vector[N] c_mom_iq = mom_iq - mean(mom_iq);
  vector[N] inter = c_mom_hs .* c_mom_iq;
  matrix[N,3] x = [c_mom_hs', c_mom_iq', inter']';
}
parameters {
  real alpha;
  vector[3] beta;
  real<lower=0> sigma;
}
model {
  kid_score ~ normal_id_glm(x, alpha, beta, sigma);
}
