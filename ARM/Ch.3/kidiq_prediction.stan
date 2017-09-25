data {
  int<lower=0> N;
  vector<lower=0, upper=200>[N] kid_score;
  vector<lower=0, upper=200>[N] mom_iq;
  vector<lower=0, upper=1>[N] mom_hs;
  real<lower=0, upper=1> mom_hs_new;           // for prediction
  real<lower=0, upper=200> mom_iq_new;
}
parameters {
  vector[3] beta;
  real<lower=0> sigma;
}
model {
  sigma ~ cauchy(0, 2.5);
  kid_score ~ normal(beta[1] + beta[2] * mom_hs + beta[3] * mom_iq, sigma);
}
generated quantities {       // prediction
  real kid_score_pred;
  kid_score_pred = normal_rng(beta[1] + beta[2] * mom_hs_new
                               + beta[3] * mom_iq_new, sigma);
}
