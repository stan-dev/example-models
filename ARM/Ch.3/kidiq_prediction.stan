data {
  int<lower=0> N;
  vector<lower=0, upper=200>[N] kid_score;
  vector<lower=0, upper=200>[N] mom_iq;
  vector<lower=0, upper=1>[N] mom_hs;
  real<lower=0, upper=1> mom_hs_new;           // for prediction
  real<lower=0, upper=200> mom_iq_new;
}
transformed data {
  matrix[N, 2] x = [mom_hs', mom_iq']';
  row_vector[2] x_new = [mom_hs_new, mom_iq_new];
}
parameters {
  real alpha;
  vector[2] beta;
  real<lower=0> sigma;
}
model {
  sigma ~ cauchy(0, 2.5);
  kid_score ~ normal_id_glm(x, alpha, beta, sigma);
}
generated quantities {       // prediction
  real kid_score_pred = normal_rng(alpha + x_new * beta, sigma);
}
