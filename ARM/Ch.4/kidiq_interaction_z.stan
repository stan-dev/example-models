data {
  int<lower=0> N;
  vector[N] kid_score;
  vector[N] mom_hs;
  vector[N] mom_iq;
}
transformed data {           // standardizing
  vector[N] z_mom_hs = (mom_hs - mean(mom_hs)) / (2 * sd(mom_hs));
  vector[N] z_mom_iq = (mom_iq - mean(mom_iq)) / (2 * sd(mom_iq));
  vector[N] inter= z_mom_hs .* z_mom_iq;
  matrix[N,3] x = [z_mom_hs', z_mom_iq', inter']';
}
parameters {
  real alpha;
  vector[3] beta;
  real<lower=0> sigma;
}
model {
  kid_score ~ normal_id_glm(x, alpha, beta, sigma);
}
