data {
  int<lower=0> N;
  vector[N] kid_score;
  vector[N] mom_hs;
  vector[N] mom_iq;
}
transformed data {           // centering on reference points
  vector[N] c2_mom_hs = mom_hs - 0.5;
  vector[N] c2_mom_iq = mom_iq - 100;
  vector[N] inter = c2_mom_hs .* c2_mom_iq;
  matrix[N,3] x = [c2_mom_hs', c2_mom_iq', inter']';
}
parameters {
  real alpha;
  vector[3] beta;
  real<lower=0> sigma;
}
model {
  kid_score ~ normal_id_glm(x, alpha, beta, sigma);
}
